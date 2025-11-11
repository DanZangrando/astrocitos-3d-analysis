import sys
import os
import json
import subprocess
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
import pandas as pd
import altair as alt
from ui.sidebar import render_sidebar
from ui import pipeline, runner
from ui.utils import detect_group

st.title("📐 Esqueletización y Análisis Sholl 2D")
st.markdown("""
Este paso realiza un **pipeline unificado 2D**:
1. **Proyección** de máscaras 3D y señal GFAP (max projection)
2. **Territorios Voronoi** con zona de exclusión para separar astrocitos cercanos
3. **Esqueletización 2D** por territorio con conexión opcional de fragmentos
4. **Análisis de Sholl 2D** nativo integrado usando SKAN

**Ventajas del enfoque 2D:**
- Resolución XY completa (0.38 µm)
- Sholl 2D más preciso y eficiente
- Territorios astrocitarios bien definidos
""")

render_sidebar(show_calibration=True)

root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
calib_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer_2d.py"

# --- Utilidades ---
def _read_global_calibration():
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())
        except Exception:
            pass
    return {}

def _save_global_calibration(cal_data: dict):
    calib_path.parent.mkdir(parents=True, exist_ok=True)
    calib_path.write_text(json.dumps(cal_data, indent=2))

def get_output_dir_for_image(img_path: Path) -> Path:
    base_name = img_path.stem
    out_dir = root / "data" / "processed" / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir

# Lectura global de parámetros
glob_calib = _read_global_calibration()

# --- Selección de imagen ---
st.markdown("---")
st.subheader("1️⃣ Seleccionar Preparado")

files = list(runner.list_raw_images(root))
if not files:
    st.warning("No se encontraron imágenes .tif/.tiff en `data/raw/`")
    st.stop()

# Agrupación por grupo (CTL vs Hipoxia)
file_groups = {}
for f in files:
    g = detect_group(f, root)
    file_groups.setdefault(g, []).append(f)

selected_group = st.selectbox(
    "Filtrar por grupo",
    options=["Todos"] + sorted(file_groups.keys()),
    index=0
)

if selected_group == "Todos":
    display_files = files
else:
    display_files = file_groups.get(selected_group, [])

file_names = [f.stem for f in display_files]
selected_name = st.selectbox("Preparado", options=file_names, index=0)
img_path = next((f for f in display_files if f.stem == selected_name), None)

if img_path is None:
    st.error("No se pudo encontrar la imagen seleccionada")
    st.stop()

st.info(f"📁 **Preparado:** `{img_path.stem}`  \n**Grupo:** {detect_group(img_path, root)}")

out_dir = get_output_dir_for_image(img_path)

# Verificar prerequisitos
mask_path = out_dir / "04_final_astrocytes_mask.tif"
if not mask_path.exists():
    st.error(f"❌ Falta `04_final_astrocytes_mask.tif`. Completá el **Paso 03** primero.")
    st.stop()

# --- Parámetros del Pipeline 2D ---
st.markdown("---")
st.subheader("2️⃣ Parámetros del Pipeline 2D")

with st.expander("⚙️ Configuración de Proyección y Territorios", expanded=True):
    col1, col2 = st.columns(2)
    
    with col1:
        proj_method = st.selectbox(
            "Método de proyección 3D→2D",
            options=["max", "mean", "sum"],
            index=0,
            help="Max projection conserva las estructuras más brillantes"
        )
        
        exclusion_gap = st.number_input(
            "Gap de exclusión Voronoi (µm)",
            value=float(glob_calib.get("TERRITORY_EXCLUSION_UM", 1.0)),
            min_value=0.0,
            max_value=10.0,
            step=0.5,
            help="Distancia de separación entre territorios adyacentes"
        )
    
    with col2:
        connect_fragments = st.checkbox(
            "Conectar fragmentos cercanos",
            value=bool(glob_calib.get("CONNECT_SKELETON_FRAGMENTS", True)),
            help="Une procesos astrocitarios fragmentados"
        )
        
        if connect_fragments:
            connection_radius = st.number_input(
                "Radio de conexión (µm)",
                value=float(glob_calib.get("CONNECTION_RADIUS_UM", 0.5)),
                min_value=0.1,
                max_value=5.0,
                step=0.1,
                help="Distancia máxima para conectar fragmentos"
            )
        else:
            connection_radius = 0.5

with st.expander("🎯 Dominio de Análisis", expanded=True):
    st.markdown("""
    **Radio máximo desde núcleo:** Define hasta qué distancia del soma se considera 
    que un proceso GFAP pertenece al astrocito. Esto elimina señal de fondo lejana 
    y mejora la precisión del análisis.
    """)
    
    max_radius_nucleus = st.number_input(
        "Radio máximo desde núcleo (µm)",
        value=float(glob_calib.get("MAX_RADIUS_FROM_NUCLEUS_UM", 100.0)),
        min_value=10.0,
        max_value=200.0,
        step=10.0,
        help="Procesos más allá de este radio desde el centroide nuclear son ignorados. "
             "Valores típicos: 80-120 µm para astrocitos hipocampales."
    )
    
    st.caption("💡 **Tip:** Si los astrocitos tienen procesos muy largos (>100 µm), "
               "aumentar este valor. Si hay mucho ruido de fondo, reducirlo.")

with st.expander("📏 Parámetros de Sholl", expanded=False):
    col1, col2, col3 = st.columns(3)
    
    with col1:
        sholl_min = st.number_input(
            "Radio mínimo (µm)",
            value=float(glob_calib.get("SHOLL_MIN_RADIUS_UM", 5.0)),
            min_value=0.0,
            step=1.0
        )
    
    with col2:
        sholl_max = st.number_input(
            "Radio máximo (µm)",
            value=float(glob_calib.get("SHOLL_MAX_RADIUS_UM", 100.0)),
            min_value=10.0,
            step=10.0
        )
    
    with col3:
        sholl_step = st.number_input(
            "Paso (µm)",
            value=float(glob_calib.get("SHOLL_STEP_UM", 2.0)),
            min_value=0.5,
            step=0.5
        )

# Botón para guardar configuración global
if st.button("💾 Guardar configuración global"):
    glob_calib.update({
        "PROJECTION_2D_METHOD": proj_method,
        "TERRITORY_EXCLUSION_UM": exclusion_gap,
        "CONNECT_SKELETON_FRAGMENTS": connect_fragments,
        "CONNECTION_RADIUS_UM": connection_radius,
        "MAX_RADIUS_FROM_NUCLEUS_UM": max_radius_nucleus,
        "SHOLL_MIN_RADIUS_UM": sholl_min,
        "SHOLL_MAX_RADIUS_UM": sholl_max,
        "SHOLL_STEP_UM": sholl_step
    })
    _save_global_calibration(glob_calib)
    st.success("✅ Configuración guardada globalmente")

# --- Ejecutar Pipeline ---
st.markdown("---")
st.subheader("3️⃣ Ejecutar Pipeline Unificado")

if st.button("🚀 Ejecutar Esqueletización + Sholl 2D", type="primary"):
    with st.spinner("Ejecutando pipeline unificado..."):
        try:
            # Configurar parámetros para esta ejecución
            cal = glob_calib.copy()
            cal.update({
                "PROJECTION_2D_METHOD": proj_method,
                "TERRITORY_EXCLUSION_UM": exclusion_gap,
                "CONNECT_SKELETON_FRAGMENTS": connect_fragments,
                "CONNECTION_RADIUS_UM": connection_radius,
                "MAX_RADIUS_FROM_NUCLEUS_UM": max_radius_nucleus,
                "SHOLL_MIN_RADIUS_UM": sholl_min,
                "SHOLL_MAX_RADIUS_UM": sholl_max,
                "SHOLL_STEP_UM": sholl_step
            })
            
            # Ejecutar pipeline unificado
            df_skeleton, df_sholl = pipeline.run_unified_2d_skeleton_and_sholl(
                img_path=img_path,
                mask_path=mask_path,
                out_dir=out_dir,
                cal=cal
            )
            
            st.success(f"✅ Pipeline completado: {len(df_skeleton)} astrocitos procesados")
            
            # Mostrar métricas rápidas
            if not df_skeleton.empty:
                st.metric("Astrocitos analizados", len(df_skeleton))
                col1, col2, col3 = st.columns(3)
                col1.metric("Long. esqueleto promedio", f"{df_skeleton['skeleton_length_um'].mean():.1f} µm")
                col2.metric("Endpoints promedio", f"{df_skeleton['n_endpoints'].mean():.1f}")
                col3.metric("Junctions promedio", f"{df_skeleton['n_junctions'].mean():.1f}")
            
        except Exception as e:
            st.error(f"❌ Error: {e}")
            import traceback
            st.code(traceback.format_exc())

# --- Visualización de Resultados ---
st.markdown("---")
st.subheader("4️⃣ Resultados")

skeleton_summary_path = out_dir / "skeletons" / "summary.csv"
sholl_path = out_dir / "sholl_2d_native.csv"

if skeleton_summary_path.exists():
    df_skel = pd.read_csv(skeleton_summary_path)
    
    st.markdown("### 📊 Métricas de Esqueletización")
    st.dataframe(df_skel, use_container_width=True, height=300)
    
    # Métricas resumidas con tabs
    tab1, tab2, tab3 = st.tabs(["📏 Básicas", "🌀 Tortuosidad", "🔀 Complejidad"])
    
    with tab1:
        col1, col2, col3, col4 = st.columns(4)
        if 'total_branch_length_um' in df_skel.columns:
            col1.metric("Long. total promedio", f"{df_skel['total_branch_length_um'].mean():.1f} µm")
        if 'n_endpoints' in df_skel.columns:
            col2.metric("Endpoints promedio", f"{df_skel['n_endpoints'].mean():.1f}")
        if 'n_junctions' in df_skel.columns:
            col3.metric("Junctions promedio", f"{df_skel['n_junctions'].mean():.1f}")
        if 'n_branches' in df_skel.columns:
            col4.metric("Ramas promedio", f"{df_skel['n_branches'].mean():.1f}")
    
    with tab2:
        col1, col2, col3 = st.columns(3)
        if 'tortuosity_mean' in df_skel.columns:
            col1.metric("Tortuosidad media", f"{df_skel['tortuosity_mean'].mean():.3f}", 
                       help="1.0 = recto, >1.0 = sinuoso")
            col2.metric("Tortuosidad máxima", f"{df_skel['tortuosity_max'].mean():.3f}")
            col3.metric("Desv. std tortuosidad", f"{df_skel['tortuosity_std'].mean():.3f}",
                       help="Mayor valor = mayor heterogeneidad")
        
        # Gráfico de tortuosidad
        if 'tortuosity_mean' in df_skel.columns:
            chart_tort = alt.Chart(df_skel).mark_bar().encode(
                x=alt.X('tortuosity_mean:Q', bin=alt.Bin(maxbins=20), title='Tortuosidad Media'),
                y=alt.Y('count()', title='Frecuencia')
            ).properties(width=600, height=250, title='Distribución de Tortuosidad (1.0 = recto)')
            st.altair_chart(chart_tort, use_container_width=True)
    
    with tab3:
        col1, col2, col3 = st.columns(3)
        if 'ramification_index' in df_skel.columns:
            col1.metric("Índice de Ramificación", f"{df_skel['ramification_index'].mean():.2f}",
                       help="Ramas / Junctions - Mayor = más complejo")
        if 'termination_index' in df_skel.columns:
            mean_term = df_skel['termination_index'].mean()
            col2.metric("Índice de Terminalización", f"{mean_term:.2f}" if not np.isnan(mean_term) else "N/A",
                       help="Endpoints / Junctions - Balance ramificación/extensión")
        if 'branch_length_cv' in df_skel.columns:
            col3.metric("CV Longitud Ramas", f"{df_skel['branch_length_cv'].mean():.2f}",
                       help="Coef. variación - Mayor = más heterogeneidad")
        
        # Gráfico de complejidad
        if 'ramification_index' in df_skel.columns and 'tortuosity_mean' in df_skel.columns:
            chart_complex = alt.Chart(df_skel).mark_circle(size=100).encode(
                x=alt.X('ramification_index:Q', title='Índice de Ramificación'),
                y=alt.Y('tortuosity_mean:Q', title='Tortuosidad Media'),
                color=alt.Color('n_junctions:Q', title='N° Junctions', scale=alt.Scale(scheme='viridis')),
                tooltip=['label', 'ramification_index', 'tortuosity_mean', 'n_junctions']
            ).properties(width=600, height=300, title='Complejidad vs Tortuosidad')
            st.altair_chart(chart_complex, use_container_width=True)
    
    # Gráfico de distribución de longitudes (movido abajo)
    if 'skeleton_length_um' in df_skel.columns:
        st.markdown("#### Distribución de Longitudes Totales")
        chart = alt.Chart(df_skel).mark_bar().encode(
            x=alt.X('skeleton_length_um:Q', bin=True, title='Longitud de esqueleto (µm)'),
            y=alt.Y('count()', title='Frecuencia')
        ).properties(width=600, height=250, title='Longitud Total del Esqueleto')
        st.altair_chart(chart, use_container_width=True)
else:
    st.info("⏳ No hay métricas de esqueleto disponibles. Ejecutá el pipeline primero.")

if sholl_path.exists():
    df_sholl = pd.read_csv(sholl_path)
    
    st.markdown("### 📈 Análisis de Sholl")
    
    # Gráfico de Sholl por célula
    chart = alt.Chart(df_sholl).mark_line(point=True).encode(
        x=alt.X('radius_um:Q', title='Radio (µm)'),
        y=alt.Y('intersections:Q', title='Intersecciones'),
        color='label:N',
        tooltip=['label', 'radius_um', 'intersections']
    ).properties(width=700, height=400, title='Perfiles de Sholl por Astrocito')
    
    st.altair_chart(chart, use_container_width=True)
    
    # Tabla de resumen
    st.markdown("**Datos crudos de Sholl:**")
    st.dataframe(df_sholl, use_container_width=True, height=200)
else:
    st.info("⏳ No hay datos de Sholl disponibles. Ejecutá el pipeline primero.")

# --- Visualización en Napari ---
st.markdown("---")
st.subheader("5️⃣ Visualización en Napari")

open_napari_2d = st.button("👁️ Abrir visualizador 2D con anillos de Sholl")

if open_napari_2d:
    rings_path = out_dir / "sholl_rings_2d_native.json"
    if not rings_path.exists():
        st.warning("⚠️ No se encontró archivo de anillos. Ejecutá el pipeline primero.")
    else:
        try:
            y, x = float(glob_calib.get('y', 0.3)), float(glob_calib.get('x', 0.3))
            dapi_idx = int(glob_calib.get("DAPI_CHANNEL_INDEX", 0))
            gfap_idx = int(glob_calib.get("GFAP_CHANNEL_INDEX", 1))
            
            cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--y", str(y), "--x", str(x),
                   "--dapi_idx", str(dapi_idx), "--gfap_idx", str(gfap_idx)]
            
            # Archivos 2D
            gfap_proj = out_dir / "gfap_projection_2d.tif"
            gfap_2d = out_dir / "03_gfap_filtered_mask_2d.tif"
            final_2d = out_dir / "04_final_astrocytes_mask_2d.tif"
            skel_2d = out_dir / "05_skeleton_labels_2d.tif"
            
            if gfap_proj.exists():
                cmd += ["--gfap_proj", str(gfap_proj)]
            if gfap_2d.exists():
                cmd += ["--gfap_2d", str(gfap_2d)]
            if final_2d.exists():
                cmd += ["--final_2d", str(final_2d)]
            if skel_2d.exists():
                cmd += ["--skeleton_2d", str(skel_2d)]
            cmd += ["--rings", str(rings_path)]
            
            subprocess.Popen(cmd, env=os.environ.copy())
            st.info("✅ Visualizador 2D lanzado con anillos de Sholl")
        except Exception as e:
            st.error(f"❌ Error al lanzar Napari: {e}")

st.markdown("---")
st.caption("💡 **Tip:** Los resultados se guardan automáticamente en `data/processed/<preparado>/`")
