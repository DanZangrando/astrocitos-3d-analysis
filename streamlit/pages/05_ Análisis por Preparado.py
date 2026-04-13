import json
from pathlib import Path
import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile
from scipy import stats
from ui.sidebar import render_sidebar
from ui.plots import GROUP_SCALE
from ui import pipeline
import subprocess
import sys
import os

st.title("📊 Análisis Detallado por Preparado")
render_sidebar(show_calibration=True)

root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
napari_script = root / "streamlit" / "napari_viewer.py"

# --- Utilidades ---
def _detect_group(p: Path, root: Path) -> str:
    try:
        rel = str(p.relative_to(root)).lower()
    except Exception:
        rel = str(p).lower()
    if "/hip/" in rel:
        return "Hypoxia"
    if "/ctl/" in rel:
        return "CTL"
    return "CTL"

def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hypoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"

def _read_global_calibration():
    calib_path = root / "streamlit" / "calibration.json"
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())
        except Exception:
            pass
    return {}

files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

# --- 1. Selección de Imagen ---
st.markdown("### 🔍 Selección de Preparado")
labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files))), format_func=lambda i: labels[i], key="prep_selector")
img_path = files[idx]
out_dir = root / "data" / "processed" / img_path.stem
group = _detect_group(img_path, root)

col_header1, col_header2 = st.columns([3, 1])
with col_header1:
    st.markdown(_group_badge_html(group), unsafe_allow_html=True)
    st.caption(f"📁 Directorio: `{out_dir.relative_to(root)}`")
with col_header2:
    st.metric("Grupo", group)

# --- 2. Resumen del Pipeline (Conteos) ---
st.markdown("---")
st.markdown("### 📈 Resumen del Pipeline")

p_cp = out_dir / "02_cellpose_mask.tif"
p_nuc_metrics = out_dir / "03_nucleus_metrics.csv"
p_final = out_dir / "04_final_astrocytes_mask_2d.tif"  # Archivo 2D
p_skel_sum = out_dir / "skeletons" / "summary.csv"  # ¡Archivo correcto!
p_sholl_sum = out_dir / "sholl_summary.csv"

n_cellpose, n_gfap, n_final, n_skel, n_sholl = 0, 0, 0, 0, 0

# Contar núcleos
if p_nuc_metrics.exists():
    try:
        df_nuc_temp = pd.read_csv(p_nuc_metrics)
        n_cellpose = df_nuc_temp.shape[0]
        n_gfap = int(df_nuc_temp["is_astrocyte_candidate"].sum())
    except Exception: pass
elif p_cp.exists():
    try: n_cellpose = int(tifffile.imread(p_cp).max())
    except Exception: pass

# Contar astrocitos finales
if p_final.exists():
    try: n_final = int((np.unique(tifffile.imread(p_final)) > 0).sum())
    except Exception: pass

# Contar esqueletos analizados
if p_skel_sum.exists():
    try: n_skel = pd.read_csv(p_skel_sum).shape[0]
    except Exception: pass

# Contar análisis Sholl
if p_sholl_sum.exists():
    try: n_sholl = pd.read_csv(p_sholl_sum).shape[0]
    except Exception: pass

c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("🔵 Núcleos Totales", n_cellpose, help="Detectados por Cellpose (paso 02)")
c2.metric("🟢 Candidatos GFAP+", n_gfap, help="Núcleos con señal GFAP (paso 03)")
c3.metric("⭐ Astrocitos Finales", n_final, help="Tras filtros espaciales (paso 04)")
c4.metric("🌿 Esqueletos", n_skel, help="Con topología analizada (paso 04)")
c5.metric("🎯 Con Sholl", n_sholl, help="Con análisis de Sholl completo (paso 04)")



# --- 4. Métricas Topológicas de Esqueleto (Paso 04) ---
st.markdown("---")
st.markdown("### 🌿 Métricas Topológicas Clave (Paso 04)")

if not p_skel_sum.exists():
    st.info("No se encontró `skeletons/summary.csv`. Ejecutá el paso 04 primero.")
else:
    try:
        df_skel = pd.read_csv(p_skel_sum)
        
        st.markdown("#### Estadísticas Globales")
        col_s1, col_s2, col_s3 = st.columns(3)
        col_s1.metric("Longitud Total Media", f"{df_skel['total_branch_length_um'].mean():.1f} µm", help="Suma de longitudes de todas las ramas")
        col_s2.metric("Índice Ramificación Medio", f"{df_skel['ramification_index'].mean():.2f}", help="Ramas / Uniones")
        col_s3.metric("Terminaciones Promedio", f"{df_skel['n_endpoints'].mean():.1f}", help="Número de puntas finales")
        
        # Tabla filtrada
        st.markdown("#### Tabla de Métricas Clave")
        cols_to_show = ['label', 'total_branch_length_um', 'ramification_index', 'n_endpoints', 'n_branches', 'n_junctions']
        # Asegurar que existan
        cols_to_show = [c for c in cols_to_show if c in df_skel.columns]
        
        st.dataframe(
            df_skel[cols_to_show].round(3),
            use_container_width=True,
            height=300,
            column_config={
                "label": st.column_config.NumberColumn("ID", format="%d"),
                "total_branch_length_um": st.column_config.NumberColumn("Longitud (µm)", format="%.1f"),
                "ramification_index": st.column_config.NumberColumn("Ramificación", format="%.2f"),
                "n_endpoints": st.column_config.NumberColumn("Terminaciones", format="%d"),
                "n_branches": st.column_config.NumberColumn("N° Ramas", format="%d"),
                "n_junctions": st.column_config.NumberColumn("N° Uniones", format="%d")
            }
        )
        
        # Gráficos Simplificados
        st.markdown("#### Visualizaciones")
        col_g1, col_g2 = st.columns(2)
        
        with col_g1:
            chart_len = alt.Chart(df_skel).mark_bar().encode(
                x=alt.X('total_branch_length_um:Q', bin=alt.Bin(maxbins=20), title='Total Length (µm)'),
                y=alt.Y('count()', title='Frequency')
            ).properties(height=250, title="Distribution: Total Length")
            st.altair_chart(chart_len, use_container_width=True)
            
        with col_g2:
            chart_ram = alt.Chart(df_skel).mark_bar().encode(
                x=alt.X('ramification_index:Q', bin=alt.Bin(maxbins=20), title='Ramification Index'),
                y=alt.Y('count()', title='Frequency')
            ).properties(height=250, title="Distribution: Ramification Index")
            st.altair_chart(chart_ram, use_container_width=True)

    except Exception as e:
        st.error(f"Error leyendo métricas topológicas: {e}")

# --- 5. Métricas de Sholl (Paso 04) ---
st.markdown("---")
st.markdown("### 🎯 Análisis de Sholl (Paso 04)")

if not p_sholl_sum.exists():
    st.info("No se encontró `sholl_summary.csv`. Ejecutá el paso 04 primero.")
else:
    try:
        df_sholl_sum = pd.read_csv(p_sholl_sum)
        
        # Métricas resumen (Solo Radio Crítico)
        st.markdown("#### Estadísticas Globales de Sholl")
        st.metric("Radio Crítico Mediano", f"{df_sholl_sum['critical_radius_um'].median():.1f} µm", help="Distancia del núcleo donde ocurre la máxima ramificación")
        
        # Tabla de resumen
        st.markdown("#### Tabla de Métricas de Sholl")
        st.dataframe(
            df_sholl_sum.round(3),
            use_container_width=True,
            height=300,
            column_config={
                "label": st.column_config.NumberColumn("ID", format="%d"),
                "auc": st.column_config.NumberColumn("AUC", format="%.1f"),
                "peak_intersections": st.column_config.NumberColumn("Pico", format="%.1f"),
                "critical_radius_um": st.column_config.NumberColumn("Radio Crítico (µm)", format="%.1f")
            }
        )
        
        # Perfiles de Sholl por célula
        p_sholl_2d = out_dir / "sholl_2d.csv"
        if p_sholl_2d.exists():
            st.markdown("#### Perfiles de Sholl por Célula")
            df_sholl_curves = pd.read_csv(p_sholl_2d)
            
            chart_sholl = alt.Chart(df_sholl_curves).mark_line(point=False).encode(
                x=alt.X('radius_um:Q', title='Radius from Nucleus (µm)'),
                y=alt.Y('intersections:Q', title='Intersections'),
                color=alt.Color('label:N', title='Astrocyte', legend=None),
                opacity=alt.value(0.6)
            ).properties(height=350, title="Sholl Curves - All Cells").interactive()
            
            st.altair_chart(chart_sholl, use_container_width=True)
            
            # Promedio con banda de confianza
            df_sholl_avg = df_sholl_curves.groupby('radius_um')['intersections'].agg(['mean', 'std']).reset_index()
            df_sholl_avg['lower'] = df_sholl_avg['mean'] - df_sholl_avg['std']
            df_sholl_avg['upper'] = df_sholl_avg['mean'] + df_sholl_avg['std']
            
            chart_avg = alt.Chart(df_sholl_avg).mark_line(color='black', size=3).encode(
                x=alt.X('radius_um:Q', title='Radius (µm)'),
                y=alt.Y('mean:Q', title='Average Intersections')
            ).properties(height=300, title="Average Sholl Profile ± 1 SD")
            
            chart_band = alt.Chart(df_sholl_avg).mark_area(opacity=0.3, color='gray').encode(
                x='radius_um:Q',
                y='lower:Q',
                y2='upper:Q'
            )
            
            st.altair_chart((chart_band + chart_avg), use_container_width=True)
        
    except Exception as e:
        st.error(f"Error leyendo métricas de Sholl: {e}")

# --- 6. Visualización Interactiva ---
st.markdown("---")
st.markdown("### 👁️ Visualización en Napari")

col_view_1, col_view_2 = st.columns(2)
with col_view_1:
    open_napari_3d = st.button("👁️ Visualización 3D Completa", use_container_width=True, help="Ver máscaras 3D y canales originales")
with col_view_2:
    open_napari_2d = st.button("📊 Visualización 2D (Sholl/Esqueletos)", use_container_width=True, help="Ver proyección, esqueletos y anillos de Sholl")

# Visualizador 3D
if open_napari_3d:
    try:
        cal = _read_global_calibration()
        z, y, x = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        dapi_idx = int(cal.get("DAPI_CHANNEL_INDEX", 0))
        gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1))
        cmd = [sys.executable, str(napari_script), 
               "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x),
               "--dapi_idx", str(dapi_idx), "--gfap_idx", str(gfap_idx)]
        
        # Máscaras 3D
        if p_final.exists():
            cmd += ["--final", str(p_final)]
        
        subprocess.Popen(cmd, env=os.environ.copy())
        st.info("✅ Visualizador 3D lanzado.")
    except Exception as e:
        st.error(f"❌ Error al lanzar Napari 3D: {e}")

# Visualizador 2D
if open_napari_2d:
    try:
        cal = _read_global_calibration()
        y, x = float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        dapi_idx = int(cal.get("DAPI_CHANNEL_INDEX", 0))
        gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1))
        napari_2d_script = root / "streamlit" / "napari_viewer_2d.py"
        cmd = [sys.executable, str(napari_2d_script), 
               "--path", str(img_path), "--y", str(y), "--x", str(x),
               "--dapi_idx", str(dapi_idx), "--gfap_idx", str(gfap_idx)]
        
        # Pasando archivos generados
        # Nota: Usamos nombres estandarizados
        gfap_candidates = out_dir / "03_gfap_filtered_mask_2d.tif"
        final_2d = out_dir / "04_final_astrocytes_mask_2d.tif" 
        skel_2d = out_dir / "05_skeleton_labels_2d.tif"
        rings = out_dir / "sholl_rings_2d_native.json"
        
        if gfap_candidates.exists(): cmd += ["--gfap_2d", str(gfap_candidates)]
        if final_2d.exists(): cmd += ["--final_2d", str(final_2d)]
        if skel_2d.exists(): cmd += ["--skeleton_2d", str(skel_2d)]
        if rings.exists(): cmd += ["--rings", str(rings)]
        
        subprocess.Popen(cmd, env=os.environ.copy())
        st.info("✅ Visualizador 2D lanzado.")
    except Exception as e:
        st.error(f"❌ Error al lanzar Napari 2D: {e}")