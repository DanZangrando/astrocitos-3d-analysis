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
        return "Hipoxia"
    if "/ctl/" in rel:
        return "CTL"
    return "CTL"

def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
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

# --- 3. Métricas de Núcleo (Paso 03) ---
st.markdown("---")
st.markdown("### 🔬 Métricas Nucleares (Paso 03)")

if not p_nuc_metrics.exists():
    st.info("No se encontró `03_nucleus_metrics.csv`. Ejecutá el paso 03 primero.")
else:
    try:
        df_nuc = pd.read_csv(p_nuc_metrics)
        
        # Métricas resumen
        st.markdown("#### Estadísticas Globales")
        col_n1, col_n2, col_n3, col_n4 = st.columns(4)
        col_n1.metric("Volumen Medio", f"{df_nuc['nucleus_volume_um3'].mean():.1f} µm³")
        col_n2.metric("Circularidad Media (2D)", f"{df_nuc['nucleus_sphericity'].mean():.3f}")
        col_n3.metric("% Candidatos", f"{100 * n_gfap / max(n_cellpose, 1):.1f}%")
        col_n4.metric("Total Células", f"{df_nuc.shape[0]}")
        
        # Tabla completa con mejor formato
        st.markdown("#### Tabla Completa de Métricas Nucleares")
        st.dataframe(
            df_nuc.round(3),
            use_container_width=True,
            height=350,
            column_config={
                "label": st.column_config.NumberColumn("ID", format="%d"),
                "nucleus_volume_um3": st.column_config.NumberColumn("Volumen (µm³)", format="%.1f"),
                "nucleus_sphericity": st.column_config.NumberColumn("Circularidad 2D", format="%.3f", 
                                                                     help="4π·área/perímetro² en proyección MIP"),
                "is_astrocyte_candidate": st.column_config.CheckboxColumn("GFAP+")
            }
        )
        
        # Gráficos de distribución
        st.markdown("#### Distribuciones")
        tab_vol, tab_circ = st.tabs(["📊 Volumen Nuclear", "📊 Circularidad"])
        
        with tab_vol:
            chart_vol = alt.Chart(df_nuc).mark_bar().encode(
                x=alt.X("nucleus_volume_um3:Q", bin=alt.Bin(maxbins=30), title="Volumen Núcleo (µm³)"),
                y=alt.Y("count()", title="Frecuencia"),
                color=alt.Color("is_astrocyte_candidate:N", title="Candidato", 
                               scale=alt.Scale(domain=[False, True], range=['#d62728', '#2ca02c']))
            ).properties(height=300, title="Distribución de Volumen Nuclear")
            st.altair_chart(chart_vol, use_container_width=True)
        
        with tab_circ:
            chart_circ = alt.Chart(df_nuc.dropna(subset=['nucleus_sphericity'])).mark_bar().encode(
                x=alt.X("nucleus_sphericity:Q", bin=alt.Bin(maxbins=30), title="Circularidad en Proyección 2D (0-1)"),
                y=alt.Y("count()", title="Frecuencia"),
                color=alt.Color("is_astrocyte_candidate:N", title="Candidato",
                               scale=alt.Scale(domain=[False, True], range=['#d62728', '#2ca02c']))
            ).properties(height=300, title="Distribución de Circularidad Nuclear (2D)")
            st.altair_chart(chart_circ, use_container_width=True)
            
    except Exception as e:
        st.error(f"Error leyendo métricas nucleares: {e}")

# --- 4. Métricas Topológicas de Esqueleto (Paso 04) ---
st.markdown("---")
st.markdown("### 🌿 Métricas Topológicas de Esqueleto (Paso 04)")

if not p_skel_sum.exists():
    st.info("No se encontró `skeletons/summary.csv`. Ejecutá el paso 04 primero.")
else:
    try:
        df_skel = pd.read_csv(p_skel_sum)
        
        # Métricas resumen en tarjetas
        st.markdown("#### Estadísticas Globales")
        col_s1, col_s2, col_s3, col_s4, col_s5 = st.columns(5)
        col_s1.metric("Longitud Media", f"{df_skel['total_branch_length_um'].mean():.1f} µm")
        col_s2.metric("Ramas Medias", f"{df_skel['n_branches'].mean():.1f}")
        col_s3.metric("Tortuosidad Media", f"{df_skel['tortuosity_mean'].mean():.3f}")
        col_s4.metric("Índice Ramificación", f"{df_skel['ramification_index'].mean():.2f}")
        col_s5.metric("Componentes", f"{df_skel['n_connected_components'].median():.0f}")
        
        # Tabla completa con mejor formato
        st.markdown("#### Tabla Completa de Métricas Topológicas")
        st.dataframe(
            df_skel.round(3),
            use_container_width=True,
            height=400,
            column_config={
                "label": st.column_config.NumberColumn("ID", format="%d"),
                "total_branch_length_um": st.column_config.NumberColumn("Longitud (µm)", format="%.1f"),
                "n_branches": st.column_config.NumberColumn("Ramas", format="%d"),
                "n_endpoints": st.column_config.NumberColumn("Terminaciones", format="%d"),
                "n_junctions": st.column_config.NumberColumn("Bifurcaciones", format="%d"),
                "tortuosity_mean": st.column_config.NumberColumn("Tortuosidad μ", format="%.3f"),
                "ramification_index": st.column_config.NumberColumn("Ramificación", format="%.2f"),
                "termination_index": st.column_config.NumberColumn("Terminación", format="%.2f"),
                "branch_length_cv": st.column_config.NumberColumn("CV Longitud", format="%.3f"),
                "n_connected_components": st.column_config.NumberColumn("Componentes", format="%d")
            }
        )
        
        # Gráficos de distribución y correlación
        st.markdown("#### Visualizaciones")
        tab_dist, tab_tort, tab_corr = st.tabs(["📊 Distribuciones", "🌀 Tortuosidad", "🔗 Correlaciones"])
        
        with tab_dist:
            col_d1, col_d2 = st.columns(2)
            with col_d1:
                chart_length = alt.Chart(df_skel).mark_bar().encode(
                    x=alt.X('total_branch_length_um:Q', bin=alt.Bin(maxbins=25), title='Longitud Total (µm)'),
                    y=alt.Y('count()', title='Frecuencia')
                ).properties(height=250, title="Distribución de Longitud")
                st.altair_chart(chart_length, use_container_width=True)
            
            with col_d2:
                chart_branches = alt.Chart(df_skel).mark_bar().encode(
                    x=alt.X('n_branches:Q', bin=alt.Bin(maxbins=20), title='Número de Ramas'),
                    y=alt.Y('count()', title='Frecuencia')
                ).properties(height=250, title="Distribución de Ramas")
                st.altair_chart(chart_branches, use_container_width=True)
        
        with tab_tort:
            chart_tort_hist = alt.Chart(df_skel).mark_bar().encode(
                x=alt.X('tortuosity_mean:Q', bin=alt.Bin(maxbins=30), title='Tortuosidad Media'),
                y=alt.Y('count()', title='Frecuencia')
            ).properties(height=250, title="Distribución de Tortuosidad")
            
            chart_tort_scatter = alt.Chart(df_skel).mark_circle(size=80).encode(
                x=alt.X('n_junctions:Q', title='Bifurcaciones'),
                y=alt.Y('tortuosity_mean:Q', title='Tortuosidad Media'),
                color=alt.Color('ramification_index:Q', title='Ramificación', scale=alt.Scale(scheme='viridis')),
                tooltip=['label', 'n_branches', 'tortuosity_mean', 'ramification_index']
            ).properties(height=250, title="Tortuosidad vs Complejidad").interactive()
            
            st.altair_chart(chart_tort_hist, use_container_width=True)
            st.altair_chart(chart_tort_scatter, use_container_width=True)
        
        with tab_corr:
            # Matriz de correlación para métricas clave
            metrics_corr = ['total_branch_length_um', 'n_branches', 'tortuosity_mean', 
                           'ramification_index', 'termination_index', 'branch_length_cv']
            df_corr = df_skel[metrics_corr].corr()
            
            st.markdown("**Matriz de Correlación (Métricas Clave)**")
            st.dataframe(df_corr.round(3), use_container_width=True)
            
            # Scatter plot: Ramificación vs Longitud
            chart_scatter = alt.Chart(df_skel).mark_circle(size=100, opacity=0.7).encode(
                x=alt.X('total_branch_length_um:Q', title='Longitud Total (µm)'),
                y=alt.Y('ramification_index:Q', title='Índice de Ramificación'),
                size=alt.Size('n_branches:Q', title='N° Ramas'),
                color=alt.Color('tortuosity_mean:Q', title='Tortuosidad', scale=alt.Scale(scheme='plasma')),
                tooltip=['label', 'total_branch_length_um', 'ramification_index', 'n_branches', 'tortuosity_mean']
            ).properties(height=400, title="Longitud vs Ramificación (tamaño=ramas, color=tortuosidad)").interactive()
            st.altair_chart(chart_scatter, use_container_width=True)
            
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
        
        # Métricas resumen
        st.markdown("#### Estadísticas Globales de Sholl")
        col_sh1, col_sh2, col_sh3 = st.columns(3)
        col_sh1.metric("AUC Mediana", f"{df_sholl_sum['auc'].median():.1f}", help="Área bajo la curva")
        col_sh2.metric("Pico Mediano", f"{df_sholl_sum['peak_intersections'].median():.1f}", help="Máximo de intersecciones")
        col_sh3.metric("Radio Crítico", f"{df_sholl_sum['critical_radius_um'].median():.1f} µm", help="Radio del pico")
        
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
                x=alt.X('radius_um:Q', title='Radio desde Núcleo (µm)'),
                y=alt.Y('intersections:Q', title='Intersecciones'),
                color=alt.Color('label:N', title='Astrocito', legend=None),
                opacity=alt.value(0.6)
            ).properties(height=350, title="Curvas de Sholl - Todas las Células").interactive()
            
            st.altair_chart(chart_sholl, use_container_width=True)
            
            # Promedio con banda de confianza
            df_sholl_avg = df_sholl_curves.groupby('radius_um')['intersections'].agg(['mean', 'std']).reset_index()
            df_sholl_avg['lower'] = df_sholl_avg['mean'] - df_sholl_avg['std']
            df_sholl_avg['upper'] = df_sholl_avg['mean'] + df_sholl_avg['std']
            
            chart_avg = alt.Chart(df_sholl_avg).mark_line(color='black', size=3).encode(
                x=alt.X('radius_um:Q', title='Radio (µm)'),
                y=alt.Y('mean:Q', title='Intersecciones Promedio')
            ).properties(height=300, title="Perfil de Sholl Promedio ± 1 SD")
            
            chart_band = alt.Chart(df_sholl_avg).mark_area(opacity=0.3, color='gray').encode(
                x='radius_um:Q',
                y='lower:Q',
                y2='upper:Q'
            )
            
            st.altair_chart((chart_band + chart_avg), use_container_width=True)
        
    except Exception as e:
        st.error(f"Error leyendo métricas de Sholl: {e}")

# --- 6. Visualización en Napari ---
st.markdown("---")
st.markdown("### 👁️ Visualización Interactiva")

open_napari_2d = st.button("� Visualización 2D Completa", use_container_width=True, 
                            help="Proyección 2D + máscaras candidatos/finales + esqueletos + anillos Sholl", 
                            key="napari_2d")

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
        
        # Archivos 2D (máscaras, esqueletos, anillos)
        gfap_candidates = out_dir / "03_gfap_filtered_mask_2d.tif"  # Candidatos GFAP+
        final_2d = out_dir / "04_final_astrocytes_mask_2d.tif"  # Astrocitos finales
        skel_2d = out_dir / "05_skeleton_labels_2d.tif"  # Esqueletos
        rings = out_dir / "sholl_rings_2d_native.json"  # Anillos Sholl
        
        if gfap_candidates.exists():
            cmd += ["--gfap_2d", str(gfap_candidates)]
        if final_2d.exists():
            cmd += ["--final_2d", str(final_2d)]
        if skel_2d.exists():
            cmd += ["--skeleton_2d", str(skel_2d)]
        if rings.exists():
            cmd += ["--rings", str(rings)]
        
        subprocess.Popen(cmd, env=os.environ.copy())
        st.success("✅ Visualizador 2D de Napari lanzado con todas las capas.")
    except Exception as e:
        st.error(f"❌ No se pudo lanzar Napari 2D: {e}")

# --- 3. Métricas de Núcleo (Paso 03) ---
if p_nuc_metrics.exists():
    st.markdown("### Métricas de Núcleo (de `03_nucleus_metrics.csv`)")
    try:
        df_nuc = pd.read_csv(p_nuc_metrics)
        st.dataframe(df_nuc.round(3), use_container_width=True)
        
        charts_nuc = []
        charts_nuc.append(alt.Chart(df_nuc).mark_bar().encode(
            x=alt.X("nucleus_volume_um3:Q", bin=alt.Bin(maxbins=30), title="Volumen Núcleo (µm³)"),
            y="count()",
            color=alt.Color("is_astrocyte_candidate:N", title="Candidato")
        ).properties(height=200))
        
        charts_nuc.append(alt.Chart(df_nuc.dropna(subset=['nucleus_sphericity'])).mark_bar().encode(
            x=alt.X("nucleus_sphericity:Q", bin=alt.Bin(maxbins=30), title="Esfericidad Núcleo"),
            y="count()",
            color=alt.Color("is_astrocyte_candidate:N", title="Candidato")
        ).properties(height=200))
        
        st.altair_chart(alt.vconcat(*charts_nuc), use_container_width=True)
    except Exception as e:
        st.error(f"No se pudo leer 03_nucleus_metrics.csv: {e}")

# --- 4. Métricas de Sholl (Paso 04) ---
if p_sholl_sum.exists():
    st.markdown("### Métricas de Sholl 2D Nativo (de `sholl_summary.csv`)")
    try:
        df_sholl_sum = pd.read_csv(p_sholl_sum)
        
        s1, s2, s3 = st.columns(3)
        s1.metric("AUC Mediana", f"{df_sholl_sum['auc'].median():.1f}")
        s2.metric("Pico Mediano", f"{df_sholl_sum['peak_intersections'].median():.1f}")
        s3.metric("Radio Crítico Mediano (µm)", f"{df_sholl_sum['critical_radius_um'].median():.1f}")
        
        st.dataframe(df_sholl_sum.round(3), use_container_width=True)

        # Perfiles de Sholl por célula
        if p_sholl_2d.exists():
            df_sholl_curves = pd.read_csv(p_sholl_2d)
            st.altair_chart(alt.Chart(df_sholl_curves).mark_line().encode(
                x=alt.X('radius_um:Q', title='Radio (µm)'), 
                y=alt.Y('intersections:Q', title='Intersecciones'), 
                color='label:N'
            ).properties(height=260, title="Perfiles de Sholl por Célula").interactive(), use_container_width=True)

    except Exception as e:
        st.error(f"No se pudo leer sholl_2d_native.csv: {e}")

# --- 6. Visualización 3D ---
st.markdown("---")
st.markdown("### Ver preparado completo en Napari")

# Botones separados para visualización 3D y 2D
col_view_1, col_view_2 = st.columns(2)
with col_view_1:
    open_napari_3d = st.button("👁️ Visualización 3D completa", use_container_width=True, help="Máscaras 3D en Napari")
with col_view_2:
    open_napari_2d = st.button("📊 Visualización 2D completa", use_container_width=True, help="Proyección 2D + esqueletos + anillos Sholl")
    
# Visualizador 3D completo
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
        st.info("Visualizador 3D de Napari lanzado.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari 3D: {e}")

# Visualizador 2D completo  
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
        
        # Archivos 2D
        skel_2d = out_dir / "05_skeleton_labels_2d.tif"
        rings = out_dir / "sholl_rings_2d_native.json"
        
        if skel_2d.exists():
            cmd += ["--skeleton_2d", str(skel_2d)]
        if rings.exists():
            cmd += ["--rings", str(rings)]
        
        subprocess.Popen(cmd, env=os.environ.copy())
        st.info("Visualizador 2D de Napari lanzado.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari 2D: {e}")