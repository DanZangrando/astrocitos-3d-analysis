import json
from pathlib import Path
import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile
from scipy import stats
from ui.sidebar import render_sidebar
from ui.plots import GROUP_SCALE # Importar la paleta
from ui import pipeline # Para cargar utils
import subprocess # <--- CORRECCIÓN
import sys # <--- CORRECCIÓN
import os # <--- CORRECCIÓN

st.title("Análisis por Preparado")
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
labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files))), format_func=lambda i: labels[i])
img_path = files[idx]
out_dir = root / "data" / "processed" / img_path.stem
group = _detect_group(img_path, root)
st.markdown(_group_badge_html(group), unsafe_allow_html=True)
st.caption(f"Directorio de procesados: {out_dir.relative_to(root)}")

# --- 2. Resumen del Pipeline (Conteos) ---
st.markdown("### Resumen de Conteo del Pipeline")
p_cp = out_dir / "02_cellpose_mask.tif"
p_nuc_metrics = out_dir / "03_nucleus_metrics.csv"
p_final = out_dir / "04_final_astrocytes_mask.tif"
p_sholl_2d = out_dir / "sholl_2d_native.csv"
p_sholl_sum = out_dir / "sholl_summary.csv"

n_cellpose, n_gfap, n_final, n_sholl = 0, 0, 0, 0
if p_nuc_metrics.exists():
    try:
        df_nuc_temp = pd.read_csv(p_nuc_metrics)
        n_cellpose = df_nuc_temp.shape[0]
        n_gfap = int(df_nuc_temp["is_astrocyte_candidate"].sum())
    except Exception: pass
elif p_cp.exists():
    try: n_cellpose = int(tifffile.imread(p_cp).max())
    except Exception: pass

if p_final.exists():
    try: n_final = int((np.unique(tifffile.imread(p_final)) > 0).sum())
    except Exception: pass
    
if p_sholl_2d.exists():
    try: n_sholl = pd.read_csv(p_sholl_2d)["label"].nunique()
    except Exception: pass

c1, c2, c3, c4 = st.columns(4)
c1.metric("Núcleos (02)", n_cellpose)
c2.metric("Candidatos (03)", n_gfap)
c3.metric("Final (04)", n_final)
c4.metric("Con Sholl (04)", n_sholl)

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