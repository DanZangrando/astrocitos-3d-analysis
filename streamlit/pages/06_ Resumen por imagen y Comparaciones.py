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

st.title("Resumen por imagen (detalle)")
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
p_skel_sum = out_dir / "skeletons" / "summary.csv"
p_sholl_sum = out_dir / "sholl_summary.csv"

n_cellpose, n_gfap, n_final, n_skel = 0, 0, 0, 0
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
    
if p_skel_sum.exists():
    try: n_skel = pd.read_csv(p_skel_sum).shape[0]
    except Exception: pass

c1, c2, c3, c4 = st.columns(4)
c1.metric("Núcleos (02)", n_cellpose)
c2.metric("Candidatos (03)", n_gfap)
c3.metric("Final (04)", n_final)
c4.metric("Esqueletos (05)", n_skel)

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

# --- 4. Métricas de Esqueleto (Paso 04) ---
if p_skel_sum.exists():
    st.markdown("### Métricas de Esqueleto (de `skeletons/summary.csv`)")
    try:
        dfs = pd.read_csv(p_skel_sum)
        
        med_len = float(dfs["total_length_um"].median()) if "total_length_um" in dfs.columns and dfs["total_length_um"].notna().any() else 0.0
        med_branch = float(dfs["n_branches"].median()) if "n_branches" in dfs.columns and dfs["n_branches"].notna().any() else 0.0
        med_tort = float(dfs["mean_tortuosity"].median()) if "mean_tortuosity" in dfs.columns and dfs["mean_tortuosity"].notna().any() else 0.0
        med_dom = float(dfs["domain_volume_um3"].median()) if "domain_volume_um3" in dfs.columns and dfs["domain_volume_um3"].notna().any() else 0.0
        med_tube_int = float(dfs["tube_mean_intensity"].median()) if "tube_mean_intensity" in dfs.columns and dfs["tube_mean_intensity"].notna().any() else 0.0

        d1,d2,d3,d4,d5 = st.columns(5)
        d1.metric("Longitud (µm)", f"{med_len:.1f}")
        d2.metric("# Ramas", f"{med_branch:.1f}")
        d3.metric("Tortuosidad", f"{med_tort:.2f}")
        d4.metric("Vol. Dominio (µm³)", f"{med_dom:.0f}")
        d5.metric("Int. Tubo (media)", f"{med_tube_int:.1f}")
        
        st.dataframe(dfs.round(3), use_container_width=True)

        charts_skel = []
        if "total_length_um" in dfs.columns:
            charts_skel.append(alt.Chart(dfs).mark_bar().encode(x=alt.X('total_length_um:Q', bin=alt.Bin(maxbins=30), title='Longitud (µm)'), y='count()').properties(height=200))
        if "n_branches" in dfs.columns:
            charts_skel.append(alt.Chart(dfs).mark_bar().encode(x=alt.X('n_branches:Q', bin=alt.Bin(maxbins=20), title='# Ramas'), y='count()').properties(height=200))
        if "mean_tortuosity" in dfs.columns:
            charts_skel.append(alt.Chart(dfs).mark_bar().encode(x=alt.X('mean_tortuosity:Q', bin=alt.Bin(maxbins=30), title='Tortuosidad Media'), y='count()').properties(height=200))
        if "domain_volume_um3" in dfs.columns:
            charts_skel.append(alt.Chart(dfs).mark_bar().encode(x=alt.X('domain_volume_um3:Q', bin=alt.Bin(maxbins=30), title='Vol. Dominio (µm³)'), y='count()').properties(height=200))
        
        if charts_skel:
            st.altair_chart(alt.vconcat(*charts_skel).resolve_scale(x='independent'), use_container_width=True)
    except Exception as e:
        st.error(f"No se pudo leer summary.csv: {e}")

# --- 5. Métricas de Sholl (Paso 05) ---
if p_sholl_sum.exists():
    st.markdown("### Métricas de Sholl (de `sholl_summary.csv` y `sholl.csv`)")
    try:
        df_sholl_sum = pd.read_csv(p_sholl_sum)
        df_sholl_curves = pd.read_csv(out_dir / "sholl.csv")
        
        s1, s2, s3 = st.columns(3)
        s1.metric("AUC Mediana", f"{df_sholl_sum['auc'].median():.1f}")
        s2.metric("Pico Mediano", f"{df_sholl_sum['peak_intersections'].median():.1f}")
        s3.metric("Radio Crítico Mediano (µm)", f"{df_sholl_sum['critical_radius_um'].median():.1f}")
        
        st.dataframe(df_sholl_sum.round(3), use_container_width=True)

        st.altair_chart(alt.Chart(df_sholl_curves).mark_line().encode(
            x='radius_um:Q', 
            y='intersections:Q', 
            color='label:N'
        ).properties(height=260).interactive(), use_container_width=True)

    except Exception as e:
        st.error(f"No se pudo leer sholl.csv o sholl_summary.csv: {e}")

# --- 6. Visualización 3D ---
st.markdown("---")
st.markdown("### Ver preparado en 3D (Napari)")
open_napari = st.button("👁️ Abrir con máscaras finales, esqueletos y anillos de Sholl")
if open_napari:
    try:
        cal = _read_global_calibration()
        z, y, x = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        cmd = [sys.executable, str(napari_script), 
               "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
        
        for pth, flag in [
            (p_final, "--final"),
            (out_dir / "05_skeleton_labels.tif", "--skeleton"),
            (out_dir / "sholl_rings.json", "--rings"),
        ]:
            if pth.exists():
                cmd += [flag, str(pth)]
        
        subprocess.Popen(cmd, env=os.environ.copy())
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")