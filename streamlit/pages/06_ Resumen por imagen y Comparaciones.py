import json
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile
from scipy import stats

from ui.sidebar import render_sidebar

st.title("Resumen por imagen (detalle)")
render_sidebar(show_calibration=True)

root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"

# Consistent palette across charts
PALETTE = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}
GROUP_SCALE = alt.Scale(domain=["CTL","Hipoxia"], range=[PALETTE["CTL"], PALETTE["Hipoxia"]])

def box_violin(df: pd.DataFrame, value_col: str, group_col: str = "group", title_x: str = "", height: int = 260):
    box = alt.Chart(df).mark_boxplot().encode(
        y=alt.Y(f"{group_col}:N", title="Grupo"),
        x=alt.X(f"{value_col}:Q", title=title_x),
        color=alt.Color(f"{group_col}:N", scale=GROUP_SCALE, legend=None),
    ).properties(height=height)
    violin = alt.Chart(df).transform_density(
        value_col, as_=[value_col, 'density'], groupby=[group_col]
    ).mark_area(orient='horizontal', opacity=0.35).encode(
        y=alt.Y(f"{value_col}:Q", title=title_x),
        x=alt.X('density:Q', axis=None),
        color=alt.Color(f"{group_col}:N", scale=GROUP_SCALE, legend=None),
    ).properties(height=height)
    return box | violin

files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

# Selección
labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files))), format_func=lambda i: labels[i])
img_path = files[idx]
out_dir = root / "data" / "processed" / img_path.stem

st.markdown("### Resumen detallado de la imagen seleccionada")
# Paths
p_otsu = out_dir / "01_otsu_mask.tif"
p_cp = out_dir / "02_cellpose_mask.tif"
p_gfap = out_dir / "03_gfap_microglia_filtered_mask.tif"
p_final = out_dir / "04_final_astrocytes_mask.tif"
p_skel = out_dir / "05_skeleton_labels.tif"
p_sum = out_dir / "skeletons" / "summary.csv"
p_sholl = out_dir / "sholl.csv"

n_cellpose = int(tifffile.imread(p_cp).max()) if p_cp.exists() else 0
n_gfap = int((np.unique(tifffile.imread(p_gfap))>0).sum()) if p_gfap.exists() else 0
n_final = int((np.unique(tifffile.imread(p_final))>0).sum()) if p_final.exists() else 0
n_skel = int((np.unique(tifffile.imread(p_skel))>0).sum()) if p_skel.exists() else 0

c1,c2,c3,c4 = st.columns(4)
c1.metric("Cellpose", n_cellpose)
c2.metric("GFAP/Microglía", n_gfap)
c3.metric("Astrocitos finales", n_final)
c4.metric("Esqueletos", n_skel)

# Summary.csv y gráficos + tarjetas adicionales
if p_sum.exists():
    try:
        dfs = pd.read_csv(p_sum)
        # Tarjetas con medianas (no medias) y sumas por preparado
        med_len = float(dfs["total_length_um"].median()) if "total_length_um" in dfs.columns else None
        med_dom = float(dfs["domain_volume_um3"].median()) if "domain_volume_um3" in dfs.columns and dfs["domain_volume_um3"].notna().any() else None
        med_tube = float(dfs["tube_volume_um3"].median()) if "tube_volume_um3" in dfs.columns and dfs["tube_volume_um3"].notna().any() else None
        sum_tube = float(dfs["tube_volume_um3"].sum()) if "tube_volume_um3" in dfs.columns and dfs["tube_volume_um3"].notna().any() else None
        d1,d2,d3,d4 = st.columns(4)
        d1.metric("Mediana longitud (µm)", f"{med_len:.1f}" if isinstance(med_len, float) else "—")
        d2.metric("Mediana volumen dominio (µm³)", f"{med_dom:.0f}" if isinstance(med_dom, float) else "—")
        d3.metric("Mediana volumen tubo (µm³)", f"{med_tube:.0f}" if isinstance(med_tube, float) else "—")
        d4.metric("Suma volumen tubo (µm³)", f"{sum_tube:.0f}" if isinstance(sum_tube, float) else "—")

        st.dataframe(dfs, use_container_width=True)

        # Violin plots modernos
        def violin_chart(df, col, title):
            # Violin por valor (sin grupos locales)
            base = alt.Chart(df).transform_density(
                col,
                as_=[col, 'density'],
            )
            return base.mark_area(opacity=0.5).encode(x=alt.X('density:Q', title='Densidad'), y=alt.Y(f'{col}:Q', title=title))

        charts = []
        if "total_length_um" in dfs.columns:
            charts.append(alt.Chart(dfs).mark_bar().encode(x=alt.X('total_length_um:Q', bin=alt.Bin(maxbins=30), title='Longitud (µm)'), y='count()').properties(height=200))
        if "domain_volume_um3" in dfs.columns:
            charts.append(alt.Chart(dfs).mark_bar().encode(x=alt.X('domain_volume_um3:Q', bin=alt.Bin(maxbins=30), title='Volumen dominio (µm³)'), y='count()').properties(height=200))
        if charts:
            st.altair_chart(alt.vconcat(*charts).resolve_scale(x='independent'), use_container_width=True)
    except Exception as e:
        st.error(f"No se pudo leer summary.csv: {e}")

# Sholl (mostrar solo si existe)
if p_sholl.exists():
    st.markdown("### Sholl")
    try:
        dfsh = pd.read_csv(p_sholl)
        st.dataframe(dfsh.head(200), use_container_width=True)
        st.altair_chart(alt.Chart(dfsh).mark_line().encode(x='radius_um:Q', y='intersections:Q', color='label:N').properties(height=260), use_container_width=True)
        peak = dfsh.sort_values(['label','intersections'], ascending=[True, False]).groupby('label', as_index=False).first()
        st.altair_chart(alt.Chart(peak).mark_bar().encode(x=alt.X('intersections:Q', title='Pico de intersecciones'), y='label:N').properties(height=240), use_container_width=True)
    except Exception as e:
        st.error(f"No se pudo leer sholl.csv: {e}")

st.markdown("---")
st.markdown("### Ver preparado en 3D (Napari)")
open_napari = st.button("👁️ Abrir con máscaras finales, esqueletos y anillos de Sholl")
if open_napari:
    try:
        calib_path = root / "streamlit" / "calibration.json"
        cal = json.loads(calib_path.read_text()) if calib_path.exists() else {}
        z = float(cal.get('z', 1.0)); y = float(cal.get('y', 1.0)); x = float(cal.get('x', 1.0))
        napari_script = root / "streamlit" / "napari_viewer.py"
        import os, sys, subprocess
        cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
        # añadir los artefactos si existen
        for pth, flag in [
            (p_final, "--final"),
            (p_skel, "--skeleton"),
            (out_dir/"sholl_rings.json", "--rings"),
        ]:
            if Path(pth).exists():
                cmd += [flag, str(pth)]
        env = os.environ.copy()
        env["NAPARI_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        subprocess.Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")
