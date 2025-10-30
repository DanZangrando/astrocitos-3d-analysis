import json
import os
import subprocess
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile

from ui.sidebar import render_sidebar
from ui.sholl import ShollParams, compute_sholl

st.title("Test de Sholl y Parámetros")
render_sidebar(show_calibration=True)

# Rutas
root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
calib_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer.py"


def _read_global_calibration():
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())
        except Exception:
            pass
    return {}


def _get_out_dir(img_path: Path) -> Path:
    out_dir = root / "data" / "processed" / img_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


# 1) Selección de imagen (con filtro global de grupo)
files_all = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files_all:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

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

group_filter = st.session_state.get("group_filter", "Todos")
if group_filter == "Todos":
    files = files_all
else:
    files = [p for p in files_all if _detect_group(p, root) == group_filter]
if not files:
    st.info("No hay preparados para el grupo seleccionado.")
    st.stop()
labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files))), format_func=lambda i: labels[i])
img_path = files[idx]
out_dir = _get_out_dir(img_path)

# --------- Recalcular por ámbito desde este paso ---------
with st.expander("Recalcular por ámbito desde este paso", expanded=False):
    scope = st.radio("Ámbito", options=["Preparado seleccionado", "Grupo", "Todos"], horizontal=True, key="p05_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p05_scope_group")
    if st.button("▶️ Recalcular (desde 05)", key="p05_recalc"):
        try:
            from ui.runner import run_scope, read_calibration
            cal = read_calibration(root/"streamlit"/"calibration.json")
            sc = "selected" if scope=="Preparado seleccionado" else ("group" if scope=="Grupo" else "all")
            sel = img_path if sc=="selected" else None
            res = run_scope(root, scope=sc, start_step="05", cal=cal, selected=sel, group=scope_group, overwrite_from_step=True)
            ok = sum(1 for _, stt in res if not stt.get("error"))
            st.success(f"Listo: {ok}/{len(res)} preparados procesados desde 05.")
        except Exception as e:
            st.error(f"Error al ejecutar: {e}")

# 2) Seteo de parámetros Sholl
st.markdown("### Parámetros de Sholl (µm)")
cal = _read_global_calibration()
z_um = float(cal.get('z', 1.0)); y_um = float(cal.get('y', 0.3)); x_um = float(cal.get('x', 0.3))
# Sugerir máximos según radio máximo del skeleton si existe resumen
try:
    max_suggest = None
    p = out_dir / "skeletons" / "summary.csv"
    if p.exists():
        df = pd.read_csv(p)
        if "max_radius_um" in df.columns and df["max_radius_um"].notna().any():
            max_suggest = float(np.nanpercentile(df["max_radius_um"].to_numpy(), 90))
except Exception:
    max_suggest = None
# Defaults desde calibration.json (si existen); si no, usar sugerencias o valores seguros
min_def = float(cal.get("SHOLL_MIN_RADIUS_UM", 0.0))
max_def = float(cal.get("SHOLL_MAX_RADIUS_UM", max_suggest if max_suggest is not None else 50.0))
step_def = float(cal.get("SHOLL_STEP_UM", 20.0))

min_r = st.number_input("Radio mínimo", value=float(min_def), min_value=0.0, step=0.5)
max_r = st.number_input("Radio máximo", value=float(max_def), min_value=0.0, step=1.0)
step_r = st.number_input("Paso entre anillos", value=float(step_def), min_value=0.1, step=0.1)

save_params = st.button("💾 Guardar parámetros de Sholl (global)")
run_sholl = st.button("🧮 Ejecutar Sholl y Guardar resultados")
open_napari_rings = st.button("👁️ Ver anillos en Napari")

if save_params:
    # Guardar SOLO a nivel global
    try:
        g = _read_global_calibration() or {}
        g.update({
            "SHOLL_MIN_RADIUS_UM": float(min_r),
            "SHOLL_MAX_RADIUS_UM": float(max_r),
            "SHOLL_STEP_UM": float(step_r),
        })
        calib_path.write_text(json.dumps(g, indent=2))
        st.success(f"Parámetros Sholl guardados globalmente en calibration.json")
    except Exception as e:
        st.error(f"No se pudieron guardar parámetros globales: {e}")


# 3) Resultados
st.markdown("---")
st.markdown("### Resultados del test de Sholl")

params = ShollParams(min_radius_um=float(min_r), max_radius_um=float(max_r), step_um=float(step_r))
spacing = (z_um, y_um, x_um)

if run_sholl:
    try:
        with st.spinner("Calculando Sholl…"):
            df = compute_sholl(out_dir, spacing, params, restrict_to_final=True, save_rings_json=True)
        if df.empty:
            st.info("No se generaron resultados (¿faltan máscaras finales/skeleton?)")
        else:
            st.success("Resultados guardados en sholl.csv")
    except Exception as e:
        st.error(f"Error en Sholl: {e}")

# Cargar si existe
try:
    df_sh = pd.read_csv(out_dir / "sholl.csv") if (out_dir / "sholl.csv").exists() else pd.DataFrame()
except Exception:
    df_sh = pd.DataFrame()

if not df_sh.empty:
    st.dataframe(df_sh.head(200), use_container_width=True)
    try:
        base = alt.Chart(df_sh)
        st.altair_chart(
            base.mark_line().encode(x=alt.X('radius_um:Q'), y=alt.Y('intersections:Q'), color='label:N').properties(height=260),
            use_container_width=True,
        )
        # Pico de intersecciones por célula
        peak = df_sh.sort_values(['label','intersections'], ascending=[True, False]).groupby('label', as_index=False).first()
        st.altair_chart(
            alt.Chart(peak).mark_bar().encode(x=alt.X('intersections:Q', title='Pico de intersecciones'), y=alt.Y('label:N')).properties(height=260),
            use_container_width=True,
        )
    except Exception:
        pass

if open_napari_rings:
    try:
        rings_path = out_dir / "sholl_rings.json"
        if not rings_path.exists():
            st.warning("No se encontró sholl_rings.json. Ejecutá Sholl primero.")
        else:
            cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z_um), "--y", str(y_um), "--x", str(x_um)]
            # Cargar máscaras si existen
            for p, flag in [
                (out_dir / "01_otsu_mask.tif", "--otsu"),
                (out_dir / "02_cellpose_mask.tif", "--cellpose"),
                (out_dir / "03_gfap_microglia_filtered_mask.tif", "--gfap"),
                (out_dir / "04_final_astrocytes_mask.tif", "--final"),
                (out_dir / "05_skeleton_labels.tif", "--skeleton"),
            ]:
                if p.exists():
                    cmd += [flag, str(p)]
            cmd += ["--rings", str(rings_path)]
            env = os.environ.copy()
            env["NAPARI_DISABLE_PLUGIN_AUTOLOAD"] = "1"
            subprocess.Popen(cmd, env=env)
            st.info("Napari lanzado con anillos de Sholl.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")
