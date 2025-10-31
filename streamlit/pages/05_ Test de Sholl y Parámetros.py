import json
import os
import subprocess # <--- CORRECCIÓN
import sys
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import streamlit as st
import tifffile
from ui.sidebar import render_sidebar
# --- Importar la lógica unificada ---
from ui import pipeline, runner

st.title("Test de Sholl y Parámetros")
render_sidebar(show_calibration=True)

# Rutas
root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
calib_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer.py"

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

def _get_out_dir(img_path: Path) -> Path:
    out_dir = root / "data" / "processed" / img_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir

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

# 1) Selección de imagen
files_all = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files_all:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

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
glob_calib = _read_global_calibration()

# Verificar archivos de entrada
if not (out_dir / "04_final_astrocytes_mask.tif").exists() or not (out_dir / "05_skeleton_labels.tif").exists():
    st.error("Faltan '04_final_astrocytes_mask.tif' o '05_skeleton_labels.tif'. Ejecutá los pasos 03 y 04 primero.")
    st.stop()

# --------- Recalcular por ámbito (Batch) ---------
with st.expander("Recalcular por ámbito (Batch)", expanded=False):
    scope = st.radio("Ámbito (Batch)", options=["Grupo", "Todos"], horizontal=True, key="p05_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p05_scope_group")
    if st.button("▶️ Recalcular (desde 05)", key="p05_recalc"):
        cal = _read_global_calibration()
        sc = "group" if scope=="Grupo" else "all"
        res = runner.run_scope(root, scope=sc, start_step="05", cal=cal, selected=None, group=scope_group, overwrite_from_step=True)
        ok = sum(1 for _, stt in res if not stt.get("error"))
        st.success(f"Listo: {ok}/{len(res)} preparados procesados desde 05.")

# 2) Seteo de parámetros Sholl
st.markdown("### Parámetros de Sholl (µm)")
min_def = float(glob_calib.get("SHOLL_MIN_RADIUS_UM", 0.0))
max_def = float(glob_calib.get("SHOLL_MAX_RADIUS_UM", 50.0))
step_def = float(glob_calib.get("SHOLL_STEP_UM", 2.0))

min_r = st.number_input("Radio mínimo", value=float(min_def), min_value=0.0, step=0.5)
max_r = st.number_input("Radio máximo", value=float(max_def), min_value=0.0, step=1.0)
step_r = st.number_input("Paso entre anillos", value=float(step_def), min_value=0.1, step=0.1)

save_params = st.button("💾 Guardar parámetros de Sholl (global)")
run_sholl = st.button("🧮 Ejecutar Sholl y Guardar (Preparado actual)")
open_napari_rings = st.button("👁️ Ver anillos en Napari")

if save_params:
    try:
        g = _read_global_calibration() or {}
        g.update({
            "SHOLL_MIN_RADIUS_UM": float(min_r),
            "SHOLL_MAX_RADIUS_UM": float(max_r),
            "SHOLL_STEP_UM": float(step_r),
        })
        _save_global_calibration(g)
        st.success(f"Parámetros Sholl guardados globalmente en calibration.json")
    except Exception as e:
        st.error(f"No se pudieron guardar parámetros globales: {e}")

# 3) Ejecución (REFACTORIZADA)
if run_sholl:
    try:
        # Guardar parámetros antes de correr
        cal = _read_global_calibration()
        cal.update({
            "SHOLL_MIN_RADIUS_UM": float(min_r),
            "SHOLL_MAX_RADIUS_UM": float(max_r),
            "SHOLL_STEP_UM": float(step_r),
        })
        _save_global_calibration(cal)
        st.info("Parámetros de Sholl guardados en calibration.json")
        
        with st.spinner("Calculando Sholl..."):
            # Ejecutar SOLO el Paso 05
            pipeline.run_sholl_and_save(
                out_dir=out_dir,
                cal=cal,
                restrict_to_final=True,
                save_rings_json=True
            )
        st.success("Resultados de Sholl guardados.")
        st.session_state["__refresh_metrics_sholl"] = True
    except Exception as e:
        st.error(f"Error en Sholl: {e}")
        st.exception(e)

# 4) Resultados (REFACTORIZADA)
st.markdown("---")
st.markdown("### Resultados del test de Sholl")

def render_sholl_results():
    sholl_csv_path = out_dir / "sholl.csv"
    sholl_summary_path = out_dir / "sholl_summary.csv"
    
    if not sholl_csv_path.exists() or not sholl_summary_path.exists():
        st.info("No se han generado resultados de Sholl. Ejecutá el análisis primero.")
        return

    try:
        df_sh = pd.read_csv(sholl_csv_path)
        df_sum = pd.read_csv(sholl_summary_path)
    except Exception as e:
        st.error(f"No se pudieron cargar los archivos de Sholl: {e}")
        return

    st.markdown("#### Métricas Escalares (por astrocito)")
    st.dataframe(df_sum.round(3), use_container_width=True)

    c1, c2, c3 = st.columns(3)
    c1.metric("AUC Mediana", f"{df_sum['auc'].median():.1f}")
    c2.metric("Pico Mediano", f"{df_sum['peak_intersections'].median():.1f}")
    c3.metric("Radio Crítico Mediano (µm)", f"{df_sum['critical_radius_um'].median():.1f}")
    
    st.markdown("#### Curvas de Sholl (por astrocito)")
    st.altair_chart(
        alt.Chart(df_sh).mark_line().encode(
            x=alt.X('radius_um:Q', title="Radio (µm)"), 
            y=alt.Y('intersections:Q', title="Intersecciones"), 
            color=alt.Color('label:N', title="Célula")
        ).properties(height=260).interactive(),
        use_container_width=True,
    )

if st.button("🔄 Actualizar resultados") or st.session_state.get("__refresh_metrics_sholl", False):
    render_sholl_results()
    st.session_state["__refresh_metrics_sholl"] = False

if open_napari_rings:
    try:
        rings_path = out_dir / "sholl_rings.json"
        if not rings_path.exists():
            st.warning("No se encontró sholl_rings.json. Ejecutá Sholl primero.")
        else:
            z, y, x = float(glob_calib.get('z', 1.0)), float(glob_calib.get('y', 0.3)), float(glob_calib.get('x', 0.3))
            cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
            for p, flag in [
                (out_dir / "04_final_astrocytes_mask.tif", "--final"),
                (out_dir / "05_skeleton_labels.tif", "--skeleton"),
            ]:
                if p.exists(): cmd += [flag, str(p)]
            cmd += ["--rings", str(rings_path)]
            subprocess.Popen(cmd, env=os.environ.copy())
            st.info("Napari lanzado con anillos de Sholl.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")