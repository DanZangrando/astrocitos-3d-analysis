import sys
import os
import json
import subprocess # <--- CORRECCIÓN
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
from ui.sidebar import render_sidebar
# --- Importar la lógica unificada ---
from ui import pipeline, runner

# --- Configuración de página y UI ---
st.title("Umbral de Otsu y Segmentación con Cellpose")
render_sidebar(show_calibration=True)

# --- Rutas base ---
root = Path(__file__).resolve().parents[2]  # /<repo>
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

def get_output_dir_for_image(img_path: Path) -> Path:
    base_name = img_path.stem
    out_dir = root / "data" / "processed" / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir

def existing_results(out_dir: Path) -> dict:
    return {
        "otsu": (out_dir / "01_otsu_mask.tif").exists(),
        "cellpose": (out_dir / "02_cellpose_mask.tif").exists(),
    }

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

# --- 1) Selección de imagen ---
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

group_filter = st.session_state.get("group_filter", "Todos")
if group_filter == "Todos":
    files_avail = files
else:
    files_avail = [p for p in files if _detect_group(p, root) == group_filter]
if not files_avail:
    st.info("No hay preparados para el grupo seleccionado.")
    st.stop()

labels = [str(p.relative_to(root)) for p in files_avail]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files_avail))), format_func=lambda i: labels[i])
img_path = files_avail[idx]
out_dir = get_output_dir_for_image(img_path)
status = existing_results(out_dir)
group = _detect_group(img_path, root)
glob_calib = _read_global_calibration()

# --- Resumen Rápido ---
st.markdown("### Resumen rápido")
otsu_frac, cp_count = None, None
try:
    if status["otsu"]:
        otsu_mask = tifffile.imread(out_dir / "01_otsu_mask.tif").astype(bool)
        otsu_frac = float(otsu_mask.mean()) if otsu_mask.size else 0.0
    if status["cellpose"]:
        cp_mask = tifffile.imread(out_dir / "02_cellpose_mask.tif")
        cp_count = int(cp_mask.max())
except Exception:
    pass

m1, m2 = st.columns(2)
m1.metric("Otsu calculado", "✅ Sí" if status["otsu"] else "❌ No")
m2.metric("Cellpose calculado", "✅ Sí" if status["cellpose"] else "❌ No")
s1, s2 = st.columns(2)
s1.metric("Fracción voxeles Otsu", f"{100*otsu_frac:.1f}%" if isinstance(otsu_frac, (int,float)) else "—")
s2.metric("Núcleos Cellpose", cp_count if isinstance(cp_count, int) else "—")

def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"
st.markdown(_group_badge_html(group) + f"&nbsp;· Directorio: {out_dir.relative_to(root)}", unsafe_allow_html=True)

# --- Recalcular por ámbito (Usa runner) ---
with st.expander("Recalcular por ámbito (Batch)", expanded=False):
    scope = st.radio("Ámbito", options=["Grupo", "Todos"], horizontal=True, key="p02_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p02_scope_group")
    include_otsu = st.checkbox("Incluir Otsu (recalcular desde 01)", value=True, key="p02_inc_otsu")
    if st.button("▶️ Recalcular", key="p02_recalc"):
        start = "01" if include_otsu else "02"
        sc = "group" if scope=="Grupo" else "all"
        res = runner.run_scope(root, scope=sc, start_step=start, cal=glob_calib, selected=None, group=scope_group, overwrite_from_step=True)
        ok = sum(1 for _, stt in res if not stt.get("error"))
        st.success(f"Listo: {ok}/{len(res)} preparados procesados desde {start}.")

# --- 3) Configuración de segmentación ---
st.markdown("### Parámetros de Segmentación")
n_channels = 1
try:
    arr_preview, axes_prev = pipeline.load_image_any(img_path)
    vol_prev = pipeline.reorder_to_zcyx(arr_preview, axes_prev)
    n_channels = vol_prev.shape[1]
except Exception as e:
    st.warning(f"No se pudo precargar la imagen: {e}")

colp1, colp2, colp3 = st.columns(3)
dapi_default = int(glob_calib.get("DAPI_CHANNEL_INDEX", 0))
dapi_default = min(max(0, dapi_default), max(0, n_channels-1))
with colp1:
    dapi_channel_index = st.number_input("Índice de canal DAPI (0=primero)", value=dapi_default, min_value=0, max_value=max(0, n_channels-1), step=1)
with colp2:
    nucleus_diameter = st.number_input("Diámetro de núcleo (px)", value=int(glob_calib.get("NUCLEUS_DIAMETER", 30)), min_value=5, step=1)
with colp3:
    use_gpu = st.checkbox("Usar GPU si está disponible", value=bool(glob_calib.get("CELLPOSE_USE_GPU", True)))

# --- 4) Acciones (REFACTORIZADAS) ---
st.markdown("### Acciones (Para el preparado actual)")
recompute_otsu = st.checkbox("Forzar recalcular Otsu (sobrescribir si existe)", value=False)
recompute_cellpose = st.checkbox("Forzar recalcular Cellpose (sobrescribir si existe)", value=False)

run_otsu = st.button("🧮 Ejecutar Otsu y Guardar")
run_cellpose = st.button("🧠 Ejecutar Cellpose y Guardar")
run_both = st.button("🚀 Ejecutar Otsu + Cellpose")

def update_calib_and_run(start_step: runner.Step, overwrite: bool):
    """Guarda los parámetros de la UI y luego ejecuta el runner."""
    glob_calib["DAPI_CHANNEL_INDEX"] = int(dapi_channel_index)
    glob_calib["NUCLEUS_DIAMETER"] = int(nucleus_diameter)
    glob_calib["CELLPOSE_USE_GPU"] = bool(use_gpu)
    _save_global_calibration(glob_calib)
    st.info("Parámetros de UI guardados en calibration.json")
    
    try:
        with st.spinner(f"Ejecutando pipeline desde el paso {start_step} para {img_path.stem}..."):
            runner.run_pipeline_for(
                root=root,
                img_path=img_path,
                cal=glob_calib,
                start_step=start_step,
                overwrite_from_step=overwrite
            )
        st.success(f"Paso {start_step} completado.")
        st.rerun() # Recargar para mostrar nuevo estado
    except Exception as e:
        st.error(f"Error en pipeline: {e}")
        st.exception(e)

if run_otsu:
    update_calib_and_run(start_step="01", overwrite=recompute_otsu)

if run_cellpose:
    update_calib_and_run(start_step="02", overwrite=recompute_cellpose)

if run_both:
    # "Run both" siempre sobrescribe 01, y luego 02 si es necesario
    update_calib_and_run(start_step="01", overwrite=True)


# --- 5) Ver en Napari (Corregido) ---
st.markdown("---")
open_napari_image = st.button("👁️ Abrir en Napari (solo imagen)")
open_napari_with_masks = st.button("🧪 Abrir en Napari con máscaras disponibles")

def _launch_napari(include_masks: bool):
    env = os.environ.copy()
    z = float(glob_calib.get('z', 1.0)); y = float(glob_calib.get('y', 0.3)); x = float(glob_calib.get('x', 0.3))
    cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
    if include_masks:
        otsu_path = out_dir / "01_otsu_mask.tif"
        cellpose_path = out_dir / "02_cellpose_mask.tif"
        if otsu_path.exists(): cmd += ["--otsu", str(otsu_path)]
        if cellpose_path.exists(): cmd += ["--cellpose", str(cellpose_path)]
    try:
        subprocess.Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")

if open_napari_image:
    _launch_napari(include_masks=False)
if open_napari_with_masks:
    _launch_napari(include_masks=True)