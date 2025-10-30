import sys
import os
import json
import subprocess
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
from skimage.filters import threshold_otsu
from ui.sidebar import render_sidebar


# --------- Configuración de página y UI ---------
st.title("Umbral de Otsu y Segmentación con Cellpose")
render_sidebar(show_calibration=True)


# --------- Rutas base ---------
root = Path(__file__).resolve().parents[2]  # /<repo>
raw_dir = root / "data" / "raw"
calib_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer.py"


# --------- Utilidades ---------
def _read_global_calibration():
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())  # Unificado: metadatos físicos + parámetros globales
        except Exception:
            pass
    return {}


def reorder_to_zcyx(arr: np.ndarray, axes: str | None):
    """Reordena un array con ejes reportados a (Z, C, Y, X)."""
    if axes is None:
        if arr.ndim != 4:
            raise ValueError(f"Forma inesperada sin ejes: {arr.shape}")
        chan_axis = int(np.argmin(arr.shape))
        if chan_axis != 1:
            arr = np.moveaxis(arr, chan_axis, 1)
        return arr
    ax_list = list(axes)
    # Seleccionar T=0 si existe
    if 'T' in ax_list:
        t_idx = ax_list.index('T')
        arr = np.take(arr, indices=0, axis=t_idx)
        ax_list.pop(t_idx)
    if 'C' not in ax_list:
        arr = np.expand_dims(arr, axis=0)
        ax_list = ['C'] + ax_list
    needed = ['Z', 'C', 'Y', 'X']
    if not all(a in ax_list for a in needed):
        raise ValueError(f"Ejes insuficientes tras seleccionar T=0: {ax_list} (se requieren Z,C,Y,X)")
    src_order = [ax_list.index(a) for a in needed]
    return np.transpose(arr, axes=src_order)


def load_image_any(path: Path):
    """Carga .tif/.tiff con tifffile y devuelve (array, axes str or None)."""
    suffix = path.suffix.lower()
    if suffix in ('.tif', '.tiff'):
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            axes = getattr(series, 'axes', None)
            arr = series.asarray()
        return arr, axes
    else:
        raise ValueError(f"Extensión no soportada: {suffix}")


def get_output_dir_for_image(img_path: Path) -> Path:
    base_name = img_path.stem  # p.ej.: "Inmuno 26-07-23.lif - CTL 1-2 a"
    out_dir = root / "data" / "processed" / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def existing_results(out_dir: Path) -> dict:
    return {
        "otsu": (out_dir / "01_otsu_mask.tif").exists(),
        "cellpose": (out_dir / "02_cellpose_mask.tif").exists(),
    }


# --------- 1) Selección de imagen ---------
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files:
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

# Filtro por grupo (unificado desde la barra lateral)
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

st.markdown("### Resumen rápido")
# Cargar métricas si existen
otsu_frac = None
cp_count = None
try:
    if (out_dir / "01_otsu_mask.tif").exists():
        otsu_mask = tifffile.imread(out_dir / "01_otsu_mask.tif").astype(bool)
        otsu_frac = float(otsu_mask.mean()) if otsu_mask.size else 0.0
    if (out_dir / "02_cellpose_mask.tif").exists():
        cp_mask = tifffile.imread(out_dir / "02_cellpose_mask.tif")
        cp_count = int(cp_mask.max())
except Exception:
    pass

m1, m2 = st.columns(2)
m1.metric("Otsu calculado", "Sí" if status["otsu"] else "No")
m2.metric("Cellpose calculado", "Sí" if status["cellpose"] else "No")

if status["otsu"] or status["cellpose"]:
    s1, s2 = st.columns(2)
    s1.metric("Fracción voxeles Otsu", f"{100*otsu_frac:.1f}%" if isinstance(otsu_frac, (int,float)) else "—")
    s2.metric("Núcleos Cellpose", cp_count if isinstance(cp_count, int) else "—")

def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"

st.markdown(_group_badge_html(group) + f"&nbsp;· Directorio: {out_dir.relative_to(root)}", unsafe_allow_html=True)

# --------- Recalcular por ámbito desde este paso ---------
with st.expander("Recalcular por ámbito desde este paso", expanded=False):
    scope = st.radio("Ámbito", options=["Preparado seleccionado", "Grupo", "Todos"], horizontal=True, key="p02_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p02_scope_group")
    include_otsu = st.checkbox("Incluir Otsu (recalcular desde 01)", value=True, key="p02_inc_otsu")
    if st.button("▶️ Recalcular", key="p02_recalc"):
        try:
            from ui.runner import run_scope, read_calibration
            cal = read_calibration(root/"streamlit"/"calibration.json")
            start = "01" if include_otsu else "02"
            sc = "selected" if scope=="Preparado seleccionado" else ("group" if scope=="Grupo" else "all")
            sel = img_path if sc=="selected" else None
            res = run_scope(root, scope=sc, start_step=start, cal=cal, selected=sel, group=scope_group, overwrite_from_step=True)
            ok = sum(1 for _, stt in res if not stt.get("error"))
            st.success(f"Listo: {ok}/{len(res)} preparados procesados desde {start}.")
        except Exception as e:
            st.error(f"Error al ejecutar: {e}")


# --------- 2) Calibración global (solo informativo) ---------
glob_calib = _read_global_calibration()
if glob_calib:
    st.caption(f"Calibración global activa (µm): Z={glob_calib.get('z')}, Y={glob_calib.get('y')}, X={glob_calib.get('x')}")
else:
    st.warning("No hay calibración global guardada en streamlit/calibration.json. Continuamos igualmente para Otsu/Cellpose.")


# --------- 3) Configuración de segmentación ---------
st.markdown("### Parámetros de Segmentación")
colp1, colp2, colp3 = st.columns(3)
with colp1:
    nucleus_diameter = st.number_input("Diámetro de núcleo (px)", value=30, min_value=5, step=1)
with colp2:
    dapi_channel_index = st.number_input("Índice de canal DAPI (0=primero)", value=0, min_value=0, step=1)
with colp3:
    use_gpu = st.checkbox("Usar GPU si está disponible", value=True)


# --------- 4) Acciones ---------
st.markdown("### Acciones")
recompute_otsu = st.checkbox("Forzar recalcular Otsu (sobrescribir si existe)", value=False)
recompute_cellpose = st.checkbox("Forzar recalcular Cellpose (sobrescribir si existe)", value=False)

run_otsu = st.button("🧮 Ejecutar Otsu y Guardar")
run_cellpose = st.button("🧠 Ejecutar Cellpose y Guardar")
run_both = st.button("🚀 Ejecutar Otsu + Cellpose")
open_napari_image = st.button("👁️ Abrir en Napari (solo imagen)")
open_napari_with_masks = st.button("🧪 Abrir en Napari con máscaras disponibles")


def compute_and_save_otsu(dapi_vol: np.ndarray, out_dir: Path):
    st.write("Calculando umbral de Otsu...")
    thr = threshold_otsu(dapi_vol)
    otsu_mask = (dapi_vol > thr)
    tifffile.imwrite(out_dir / "01_otsu_mask.tif", otsu_mask.astype(np.uint8))
    st.success("Máscara de Otsu guardada (01_otsu_mask.tif)")
    st.info(f"Umbral Otsu: {thr:.1f} | Fracción en máscara: {100*otsu_mask.mean():.1f}%")
    return otsu_mask


def compute_and_save_cellpose(dapi_vol_clean: np.ndarray, out_dir: Path):
    from cellpose import models
    st.write("Cargando modelo Cellpose...")
    try:
        model = models.CellposeModel(gpu=use_gpu)
    except Exception:
        st.warning("No se pudo inicializar GPU; se usará CPU.")
        model = models.CellposeModel(gpu=False)

    st.write("Ejecutando segmentación 3D con Cellpose...")
    masks, _, _ = model.eval(
        dapi_vol_clean,
        diameter=int(nucleus_diameter),
        z_axis=0,
        do_3D=True,
    )
    tifffile.imwrite(out_dir / "02_cellpose_mask.tif", masks.astype(np.uint16))
    st.success("Máscara de Cellpose guardada (02_cellpose_mask.tif)")
    return masks


def load_image_and_get_dapi():
    arr, axes = load_image_any(img_path)
    vol = reorder_to_zcyx(arr, axes)  # Z,C,Y,X
    if dapi_channel_index >= vol.shape[1]:
        raise IndexError(f"Índice de canal DAPI={dapi_channel_index} fuera de rango (n_canales={vol.shape[1]})")
    dapi = vol[:, int(dapi_channel_index), :, :]
    return dapi


if run_otsu or run_both:
    try:
        dapi_vol = load_image_and_get_dapi()
        if status["otsu"] and not recompute_otsu:
            st.info("Ya existe 01_otsu_mask.tif. Marca 'Forzar recalcular' para sobrescribir.")
        else:
            with st.spinner("Calculando Otsu..."):
                compute_and_save_otsu(dapi_vol, out_dir)
    except Exception as e:
        st.error(f"Error en Otsu: {e}")


if run_cellpose or run_both:
    try:
        dapi_vol = load_image_and_get_dapi()
        otsu_path = out_dir / "01_otsu_mask.tif"
        if otsu_path.exists() and not (run_both and recompute_otsu):
            otsu_mask = tifffile.imread(otsu_path).astype(bool)
        else:
            st.info("No se encontró Otsu previo o se forzó recalcular; generando uno temporal…")
            otsu_mask = compute_and_save_otsu(dapi_vol, out_dir)

        dapi_clean = np.where(otsu_mask, dapi_vol, 0)
        if status["cellpose"] and not recompute_cellpose:
            st.info("Ya existe 02_cellpose_mask.tif. Marca 'Forzar recalcular' para sobrescribir.")
        else:
            with st.spinner("Corriendo Cellpose (puede tardar)…"):
                compute_and_save_cellpose(dapi_clean, out_dir)
    except Exception as e:
        st.error(f"Error en Cellpose: {e}")


# --------- 5) Ver en Napari ---------
def _launch_napari(include_masks: bool):
    env = os.environ.copy()
    # Calibración para escala
    cal = _read_global_calibration()
    z = float(cal.get('z', 1.0))
    y = float(cal.get('y', 1.0))
    x = float(cal.get('x', 1.0))
    cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
    if include_masks:
        otsu_path = out_dir / "01_otsu_mask.tif"
        cellpose_path = out_dir / "02_cellpose_mask.tif"
        if otsu_path.exists():
            cmd += ["--otsu", str(otsu_path)]
        if cellpose_path.exists():
            cmd += ["--cellpose", str(cellpose_path)]
    try:
        env["NAPARI_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        subprocess.Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")

if open_napari_image:
    _launch_napari(include_masks=False)

if open_napari_with_masks:
    _launch_napari(include_masks=True)
