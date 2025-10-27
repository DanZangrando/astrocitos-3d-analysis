import sys
import os
import json
import subprocess
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
from skimage.measure import regionprops
from skimage.morphology import binary_dilation
import concurrent.futures
from ui.sidebar import render_sidebar
import pandas as pd
import altair as alt


st.title("Filtrado de Núcleos de Astrocitos")
render_sidebar(show_calibration=True)

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


def reorder_to_zcyx(arr: np.ndarray, axes: str | None):
    if axes is None:
        if arr.ndim != 4:
            raise ValueError(f"Forma inesperada sin ejes: {arr.shape}")
        chan_axis = int(np.argmin(arr.shape))
        if chan_axis != 1:
            arr = np.moveaxis(arr, chan_axis, 1)
        return arr
    ax_list = list(axes)
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
    suffix = path.suffix.lower()
    if suffix in ('.tif', '.tiff'):
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            axes = getattr(series, 'axes', None)
            arr = series.asarray()
        return arr, axes
    elif suffix == '.lif':
        import readlif
        rdr = readlif.Reader(str(path))
        try:
            img0 = rdr.get_image(0)
            arr = np.asarray(img0)
        except Exception:
            series_list = rdr.getSeries() if hasattr(rdr, 'getSeries') else []
            if series_list:
                s0 = series_list[0]
                try:
                    arr = np.asarray(s0)
                except Exception as e:
                    raise RuntimeError(f"No se pudo leer la primera serie del LIF: {e}")
            else:
                raise RuntimeError("No se encontraron series en el LIF")

        if arr.ndim == 3:  # Z,Y,X
            arr = arr[:, None, :, :]
            axes = 'ZCYX'
        elif arr.ndim == 4:
            chan_axis = int(np.argmin(arr.shape))
            if chan_axis != 1:
                arr = np.moveaxis(arr, chan_axis, 1)
            axes = None
        elif arr.ndim == 5:
            arr = arr[0]
            if arr.ndim == 4:
                chan_axis = int(np.argmin(arr.shape))
                if chan_axis != 1:
                    arr = np.moveaxis(arr, chan_axis, 1)
                axes = None
            else:
                raise ValueError(f"Forma LIF no soportada: {arr.shape}")
        else:
            raise ValueError(f"Forma LIF no soportada: {arr.shape}")
        return arr, axes
    else:
        raise ValueError(f"Extensión no soportada: {suffix}")


def get_output_dir_for_image(img_path: Path) -> Path:
    base_name = img_path.stem
    out_dir = root / "data" / "processed" / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def existing_results(out_dir: Path) -> dict:
    return {
        "cellpose": (out_dir / "02_cellpose_mask.tif").exists(),
        "gfap_filtered": (out_dir / "03_gfap_microglia_filtered_mask.tif").exists(),
        "params": (out_dir / "params.json").exists(),
    }


def save_params(out_dir: Path, updates: dict):
    p = out_dir / "params.json"
    cur = {}
    if p.exists():
        try:
            cur = json.loads(p.read_text())
        except Exception:
            cur = {}
    cur.update(updates)
    p.write_text(json.dumps(cur, indent=2))


# --------- 1) Selección de imagen ---------
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.lif")])
if not files:
    st.warning("No se encontraron archivos .tif/.lif en data/raw.")
    st.stop()

labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files))), format_func=lambda i: labels[i])
img_path = files[idx]
out_dir = get_output_dir_for_image(img_path)
status = existing_results(out_dir)

st.markdown("### Estado de resultados guardados")
st.write({
    "cellpose_mask": status["cellpose"],
    "gfap_filtered_mask": status["gfap_filtered"],
    "params_json": status["params"],
    "output_dir": str(out_dir.relative_to(root)),
})


# --------- 2) Canales y parámetros ---------
arr_preview, axes_prev = load_image_any(img_path)
vol_prev = reorder_to_zcyx(arr_preview, axes_prev)  # Z,C,Y,X
n_channels = vol_prev.shape[1]

st.markdown("### Selección de canales")
cc1, cc2 = st.columns(2)
with cc1:
    gfap_idx = st.number_input("Índice de canal GFAP", value=1 if n_channels > 1 else 0, min_value=0, max_value=max(0, n_channels-1), step=1)
with cc2:
    micro_idx = st.number_input("Índice de canal Microglía", value=2 if n_channels > 2 else min(1, n_channels-1), min_value=0, max_value=max(0, n_channels-1), step=1)

st.markdown("### Parámetros de filtrado (GFAP / Microglía)")
fcol1, fcol2, fcol3 = st.columns(3)
with fcol1:
    max_dilation_iterations = st.number_input("Iteraciones máximo del anillo", value=10, min_value=1, step=1)
with fcol2:
    gfap_intensity_threshold = st.number_input("Umbral GFAP (intensidad)", value=40, min_value=0, step=1)
with fcol3:
    microglia_intensity_threshold = st.number_input("Umbral Microglía (intensidad)", value=200, min_value=0, step=1)

# --- Parámetros de tamaño (volumen) ---
st.markdown("### Filtro por tamaño (volumen físico)")
vcol1, vcol2 = st.columns(2)
with vcol1:
    min_volume_um3 = st.number_input("Volumen mínimo (µm³)", value=75, min_value=0, step=1)
with vcol2:
    apply_size_after_filter = st.checkbox("Aplicar filtro de volumen tras GFAP/Microglía", value=True)

save_experiment_params = st.button("💾 Guardar parámetros del experimento (sidebar)")
if save_experiment_params:
    exp_params_path = root / "streamlit" / "experiment_params.json"
    exp = {
        "MAX_DILATION_ITERATIONS": int(max_dilation_iterations),
        "GFAP_INTENSITY_THRESHOLD": int(gfap_intensity_threshold),
        "MICROGLIA_INTENSITY_THRESHOLD": int(microglia_intensity_threshold),
        "MIN_VOLUME_UM3": int(min_volume_um3),
    }
    exp_params_path.parent.mkdir(parents=True, exist_ok=True)
    exp_params_path.write_text(json.dumps(exp, indent=2))
    st.success(f"Parámetros del experimento guardados en {exp_params_path.relative_to(root)}")


# --------- 3) Lógica de filtrado ---------
def compute_and_save_filtering(cellpose_masks: np.ndarray, gfap_channel: np.ndarray, microglia_channel: np.ndarray, out_dir: Path):
    st.write("Iniciando filtrado combinado por GFAP / Microglía...")
    if cellpose_masks is None or cellpose_masks.max() == 0:
        raise RuntimeError("No hay máscaras de Cellpose válidas para filtrar. Generá 02_cellpose_mask.tif primero.")

    nuclei_props = regionprops(cellpose_masks)
    astrocyte_labels_candidate = []

    def _proc(nucleus):
        nucleus_mask = (cellpose_masks == nucleus.label)
        current_mask = nucleus_mask
        for _ in range(int(max_dilation_iterations)):
            dilated_mask = binary_dilation(current_mask)
            shell_mask = dilated_mask & ~current_mask
            if not np.any(shell_mask):
                break
            shell_gfap_intensity = float(gfap_channel[shell_mask].mean()) if np.any(shell_mask) else 0.0
            shell_microglia_intensity = float(microglia_channel[shell_mask].mean()) if np.any(shell_mask) else 0.0
            if shell_microglia_intensity > float(microglia_intensity_threshold):
                return None
            if shell_gfap_intensity > float(gfap_intensity_threshold):
                return nucleus.label
            current_mask = dilated_mask
        return None

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(_proc, n) for n in nuclei_props]
        for f in concurrent.futures.as_completed(futures):
            res = f.result()
            if res is not None:
                astrocyte_labels_candidate.append(res)

    gfap_filtered_mask = np.where(np.isin(cellpose_masks, astrocyte_labels_candidate), cellpose_masks, 0)
    tifffile.imwrite(out_dir / "03_gfap_microglia_filtered_mask.tif", gfap_filtered_mask.astype(np.uint16))
    save_params(out_dir, {
        "MAX_DILATION_ITERATIONS": int(max_dilation_iterations),
        "GFAP_INTENSITY_THRESHOLD": int(gfap_intensity_threshold),
        "MICROGLIA_INTENSITY_THRESHOLD": int(microglia_intensity_threshold),
        "candidates_after_filter": len(astrocyte_labels_candidate)
    })
    st.success(f"Máscara filtrada guardada (03_gfap_microglia_filtered_mask.tif). {len(astrocyte_labels_candidate)} candidatos retenidos.")
    return gfap_filtered_mask


def apply_size_filter_and_save(gfap_filtered_mask: np.ndarray, min_volume_um3: float, out_dir: Path):
    """Aplica filtro por tamaño físico y guarda 04_final_astrocytes_mask.tif"""
    cal = _read_global_calibration()
    z_um = float(cal.get('z', 1.0))
    y_um = float(cal.get('y', 1.0))
    x_um = float(cal.get('x', 1.0))
    voxel_vol_um3 = z_um * y_um * x_um
    if voxel_vol_um3 <= 0:
        st.warning("Calibración inválida; usando 1 µm³ por vóxel como fallback.")
        voxel_vol_um3 = 1.0

    min_voxels = int(np.ceil(float(min_volume_um3) / voxel_vol_um3)) if min_volume_um3 > 0 else 0
    if min_voxels <= 1:
        st.info("Umbral de volumen mínimo es <= 1 vóxel; el filtro de tamaño no eliminará objetos.")

    st.write(f"Aplicando filtro de tamaño: {min_volume_um3} µm³ → {min_voxels} vóxeles (voxel={voxel_vol_um3:.3f} µm³)")
    props = regionprops(gfap_filtered_mask)
    kept = [p.label for p in props if p.area >= min_voxels]
    final_mask = np.where(np.isin(gfap_filtered_mask, kept), gfap_filtered_mask, 0)
    tifffile.imwrite(out_dir / "04_final_astrocytes_mask.tif", final_mask.astype(np.uint16))
    save_params(out_dir, {
        "MIN_VOLUME_UM3": float(min_volume_um3),
        "MIN_VOLUME_VOXELS": int(min_voxels),
        "final_astrocyte_count": int(len(kept)),
    })
    st.success(f"Máscara final guardada (04_final_astrocytes_mask.tif). Objetos finales: {len(kept)}")
    return final_mask


st.markdown("### Ejecutar Filtrado")
run_filter = st.button("🔬 Ejecutar filtrado GFAP/Microglía y Guardar")
if run_filter:
    try:
        # Volumen completo para extraer canales
        arr, axes = load_image_any(img_path)
        vol = reorder_to_zcyx(arr, axes)
        gfap_channel = vol[:, int(gfap_idx), :, :]
        microglia_channel = vol[:, int(micro_idx), :, :]

        cellpose_path = out_dir / "02_cellpose_mask.tif"
        if not cellpose_path.exists():
            st.error("No se encontró 02_cellpose_mask.tif. Ejecutá primero la página de Otsu/Cellpose.")
        else:
            cellpose_masks = tifffile.imread(cellpose_path)
            with st.spinner("Ejecutando filtrado combinado (puede tardar)…"):
                gfap_filtered_mask = compute_and_save_filtering(cellpose_masks, gfap_channel, microglia_channel, out_dir)
            if apply_size_after_filter:
                with st.spinner("Aplicando filtro por tamaño…"):
                    apply_size_filter_and_save(gfap_filtered_mask, min_volume_um3, out_dir)
        # Mostrar métricas luego de procesar
        st.session_state["__refresh_metrics"] = True
    except Exception as e:
        st.error(f"Error en filtrado combinado: {e}")


# Botón adicional para aplicar solo el filtro de tamaño cuando ya existe 03_
apply_size_only = st.button("📦 Aplicar solo filtro por tamaño y Guardar")
if apply_size_only:
    try:
        gfap_path = out_dir / "03_gfap_microglia_filtered_mask.tif"
        if not gfap_path.exists():
            st.error("No se encontró 03_gfap_microglia_filtered_mask.tif. Ejecutá primero el filtrado GFAP/Microglía.")
        else:
            gfap_filtered_mask = tifffile.imread(gfap_path)
            with st.spinner("Aplicando filtro por tamaño…"):
                apply_size_filter_and_save(gfap_filtered_mask, min_volume_um3, out_dir)
    except Exception as e:
        st.error(f"Error en filtro por tamaño: {e}")
    else:
        # Mostrar métricas luego de aplicar tamaño
        st.session_state["__refresh_metrics"] = True


# --------- 4) Ver en Napari ---------
open_napari_image = st.button("👁️ Abrir en Napari (solo imagen)")
open_napari_with_masks = st.button("🧪 Abrir en Napari con máscaras disponibles")


def _launch_napari(include_masks: bool):
    env = os.environ.copy()
    cal = _read_global_calibration()
    z = float(cal.get('z', 1.0))
    y = float(cal.get('y', 1.0))
    x = float(cal.get('x', 1.0))
    cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
    if include_masks:
        otsu_path = out_dir / "01_otsu_mask.tif"
        cellpose_path = out_dir / "02_cellpose_mask.tif"
        gfap_path = out_dir / "03_gfap_microglia_filtered_mask.tif"
        final_path = out_dir / "04_final_astrocytes_mask.tif"
        if otsu_path.exists():
            cmd += ["--otsu", str(otsu_path)]
        if cellpose_path.exists():
            cmd += ["--cellpose", str(cellpose_path)]
        if gfap_path.exists():
            cmd += ["--gfap", str(gfap_path)]
        if final_path.exists():
            cmd += ["--final", str(final_path)]
    try:
        subprocess.Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")


if open_napari_image:
    _launch_napari(include_masks=False)

if open_napari_with_masks:
    _launch_napari(include_masks=True)


# --------- 5) Métricas y resumen ---------
st.markdown("### 📈 Métricas del pipeline y resumen")

def _load_masks_for_metrics(out_dir: Path):
    cp = out_dir / "02_cellpose_mask.tif"
    gf = out_dir / "03_gfap_microglia_filtered_mask.tif"
    fn = out_dir / "04_final_astrocytes_mask.tif"
    cellpose_masks = tifffile.imread(cp) if cp.exists() else None
    gfap_mask = tifffile.imread(gf) if gf.exists() else None
    final_mask = tifffile.imread(fn) if fn.exists() else None
    return cellpose_masks, gfap_mask, final_mask


def _compute_counts(cellpose_masks, gfap_mask, final_mask):
    n_cellpose = int(np.max(cellpose_masks)) if cellpose_masks is not None else 0
    # Contar etiquetas presentes (>0) de forma robusta
    def count_labels(lbl):
        if lbl is None:
            return 0
        unique = np.unique(lbl)
        return int((unique > 0).sum())
    n_gfap = count_labels(gfap_mask)
    n_final = count_labels(final_mask)
    return n_cellpose, n_gfap, n_final


def _compute_table(cellpose_masks, gfap_mask, final_mask):
    if cellpose_masks is None:
        return pd.DataFrame()
    # Volumen por etiqueta en voxeles usando bincount
    flat = cellpose_masks.ravel()
    max_label = int(flat.max())
    counts = np.bincount(flat, minlength=max_label + 1)
    # Conjuntos de etiquetas
    labels_all = np.arange(1, max_label + 1)
    kept_gfap = set(np.unique(gfap_mask)) - {0} if gfap_mask is not None else set()
    kept_final = set(np.unique(final_mask)) - {0} if final_mask is not None else set()

    # Calibración para volumen físico
    cal = _read_global_calibration()
    z_um = float(cal.get('z', 1.0))
    y_um = float(cal.get('y', 1.0))
    x_um = float(cal.get('x', 1.0))
    voxel_vol_um3 = z_um * y_um * x_um

    data = []
    for label in labels_all:
        vox = int(counts[label]) if label < len(counts) else 0
        um3 = float(vox) * voxel_vol_um3
        data.append({
            "label": int(label),
            "volume_voxels": int(vox),
            "volume_um3": float(um3),
            "kept_gfap": bool(label in kept_gfap) if gfap_mask is not None else None,
            "kept_final": bool(label in kept_final) if final_mask is not None else None,
            "filtered_final": (False if final_mask is None else (label not in kept_final)),
        })
    df = pd.DataFrame(data)
    return df


def render_metrics_section():
    cellpose_masks, gfap_mask, final_mask = _load_masks_for_metrics(out_dir)
    if cellpose_masks is None:
        st.info("No hay 02_cellpose_mask.tif aún. Ejecutá los pasos previos para ver métricas.")
        return

    n_cellpose, n_gfap, n_final = _compute_counts(cellpose_masks, gfap_mask, final_mask)

    # Métricas numéricas
    mcol1, mcol2, mcol3 = st.columns(3)
    mcol1.metric("Núcleos tras Cellpose", n_cellpose)
    mcol2.metric("Candidatos GFAP/Microglía", n_gfap)
    mcol3.metric("Astrocitos finales", n_final)

    # Gráfico de barras del pipeline
    chart_df = pd.DataFrame({
        "Etapa": ["Cellpose", "GFAP/Microglía", "Final (Volumen)"],
        "Cantidad": [n_cellpose, n_gfap, n_final],
    })
    bar = alt.Chart(chart_df).mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3).encode(
        x=alt.X("Etapa:N", sort=None),
        y=alt.Y("Cantidad:Q"),
        color=alt.Color("Etapa:N", legend=None),
        tooltip=["Etapa", "Cantidad"],
    ).properties(height=220)
    st.altair_chart(bar, use_container_width=True)

    # Tabla por núcleo
    df = _compute_table(cellpose_masks, gfap_mask, final_mask)
    if not df.empty:
        st.markdown("#### Detalle por núcleo")
        st.dataframe(df, use_container_width=True)
    else:
        st.info("Aún no hay datos de núcleos para tabular.")


# Botón manual para refrescar métricas
refresh = st.button("🔄 Actualizar métricas")
if refresh or st.session_state.get("__refresh_metrics", False):
    render_metrics_section()
    st.session_state["__refresh_metrics"] = False


st.markdown("---")
st.caption("Esta página genera 03_gfap_microglia_filtered_mask.tif en data/processed/<preparado>/. Los parámetros también pueden guardarse globalmente para mostrarlos en el sidebar.")
