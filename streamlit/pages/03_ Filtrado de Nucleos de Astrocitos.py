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

# Filtro por grupo (unificado desde sidebar)
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

st.markdown("### Estado de resultados guardados")
st.write({
    "cellpose_mask": status["cellpose"],
    "gfap_filtered_mask": status["gfap_filtered"],
    "output_dir": str(out_dir.relative_to(root)),
})
def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"
st.markdown(_group_badge_html(group), unsafe_allow_html=True)

# --------- Recalcular por ámbito desde este paso ---------
with st.expander("Recalcular por ámbito desde este paso", expanded=False):
    scope = st.radio("Ámbito", options=["Preparado seleccionado", "Grupo", "Todos"], horizontal=True, key="p03_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p03_scope_group")
    if st.button("▶️ Recalcular (desde 03)", key="p03_recalc"):
        try:
            from ui.runner import run_scope, read_calibration
            cal = read_calibration(root/"streamlit"/"calibration.json")
            sc = "selected" if scope=="Preparado seleccionado" else ("group" if scope=="Grupo" else "all")
            sel = img_path if sc=="selected" else None
            res = run_scope(root, scope=sc, start_step="03", cal=cal, selected=sel, group=scope_group, overwrite_from_step=True)
            ok = sum(1 for _, stt in res if not stt.get("error"))
            st.success(f"Listo: {ok}/{len(res)} preparados procesados desde 03.")
        except Exception as e:
            st.error(f"Error al ejecutar: {e}")


# --------- 2) Canales y parámetros ---------
arr_preview, axes_prev = load_image_any(img_path)
vol_prev = reorder_to_zcyx(arr_preview, axes_prev)  # Z,C,Y,X
n_channels = vol_prev.shape[1]

st.markdown("### Selección de canales")
cc1, cc2 = st.columns(2)
cal = _read_global_calibration()
with cc1:
    gfap_default = int(cal.get("GFAP_CHANNEL_INDEX", 1 if n_channels > 1 else 0))
    gfap_default = min(max(0, gfap_default), max(0, n_channels-1))
    gfap_idx = st.number_input("Índice de canal GFAP", value=int(gfap_default), min_value=0, max_value=max(0, n_channels-1), step=1)
with cc2:
    micro_default = int(cal.get("MICROGLIA_CHANNEL_INDEX", 2 if n_channels > 2 else min(1, n_channels-1)))
    micro_default = min(max(0, micro_default), max(0, n_channels-1))
    micro_idx = st.number_input("Índice de canal Microglía", value=int(micro_default), min_value=0, max_value=max(0, n_channels-1), step=1)

st.markdown("### Parámetros de filtrado (GFAP / Microglía)")
fcol1, fcol2, fcol3 = st.columns(3)
with fcol1:
    max_dilation_iterations = st.number_input("Iteraciones máximo del anillo", value=int(cal.get("MAX_DILATION_ITERATIONS", 10)), min_value=1, step=1)
with fcol2:
    gfap_intensity_threshold = st.number_input("Umbral GFAP (intensidad)", value=int(cal.get("GFAP_INTENSITY_THRESHOLD", 40)), min_value=0, step=1)
with fcol3:
    microglia_intensity_threshold = st.number_input("Umbral Microglía (intensidad)", value=int(cal.get("MICROGLIA_INTENSITY_THRESHOLD", 200)), min_value=0, step=1)

# --- Parámetros de tamaño (volumen) ---
st.markdown("### Filtro por tamaño (volumen físico)")
vcol1, vcol2 = st.columns(2)
with vcol1:
    min_volume_um3 = st.number_input("Volumen mínimo (µm³)", value=int(cal.get("MIN_VOLUME_UM3", 75)), min_value=0, step=1)
with vcol2:
    apply_size_after_filter = st.checkbox("Aplicar filtro de volumen tras GFAP/Microglía", value=True)

save_experiment_params = st.button("💾 Guardar parámetros del experimento (sidebar)")
if save_experiment_params:
    # Unificar parámetros globales en streamlit/calibration.json
    calib_all_path = root / "streamlit" / "calibration.json"
    cur = {}
    if calib_all_path.exists():
        try:
            cur = json.loads(calib_all_path.read_text())
        except Exception:
            cur = {}
    cur.update({
        "MAX_DILATION_ITERATIONS": int(max_dilation_iterations),
        "GFAP_INTENSITY_THRESHOLD": int(gfap_intensity_threshold),
        "MICROGLIA_INTENSITY_THRESHOLD": int(microglia_intensity_threshold),
        "MIN_VOLUME_UM3": int(min_volume_um3),
    })
    calib_all_path.parent.mkdir(parents=True, exist_ok=True)
    calib_all_path.write_text(json.dumps(cur, indent=2))
    st.success(f"Parámetros del experimento guardados en {calib_all_path.relative_to(root)}")


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

    # Tasas de retención
    ret_gfap = (n_gfap / n_cellpose * 100.0) if n_cellpose else 0.0
    ret_final = (n_final / n_cellpose * 100.0) if n_cellpose else 0.0
    drop_final = n_cellpose - n_final

    # Métricas numéricas con porcentajes
    mcol1, mcol2, mcol3 = st.columns(3)
    mcol1.metric("Núcleos tras Cellpose", n_cellpose)
    mcol2.metric("Candidatos GFAP/Microglía", n_gfap, delta=f"{ret_gfap:.1f}% retenidos")
    mcol3.metric("Astrocitos finales", n_final, delta=f"−{drop_final}" if drop_final>0 else None)

    # Gráfico de barras del pipeline
    chart_df = pd.DataFrame({
        "Etapa": ["Cellpose", "GFAP/Microglía", "Final (Volumen)"],
        "Cantidad": [n_cellpose, n_gfap, n_final],
        "Retención %": [100.0, ret_gfap, ret_final],
    })
    bar = alt.Chart(chart_df).mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3).encode(
        x=alt.X("Etapa:N", sort=None),
        y=alt.Y("Cantidad:Q"),
        color=alt.Color("Etapa:N", legend=None),
        tooltip=["Etapa", "Cantidad"],
    ).properties(height=220)
    st.altair_chart(bar, use_container_width=True)

    # Línea de retención (%)
    line = alt.Chart(chart_df).mark_line(point=alt.OverlayMarkDef(size=80, filled=True)).encode(
        x=alt.X("Etapa:N", sort=None),
        y=alt.Y("Retención %:Q", scale=alt.Scale(domain=[0,100])),
        tooltip=["Etapa", alt.Tooltip("Retención %:Q", format=".1f")],
    ).properties(height=220)
    st.altair_chart(line, use_container_width=True)

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
st.caption("Esta página genera 03_gfap_microglia_filtered_mask.tif en data/processed/<preparado>/. Los parámetros se gestionan SIEMPRE en calibration.json para asegurar reproducibilidad.")
