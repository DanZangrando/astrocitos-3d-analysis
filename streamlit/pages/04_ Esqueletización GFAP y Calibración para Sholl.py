import sys
import os
import json
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
from skimage.morphology import binary_closing, binary_dilation, ball, skeletonize as _skeletonize_2d
from skimage.measure import regionprops, label as cc_label
from scipy.ndimage import distance_transform_edt
from skimage.transform import resize
from ui.sidebar import render_sidebar


st.title("Esqueletización GFAP y Calibración para Sholl")
render_sidebar(show_calibration=True)

root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
calib_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer.py"
# Intentar usar skeletonize_3d real; si no está disponible, caer a slice-wise 2D
try:
    from skimage.morphology import skeletonize_3d as _skeletonize_3d

    def skeletonize3d(arr: np.ndarray) -> np.ndarray:
        return _skeletonize_3d(arr)
except Exception:
    # Fallback: aplicar skeletonize 2D por slice (aproximación)
    def skeletonize3d(arr: np.ndarray) -> np.ndarray:
        out = np.zeros_like(arr, dtype=bool)
        for z in range(arr.shape[0]):
            out[z] = _skeletonize_2d(arr[z].astype(bool))
        return out


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


# 1) Selección de imagen y opciones de máscara
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
group = _detect_group(img_path, root)

mask_final = out_dir / "04_final_astrocytes_mask.tif"
mask_gfap = out_dir / "03_gfap_microglia_filtered_mask.tif"

mask_choice = st.radio(
    "Máscara a usar",
    options=["Final (04)", "GFAP/Microglía (03)"],
    horizontal=True,
)
mask_path = mask_final if (mask_choice == "Final (04)" and mask_final.exists()) else mask_gfap
if not mask_path.exists():
    st.error("No se encontró la máscara seleccionada. Ejecutá antes las etapas previas.")
    st.stop()
def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"
st.markdown(_group_badge_html(group), unsafe_allow_html=True)

# --------- Recalcular por ámbito desde este paso ---------
with st.expander("Recalcular por ámbito desde este paso", expanded=False):
    scope = st.radio("Ámbito", options=["Preparado seleccionado", "Grupo", "Todos"], horizontal=True, key="p04_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p04_scope_group")
    if st.button("▶️ Recalcular (desde 04)", key="p04_recalc"):
        try:
            from ui.runner import run_scope, read_calibration
            cal = read_calibration(root/"streamlit"/"calibration.json")
            sc = "selected" if scope=="Preparado seleccionado" else ("group" if scope=="Grupo" else "all")
            sel = img_path if sc=="selected" else None
            res = run_scope(root, scope=sc, start_step="04", cal=cal, selected=sel, group=scope_group, overwrite_from_step=True)
            ok = sum(1 for _, stt in res if not stt.get("error"))
            st.success(f"Listo: {ok}/{len(res)} preparados procesados desde 04.")
        except Exception as e:
            st.error(f"Error al ejecutar: {e}")

# 2) Parámetros de esqueletización calibrada
arr_preview, axes_prev = load_image_any(img_path)
vol_prev = reorder_to_zcyx(arr_preview, axes_prev)  # Z,C,Y,X
n_channels = vol_prev.shape[1]

gcol1, gcol2, gcol3 = st.columns(3)
with gcol1:
    # Prefill desde calibración global si existe
    gfap_default = int(_read_global_calibration().get("GFAP_CHANNEL_INDEX", 1 if n_channels > 1 else 0))
    gfap_default = min(max(0, gfap_default), max(0, n_channels-1))
    gfap_idx = st.number_input("Índice de canal GFAP", value=int(gfap_default), min_value=0, max_value=max(0, n_channels-1), step=1)
cal = _read_global_calibration()
z_um = float(cal.get('z', 1.0))
y_um = float(cal.get('y', 1.0))
x_um = float(cal.get('x', 1.0))
def_iso = float(min(z_um, y_um, x_um))

# Cargar parámetros globales para autocompletar UI (SIEMPRE globales)
try:
    params_global = json.loads((root / "streamlit" / "calibration.json").read_text()) if (root / "streamlit" / "calibration.json").exists() else {}
except Exception:
    params_global = {}

def _get_param(name: str, fallback):
    return params_global.get(name, fallback)

# Defaults desde JSON si existen
default_target_iso_um = float(_get_param("SKELETON_TARGET_ISO_UM", def_iso))
default_padding_um = float(_get_param("SKELETON_PADDING_UM", 2.0))
default_seed_dilate_um = float(_get_param("SKELETON_SEED_DILATE_UM", 2.0))
default_connectivity = int(_get_param("SKELETON_CONNECTIVITY", 26))
default_closing_um = float(_get_param("SKELETON_CLOSING_UM", 0.8))
default_territory_voronoi = bool(_get_param("SKELETON_TERRITORY_VORONOI", False))
default_territory_excl_um = float(_get_param("SKELETON_TERRITORY_EXCLUSION_UM", 1.0))
# Nuevos defaults: fuente para volumen de dominio y pruning topológico
default_domain_source = str(_get_param("SKELETON_DOMAIN_VOLUME_SOURCE", "gfap")).lower()
default_prune_enable = bool(_get_param("SKELETON_PRUNE_ENABLE", False))
default_prune_min_len_um = float(_get_param("SKELETON_PRUNE_MIN_LEN_UM", max(2.0, 2.0 * float(def_iso))))
default_tube_radius_um = float(_get_param("SKELETON_TUBE_RADIUS_UM", 1.5))

# Calcular default de radio máximo: preferir SKELETON_MAX_RADIUS_UM; si no, derivar de MAX_DILATION_ITERATIONS
mdi = params_global.get("MAX_DILATION_ITERATIONS")
if "SKELETON_MAX_RADIUS_UM" in params_global:
    default_max_radius_um = float(_get_param("SKELETON_MAX_RADIUS_UM", 10.0 * default_target_iso_um))
elif mdi is not None:
    try:
        default_max_radius_um = float(mdi) * float(default_target_iso_um)
    except Exception:
        default_max_radius_um = float(10.0) * float(default_target_iso_um)
else:
    default_max_radius_um = float(10.0) * float(default_target_iso_um)
with gcol2:
    target_iso_um = st.number_input("Vóxel isotrópico objetivo (µm)", value=float(default_target_iso_um), min_value=0.05, step=0.01, format="%.2f")
with gcol3:
    padding_um = st.number_input("Padding ROI (µm)", value=float(default_padding_um), min_value=0.0, step=0.5)

mcol1, mcol2, mcol3 = st.columns(3)
thresh_options = ["Otsu (ROI)", "Manual"]
_saved_thresh_mode = str(_get_param("SKELETON_THRESHOLD_MODE", thresh_options[0]))
_saved_thresh_mode_l = _saved_thresh_mode.lower()
_def_idx = 0 if _saved_thresh_mode_l.startswith("otsu") else 1 if _saved_thresh_mode_l.startswith("manual") else 0
with mcol1:
    thresh_mode = st.selectbox("Umbral GFAP", options=thresh_options, index=_def_idx)
with mcol2:
    default_manual_thr = float(_get_param("SKELETON_MANUAL_THRESHOLD", 50))
    manual_thr = st.number_input("Umbral manual", value=float(default_manual_thr), min_value=0.0, step=1.0, format="%.0f")
with mcol3:
    closing_um = st.number_input("Cierre morfológico (µm)", value=float(default_closing_um), min_value=0.0, step=0.2, format="%.2f")

# Parámetro para capturar GFAP conectado más allá del núcleo
seedcol1, seedcol2 = st.columns(2)
with seedcol1:
    seed_dilate_um = st.number_input("Dilatación de semilla (µm)", value=float(default_seed_dilate_um), min_value=0.0, step=0.5, help="Se dilata la máscara del núcleo en el espacio isotrópico para seleccionar el GFAP conectado a ese núcleo.")
with seedcol2:
    connectivity = st.selectbox("Conectividad 3D", options=[6, 26], index=(1 if default_connectivity == 26 else 0), help="Usada para seleccionar el componente de GFAP conectado a la semilla.")

# Radio máximo desde núcleo (usar defaults ya calculados desde JSON/MDI)
radiuscol = st.columns(1)[0]
with radiuscol:
    max_radius_um = st.number_input("Radio máximo desde núcleo (µm)", value=float(round(default_max_radius_um, 2)), min_value=0.0, step=0.5, help="Limita el GFAP a una distancia radial equivalente al cascarón usado en la etapa de filtrado.")

# Territorios de exclusión entre células
tcol1, tcol2 = st.columns(2)
with tcol1:
    territory_voronoi = st.checkbox("Territorios por proximidad al núcleo (Voronoi)", value=bool(default_territory_voronoi), help="Restringe el GFAP de cada célula a su territorio más cercano al núcleo.")
with tcol2:
    territory_excl_um = st.number_input("Zona de exclusión en frontera (µm)", value=float(default_territory_excl_um), min_value=0.0, step=0.2, format="%.2f", help="Crea un ‘gap’ alrededor de los límites entre territorios para evitar entrelazados ambíguos.")

# Fuente para volumen del dominio y opciones de pruning (avanzado)
vcol1, vcol2 = st.columns(2)
with vcol1:
    domain_volume_source = st.selectbox(
        "Volumen de dominio por célula",
        options=["GFAP conectado", "Territorio Voronoi"],
        index=(0 if default_domain_source.startswith("gfap") else 1),
        help="Define de dónde calcular el volumen del dominio astrocitario (µm³)."
    )
with vcol2:
    with st.expander("Opciones avanzadas: pruning topológico", expanded=False):
        prune_enable = st.checkbox("Pruning: eliminar espículas cortas", value=bool(default_prune_enable))
        prune_min_len_um = st.number_input(
            "Longitud mínima de espícula a eliminar (µm)",
            value=float(default_prune_min_len_um), min_value=0.0, step=0.5,
            help="Se eliminan iterativamente ramas terminales con longitud menor a este umbral antes del resumen."
        )

# Radio del tubo alrededor del esqueleto (para muestrear intensidad GFAP)
tube_radius_um = st.number_input(
    "Radio del tubo alrededor del esqueleto (µm)",
    value=float(default_tube_radius_um), min_value=0.0, step=0.2, format="%.2f",
    help="Se mide la señal de GFAP dentro de un tubo centrado en el esqueleto con este radio (en isotrópico)."
)

save_skel_params = st.button("💾 Guardar parámetros del skeleton (sidebar)")
if save_skel_params:
    exp_params_path = root / "streamlit" / "calibration.json"
    exp = params_global.copy()
    _mode_str = "manual" if str(thresh_mode).lower().startswith("manual") else "otsu"
    exp.update({
        "SKELETON_TARGET_ISO_UM": float(target_iso_um),
        "SKELETON_PADDING_UM": float(padding_um),
        "SKELETON_SEED_DILATE_UM": float(seed_dilate_um),
        "SKELETON_CONNECTIVITY": int(connectivity),
        "SKELETON_CLOSING_UM": float(closing_um),
        "SKELETON_MAX_RADIUS_UM": float(max_radius_um),
        "SKELETON_THRESHOLD_MODE": _mode_str,
        "SKELETON_MANUAL_THRESHOLD": float(manual_thr),
        "SKELETON_TERRITORY_VORONOI": bool(territory_voronoi),
        "SKELETON_TERRITORY_EXCLUSION_UM": float(territory_excl_um),
        "SKELETON_DOMAIN_VOLUME_SOURCE": ("gfap" if domain_volume_source.lower().startswith("gfap") else "voronoi"),
        "SKELETON_PRUNE_ENABLE": bool(prune_enable),
        "SKELETON_PRUNE_MIN_LEN_UM": float(prune_min_len_um),
        "SKELETON_TUBE_RADIUS_UM": float(tube_radius_um),
    })
    exp_params_path.parent.mkdir(parents=True, exist_ok=True)
    exp_params_path.write_text(json.dumps(exp, indent=2))
    st.success(f"Parámetros del skeleton guardados en {exp_params_path.relative_to(root)}")

conflict_resolve = st.checkbox(
    "Resolver solapamientos por cercanía al núcleo (µm)",
    value=True,
    help="Si dos esqueletos se superponen, se asigna el vóxel al astro más cercano al núcleo en distancia física.")

run_skel = st.button("🕸️ Esqueletizar y Guardar")
open_napari = st.button("👁️ Abrir en Napari (con skeleton)")


# 3) Esqueletización por etiqueta con re-muestreo isotrópico

def _zoom_factors(z_um, y_um, x_um, target_iso_um):
    return (
        max(z_um / target_iso_um, 1e-6),
        max(y_um / target_iso_um, 1e-6),
        max(x_um / target_iso_um, 1e-6),
    )


def _resize_to_isotropic(vol, factors, order=1, preserve_range=True):
    zf, yf, xf = factors
    new_shape = (
        int(round(vol.shape[0] * zf)),
        int(round(vol.shape[1] * yf)),
        int(round(vol.shape[2] * xf)),
    )
    if any(s <= 0 for s in new_shape):
        raise ValueError("Forma objetivo inválida para remuestreo")
    return resize(
        vol,
        new_shape,
        order=order,
        preserve_range=preserve_range,
        mode="edge",
        anti_aliasing=(order != 0),
    )


def _resize_to_original(vol_iso, orig_shape, order=0):
    return resize(
        vol_iso,
        orig_shape,
        order=order,
        preserve_range=True,
        mode="edge",
        anti_aliasing=(order != 0),
    )


def _um_to_vox(um_val, voxel_um):
    return int(round(max(0.0, float(um_val)) / max(voxel_um, 1e-6)))


def run_skeletonization():
    # Datos completos
    arr, axes = load_image_any(img_path)
    vol = reorder_to_zcyx(arr, axes)
    gfap = vol[:, int(gfap_idx), :, :]
    mask = tifffile.imread(mask_path)
    if mask.shape != gfap.shape:
        st.warning("La máscara y el volumen no coinciden en forma; se intentará adaptar si es posible.")
        # Intento forzado si canal GFAP tiene misma Z,Y,X
        if mask.ndim == 3 and gfap.shape == mask.shape:
            pass
        else:
            raise RuntimeError("Dimensiones incompatibles entre GFAP y máscara")

    # Preparar salidas
    skel_dir = out_dir / "skeletons"
    skel_dir.mkdir(parents=True, exist_ok=True)
    combined_labels = np.zeros_like(mask, dtype=np.uint16)
    combined_dist_um = np.full(mask.shape, np.inf, dtype=np.float32) if conflict_resolve else None

    # Centroides de núcleos en coordenadas físicas (µm) para resolver solapamientos
    label_centroid_um = {}
    if conflict_resolve:
        for p in regionprops(mask.astype(np.int32)):
            lab_id = int(p.label)
            cz, cy, cx = p.centroid  # en índices de vóxel (Z,Y,X)
            label_centroid_um[lab_id] = (
                float(cz) * float(z_um),
                float(cy) * float(y_um),
                float(cx) * float(x_um),
            )

    # Parámetros calibrados
    zf, yf, xf = _zoom_factors(z_um, y_um, x_um, target_iso_um)
    closing_r_vox = _um_to_vox(closing_um, target_iso_um) if closing_um > 0 else 0
    seed_r_vox = _um_to_vox(seed_dilate_um, target_iso_um) if seed_dilate_um > 0 else 0
    max_r_vox = _um_to_vox(max_radius_um, target_iso_um) if max_radius_um > 0 else 0
    excl_r_vox = _um_to_vox(territory_excl_um, target_iso_um) if 'territory_excl_um' in globals() else 0
    tube_r_vox = _um_to_vox(tube_radius_um, target_iso_um) if 'tube_radius_um' in globals() else 0

    # Padding en voxeles originales (anisotrópicos) por eje
    # Importante: incluir el radio máximo solicitado, de lo contrario el recorte de la ROI
    # puede truncar procesos largos aunque max_radius_um sea alto.
    pad_z = _um_to_vox(padding_um + max_radius_um, z_um)
    pad_y = _um_to_vox(padding_um + max_radius_um, y_um)
    pad_x = _um_to_vox(padding_um + max_radius_um, x_um)

    labels = np.unique(mask)
    labels = labels[labels > 0]

    metrics = []

    def _prune_short_spurs(skel: np.ndarray, min_len_um: float, vox_um: float) -> np.ndarray:
        """Elimina iterativamente ramas terminales más cortas que min_len_um.
        skel: bool 3D skeleton en vóxeles isotrópicos.
        vox_um: tamaño de vóxel isotrópico en µm.
        """
        from scipy.ndimage import convolve
        sk = skel.copy()
        if min_len_um <= 0:
            return sk
        # 26 vecinos
        kernel = np.ones((3, 3, 3), dtype=np.int16)
        kernel[1, 1, 1] = 0
        # Offsets y sus distancias euclidianas en unidades de vóxel
        offsets = []
        dists = []
        for dz in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if dz == 0 and dy == 0 and dx == 0:
                        continue
                    offsets.append((dz, dy, dx))
                    dists.append(float(np.sqrt(dz*dz + dy*dy + dx*dx)))
        offsets = np.array(offsets, dtype=np.int8)
        dists = np.array(dists, dtype=np.float32) * float(vox_um)

        def neighbors_of(p):
            z, y, x = p
            nbrs = []
            for (dz, dy, dx) in offsets:
                zz, yy, xx = z + int(dz), y + int(dy), x + int(dx)
                if 0 <= zz < sk.shape[0] and 0 <= yy < sk.shape[1] and 0 <= xx < sk.shape[2] and sk[zz, yy, xx]:
                    nbrs.append((zz, yy, xx))
            return nbrs

        while True:
            deg = convolve(sk.astype(np.uint8), kernel, mode='constant', cval=0)
            endpoints = np.argwhere((sk > 0) & (deg == 1))
            if endpoints.size == 0:
                break
            to_remove = []
            # Evaluar cada endpoint de forma independiente
            for pz, py, px in endpoints:
                path = [(int(pz), int(py), int(px))]
                length_um = 0.0
                prev = None
                curr = (int(pz), int(py), int(px))
                while True:
                    nbrs = neighbors_of(curr)
                    if prev is not None and prev in nbrs:
                        nbrs.remove(prev)
                    deg_curr = int(deg[curr])
                    # Si alcanzamos otro endpoint sin pasar por unión, salimos
                    if deg_curr <= 1:
                        break
                    # Si alcanzamos unión (>=3), detener antes de incluir la unión
                    if deg_curr >= 3:
                        break
                    # Continuar a lo largo del único vecino disponible (deg==2)
                    if len(nbrs) == 0:
                        break
                    nxt = nbrs[0]
                    # distancia paso actual
                    dz = abs(nxt[0] - curr[0]); dy = abs(nxt[1] - curr[1]); dx = abs(nxt[2] - curr[2])
                    step_um = float(np.sqrt(dz*dz + dy*dy + dx*dx)) * float(vox_um)
                    length_um += step_um
                    path.append(nxt)
                    prev, curr = curr, nxt
                if length_um < float(min_len_um):
                    # remover todo el path excepto la unión final (curr) si es unión
                    # detectar si curr es unión
                    if int(deg[curr]) >= 3 and len(path) > 0:
                        to_remove.extend(path[:-1])
                    else:
                        to_remove.extend(path)
            if not to_remove:
                break
            rr = tuple(np.array(to_remove).T)
            sk[rr] = False
        return sk
    for lab in labels:
        terr_lab = None  # reiniciar territorio por ROI/etiqueta
        reg_mask = (mask == lab)
        if not np.any(reg_mask):
            continue
        props = regionprops(reg_mask.astype(np.uint8))
        if not props:
            continue
        # Usamos la bbox de la región
        z0, y0, x0, z1, y1, x1 = (*props[0].bbox[0:3], *props[0].bbox[3:6])
        # Expandir con padding y clip a límites
        z0 = max(0, z0 - pad_z); y0 = max(0, y0 - pad_y); x0 = max(0, x0 - pad_x)
        z1 = min(mask.shape[0], z1 + pad_z); y1 = min(mask.shape[1], y1 + pad_y); x1 = min(mask.shape[2], x1 + pad_x)

        roi_mask = reg_mask[z0:z1, y0:y1, x0:x1]
        roi_gfap = gfap[z0:z1, y0:y1, x0:x1].astype(np.float32)
        roi_all_labels = mask[z0:z1, y0:y1, x0:x1]

        # Remuestreo a isotrópico
        roi_mask_iso = _resize_to_isotropic(roi_mask.astype(np.uint8), (zf, yf, xf), order=0) > 0.5
        roi_gfap_iso = _resize_to_isotropic(roi_gfap, (zf, yf, xf), order=1)
        roi_all_labels_iso = _resize_to_isotropic(roi_all_labels.astype(np.int32), (zf, yf, xf), order=0).astype(np.int32)

        # Umbral en ROI
        if thresh_mode.startswith("Otsu"):
            from skimage.filters import threshold_otsu
            thr = threshold_otsu(roi_gfap_iso[roi_mask_iso]) if np.any(roi_mask_iso) else threshold_otsu(roi_gfap_iso)
        else:
            thr = float(manual_thr)
        # 1) Umbral GFAP en ROI sin limitarlo a la máscara de núcleo
        bin_gfap_iso = (roi_gfap_iso >= thr)

        # 2) Cierre morfológico opcional para conectar hebras
        if closing_r_vox > 0:
            try:
                bin_gfap_iso = binary_closing(bin_gfap_iso, ball(int(closing_r_vox)))
            except Exception:
                pass

        # 3) Seleccionar solo el GFAP conectado a la semilla (núcleo dilatado) en espacio isotrópico
        seed_iso = roi_mask_iso.copy()
        if seed_r_vox > 0:
            try:
                seed_iso = binary_dilation(seed_iso, ball(int(seed_r_vox)))
            except Exception:
                pass
        if np.any(bin_gfap_iso):
            cc = cc_label(bin_gfap_iso.astype(np.uint8), connectivity=1 if connectivity == 6 else 3)
            keep_ids = np.unique(cc[seed_iso > 0])
            keep_ids = keep_ids[keep_ids > 0] if keep_ids.size else keep_ids
            bin_iso = np.isin(cc, keep_ids)
        else:
            bin_iso = np.zeros_like(bin_gfap_iso, dtype=bool)

        # 3.5) Territorios por proximidad al núcleo + zona de exclusión en frontera
        if 'territory_voronoi' in globals() and territory_voronoi:
            labs_here = np.unique(roi_all_labels_iso)
            labs_here = labs_here[labs_here > 0]
            if labs_here.size > 0:
                # Calcular centroides en coordenadas de vóxel isotrópicas
                props_here = regionprops(roi_all_labels_iso)
                centroids = {int(p.label): p.centroid for p in props_here}
                # Greedy Voronoi: mantener mejor distancia por vóxel
                Z, Y, X = roi_all_labels_iso.shape
                zz = np.arange(Z, dtype=np.float32)[:, None, None]
                yy = np.arange(Y, dtype=np.float32)[None, :, None]
                xx = np.arange(X, dtype=np.float32)[None, None, :]
                best_d2 = np.full((Z, Y, X), np.inf, dtype=np.float32)
                best_lab = np.zeros((Z, Y, X), dtype=np.int32)
                for l in labs_here:
                    cz, cy, cx = centroids.get(int(l), (None, None, None))
                    if cz is None:
                        continue
                    dz = zz - float(cz)
                    dy = yy - float(cy)
                    dx = xx - float(cx)
                    d2 = dz*dz + dy*dy + dx*dx
                    better = d2 < best_d2
                    best_lab[better] = int(l)
                    best_d2[better] = d2[better]
                terr_lab = (best_lab == int(lab))
                if excl_r_vox and excl_r_vox > 0:
                    # Gap desde la frontera: mantener solo interior a distancia >= excl_r_vox
                    from scipy.ndimage import distance_transform_edt as _edt
                    interior = _edt(terr_lab) >= float(excl_r_vox)
                    terr_lab = terr_lab & interior
                # Restringir GFAP al territorio (y posible gap)
                bin_iso &= terr_lab

        # 4) Limitar a radio máximo desde la semilla (distancia euclidiana en isotrópico)
        if max_r_vox > 0 and np.any(bin_iso):
            dist = distance_transform_edt(~seed_iso)
            bin_iso &= (dist <= float(max_r_vox))

        # Skeleton 3D (o fallback 2D slice-wise si no disponible)
        skel_iso = skeletonize3d(bin_iso.astype(np.uint8)) > 0
        # Pruning topológico opcional
        if 'prune_enable' in globals() and prune_enable and float(prune_min_len_um) > 0:
            try:
                skel_iso = _prune_short_spurs(skel_iso, float(prune_min_len_um), float(target_iso_um))
            except Exception:
                pass

        # Métricas simples: conteo de vóxeles de esqueleto y longitud aprox (en µm)
        skel_voxels = int(np.count_nonzero(skel_iso))
        approx_length_um = skel_voxels * float(target_iso_um)
        # Volúmenes de dominio
        gfap_connected_vox = int(np.count_nonzero(bin_iso))
        gfap_connected_um3 = float(gfap_connected_vox) * float(target_iso_um) ** 3
        voronoi_um3 = np.nan
        if 'territory_voronoi' in globals() and territory_voronoi:
            try:
                # terr_lab se define si territory_voronoi estuvo activo
                vor_vox = int(np.count_nonzero(terr_lab)) if 'terr_lab' in locals() else 0
                voronoi_um3 = float(vor_vox) * float(target_iso_um) ** 3
            except Exception:
                voronoi_um3 = np.nan
        # Selección de volumen de dominio
        if domain_volume_source.lower().startswith("gfap"):
            domain_um3 = gfap_connected_um3
        else:
            domain_um3 = voronoi_um3

        # Medición de señal GFAP alrededor del esqueleto (tubo de radio fijo)
        tube_sum_intensity = np.nan
        tube_mean_intensity = np.nan
        tube_voxels = 0
        tube_volume_um3 = np.nan
        if tube_r_vox > 0 and np.any(skel_iso):
            try:
                # Distancia a esqueleto en espacio isotrópico (en vóxeles)
                dist_to_skel = distance_transform_edt(~skel_iso)
                tube_mask = dist_to_skel <= float(tube_r_vox)
                # Medir intensidad por encima del mismo umbral usado para binarizar
                tube_sel = tube_mask & (roi_gfap_iso >= float(thr))
                tube_voxels = int(np.count_nonzero(tube_sel))
                tube_sum_intensity = float(np.sum(roi_gfap_iso[tube_sel])) if tube_voxels > 0 else 0.0
                tube_mean_intensity = float(np.mean(roi_gfap_iso[tube_sel])) if tube_voxels > 0 else 0.0
                tube_volume_um3 = float(tube_voxels) * float(target_iso_um) ** 3
            except Exception:
                pass

        # Radio local (espesor) estimado en puntos del esqueleto a partir del GFAP conectado
        local_radius_um_mean = np.nan
        local_radius_um_median = np.nan
        local_radius_um_p95 = np.nan
        if np.any(bin_iso) and np.any(skel_iso):
            try:
                dt_gfap = distance_transform_edt(bin_iso)
                loc_r_um_vals = dt_gfap[skel_iso].astype(np.float32) * float(target_iso_um)
                if loc_r_um_vals.size > 0:
                    local_radius_um_mean = float(np.mean(loc_r_um_vals))
                    local_radius_um_median = float(np.median(loc_r_um_vals))
                    local_radius_um_p95 = float(np.percentile(loc_r_um_vals, 95))
            except Exception:
                pass

        # Reescalar a forma original del ROI
        skel_roi = _resize_to_original(skel_iso.astype(np.uint8), roi_mask.shape, order=0) > 0.5

        # Pegar en canvas combinado con resolución de conflictos opcional
        target_slice = combined_labels[z0:z1, y0:y1, x0:x1]
        if conflict_resolve:
            dist_slice = combined_dist_um[z0:z1, y0:y1, x0:x1]
            # Distancia física de cada vóxel del ROI al centroide del núcleo actual
            cz_um, cy_um, cx_um = label_centroid_um.get(int(lab), (None, None, None))
            if cz_um is None:
                # Sin centroide: fallback a sobrescribir simple
                target_slice[skel_roi] = lab
            else:
                # Construir mallas 1D y broadcast para evitar alta memoria
                zz = (np.arange(z0, z1, dtype=np.float32) * float(z_um))[:, None, None]
                yy = (np.arange(y0, y1, dtype=np.float32) * float(y_um))[None, :, None]
                xx = (np.arange(x0, x1, dtype=np.float32) * float(x_um))[None, None, :]
                # Distancia euclidiana en µm
                dist_um = np.sqrt((zz - cz_um) ** 2 + (yy - cy_um) ** 2 + (xx - cx_um) ** 2)
                # Voxeles candidatos del esqueleto actual
                cand = skel_roi
                # Aceptar si estaba vacío o si estamos más cerca que lo asignado previamente
                better = cand & ((target_slice == 0) | (dist_um < dist_slice))
                target_slice[better] = lab
                dist_slice[better] = dist_um[better]
        else:
            # Asignación simple (último gana)
            target_slice[skel_roi] = lab

        # Guardar por etiqueta
        tifffile.imwrite(skel_dir / f"astro_{int(lab)}_skeleton_roi.tif", skel_roi.astype(np.uint8))
        metrics.append({
            "label": int(lab),
            "group": group,
            "skeleton_voxels_iso": skel_voxels,
            "approx_length_um": float(approx_length_um),
            "threshold": float(thr),
            "max_radius_um": float(max_radius_um),
            "gfap_connected_volume_um3": float(gfap_connected_um3),
            "voronoi_territory_volume_um3": float(voronoi_um3) if not (isinstance(voronoi_um3, float) and np.isnan(voronoi_um3)) else np.nan,
            "domain_volume_um3": float(domain_um3) if domain_um3 == domain_um3 else np.nan,
            "tube_radius_um": float(tube_radius_um),
            "tube_voxels": int(tube_voxels),
            "tube_volume_um3": float(tube_volume_um3) if tube_volume_um3 == tube_volume_um3 else np.nan,
            "tube_mean_intensity": float(tube_mean_intensity) if tube_mean_intensity == tube_mean_intensity else np.nan,
            "tube_sum_intensity": float(tube_sum_intensity) if tube_sum_intensity == tube_sum_intensity else np.nan,
            "tube_sum_intensity_per_um": (float(tube_sum_intensity) / float(approx_length_um)) if approx_length_um > 0 else np.nan,
            "local_radius_um_mean": float(local_radius_um_mean) if local_radius_um_mean == local_radius_um_mean else np.nan,
            "local_radius_um_median": float(local_radius_um_median) if local_radius_um_median == local_radius_um_median else np.nan,
            "local_radius_um_p95": float(local_radius_um_p95) if local_radius_um_p95 == local_radius_um_p95 else np.nan,
        })

    # Guardar combinado
    out_skel_labels = out_dir / "05_skeleton_labels.tif"
    tifffile.imwrite(out_skel_labels, combined_labels.astype(np.uint16))

    # Guardar métricas
    import pandas as _pd
    df = _pd.DataFrame(metrics)
    if not df.empty:
        df.to_csv(skel_dir / "summary.csv", index=False)
    st.success(f"Esqueletización completada. Guardado: {out_skel_labels.relative_to(root)} y skeletons/summary.csv")


if run_skel:
    try:
        run_skeletonization()
    except Exception as e:
        st.error(f"Error en esqueletización: {e}")


if open_napari:
    try:
        # Lanzar Napari con skeleton
        cal = _read_global_calibration()
        z = float(cal.get('z', 1.0)); y = float(cal.get('y', 1.0)); x = float(cal.get('x', 1.0))
        out_skel_labels = out_dir / "05_skeleton_labels.tif"
        cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
        if mask_path.exists():
            # pasar final si es la seleccionada
            if mask_path.name.startswith("04_"):
                cmd += ["--final", str(mask_path)]
            else:
                cmd += ["--gfap", str(mask_path)]
        if out_skel_labels.exists():
            cmd += ["--skeleton", str(out_skel_labels)]
        env = os.environ.copy()
        env["NAPARI_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        __import__("subprocess").Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")

# --------- 4) Métricas y visualización integradas ---------
st.markdown("---")
st.markdown("### 📈 Métricas del esqueleto por astrocito")
min_branch_len_um = st.number_input(
    "Longitud mínima de rama (µm) para métricas",
    value=float(max(2.0, 2.0 * float(target_iso_um))), min_value=0.0, step=0.5,
    help="Ramas más cortas que este umbral se consideran espículas/artefactos y se excluyen del conteo.")

def _compute_skeleton_metrics(out_dir: Path, mask_path: Path, calib: dict, target_iso_um: float, min_branch_len_um: float) -> tuple:
    import pandas as _pd
    z_um = float(calib.get('z', 1.0)); y_um = float(calib.get('y', 1.0)); x_um = float(calib.get('x', 1.0))
    out_skel_labels = out_dir / "05_skeleton_labels.tif"
    if not out_skel_labels.exists():
        return _pd.DataFrame(), None
    skel_lab_img = tifffile.imread(out_skel_labels)
    labels = np.unique(skel_lab_img); labels = labels[labels > 0]
    # Centroides del núcleo por label (para radio máximo)
    try:
        base_mask = tifffile.imread(mask_path)
    except Exception:
        base_mask = None
    centroids = {}
    if base_mask is not None and base_mask.shape == skel_lab_img.shape:
        try:
            for p in regionprops(base_mask.astype(np.int32)):
                cz, cy, cx = p.centroid
                centroids[int(p.label)] = (cz * z_um, cy * y_um, cx * x_um)
        except Exception:
            pass
    # Calcular métricas
    rows = []
    skan_ok = True
    try:
        from skan import Skeleton, summarize as sk_summarize
    except Exception:
        skan_ok = False
    for lab in labels:
        vox = (skel_lab_img == lab)
        # Afinar para asegurar 1 vóxel de grosor (evita sobrecuentas en Skan)
        try:
            from skimage.morphology import skeletonize_3d as _skel3d
            vox_thin = _skel3d(vox.astype(np.uint8)) > 0
        except Exception:
            vox_thin = vox
        skel_vox = int(np.count_nonzero(vox_thin))
        approx_len_um = float(skel_vox) * float(target_iso_um)
        total_len_um = approx_len_um
        n_branches = None; mean_branch_len_um = None; endpoints = None; junctions = None
        raw_n_branches = None
        if skan_ok:
            try:
                sk = Skeleton(vox_thin.astype(np.uint8), spacing=(z_um, y_um, x_um))
                dfb = sk_summarize(sk, separator="_")
                if not dfb.empty:
                    # Filtrar espículas cortas
                    dist = dfb["branch-distance"].to_numpy()
                    raw_n_branches = int(dfb.shape[0])
                    keep = dist >= float(min_branch_len_um)
                    total_len_um = float(np.nansum(dist[keep]))
                    n_branches = int(np.count_nonzero(keep))
                    mean_branch_len_um = float(np.nanmean(dist[keep])) if np.any(keep) else 0.0
                    # Aproximación: contar branches con tipo 1 (endpoint-junction)
                    if "branch-type" in dfb.columns:
                        endpoints = int(np.count_nonzero((dfb["branch-type"].to_numpy() == 1) & keep))
                    # Heurística para junctions: ramas con tipo 2 o más
                    if "branch-type" in dfb.columns:
                        junctions = int(np.count_nonzero((dfb["branch-type"].to_numpy() >= 2) & keep))
            except Exception:
                pass
        # Radio máximo desde el núcleo hacia la punta más distante del skeleton
        max_radius_um = None
        if lab in centroids:
            cz_um, cy_um, cx_um = centroids[lab]
            zz, yy, xx = np.where(vox)
            if zz.size:
                d_um = np.sqrt(((zz * z_um) - cz_um) ** 2 + ((yy * y_um) - cy_um) ** 2 + ((xx * x_um) - cx_um) ** 2)
                max_radius_um = float(np.max(d_um))
        density = (n_branches / total_len_um * 100.0) if (n_branches and total_len_um and total_len_um > 0) else np.nan
        rows.append({
            "label": int(lab),
            "skeleton_voxels": skel_vox,
            "total_length_um": float(total_len_um),
            "approx_length_um": float(approx_len_um),
            "mean_branch_len_um": mean_branch_len_um if mean_branch_len_um is not None else np.nan,
            "n_branches": n_branches if n_branches is not None else np.nan,
            "raw_n_branches": raw_n_branches if raw_n_branches is not None else np.nan,
            "endpoints": endpoints if endpoints is not None else np.nan,
            "junctions": junctions if junctions is not None else np.nan,
            "max_radius_um": max_radius_um if max_radius_um is not None else np.nan,
            "branch_density_per_100um": float(density) if not np.isnan(density) else np.nan,
        })
    dfm = _pd.DataFrame(rows).sort_values("label") if rows else _pd.DataFrame()
    return dfm, out_skel_labels


def _render_metrics_table_and_charts(out_dir: Path, mask_path: Path, target_iso_um: float):
    cal = _read_global_calibration()
    import pandas as _pd
    # Tomar summary.csv si existe y fusionar con métricas enriquecidas
    df_sum = None
    try:
        p = out_dir / "skeletons" / "summary.csv"
        if p.exists():
            df_sum = _pd.read_csv(p)
    except Exception:
        df_sum = None
    dfm, sk_path = _compute_skeleton_metrics(out_dir, mask_path, cal, target_iso_um, float(min_branch_len_um))
    if df_sum is not None and not dfm.empty:
        df = dfm.merge(df_sum, on="label", how="left", suffixes=("", "_old"))
    else:
        df = dfm if not dfm.empty else df_sum
    if df is None or df.empty:
        st.info("No hay resultados de esqueletización aún. Ejecutá 'Esqueletizar y Guardar'.")
        return
    # Resumen moderno: tarjetas con métricas clave
    try:
        n_cells = int(df.shape[0])
        mlen = float(df["total_length_um"].median()) if "total_length_um" in df.columns and not df["total_length_um"].empty else None
        mint_per_um = float(df["tube_sum_intensity_per_um"].median()) if "tube_sum_intensity_per_um" in df.columns and df["tube_sum_intensity_per_um"].notna().any() else None
        mrad = float(df["local_radius_um_median"].median()) if "local_radius_um_median" in df.columns and df["local_radius_um_median"].notna().any() else None
        mvol = float(df["domain_volume_um3"].median()) if "domain_volume_um3" in df.columns and df["domain_volume_um3"].notna().any() else None
    except Exception:
        n_cells = df.shape[0]
        mlen = mint_per_um = mrad = mvol = None
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Astrocitos (n)", n_cells)
    c2.metric("Longitud mediana (µm)", f"{mlen:.1f}" if isinstance(mlen, float) else "—")
    c3.metric("Intensidad/µm mediana", f"{mint_per_um:.2f}" if isinstance(mint_per_um, float) else "—")
    c4.metric("Volumen dominio mediano (µm³)", f"{mvol:.0f}" if isinstance(mvol, float) else "—")
    st.caption("El radio local se muestra en tablas/gráficos; usá el inspector para más detalle.")
    st.dataframe(df, use_container_width=True)
    # Gráficas
    try:
        import altair as alt
        base = alt.Chart(df)
        charts = []
        charts.append(base.mark_bar().encode(x=alt.X("total_length_um:Q", bin=alt.Bin(maxbins=30), title="Longitud total (µm)"), y="count()").properties(height=200))
        if "n_branches" in df.columns:
            charts.append(base.mark_bar().encode(x=alt.X("n_branches:Q", bin=alt.Bin(maxbins=20), title="# ramas (filtradas)"), y="count()").properties(height=200))
            charts.append(base.mark_point(opacity=0.7).encode(x=alt.X("n_branches:Q", title="# ramas"), y=alt.Y("total_length_um:Q", title="Longitud total (µm)"), tooltip=["label","total_length_um","n_branches","max_radius_um"]).properties(height=260))
        # Nuevas gráficas: volumen de dominio si está disponible
        if "domain_volume_um3" in df.columns and df["domain_volume_um3"].notna().any():
            charts.append(base.mark_bar().encode(x=alt.X("domain_volume_um3:Q", bin=alt.Bin(maxbins=30), title="Volumen de dominio (µm³)"), y="count()").properties(height=200))
            charts.append(base.mark_point(opacity=0.7).encode(x=alt.X("domain_volume_um3:Q", title="Volumen dominio (µm³)"), y=alt.Y("total_length_um:Q", title="Longitud total (µm)"), tooltip=["label","domain_volume_um3","total_length_um"]).properties(height=260))
        # Gráficas de señal en el tubo
        if "tube_sum_intensity_per_um" in df.columns and df["tube_sum_intensity_per_um"].notna().any():
            charts.append(base.mark_bar().encode(x=alt.X("tube_sum_intensity_per_um:Q", bin=alt.Bin(maxbins=30), title="Intensidad en tubo por µm"), y="count()").properties(height=200))
            charts.append(base.mark_point(opacity=0.7).encode(x=alt.X("tube_sum_intensity_per_um:Q", title="Intensidad/µm"), y=alt.Y("total_length_um:Q", title="Longitud total (µm)"), tooltip=["label","tube_sum_intensity_per_um","total_length_um"]).properties(height=260))
        st.altair_chart(alt.vconcat(*charts).resolve_scale(x='independent'), use_container_width=True)
    except Exception:
        st.bar_chart(df.set_index("label")["total_length_um"], height=240)
        if "n_branches" in df.columns:
            st.bar_chart(df.set_index("label")["n_branches"], height=240)


# Render si hay resultados o tras ejecución
_render_metrics_table_and_charts(out_dir, mask_path, target_iso_um)
