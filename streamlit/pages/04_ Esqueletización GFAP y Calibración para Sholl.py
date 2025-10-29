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
                arr = np.asarray(s0)
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


# 1) Selección de imagen y opciones de máscara
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.lif")])
if not files:
    st.warning("No se encontraron archivos .tif/.lif en data/raw.")
    st.stop()
labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files))), format_func=lambda i: labels[i])
img_path = files[idx]
out_dir = get_output_dir_for_image(img_path)

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

# 2) Parámetros de esqueletización calibrada
arr_preview, axes_prev = load_image_any(img_path)
vol_prev = reorder_to_zcyx(arr_preview, axes_prev)  # Z,C,Y,X
n_channels = vol_prev.shape[1]

gcol1, gcol2, gcol3 = st.columns(3)
with gcol1:
    gfap_idx = st.number_input("Índice de canal GFAP", value=1 if n_channels > 1 else 0, min_value=0, max_value=max(0, n_channels-1), step=1)
cal = _read_global_calibration()
z_um = float(cal.get('z', 1.0))
y_um = float(cal.get('y', 1.0))
x_um = float(cal.get('x', 1.0))
def_iso = float(min(z_um, y_um, x_um))

# Cargar parámetros previos (prioridad: locales del preparado > globales) para autocompletar UI
try:
    params_local = json.loads((out_dir / "params.json").read_text()) if (out_dir / "params.json").exists() else {}
except Exception:
    params_local = {}
try:
    params_global = json.loads((root / "streamlit" / "experiment_params.json").read_text()) if (root / "streamlit" / "experiment_params.json").exists() else {}
except Exception:
    params_global = {}

def _get_param(name: str, fallback):
    if name in params_local:
        return params_local.get(name)
    if name in params_global:
        return params_global.get(name)
    return fallback

# Defaults desde JSON si existen
default_target_iso_um = float(_get_param("SKELETON_TARGET_ISO_UM", def_iso))
default_padding_um = float(_get_param("SKELETON_PADDING_UM", 2.0))
default_seed_dilate_um = float(_get_param("SKELETON_SEED_DILATE_UM", 2.0))
default_connectivity = int(_get_param("SKELETON_CONNECTIVITY", 26))
default_closing_um = float(_get_param("SKELETON_CLOSING_UM", 0.8))
default_territory_voronoi = bool(_get_param("SKELETON_TERRITORY_VORONOI", False))
default_territory_excl_um = float(_get_param("SKELETON_TERRITORY_EXCLUSION_UM", 1.0))

# Calcular default de radio máximo: preferir SKELETON_MAX_RADIUS_UM; si no, derivar de MAX_DILATION_ITERATIONS
mdi = params_local.get("MAX_DILATION_ITERATIONS") or params_global.get("MAX_DILATION_ITERATIONS")
if "SKELETON_MAX_RADIUS_UM" in params_local or "SKELETON_MAX_RADIUS_UM" in params_global:
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

save_skel_params = st.button("� Guardar parámetros del skeleton (sidebar)")
if save_skel_params:
    exp_params_path = root / "streamlit" / "experiment_params.json"
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
    })
    exp_params_path.parent.mkdir(parents=True, exist_ok=True)
    exp_params_path.write_text(json.dumps(exp, indent=2))
    st.success(f"Parámetros del skeleton guardados en {exp_params_path.relative_to(root)}")

conflict_resolve = st.checkbox(
    "Resolver solapamientos por cercanía al núcleo (µm)",
    value=True,
    help="Si dos esqueletos se superponen, se asigna el vóxel al astro más cercano al núcleo en distancia física.")

run_skel = st.button("�🕸️ Esqueletizar y Guardar")
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

    # Padding en voxeles originales (anisotrópicos) por eje
    # Importante: incluir el radio máximo solicitado, de lo contrario el recorte de la ROI
    # puede truncar procesos largos aunque max_radius_um sea alto.
    pad_z = _um_to_vox(padding_um + max_radius_um, z_um)
    pad_y = _um_to_vox(padding_um + max_radius_um, y_um)
    pad_x = _um_to_vox(padding_um + max_radius_um, x_um)

    labels = np.unique(mask)
    labels = labels[labels > 0]

    metrics = []
    for lab in labels:
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

        # Métricas simples: conteo de vóxeles de esqueleto y longitud aprox (en µm)
        skel_voxels = int(np.count_nonzero(skel_iso))
        approx_length_um = skel_voxels * float(target_iso_um)

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
            "skeleton_voxels_iso": skel_voxels,
            "approx_length_um": float(approx_length_um),
            "threshold": float(thr),
            "max_radius_um": float(max_radius_um),
        })

    # Guardar combinado
    out_skel_labels = out_dir / "05_skeleton_labels.tif"
    tifffile.imwrite(out_skel_labels, combined_labels.astype(np.uint16))

    # Guardar métricas
    import pandas as _pd
    df = _pd.DataFrame(metrics)
    if not df.empty:
        df.to_csv(skel_dir / "summary.csv", index=False)
    # Persistir parámetros de skeleton en params por preparado
    try:
        _mode_str = "manual" if str(thresh_mode).lower().startswith("manual") else "otsu"
        run_params = {
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
        }
        params_path = out_dir / "params.json"
        current = {}
        if params_path.exists():
            try:
                current = json.loads(params_path.read_text())
            except Exception:
                current = {}
        current.update(run_params)
        params_path.write_text(json.dumps(current, indent=2))
    except Exception:
        pass
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
        os.environ.update({})
        __import__("subprocess").Popen(cmd)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")
