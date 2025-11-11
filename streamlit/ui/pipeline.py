import json
from pathlib import Path
from typing import Tuple, Literal
from concurrent.futures import ProcessPoolExecutor
from functools import partial # Necesario para pasar argumentos a los workers

import numpy as np
import tifffile
import pandas as pd
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.measure import regionprops, label as cc_label
from skimage.morphology import binary_closing, binary_dilation, ball, skeletonize as _skeletonize_2d
from skimage.transform import resize
from scipy.ndimage import distance_transform_edt

# --- Dependencias Opcionales ---
try:
    from cellpose import models as _cp_models  # type: ignore
except Exception:
    _cp_models = None

try:
    from skimage.morphology import skeletonize_3d as _skeletonize_3d
    def skeletonize3d(arr: np.ndarray) -> np.ndarray:
        return _skeletonize_3d(arr)
except Exception:
    def skeletonize3d(arr: np.ndarray) -> np.ndarray:
        out = np.zeros_like(arr, dtype=bool)
        for z in range(arr.shape[0]):
            out[z] = _skeletonize_2d(arr[z].astype(bool))
        return out

try:
    from skan import Skeleton, summarize as sk_summarize, sholl_analysis as skan_sholl
    import skan
except Exception:
    sk_summarize = None
    Skeleton = None
    skan_sholl = None
    skan = None

# --- PASO 0: Utilidades de Carga y Calibración ---

def _read_global_calibration(calib_path: Path) -> dict:
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())
        except Exception:
            pass
    return {}

def load_image_any(path: Path) -> tuple[np.ndarray, str | None]:
    """Carga .tif/.tiff y devuelve (array, axes)."""
    suffix = path.suffix.lower()
    if suffix in (".tif", ".tiff"):
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            axes = getattr(series, 'axes', None)
            arr = series.asarray()
        return arr, axes
    raise ValueError(f"Extensión no soportada: {suffix}")

def reorder_to_zcyx(arr: np.ndarray, axes: str | None) -> np.ndarray:
    """Reordena un array con ejes reportados a (Z, C, Y, X)."""
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
    needed = ['Z','C','Y','X']
    if not all(a in ax_list for a in needed):
        raise ValueError(f"Ejes insuficientes: {ax_list}")
    src_order = [ax_list.index(a) for a in needed]
    return np.transpose(arr, axes=src_order)

def get_voxel_volume(cal: dict) -> float:
    z_um = float(cal.get('z', 1.0))
    y_um = float(cal.get('y', 1.0))
    x_um = float(cal.get('x', 1.0))
    voxel_vol_um3 = z_um * y_um * x_um
    return voxel_vol_um3 if voxel_vol_um3 > 0 else 1.0

# --- FUNCIONES DE PROYECCIÓN 2D ---

def project_3d_to_2d(volume_3d: np.ndarray, method: str = "max", axis: int = 0) -> np.ndarray:
    """
    Proyecta un volumen 3D a una imagen 2D.
    
    Parameters:
    -----------
    volume_3d : np.ndarray
        Volumen 3D de entrada
    method : str
        Método de proyección: "max", "mean", "sum", "any"
    axis : int
        Eje a lo largo del cual proyectar (default: 0 para Z)
    
    Returns:
    --------
    np.ndarray : Imagen 2D proyectada
    """
    if method == "max":
        return np.max(volume_3d, axis=axis)
    elif method == "mean":
        return np.mean(volume_3d, axis=axis)
    elif method == "sum":
        return np.sum(volume_3d, axis=axis)
    elif method == "any":
        return np.any(volume_3d, axis=axis)
    else:
        raise ValueError(f"Método de proyección no soportado: {method}")

def project_mask_labels_2d(mask_labels_3d: np.ndarray, method: str = "max") -> np.ndarray:
    """
    Proyecta máscaras con labels preservando las etiquetas.
    Para máscaras con labels, usa método inteligente para preservar información.
    """
    if method in ["max", "any"]:
        # Para labels, buscar el label más frecuente por pixel en Z
        projection = np.zeros(mask_labels_3d.shape[1:], dtype=mask_labels_3d.dtype)
        
        for y in range(mask_labels_3d.shape[1]):
            for x in range(mask_labels_3d.shape[2]):
                z_vals = mask_labels_3d[:, y, x]
                non_zero = z_vals[z_vals > 0]
                if len(non_zero) > 0:
                    # Tomar el label más frecuente o el primero si hay empate
                    unique, counts = np.unique(non_zero, return_counts=True)
                    projection[y, x] = unique[np.argmax(counts)]
        
        return projection
    else:
        # Para otros métodos, usar proyección simple
        return project_3d_to_2d(mask_labels_3d, method)

# --- PASO 1: Otsu ---

def run_otsu_and_save(dapi_vol: np.ndarray, out_dir: Path) -> Tuple[np.ndarray, float]:
    """Calcula y guarda la máscara de Otsu."""
    thr = threshold_otsu(dapi_vol)
    otsu_mask = (dapi_vol > thr)
    tifffile.imwrite(out_dir / "01_otsu_mask.tif", otsu_mask.astype(np.uint8))
    return otsu_mask.astype(np.uint8), float(thr)

# --- PASO 2: Cellpose ---

def run_cellpose_and_save(dapi_vol_clean: np.ndarray, out_dir: Path, nucleus_diameter: int, use_gpu: bool) -> np.ndarray:
    """Ejecuta Cellpose y guarda la máscara de núcleos."""
    if _cp_models is None:
        raise RuntimeError("Cellpose no está instalado. Por favor, instale 'cellpose'.")
    try:
        model = _cp_models.CellposeModel(gpu=use_gpu)
    except Exception:
        print("ADVERTENCIA: No se pudo inicializar GPU; se usará CPU.")
        model = _cp_models.CellposeModel(gpu=False)

    masks, _, _ = model.eval(
        dapi_vol_clean,
        diameter=int(nucleus_diameter),
        z_axis=0,
        do_3D=True,
    )
    tifffile.imwrite(out_dir / "02_cellpose_mask.tif", masks.astype(np.uint16))
    return masks.astype(np.uint16)

# --- PASO 3: Filtrado de Astrocitos (Paralelizado) ---

def _calculate_background_stats(channel: np.ndarray, otsu_mask: np.ndarray) -> Tuple[float, float]:
    """Calcula la media y std del fondo (definido como fuera de la máscara de Otsu)."""
    background_voxels = channel[otsu_mask == 0]
    if background_voxels.size < 100:
        background_voxels = channel
    
    bg_mean = float(np.mean(background_voxels))
    bg_std = float(np.std(background_voxels))
    return bg_mean, (bg_std if bg_std > 1e-6 else 1.0) 

def _process_nucleus_filter(nucleus_props, cellpose_masks, gfap_channel, se, gfap_decision_thr, voxel_vol_um3):
    """
    Función worker para procesar un solo núcleo (usada por ProcessPoolExecutor).
    
    Filtrado basado únicamente en la intensidad de GFAP en el shell perinuclear.
    
    ESFERICIDAD: Se calcula en 2D a partir de la proyección MIP del núcleo 3D.
    Usa la fórmula de circularidad: 4π * area / perimeter²
    """
    label = nucleus_props.label
    nucleus_mask = (cellpose_masks == label)
    
    dilated_mask = binary_dilation(nucleus_mask, footprint=se)
    shell_mask = dilated_mask & ~nucleus_mask
    
    shell_gfap_intensity = 0.0
    
    if np.any(shell_mask):
        shell_gfap_intensity = float(gfap_channel[shell_mask].mean())
    
    # Decisión basada SOLO en GFAP
    is_astrocyte = shell_gfap_intensity > gfap_decision_thr
    
    volume_vox = nucleus_props.area
    volume_um3 = float(volume_vox) * voxel_vol_um3
    
    # --- CÁLCULO DE ESFERICIDAD EN 2D (Circularidad) ---
    # Proyección MIP del núcleo
    nucleus_2d = np.max(nucleus_mask, axis=0).astype(np.uint8)
    
    sphericity_2d = np.nan
    if np.any(nucleus_2d):
        # Calcular propiedades en 2D
        props_2d = regionprops(nucleus_2d)
        if len(props_2d) > 0:
            area_2d = props_2d[0].area
            perimeter_2d = props_2d[0].perimeter
            
            if perimeter_2d > 0:
                # Circularidad 2D: 4π * area / perimeter²
                # Valor = 1 para círculo perfecto, < 1 para formas irregulares
                sphericity_2d = (4.0 * np.pi * area_2d) / (perimeter_2d ** 2)
                # Limitar al rango [0, 1] para evitar artefactos de discretización
                sphericity_2d = min(sphericity_2d, 1.0)

    metric_dict = {
        "label": label,
        "nucleus_volume_um3": volume_um3,
        "nucleus_sphericity": sphericity_2d,  # Circularidad 2D
        "shell_gfap_mean": shell_gfap_intensity,
        "is_astrocyte_candidate": is_astrocyte
    }
    
    return label if is_astrocyte else None, metric_dict


def run_filter_and_save(
    cellpose_masks: np.ndarray, 
    gfap_channel: np.ndarray, 
    otsu_mask: np.ndarray,
    cal: dict,
    out_dir: Path
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Filtra núcleos de Cellpose basado en umbrales relativos (std dev) 
    en un shell alrededor del núcleo. (Ejecución PARALELIZADA)
    
    Filtrado basado únicamente en la señal de GFAP en el shell perinuclear.
    """
    
    z_um, y_um, x_um = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
    voxel_vol_um3 = get_voxel_volume(cal)
    
    shell_radius_um = float(cal.get("SHELL_RADIUS_UM", 2.0))
    shell_r_vox_z = max(1, int(round(shell_radius_um / z_um)))
    shell_r_vox_y = max(1, int(round(shell_radius_um / y_um)))
    shell_r_vox_x = max(1, int(round(shell_radius_um / x_um)))

    se = np.zeros((shell_r_vox_z*2+1, shell_r_vox_y*2+1, shell_r_vox_x*2+1), dtype=np.uint8)
    cz, cy, cx = shell_r_vox_z, shell_r_vox_y, shell_r_vox_x
    z, y, x = np.ogrid[-cz:cz+1, -cy:cy+1, -cx:cx+1]
    mask_struct = (z/shell_r_vox_z)**2 + (y/shell_r_vox_y)**2 + (x/shell_r_vox_x)**2 <= 1
    se[mask_struct] = 1

    gfap_std_thr = float(cal.get("GFAP_STD_DEV_THRESHOLD", 3.0))

    gfap_bg_mean, gfap_bg_std = _calculate_background_stats(gfap_channel, otsu_mask)
    
    gfap_decision_thr = gfap_bg_mean + (gfap_std_thr * gfap_bg_std)

    print(f"--- Filtrado Relativo (Solo GFAP) ---")
    print(f"GFAP Fondo: {gfap_bg_mean:.2f} ± {gfap_bg_std:.2f}. Umbral Decisión: > {gfap_decision_thr:.2f}")

    nuclei_props = regionprops(cellpose_masks)
    astrocyte_labels = []
    metrics_data = []

    task_func = partial(
        _process_nucleus_filter,
        cellpose_masks=cellpose_masks,
        gfap_channel=gfap_channel,
        se=se,
        gfap_decision_thr=gfap_decision_thr,
        voxel_vol_um3=voxel_vol_um3
    )
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(task_func, nuclei_props))

    for label, metric_dict in results:
        if label is not None:
            astrocyte_labels.append(label)
        metrics_data.append(metric_dict)

    gfap_filtered_mask = np.where(np.isin(cellpose_masks, astrocyte_labels), cellpose_masks, 0)
    tifffile.imwrite(out_dir / "03_gfap_filtered_mask.tif", gfap_filtered_mask.astype(np.uint16))
    
    # Generar proyección 2D si está habilitada
    proj_2d_enable = bool(cal.get("PROJECTION_2D_ENABLE", True))
    if proj_2d_enable:
        proj_method = str(cal.get("PROJECTION_2D_METHOD", "max")).lower()
        gfap_filtered_2d = project_mask_labels_2d(gfap_filtered_mask, method=proj_method)
        tifffile.imwrite(out_dir / "03_gfap_filtered_mask_2d.tif", gfap_filtered_2d.astype(np.uint16))
    
    df_metrics = pd.DataFrame(metrics_data)
    df_metrics.to_csv(out_dir / "03_nucleus_metrics.csv", index=False)
    
    print(f"Máscara filtrada guardada. {len(astrocyte_labels)} candidatos retenidos.")
    return gfap_filtered_mask.astype(np.uint16), df_metrics


def run_size_filter_and_save(gfap_filtered_mask: np.ndarray, cal: dict, out_dir: Path) -> np.ndarray:
    """Aplica filtro por tamaño físico y guarda 04_final_astrocytes_mask.tif"""
    min_volume_um3 = float(cal.get("MIN_VOLUME_UM3", 75))
    voxel_vol_um3 = get_voxel_volume(cal)
    
    min_voxels = int(np.ceil(min_volume_um3 / voxel_vol_um3)) if min_volume_um3 > 0 else 0
    
    if min_voxels <= 1:
        print("Filtro de tamaño deshabilitado (min_voxels <= 1).")
        tifffile.imwrite(out_dir / "04_final_astrocytes_mask.tif", gfap_filtered_mask.astype(np.uint16))
        return gfap_filtered_mask

    props = regionprops(gfap_filtered_mask)
    kept_labels = [p.label for p in props if p.area >= min_voxels]
    
    final_mask = np.where(np.isin(gfap_filtered_mask, kept_labels), gfap_filtered_mask, 0)
    tifffile.imwrite(out_dir / "04_final_astrocytes_mask.tif", final_mask.astype(np.uint16))
    print(f"Máscara final guardada. Objetos finales: {len(kept_labels)}")
    return final_mask.astype(np.uint16)

# --- PASO 4: Esqueletización (Paralelizado) ---

def _um_to_vox(um_val: float, voxel_um: float) -> int:
    return int(round(max(0.0, um_val) / max(voxel_um, 1e-6)))

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
        vol, new_shape, order=order, preserve_range=preserve_range, 
        mode="edge", anti_aliasing=(order != 0)
    )

def _resize_to_original(vol_iso, orig_shape, order=0):
    return resize(
        vol_iso, orig_shape, order=order, preserve_range=True, 
        mode="edge", anti_aliasing=(order != 0)
    )

def _prune_short_spurs(skel: np.ndarray, min_len_um: float, vox_um: float) -> np.ndarray:
    """Elimina iterativamente ramas terminales más cortas que min_len_um."""
    if sk_summarize is None or Skeleton is None:
        print("Advertencia: `skan` no está instalado. No se puede realizar pruning.")
        return skel
    if min_len_um <= 0:
        return skel
    
    try:
        sk = Skeleton(skel.astype(np.uint8), spacing=vox_um)
        summary = sk_summarize(sk, separator='_') 
        
        # --- ¡CORRECCIÓN 1: Comprobar si 'branch_distance' y 'branch_type' existen! ---
        if not all(col in summary.columns for col in ['branch_distance', 'branch_type']):
            print("Advertencia: Skan no devolvió 'branch_distance' o 'branch_type' (esqueleto trivial). No se puede podar.")
            return skel
        # -------------------------------------------------------------------------
        
        terminal_branches = summary[summary['branch_type'] == 1]
        short_terminal_branches = terminal_branches[terminal_branches['branch_distance'] < min_len_um]
        
        if short_terminal_branches.empty:
            return skel 

        pruned_skel = skel.copy()
        
        # --- ¡CORRECCIÓN 2: Usar 'index' en lugar de ['branch-id']! ---
        for branch_id in short_terminal_branches.index:
        # ---------------------------------------------------------
            branch_coords = sk.path_coordinates(branch_id)
            if branch_coords.shape[0] > 0:
                indices = np.round(branch_coords / vox_um).astype(int)
                indices = np.clip(indices, 0, np.array(skel.shape) - 1)
                pruned_skel[indices[:, 0], indices[:, 1], indices[:, 2]] = False
        
        return pruned_skel

    except Exception as e:
        print(f"Advertencia: Falla en pruning con Skan: {e}. Se devuelve esqueleto sin podar.")
        return skel

def _process_skeleton_label(
    lab, 
    mask: np.ndarray, 
    gfap: np.ndarray, 
    cal: dict, 
    label_centroid_um: dict
) -> tuple:
    """Función worker para procesar un solo esqueleto (usada por ProcessPoolExecutor)."""
    
    # Extraer parámetros de cal
    z_um, y_um, x_um = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
    target_iso_um = float(cal.get("SKELETON_TARGET_ISO_UM", min(z_um, y_um, x_um)))
    padding_um = float(cal.get("SKELETON_PADDING_UM", 2.0))
    seed_dilate_um = float(cal.get("SKELETON_SEED_DILATE_UM", 2.0))
    connectivity = 3 if int(cal.get("SKELETON_CONNECTIVITY", 26)) == 26 else 1
    closing_um = float(cal.get("SKELETON_CLOSING_UM", 0.8))
    max_radius_um = float(cal.get("SKELETON_MAX_RADIUS_UM", 50.0))
    thresh_mode = str(cal.get("SKELETON_THRESHOLD_MODE", "Otsu (ROI)")).lower()
    manual_thr = float(cal.get("SKELETON_MANUAL_THRESHOLD", 50))
    territory_voronoi = bool(cal.get("SKELETON_TERRITORY_VORONOI", False))
    territory_excl_um = float(cal.get("SKELETON_TERRITORY_EXCLUSION_UM", 1.0))
    domain_source = str(cal.get("SKELETON_DOMAIN_VOLUME_SOURCE", "gfap")).lower()
    prune_enable = bool(cal.get("SKELETON_PRUNE_ENABLE", False))
    prune_min_len_um = float(cal.get("SKELETON_PRUNE_MIN_LEN_UM", 2.0))
    tube_radius_um = float(cal.get("SKELETON_TUBE_RADIUS_UM", 1.5))
    conflict_resolve = bool(cal.get("SKELETON_CONFLICT_RESOLVE", True))

    # Parámetros en Vóxeles
    zf, yf, xf = _zoom_factors(z_um, y_um, x_um, target_iso_um)
    closing_r_vox = _um_to_vox(closing_um, target_iso_um)
    seed_r_vox = _um_to_vox(seed_dilate_um, target_iso_um)
    max_r_vox = _um_to_vox(max_radius_um, target_iso_um)
    excl_r_vox = _um_to_vox(territory_excl_um, target_iso_um)
    tube_r_vox = _um_to_vox(tube_radius_um, target_iso_um)
    pad_z = _um_to_vox(padding_um + max_radius_um, z_um)
    pad_y = _um_to_vox(padding_um + max_radius_um, y_um)
    pad_x = _um_to_vox(padding_um + max_radius_um, x_um)
    
    # --- 3. Procesar Label ---
    reg_mask = (mask == lab)
    if not np.any(reg_mask): return (lab, None, None, None, None) # Retornar Nones
    props = regionprops(reg_mask.astype(np.uint8))
    if not props: return (lab, None, None, None, None)
    
    z0, y0, x0, z1, y1, x1 = (*props[0].bbox[0:3], *props[0].bbox[3:6])
    z0=max(0,z0-pad_z); y0=max(0,y0-pad_y); x0=max(0,x0-pad_x)
    z1=min(mask.shape[0],z1+pad_z); y1=min(mask.shape[1],y1+pad_y); x1=min(mask.shape[2],x1+pad_x)

    roi_mask = reg_mask[z0:z1, y0:y1, x0:x1]
    roi_gfap = gfap[z0:z1, y0:y1, x0:x1]
    roi_all_labels = mask[z0:z1, y0:y1, x0:x1]

    roi_mask_iso = _resize_to_isotropic(roi_mask.astype(np.uint8), (zf, yf, xf), order=0) > 0.5
    roi_gfap_iso = _resize_to_isotropic(roi_gfap, (zf, yf, xf), order=1)
    roi_all_labels_iso = _resize_to_isotropic(roi_all_labels.astype(np.int32), (zf, yf, xf), order=0).astype(np.int32)

    if thresh_mode.startswith("otsu"):
        thr = threshold_otsu(roi_gfap_iso[roi_mask_iso]) if np.any(roi_mask_iso) else threshold_otsu(roi_gfap_iso)
    else:
        thr = manual_thr
    
    bin_gfap_iso = (roi_gfap_iso >= thr)
    if closing_r_vox > 0: bin_gfap_iso = binary_closing(bin_gfap_iso, ball(int(closing_r_vox)))

    seed_iso = roi_mask_iso.copy()
    if seed_r_vox > 0: seed_iso = binary_dilation(seed_iso, ball(int(seed_r_vox)))
    
    bin_iso = np.zeros_like(bin_gfap_iso, dtype=bool)
    if np.any(bin_gfap_iso):
        cc = cc_label(bin_gfap_iso.astype(np.uint8), connectivity=connectivity)
        keep_ids = np.unique(cc[seed_iso > 0])
        keep_ids = keep_ids[keep_ids > 0]
        if keep_ids.size > 0:
            bin_iso = np.isin(cc, keep_ids)
    
    terr_lab = None
    if territory_voronoi:
        labs_here = np.unique(roi_all_labels_iso)
        labs_here = labs_here[labs_here > 0]
        if labs_here.size > 0:
            props_here = regionprops(roi_all_labels_iso)
            centroids = {int(p.label): p.centroid for p in props_here}
            Z, Y, X = roi_all_labels_iso.shape
            zz, yy, xx = np.indices((Z, Y, X)).astype(np.float32)
            best_d2 = np.full((Z, Y, X), np.inf, dtype=np.float32)
            best_lab = np.zeros((Z, Y, X), dtype=np.int32)
            for l in labs_here:
                if int(l) not in centroids: continue
                cz, cy, cx = centroids[int(l)]
                d2 = (zz - cz)**2 + (yy - cy)**2 + (xx - cx)**2
                better = d2 < best_d2
                best_lab[better] = int(l)
                best_d2[better] = d2[better]
            terr_lab = (best_lab == int(lab))
            if excl_r_vox > 0:
                interior = distance_transform_edt(terr_lab) >= float(excl_r_vox)
                terr_lab = terr_lab & interior
            bin_iso &= terr_lab

    if max_r_vox > 0 and np.any(bin_iso):
        dist = distance_transform_edt(~seed_iso)
        bin_iso &= (dist <= float(max_r_vox))

    skel_iso = skeletonize3d(bin_iso.astype(np.uint8)) > 0
    if prune_enable and prune_min_len_um > 0:
        skel_iso = _prune_short_spurs(skel_iso, prune_min_len_um, target_iso_um)

    # --- 4. CÁLCULO DE MÉTRICAS (UNIFICADO) ---
    
    skel_voxels = int(np.count_nonzero(skel_iso))
    iso_vox_vol_um3 = float(target_iso_um)**3
    
    total_length_um, n_branches, n_junctions, n_endpoints, mean_branch_len_um, mean_tortuosity = (np.nan, 0, 0, 0, np.nan, np.nan)
    
    if sk_summarize is not None and Skeleton is not None and skel_voxels > 0:
        try:
            skel_obj = Skeleton(skel_iso.astype(np.uint8), spacing=target_iso_um)
            df_skan = sk_summarize(skel_obj, separator='_')
            
            # --- ¡CORRECCIÓN 'branch-distance' y 'Mean of empty slice'! ---
            if not df_skan.empty and 'branch_distance' in df_skan.columns:
                df_skan_filt = df_skan[df_skan['branch_distance'] > 1e-6]
                # prefer underscore named column if present
                if 'branch_distance' in df_skan_filt.columns:
                    total_length_um = float(df_skan_filt['branch_distance'].sum())
                else:
                    total_length_um = float(df_skan_filt.iloc[:, 0].sum()) if not df_skan_filt.empty else 0.0
                n_branches = int(df_skan_filt.shape[0])
                
                branch_len_slice = df_skan_filt['branch_distance'] if 'branch_distance' in df_skan_filt.columns else df_skan_filt.iloc[:, 0]
                if n_branches > 0 and branch_len_slice.notnull().any():
                    mean_branch_len_um = float(branch_len_slice.mean())
                else:
                    mean_branch_len_um = 0.0

                if 'tortuosity' in df_skan_filt.columns:
                    tort_slice = df_skan_filt['tortuosity']
                    if n_branches > 0 and tort_slice.notnull().any():
                        mean_tortuosity = float(tort_slice.mean())
                    else:
                        mean_tortuosity = 1.0
                else:
                    mean_tortuosity = 1.0

                # branch type column may be named with underscore
                bt_col = 'branch_type' if 'branch_type' in df_skan.columns else ('branch-type' if 'branch-type' in df_skan.columns else None)
                if bt_col is not None:
                    n_junctions = int(((df_skan[bt_col] == 2).sum() + (df_skan[bt_col] == 3).sum()))
                    n_endpoints = int((df_skan[bt_col] == 1).sum())
                else:
                    n_junctions = 0
                    n_endpoints = 0
            
            # Si el df de skan no tenía 'branch-distance', usar el fallback
            elif skel_voxels > 0:
                total_length_um = skel_voxels * float(target_iso_um) 

        except Exception as e:
            # Detailed logging for Skan failures: include traceback and df_skan diagnostics
            try:
                import traceback as _tb
                print(f"Error en Skan (Label {lab}): {e} ({type(e).__name__})")
                print("Traceback (most recent call last):")
                _tb.print_exc()
            except Exception:
                # fallback to simple message if traceback printing fails
                print(f"Error en Skan (Label {lab}): {e}")

            # If df_skan was created, print available columns and a small preview to help debug
            try:
                if 'df_skan' in locals() and df_skan is not None:
                    try:
                        print(f"[skan diagnostics] df_skan.columns: {list(df_skan.columns)}")
                        # print small preview
                        try:
                            print("[skan diagnostics] df_skan.head():\n" + df_skan.head().to_string())
                        except Exception:
                            pass
                    except Exception:
                        pass
            except Exception:
                pass

            # Additional runtime info
            try:
                print(f"[skan diagnostics] skel_voxels={skel_voxels}")
                if 'skel_iso' in locals():
                    print(f"[skan diagnostics] skel_iso.shape={getattr(skel_iso, 'shape', 'N/A')}")
            except Exception:
                pass

            # Fallback: estimate length from voxel count to continue processing
            if skel_voxels > 0:
                total_length_um = skel_voxels * float(target_iso_um) # Fallback

    elif skel_voxels > 0: # Si Skan no está instalado
        total_length_um = skel_voxels * float(target_iso_um)

    gfap_connected_vox = int(np.count_nonzero(bin_iso))
    gfap_connected_um3 = float(gfap_connected_vox) * iso_vox_vol_um3
    voronoi_um3 = np.nan
    if territory_voronoi and terr_lab is not None:
        voronoi_um3 = float(np.count_nonzero(terr_lab)) * iso_vox_vol_um3
    
    domain_um3 = voronoi_um3 if domain_source.startswith("voro") else gfap_connected_um3

    tube_sum_intensity, tube_mean_intensity, tube_volume_um3 = (np.nan, np.nan, np.nan)
    if tube_r_vox > 0 and np.any(skel_iso):
        dist_to_skel = distance_transform_edt(~skel_iso)
        tube_mask = dist_to_skel <= float(tube_r_vox)
        tube_sel = tube_mask & (roi_gfap_iso >= float(thr))
        tube_voxels = int(np.count_nonzero(tube_sel))
        if tube_voxels > 0:
            tube_sum_intensity = float(np.sum(roi_gfap_iso[tube_sel]))
            tube_mean_intensity = float(np.mean(roi_gfap_iso[tube_sel]))
            tube_volume_um3 = float(tube_voxels) * iso_vox_vol_um3
    
    local_radius_um_mean, local_radius_um_median = (np.nan, np.nan)
    if np.any(bin_iso) and np.any(skel_iso):
        dt_gfap = distance_transform_edt(bin_iso)
        loc_r_um_vals = dt_gfap[skel_iso].astype(np.float32) * float(target_iso_um)
        if loc_r_um_vals.size > 0:
            local_radius_um_mean = float(np.mean(loc_r_um_vals))
            local_radius_um_median = float(np.median(loc_r_um_vals))

    # --- 5. Guardar Métricas y Reescalar ---
    skel_roi = _resize_to_original(skel_iso.astype(np.uint8), roi_mask.shape, order=0) > 0.5
    
    dist_um = None
    if conflict_resolve and int(lab) in label_centroid_um:
        cz_um, cy_um, cx_um = label_centroid_um[int(lab)]
        zz = (np.arange(z0, z1, dtype=np.float32) * z_um)[:, None, None]
        yy = (np.arange(y0, y1, dtype=np.float32) * y_um)[None, :, None]
        xx = (np.arange(x0, x1, dtype=np.float32) * x_um)[None, None, :]
        dist_um = np.sqrt((zz - cz_um)**2 + (yy - cy_um)**2 + (xx - cx_um)**2)

    metric_dict = {
        "label": int(lab),
        "total_length_um": total_length_um,
        "n_branches": n_branches,
        "n_junctions": n_junctions,
        "n_endpoints": n_endpoints,
        "mean_branch_len_um": mean_branch_len_um,
        "mean_tortuosity": mean_tortuosity,
        "domain_volume_um3": domain_um3,
        "gfap_connected_volume_um3": gfap_connected_um3,
        "voronoi_territory_volume_um3": voronoi_um3,
        "tube_volume_um3": tube_volume_um3,
        "tube_mean_intensity": tube_mean_intensity,
        "tube_sum_intensity": tube_sum_intensity,
        "local_radius_um_mean": local_radius_um_mean,
        "local_radius_um_median": local_radius_um_median
    }
    
    return (lab, (z0,y0,x0,z1,y1,x1), skel_roi, dist_um, metric_dict)


def run_skeleton_2d_native_skan(
    mask: np.ndarray,
    gfap: np.ndarray, 
    cal: dict,
    out_dir: Path,
    territory_expansion_um: float = None
) -> pd.DataFrame:
    """
    Esqueletización 2D nativa usando SKAN después de proyección.
    Flujo: 3D -> 2D proyección -> SKAN nativo -> métricas oficiales
    """
    print("=== ESQUELETIZACIÓN 2D GFAP TERRITORIAL ===")
    
    # Obtener parámetro de expansión territorial
    if territory_expansion_um is None:
        territory_expansion_um = float(cal.get("TERRITORY_EXPANSION_UM", 10.0))
    
    # 1. Proyectar las máscaras 3D a 2D
    proj_method = str(cal.get("PROJECTION_2D_METHOD", "max")).lower()
    print(f"Proyectando máscaras 3D a 2D usando método: {proj_method}")
    print(f"Expansión territorial: {territory_expansion_um:.1f} µm")
    
    # Proyectar máscara de núcleos (astrocitos)
    mask_2d = project_mask_labels_2d(mask, method=proj_method)
    out_mask_2d = out_dir / "04_final_astrocytes_mask_2d.tif" 
    tifffile.imwrite(out_mask_2d, mask_2d.astype(np.uint16))
    
    # Proyectar canal GFAP
    gfap_2d = project_3d_to_2d(gfap, method=proj_method)
    out_gfap_2d = out_dir / "gfap_projection_2d.tif"
    tifffile.imwrite(out_gfap_2d, gfap_2d.astype(np.float32))
    
    # 2. Parámetros para análisis 2D
    y_um, x_um = float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
    spacing_2d = (y_um, x_um)
    
    # Calcular radio de expansión en píxeles
    expansion_radius_px_y = int(territory_expansion_um / y_um)
    expansion_radius_px_x = int(territory_expansion_um / x_um)
    
    skel_dir = out_dir / "skeletons"
    skel_dir.mkdir(parents=True, exist_ok=True)
    
    # Obtener labels únicos (núcleos)
    labels = np.unique(mask_2d)
    labels = labels[labels > 0]
    
    all_metrics = []
    combined_skeleton_2d = np.zeros_like(mask_2d, dtype=np.uint16)
    
    print(f"Procesando {len(labels)} astrocitos con territorio expandido...")
    
    from skimage.morphology import skeletonize, binary_dilation, disk
    from skimage.segmentation import find_boundaries
    import skan
    
    for label in labels:
        try:
            # 1. Obtener máscara nuclear de este astrocito
            nuclear_mask = (mask_2d == label).astype(bool)
            
            if not np.any(nuclear_mask):
                continue
            
            # 2. Expandir territorio alrededor del núcleo
            # Usar elemento estructurante elíptico para mejor expansión isotrópica
            expansion_element = disk(max(expansion_radius_px_y, expansion_radius_px_x))
            territorial_mask = binary_dilation(nuclear_mask, expansion_element)
            
            # 3. Extraer señal GFAP en territorio expandido
            gfap_territorial = gfap_2d * territorial_mask
            
            # 4. Binarizar señal GFAP usando umbralización automática
            from skimage.filters import threshold_otsu
            try:
                gfap_region = gfap_territorial[territorial_mask]
                if len(gfap_region) > 0 and np.std(gfap_region) > 0:
                    threshold_val = threshold_otsu(gfap_region)
                    gfap_binary = gfap_territorial > threshold_val
                else:
                    # Si no hay variación, usar percentil
                    threshold_val = np.percentile(gfap_region, 75) if len(gfap_region) > 0 else 0
                    gfap_binary = gfap_territorial > threshold_val
            except:
                # Fallback: usar percentil 75
                gfap_region = gfap_territorial[territorial_mask]
                threshold_val = np.percentile(gfap_region, 75) if len(gfap_region) > 0 else 0
                gfap_binary = gfap_territorial > threshold_val
            
            # 5. Restringir a territorio (quitar ruido externo)
            gfap_binary = gfap_binary & territorial_mask
            
            # 6. Esqueletizar señal GFAP binarizada
            if not np.any(gfap_binary):
                print(f"  Astrocito {label}: sin señal GFAP significativa")
                continue
                
            skeleton_binary = skeletonize(gfap_binary)
            
            if not np.any(skeleton_binary):
                print(f"  Astrocito {label}: sin esqueleto GFAP")
                continue
            
            # Convertir a SKAN Skeleton
            import skan
            skeleton_obj = skan.Skeleton(skeleton_binary, spacing=spacing_2d)
            
            # Obtener métricas nativas de SKAN con separator para evitar warnings
            summary_stats = skan.summarize(skeleton_obj, separator='_')
            
            # Agregar al esqueleto combinado
            combined_skeleton_2d[skeleton_binary] = label
            
            # 7. Extraer métricas principales calculándolas correctamente
            # SKAN summarize devuelve datos por rama, no métricas agregadas
            if len(summary_stats) > 0:
                # Longitud total: sumar todas las branch_distance
                total_length_um = float(summary_stats['branch_distance'].sum())
                
                # Número de ramas: contar filas del summary
                n_branches = int(len(summary_stats))
                
                # Número de junctions: usar grados del skeleton object
                # Junctions son nodos con grado > 2
                degrees = skeleton_obj.degrees
                n_junctions = int(np.sum(degrees > 2))
            else:
                total_length_um = 0.0
                n_branches = 0
                n_junctions = 0
            
            # 8. Métricas territoriales y de área
            nuclear_area_px = np.sum(nuclear_mask)
            nuclear_area_um2 = nuclear_area_px * y_um * x_um
            
            territorial_area_px = np.sum(territorial_mask) 
            territorial_area_um2 = territorial_area_px * y_um * x_um
            
            gfap_positive_area_px = np.sum(gfap_binary)
            gfap_positive_area_um2 = gfap_positive_area_px * y_um * x_um
            
            # Calcular densidades
            skeleton_density_territorial = total_length_um / territorial_area_um2 if territorial_area_um2 > 0 else 0.0
            skeleton_density_gfap = total_length_um / gfap_positive_area_um2 if gfap_positive_area_um2 > 0 else 0.0
            
            # Ratio de ocupación GFAP en territorio
            gfap_occupancy_ratio = gfap_positive_area_um2 / territorial_area_um2 if territorial_area_um2 > 0 else 0.0
            
            # Para métricas de tortuosidad, usar datos del skeleton object directamente
            coords = skeleton_obj.coordinates
            if len(coords) > 1:
                # Calcular tortuosidad simple como longitud total / distancia euclidiana máxima
                coords_scaled = coords * np.array([y_um, x_um])
                distances = np.linalg.norm(np.diff(coords_scaled, axis=0), axis=1)
                if len(distances) > 0:
                    path_length = np.sum(distances)
                    euclidean_dist = np.linalg.norm(coords_scaled[-1] - coords_scaled[0])
                    mean_tortuosity = path_length / max(euclidean_dist, 1e-6)
                    max_tortuosity = mean_tortuosity  # Simplificado
                else:
                    mean_tortuosity = max_tortuosity = 1.0
            else:
                mean_tortuosity = max_tortuosity = 1.0
            
            # Longitud promedio de rama
            avg_branch_length = total_length_um / max(n_branches, 1)
            
            # 9. Construir métricas completas
            metrics_row = {
                "label": int(label),
                "total_length_um": float(total_length_um),
                "nuclear_area_um2": float(nuclear_area_um2),
                "territorial_area_um2": float(territorial_area_um2),
                "gfap_positive_area_um2": float(gfap_positive_area_um2),
                "territory_area_um2": float(territorial_area_um2),  # Para compatibilidad
                "skeleton_density_territorial": float(skeleton_density_territorial),
                "skeleton_density_gfap": float(skeleton_density_gfap), 
                "skeleton_density_um_per_um2": float(skeleton_density_territorial),  # Para compatibilidad
                "gfap_occupancy_ratio": float(gfap_occupancy_ratio),
                "territory_expansion_um": float(territory_expansion_um),
                "n_branches": int(n_branches),
                "n_junctions": int(n_junctions),
                "avg_branch_length_um": float(avg_branch_length),
                "mean_tortuosity": float(mean_tortuosity),
                "max_tortuosity": float(max_tortuosity),
                "skeleton_method": "skan_gfap_territorial_2d",
                "projection_method": proj_method
            }
            
            all_metrics.append(metrics_row)
            
            # Guardar datos SKAN individuales
            summary_stats.to_csv(skel_dir / f"skan_summary_{label}.csv", index=False)
            
            # 10. Guardar imágenes territoriales para debugging/visualización
            tifffile.imwrite(skel_dir / f"nuclear_mask_{label}.tif", nuclear_mask.astype(np.uint8))
            tifffile.imwrite(skel_dir / f"territorial_mask_{label}.tif", territorial_mask.astype(np.uint8))
            tifffile.imwrite(skel_dir / f"gfap_binary_{label}.tif", gfap_binary.astype(np.uint8))
            tifffile.imwrite(skel_dir / f"skeleton_{label}.tif", skeleton_binary.astype(np.uint8))
            
            print(f"  Astrocito {label}: {total_length_um:.1f}µm GFAP, {n_branches} ramas, territorio {territorial_area_um2:.1f}µm², GFAP+ {gfap_positive_area_um2:.1f}µm²")
            
        except Exception as e:
            print(f"  Error procesando astrocito {label}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Guardar esqueleto combinado 2D
    out_skel_2d = out_dir / "05_skeleton_labels_2d.tif"
    tifffile.imwrite(out_skel_2d, combined_skeleton_2d)
    
    # Guardar métricas
    df_metrics = pd.DataFrame(all_metrics)
    if not df_metrics.empty:
        df_metrics.to_csv(skel_dir / "skan_native_2d_summary.csv", index=False)
        print(f"Esqueletización 2D completada. {len(all_metrics)} astrocitos procesados.")
    else:
        print("No se generaron métricas de esqueletización.")
    
    return df_metrics


def run_unified_2d_skeleton_and_sholl(
    img_path: Path,
    mask_path: Path,
    out_dir: Path,
    cal: dict
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Pipeline 2D unificado: esqueletización + Sholl en un solo paso.
    Reemplaza run_advanced_skeleton_and_save() + run_sholl_and_save().
    
    Returns:
        (skeleton_metrics_df, sholl_df)
    """
    from .pipeline_2d_unified import run_unified_2d_skeleton_and_sholl as _run_unified
    
    # Cargar datos 3D
    gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1))
    arr, axes = load_image_any(img_path)
    vol = reorder_to_zcyx(arr, axes)
    gfap_3d = vol[:, gfap_idx, :, :].astype(np.float32)
    mask_3d = tifffile.imread(mask_path)
    
    if mask_3d.shape != gfap_3d.shape:
        raise RuntimeError(f"Dimensiones incompatibles: mask {mask_3d.shape} vs gfap {gfap_3d.shape}")
    
    # Ejecutar pipeline unificado
    return _run_unified(mask_3d, gfap_3d, cal, out_dir)


def run_advanced_skeleton_and_save(
    img_path: Path,
    mask_path: Path,
    out_dir: Path,
    cal: dict,
    conflict_resolve: bool = True
) -> pd.DataFrame:
    """
    DEPRECATED: Usar run_unified_2d_skeleton_and_sholl() en su lugar.
    
    Esta función se mantiene temporalmente por compatibilidad,
    pero ahora redirige al pipeline 2D unificado.
    """
    print("⚠️  run_advanced_skeleton_and_save() está deprecated")
    print("   Redirigiendo a pipeline 2D unificado...")
    
    df_skeleton, _ = run_unified_2d_skeleton_and_sholl(img_path, mask_path, out_dir, cal)
    return df_skeleton
    
    z_um, y_um, x_um = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
    gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1))
    
    arr, axes = load_image_any(img_path)
    vol = reorder_to_zcyx(arr, axes)
    gfap = vol[:, gfap_idx, :, :].astype(np.float32)
    mask = tifffile.imread(mask_path)
    if mask.shape != gfap.shape:
        raise RuntimeError("Dimensiones incompatibles entre GFAP y máscara")

    skel_dir = out_dir / "skeletons"
    skel_dir.mkdir(parents=True, exist_ok=True)
    combined_labels = np.zeros_like(mask, dtype=np.uint16)
    combined_dist_um = np.full(mask.shape, np.inf, dtype=np.float32) if conflict_resolve else None

    label_centroid_um = {}
    if conflict_resolve:
        for p in regionprops(mask.astype(np.int32)):
            lab_id = int(p.label)
            cz, cy, cx = p.centroid
            label_centroid_um[lab_id] = (float(cz) * z_um, float(cy) * y_um, float(cx) * x_um)

    labels = np.unique(mask)
    labels = labels[labels > 0]
    metrics = []

    task_func = partial(
        _process_skeleton_label,
        mask=mask,
        gfap=gfap,
        cal=cal,
        label_centroid_um=label_centroid_um
    )

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(task_func, labels))
    
    for (lab, bbox, skel_roi, dist_um, metric_dict) in results:
        if bbox is None or metric_dict is None:
            print(f"Label {lab} falló en el procesamiento.")
            continue
            
        metrics.append(metric_dict)
        (z0,y0,x0,z1,y1,x1) = bbox
        
        target_slice = combined_labels[z0:z1, y0:y1, x0:x1]
        
        if conflict_resolve and dist_um is not None:
            dist_slice = combined_dist_um[z0:z1, y0:y1, x0:x1]
            better = skel_roi & ((target_slice == 0) | (dist_um < dist_slice))
            target_slice[better] = lab
            dist_slice[better] = dist_um[better]
        else:
            target_slice[skel_roi] = lab 

    out_skel_labels = out_dir / "05_skeleton_labels.tif"
    tifffile.imwrite(out_skel_labels, combined_labels.astype(np.uint16))
    
    # Generar proyecciones 2D si está habilitado
    proj_2d_enable = bool(cal.get("PROJECTION_2D_ENABLE", True))
    if proj_2d_enable:
        proj_method = str(cal.get("PROJECTION_2D_METHOD", "max")).lower()
        
        # Proyectar esqueleto 3D a 2D
        skel_2d = project_mask_labels_2d(combined_labels, method=proj_method)
        out_skel_2d = out_dir / "05_skeleton_labels_2d.tif"
        tifffile.imwrite(out_skel_2d, skel_2d.astype(np.uint16))
        
        # También proyectar las máscaras originales
        mask_2d = project_mask_labels_2d(mask, method=proj_method)
        out_mask_2d = out_dir / "04_final_astrocytes_mask_2d.tif"
        tifffile.imwrite(out_mask_2d, mask_2d.astype(np.uint16))
        
        # Proyectar GFAP también para referencia
        gfap_2d = project_3d_to_2d(gfap, method=proj_method)
        out_gfap_2d = out_dir / "gfap_projection_2d.tif"
        tifffile.imwrite(out_gfap_2d, gfap_2d.astype(np.float32))

    df = pd.DataFrame(metrics)
    df.to_csv(skel_dir / "summary.csv", index=False)
    print(f"Esqueletización y análisis Skan completados. Guardado en {skel_dir.name}/summary.csv")
    return df

# --- PASO 5: Sholl (Paralelizado) ---

@pd.api.extensions.register_dataframe_accessor("sholl")
class ShollAnalysis:
    """Accesor de Pandas para análisis de Sholl."""
    def __init__(self, pandas_obj):
        if not all(col in pandas_obj.columns for col in ['label', 'radius_um', 'intersections']):
            raise AttributeError("El DataFrame debe tener columnas 'label', 'radius_um', y 'intersections'.")
        self._df = pandas_obj

    def get_scalar_metrics(self) -> pd.DataFrame:
        """Calcula métricas escalares (AUC, radio crítico, pico) por 'label'."""
        if self._df.empty:
            return pd.DataFrame(columns=['label', 'critical_radius_um', 'peak_intersections', 'auc'])
        
        from scipy.integrate import trapezoid
        
        results = []
        for label, group_df in self._df.groupby('label'):
            radii = group_df['radius_um'].to_numpy()
            intersections = group_df['intersections'].to_numpy()
            
            if radii.size < 2:
                results.append({'label': label, 'critical_radius_um': np.nan, 'peak_intersections': np.nan, 'auc': np.nan})
                continue

            peak_idx = np.argmax(intersections)
            critical_radius_um = radii[peak_idx]
            peak_intersections = intersections[peak_idx]
            auc = trapezoid(intersections, radii)
            
            results.append({
                'label': label,
                'critical_radius_um': float(critical_radius_um),
                'peak_intersections': float(peak_intersections),
                'auc': float(auc)
            })
        
        return pd.DataFrame(results)

def _process_sholl_label(lab, centroids, skel_points, spacing, radii, save_rings_json):
    """Función worker para procesar Sholl en un solo label."""
    if lab not in skel_points or lab not in centroids:
        return ([], None)
        
    centroid_vox = centroids[lab]
    points_vox = skel_points[lab]
    
    centroid_um = np.array(centroid_vox) * spacing
    points_um = points_vox.astype(np.float32) * spacing
    
    distances = np.sqrt(np.sum((points_um - centroid_um)**2, axis=1))
    counts, _ = np.histogram(distances, bins=radii)
    
    rings_data = None
    if save_rings_json:
        rings_data = {str(lab): {
            "centroid_um": centroid_um.tolist(),
            "radii_um": radii.tolist(),
        }}

    results_list = []
    for i in range(len(radii) - 1):
        results_list.append({
            "label": int(lab),
            "radius_um": float(radii[i]), # Radio interior
            "intersections": int(counts[i]),
        })
    
    return (results_list, rings_data)


def run_sholl_2d_native_skan_simple(
    skeleton_labels: np.ndarray,
    cal: dict,
    out_dir: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Análisis de Sholl usando SKAN nativo sobre esqueletos 2D."""
    import skan
    
    y_um = float(cal.get('y', 0.3))
    x_um = float(cal.get('x', 0.3))
    spacing = (y_um, x_um)
    
    min_radius_um = float(cal.get('SHOLL_MIN_RADIUS_UM', 5.0))
    max_radius_um = float(cal.get('SHOLL_MAX_RADIUS_UM', 100.0))
    step_um = float(cal.get('SHOLL_STEP_UM', 2.0))
    
    radii_um = np.arange(min_radius_um, max_radius_um + step_um, step_um)
    
    labels = np.unique(skeleton_labels)
    labels = labels[labels > 0]
    
    all_results = []
    summary_results = []
    rings_data = {}
    
    for label in labels:
        try:
            # Crear máscara binaria del esqueleto
            skel_binary = (skeleton_labels == label).astype(bool)
            
            if not np.any(skel_binary):
                continue
            
            # Crear objeto SKAN
            skeleton_obj = skan.Skeleton(skel_binary, spacing=spacing)
            
            # Obtener centroide del esqueleto
            coords = skeleton_obj.coordinates
            if len(coords) == 0:
                continue
                
            # Usar centroide de masa de las coordenadas del esqueleto
            centroid_voxels = np.mean(coords, axis=0)  # (y, x) en voxeles
            centroid_um = centroid_voxels * np.array([y_um, x_um])  # (y, x) en µm
            
            # Realizar análisis de Sholl nativo con SKAN
            # Nota: centroid debe ser en coordenadas físicas (µm)
            sholl_results = skan.sholl_analysis(skeleton_obj, centroid_um, radii_um)
            
            # Extraer datos para cada radio
            for i, radius in enumerate(radii_um):
                intersections = sholl_results[i] if i < len(sholl_results) else 0
                
                all_results.append({
                    'label': int(label),
                    'radius_um': float(radius),
                    'intersections': int(intersections),
                    'method': 'skan_native'
                })
            
            # Calcular métricas resumen
            total_intersections = np.sum(sholl_results)
            if total_intersections > 0:
                # AUC usando integración trapezoidal
                auc = np.trapz(sholl_results, radii_um)
                
                # Radio crítico (primer radio con intersecciones)
                nonzero_indices = np.where(sholl_results > 0)[0]
                critical_radius = radii_um[nonzero_indices[0]] if len(nonzero_indices) > 0 else 0.0
                
                # Pico de intersecciones
                peak_intersections = np.max(sholl_results)
                peak_radius = radii_um[np.argmax(sholl_results)]
            else:
                auc = critical_radius = peak_intersections = peak_radius = 0.0
            
            summary_results.append({
                'label': int(label),
                'auc': float(auc),
                'critical_radius_um': float(critical_radius),
                'peak_intersections': int(peak_intersections),
                'peak_radius_um': float(peak_radius),
                'total_intersections': int(total_intersections),
                'method': 'skan_native'
            })
            
            # Guardar anillos para visualización (agregar Z=0 para compatibilidad)
            rings_data[str(label)] = {
                'centroid_um': [0.0, float(centroid_um[0]), float(centroid_um[1])],  # Z,Y,X
                'radii_um': radii_um.tolist()
            }
            
        except Exception as e:
            print(f"Error procesando Sholl para label {label}: {e}")
            continue
    
    # Crear DataFrames
    df_sholl = pd.DataFrame(all_results)
    df_summary = pd.DataFrame(summary_results)
    
    # Guardar resultados
    if not df_sholl.empty:
        df_sholl.to_csv(out_dir / "sholl_2d_native.csv", index=False)
    if not df_summary.empty:
        df_summary.to_csv(out_dir / "sholl_summary_2d_native.csv", index=False)
    
    # Guardar anillos
    if rings_data:
        rings_path = out_dir / "sholl_rings_2d_native.json"
        with open(rings_path, 'w') as f:
            import json
            json.dump(rings_data, f, indent=2)
    
    print(f"Análisis Sholl SKAN nativo completado. {len(labels)} astrocitos procesados.")
    return df_sholl, df_summary


def run_sholl_native_skan(
    skeleton_labels: np.ndarray,
    mask_labels: np.ndarray,
    cal: dict,
    out_dir: Path,
    is_2d: bool = False
) -> pd.DataFrame:
    """
    Análisis de Sholl usando skan.sholl_analysis nativo.
    Trabaja tanto con datos 2D como 3D.
    """
    if sk_summarize is None or Skeleton is None:
        print("Advertencia: skan no está disponible, usando método manual.")
        return run_sholl_manual(skeleton_labels, mask_labels, cal, out_dir, is_2d)
    
    min_r = float(cal.get("SHOLL_MIN_RADIUS_UM", 5.0))
    max_r = float(cal.get("SHOLL_MAX_RADIUS_UM", 100.0))
    step_r = float(cal.get("SHOLL_STEP_UM", 2.0))
    
    # Configurar spacing según dimensionalidad
    if is_2d:
        y_um, x_um = float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        spacing = (y_um, x_um)
    else:
        z_um, y_um, x_um = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        spacing = (z_um, y_um, x_um)
    
    # Generar radios para Sholl
    radii = np.arange(min_r, max_r + step_r, step_r)
    
    results = []
    labels = np.unique(mask_labels)
    labels = labels[labels > 0]
    
    print(f"Ejecutando Sholl {'2D' if is_2d else '3D'} nativo de SKAN para {len(labels)} astrocitos...")
    
    for lab in labels:
        try:
            # Obtener esqueleto de esta etiqueta
            skel_mask = (skeleton_labels == lab).astype(np.uint8)
            if not np.any(skel_mask):
                continue
            
            # Obtener centro desde la máscara (centroide del núcleo)
            mask_this = (mask_labels == lab)
            if not np.any(mask_this):
                continue
                
            # Calcular centroide de la máscara
            coords = np.nonzero(mask_this)
            if is_2d:
                center = (np.mean(coords[0]) * spacing[0], np.mean(coords[1]) * spacing[1])
            else:
                center = (np.mean(coords[0]) * spacing[0], 
                         np.mean(coords[1]) * spacing[1], 
                         np.mean(coords[2]) * spacing[2])
            
            # Crear objeto Skeleton de skan
            skel_obj = Skeleton(skel_mask, spacing=spacing)
            
            # Ejecutar análisis de Sholl nativo
            center_result, shell_radii, intersections = skan.sholl_analysis(
                skel_obj, 
                center=center,
                shells=radii
            )
            
            # Calcular métricas adicionales
            if len(intersections) > 0:
                peak_idx = np.argmax(intersections)
                critical_radius_um = float(shell_radii[peak_idx])
                peak_intersections = float(intersections[peak_idx])
                total_intersections = float(np.sum(intersections))
                
                # Calcular AUC (área bajo la curva)
                from scipy.integrate import trapezoid
                auc = float(trapezoid(intersections, shell_radii))
                
                results.append({
                    'label': int(lab),
                    'critical_radius_um': critical_radius_um,
                    'peak_intersections': peak_intersections,
                    'total_intersections': total_intersections,
                    'auc': auc,
                    'center_y_um': float(center_result[0] if is_2d else center_result[1]),
                    'center_x_um': float(center_result[1] if is_2d else center_result[2]),
                    'center_z_um': 0.0 if is_2d else float(center_result[0]),
                    'is_2d': is_2d
                })
                
        except Exception as e:
            print(f"Error en Sholl SKAN para label {lab}: {e}")
            continue
    
    if results:
        df = pd.DataFrame(results)
        suffix = "_2d" if is_2d else ""
        df.to_csv(out_dir / f"sholl{suffix}.csv", index=False)
        print(f"Sholl {'2D' if is_2d else '3D'} completado: {len(results)} astrocitos analizados")
        return df
    else:
        print("No se pudieron analizar astrocitos con Sholl SKAN")
        return pd.DataFrame()

def run_sholl_and_save(
    out_dir: Path,
    cal: dict,
    restrict_to_final: bool = True,
    save_rings_json: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    DEPRECATED: El análisis de Sholl ahora se ejecuta automáticamente
    junto con la esqueletización en run_unified_2d_skeleton_and_sholl().
    
    Esta función se mantiene por compatibilidad pero ya no hace nada.
    Los resultados se encuentran en: sholl_2d_native.csv
    """
    print("⚠️  run_sholl_and_save() está deprecated")
    print("   El análisis de Sholl se ejecuta automáticamente con la esqueletización")
    
    # Cargar resultados existentes si están disponibles
    sholl_path = out_dir / "sholl_2d_native.csv"
    if sholl_path.exists():
        df_sholl = pd.read_csv(sholl_path)
        return pd.DataFrame(), df_sholl  # Retornar vacío 3D, lleno 2D
    
    return pd.DataFrame(), pd.DataFrame()

def run_sholl_manual(
    skeleton_labels: np.ndarray,
    mask_labels: np.ndarray, 
    cal: dict,
    out_dir: Path,
    is_2d: bool = False
) -> pd.DataFrame:
    """
    Análisis de Sholl manual (respaldo cuando SKAN no está disponible).
    """
    min_r = float(cal.get("SHOLL_MIN_RADIUS_UM", 5.0))
    max_r = float(cal.get("SHOLL_MAX_RADIUS_UM", 100.0))
    step_r = float(cal.get("SHOLL_STEP_UM", 2.0))
    
    # Configurar spacing según dimensionalidad
    if is_2d:
        y_um, x_um = float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        spacing = np.array([y_um, x_um])
    else:
        z_um, y_um, x_um = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
        spacing = np.array([z_um, y_um, x_um])
    
    radii = np.arange(min_r, max_r + step_r, step_r)
    labels = np.unique(mask_labels)
    labels = labels[labels > 0]
    
    results = []
    
    for lab in labels:
        try:
            # Obtener centroide de la máscara
            mask_this = (mask_labels == lab)
            if not np.any(mask_this):
                continue
                
            coords = np.nonzero(mask_this)
            if is_2d:
                centroid = np.array([np.mean(coords[0]), np.mean(coords[1])])
            else:
                centroid = np.array([np.mean(coords[0]), np.mean(coords[1]), np.mean(coords[2])])
            
            # Obtener puntos del esqueleto
            skel_points = np.argwhere(skeleton_labels == lab)
            if len(skel_points) == 0:
                continue
            
            # Convertir a coordenadas físicas
            centroid_um = centroid * spacing
            points_um = skel_points.astype(np.float32) * spacing
            
            # Calcular distancias desde el centroide
            distances = np.sqrt(np.sum((points_um - centroid_um)**2, axis=1))
            
            # Contar intersecciones por anillo
            intersections, _ = np.histogram(distances, bins=np.append(radii - step_r/2, radii[-1] + step_r/2))
            
            if len(intersections) > 0 and np.max(intersections) > 0:
                peak_idx = np.argmax(intersections)
                critical_radius_um = float(radii[peak_idx])
                peak_intersections = float(intersections[peak_idx])
                total_intersections = float(np.sum(intersections))
                
                # Calcular AUC usando método trapezoidal
                auc = float(np.trapz(intersections, radii))
                
                results.append({
                    'label': int(lab),
                    'critical_radius_um': critical_radius_um,
                    'peak_intersections': peak_intersections,
                    'total_intersections': total_intersections,
                    'auc': auc,
                    'center_y_um': float(centroid_um[0] if is_2d else centroid_um[1]),
                    'center_x_um': float(centroid_um[1] if is_2d else centroid_um[2]),
                    'center_z_um': 0.0 if is_2d else float(centroid_um[0]),
                    'is_2d': is_2d
                })
                
        except Exception as e:
            print(f"Error en Sholl manual para label {lab}: {e}")
            continue
    
    if results:
        df = pd.DataFrame(results)
        suffix = "_2d" if is_2d else ""
        df.to_csv(out_dir / f"sholl{suffix}.csv", index=False)
        return df
    else:
        return pd.DataFrame()