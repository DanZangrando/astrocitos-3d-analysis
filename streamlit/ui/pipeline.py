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
    from skan import Skeleton, summarize as sk_summarize
except Exception:
    sk_summarize = None
    Skeleton = None

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

def _process_nucleus_filter(nucleus_props, cellpose_masks, gfap_channel, microglia_channel, se, gfap_decision_thr, micro_decision_thr, voxel_vol_um3):
    """Función worker para procesar un solo núcleo (usada por ProcessPoolExecutor)."""
    label = nucleus_props.label
    nucleus_mask = (cellpose_masks == label)
    
    dilated_mask = binary_dilation(nucleus_mask, footprint=se)
    shell_mask = dilated_mask & ~nucleus_mask
    
    shell_gfap_intensity = 0.0
    shell_microglia_intensity = 0.0
    
    if np.any(shell_mask):
        shell_gfap_intensity = float(gfap_channel[shell_mask].mean())
        shell_microglia_intensity = float(microglia_channel[shell_mask].mean())
    
    is_astrocyte = (
        shell_gfap_intensity > gfap_decision_thr and 
        shell_microglia_intensity < micro_decision_thr
    )
    
    volume_vox = nucleus_props.area
    volume_um3 = float(volume_vox) * voxel_vol_um3
    sphericity = np.nan
    if hasattr(nucleus_props, 'surface_area') and nucleus_props.surface_area > 0:
        try:
            sphericity = (np.pi**(1/3) * (6 * volume_vox)**(2/3)) / nucleus_props.surface_area
        except (ValueError, ZeroDivisionError):
            pass

    metric_dict = {
        "label": label,
        "nucleus_volume_um3": volume_um3,
        "nucleus_sphericity": sphericity,
        "shell_gfap_mean": shell_gfap_intensity,
        "shell_microglia_mean": shell_microglia_intensity,
        "is_astrocyte_candidate": is_astrocyte
    }
    
    return label if is_astrocyte else None, metric_dict


def run_filter_and_save(
    cellpose_masks: np.ndarray, 
    gfap_channel: np.ndarray, 
    microglia_channel: np.ndarray,
    otsu_mask: np.ndarray,
    cal: dict,
    out_dir: Path
) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Filtra núcleos de Cellpose basado en umbrales relativos (std dev) 
    en un shell alrededor del núcleo. (Ejecución PARALELIZADA)
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
    micro_std_thr = float(cal.get("MICROGLIA_STD_DEV_THRESHOLD", 5.0))

    gfap_bg_mean, gfap_bg_std = _calculate_background_stats(gfap_channel, otsu_mask)
    micro_bg_mean, micro_bg_std = _calculate_background_stats(microglia_channel, otsu_mask)
    
    gfap_decision_thr = gfap_bg_mean + (gfap_std_thr * gfap_bg_std)
    micro_decision_thr = micro_bg_mean + (micro_std_thr * micro_bg_std)

    print(f"--- Filtrado Relativo ---")
    print(f"GFAP Fondo: {gfap_bg_mean:.2f} ± {gfap_bg_std:.2f}. Umbral Decisión: > {gfap_decision_thr:.2f}")
    print(f"Micro Fondo: {micro_bg_mean:.2f} ± {micro_bg_std:.2f}. Umbral Decisión: < {micro_decision_thr:.2f}")

    nuclei_props = regionprops(cellpose_masks)
    astrocyte_labels = []
    metrics_data = []

    task_func = partial(
        _process_nucleus_filter,
        cellpose_masks=cellpose_masks,
        gfap_channel=gfap_channel,
        microglia_channel=microglia_channel,
        se=se,
        gfap_decision_thr=gfap_decision_thr,
        micro_decision_thr=micro_decision_thr,
        voxel_vol_um3=voxel_vol_um3
    )
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(task_func, nuclei_props))

    for label, metric_dict in results:
        if label is not None:
            astrocyte_labels.append(label)
        metrics_data.append(metric_dict)

    gfap_filtered_mask = np.where(np.isin(cellpose_masks, astrocyte_labels), cellpose_masks, 0)
    tifffile.imwrite(out_dir / "03_gfap_microglia_filtered_mask.tif", gfap_filtered_mask.astype(np.uint16))
    
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


def run_advanced_skeleton_and_save(
    img_path: Path,
    mask_path: Path,
    out_dir: Path,
    cal: dict,
    conflict_resolve: bool = True
) -> pd.DataFrame:
    """
    Lógica completa de esqueletización (de 04_...) y análisis con Skan.
    (Ejecución PARALELIZADA)
    """
    
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


def run_sholl_and_save(
    out_dir: Path,
    cal: dict,
    restrict_to_final: bool = True,
    save_rings_json: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Ejecuta el análisis de Sholl y guarda los resultados.
    (Ejecución PARALELIZADA)
    """
    
    min_r = float(cal.get("SHOLL_MIN_RADIUS_UM", 0.0))
    max_r = float(cal.get("SHOLL_MAX_RADIUS_UM", 50.0))
    step_r = float(cal.get("SHOLL_STEP_UM", 20.0))
    if step_r <= 0:
        raise ValueError("El paso de Sholl (SHOLL_STEP_UM) debe ser positivo.")

    z_um, y_um, x_um = float(cal.get('z', 1.0)), float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
    spacing = (z_um, y_um, x_um)

    skeleton_path = out_dir / "05_skeleton_labels.tif"
    mask_path = out_dir / "04_final_astrocytes_mask.tif" if restrict_to_final else (out_dir / "02_cellpose_mask.tif")
    
    if not skeleton_path.exists() or not mask_path.exists():
        print(f"Advertencia: Faltan archivos para Sholl (Skeleton o Mask) en {out_dir}")
        return pd.DataFrame(), pd.DataFrame()

    skel_labels = tifffile.imread(skeleton_path)
    mask_labels = tifffile.imread(mask_path)
    
    labels = np.unique(mask_labels)
    labels = labels[labels > 0]
    
    radii = np.arange(min_r, max_r + step_r, step_r)
    all_results_rows = []
    all_rings_data = {}

    props = regionprops(mask_labels)
    centroids = {p.label: p.centroid for p in props}

    skel_points = {}
    for lab in labels:
        if lab not in centroids: continue
        points = np.argwhere(skel_labels == lab)
        if points.shape[0] > 0:
            skel_points[lab] = points

    labels_to_process = [lab for lab in labels if lab in skel_points]
    if not labels_to_process:
        print(f"No se encontraron puntos de esqueleto para las etiquetas en {out_dir}")
        return pd.DataFrame(), pd.DataFrame()

    task_func = partial(
        _process_sholl_label,
        centroids=centroids,
        skel_points=skel_points,
        spacing=spacing,
        radii=radii,
        save_rings_json=save_rings_json
    )
    
    with ProcessPoolExecutor() as executor:
        results = list(executor.map(task_func, labels_to_process))

    for (results_list, rings_data) in results:
        all_results_rows.extend(results_list)
        if rings_data:
            all_rings_data.update(rings_data)

    df_sholl = pd.DataFrame(all_results_rows)
    
    df_summary = pd.DataFrame()
    if not df_sholl.empty:
        df_summary = df_sholl.sholl.get_scalar_metrics()

    # Guardar resultados
    df_sholl.to_csv(out_dir / "sholl.csv", index=False)
    df_summary.to_csv(out_dir / "sholl_summary.csv", index=False)
    if save_rings_json:
        (out_dir / "sholl_rings.json").write_text(json.dumps(all_rings_data, indent=2))
        
    return df_sholl, df_summary