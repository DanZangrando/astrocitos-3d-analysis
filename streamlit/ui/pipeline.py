from __future__ import annotations
from pathlib import Path
import json
from typing import Tuple

import numpy as np
import tifffile
from skimage.filters import threshold_otsu
from scipy import ndimage as ndi

# Optional deps
try:
    from cellpose import models as _cp_models  # type: ignore
except Exception:  # pragma: no cover
    _cp_models = None
try:
    from skimage.morphology import skeletonize_3d as _skeletonize_3d  # type: ignore
except Exception:  # pragma: no cover
    _skeletonize_3d = None
try:
    from skimage.morphology import skeletonize as _skeletonize_2d  # type: ignore
except Exception:  # pragma: no cover
    _skeletonize_2d = None


def _read_calibration(cal_path: Path) -> dict:
    try:
        return json.loads(cal_path.read_text()) if cal_path.exists() else {}
    except Exception:
        return {}


def load_image_any(path: Path) -> tuple[np.ndarray, str | None]:
    """Load .tif/.tiff and return (array, axes)."""
    suffix = path.suffix.lower()
    if suffix in (".tif", ".tiff"):
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            axes = getattr(series, 'axes', None)
            arr = series.asarray()
        return arr, axes
    raise ValueError(f"Extensión no soportada: {suffix}")


def reorder_to_zcyx(arr: np.ndarray, axes: str | None) -> np.ndarray:
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


def run_otsu_and_save(dapi: np.ndarray, out_dir: Path) -> float:
    thr = float(threshold_otsu(dapi))
    otsu_mask = (dapi > thr).astype(np.uint8)
    tifffile.imwrite(out_dir / "01_otsu_mask.tif", otsu_mask)
    # Persistencia local deshabilitada para reproducibilidad (parámetros son globales)
    return thr


def run_cellpose_and_save(dapi_vol: np.ndarray, out_dir: Path, nucleus_diameter: int = 30, use_gpu: bool = True) -> np.ndarray:
    if _cp_models is None:
        raise RuntimeError("Cellpose no disponible. Instalá cellpose para este paso.")
    try:
        model = _cp_models.CellposeModel(gpu=bool(use_gpu))
    except Exception:
        model = _cp_models.CellposeModel(gpu=False)
    masks, _, _ = model.eval(dapi_vol, diameter=int(nucleus_diameter), z_axis=0, do_3D=True)
    tifffile.imwrite(out_dir / "02_cellpose_mask.tif", masks.astype(np.uint16))
    # Persistencia local deshabilitada (usar calibration.json global)
    return masks


def run_filter_and_save(cellpose_masks: np.ndarray, gfap: np.ndarray, microglia: np.ndarray,
                        gfap_thr: float, microglia_thr: float, shell_um: float,
                        spacing_um: Tuple[float,float,float], out_dir: Path) -> np.ndarray:
    """Simple GFAP/Microglia shell filtering around each nucleus label."""
    z_um, y_um, x_um = spacing_um
    # Cascarón ~ 2 vóxeles en Y/X por defecto si shell_um no definido
    r_vox = max(1, int(round((shell_um or (2.0 * y_um)) / y_um)))
    labels = [int(v) for v in np.unique(cellpose_masks) if v > 0]
    keep = []
    for lab in labels:
        sel = (cellpose_masks == lab)
        dil = ndi.binary_dilation(sel, iterations=r_vox)
        shell = np.logical_and(dil, np.logical_not(sel))
        mic_ok = float(microglia[shell].mean()) if np.any(shell) else 0.0
        gf_ok = float(gfap[shell].mean()) if np.any(shell) else 0.0
        if mic_ok <= float(microglia_thr) and gf_ok >= float(gfap_thr):
            keep.append(lab)
    gfap_filtered = np.where(np.isin(cellpose_masks, keep), cellpose_masks, 0).astype(np.uint16)
    tifffile.imwrite(out_dir / "03_gfap_microglia_filtered_mask.tif", gfap_filtered)
    # Persistencia local deshabilitada (usar calibration.json global)
    return gfap_filtered


def run_size_filter_and_save(gfap_filtered: np.ndarray, min_volume_um3: float, spacing_um: Tuple[float,float,float], out_dir: Path) -> np.ndarray:
    z_um, y_um, x_um = spacing_um
    vox_um3 = float(z_um) * float(y_um) * float(x_um)
    kept = []
    for lab in [int(v) for v in np.unique(gfap_filtered) if v > 0]:
        v = int((gfap_filtered == lab).sum()) * vox_um3
        if v >= float(min_volume_um3):
            kept.append(lab)
    final_mask = np.where(np.isin(gfap_filtered, kept), gfap_filtered, 0).astype(np.uint16)
    tifffile.imwrite(out_dir / "04_final_astrocytes_mask.tif", final_mask)
    # Persistencia local deshabilitada (usar calibration.json global)
    return final_mask


def run_minimal_skeleton_and_save(final_mask: np.ndarray, gfap: np.ndarray, spacing_um: Tuple[float,float,float], out_dir: Path) -> Path:
    """Very simplified skeletonization: threshold GFAP per cell ROI and skeletonize3d.
    Saves 05_skeleton_labels.tif
    """
    labels = [int(v) for v in np.unique(final_mask) if v > 0]
    skel_labels = np.zeros_like(final_mask, dtype=np.uint16)
    # métricas mínimas por etiqueta
    metrics = []
    for lab in labels:
        sel = (final_mask == lab)
        # Simple local otsu on ROI
        if not np.any(sel):
            continue
        z0,z1 = np.where(sel)[0].min(), np.where(sel)[0].max()+1
        y0,y1 = np.where(sel)[1].min(), np.where(sel)[1].max()+1
        x0,x1 = np.where(sel)[2].min(), np.where(sel)[2].max()+1
        roi_gf = gfap[z0:z1, y0:y1, x0:x1]
        thr = threshold_otsu(roi_gf)
        bin_roi = (roi_gf >= thr)
        # skeletonize con fallback 2D si 3D no está disponible
        if _skeletonize_3d is not None:
            sk = _skeletonize_3d(bin_roi.astype(np.uint8)) > 0
        elif _skeletonize_2d is not None:
            sk = np.zeros_like(bin_roi, dtype=bool)
            for zz in range(bin_roi.shape[0]):
                try:
                    sk[zz] = _skeletonize_2d(bin_roi[zz].astype(bool))
                except Exception:
                    sk[zz] = bin_roi[zz] > 0
        else:
            # último recurso: usar binario original como "skeleton" (degradado)
            sk = bin_roi > 0
        # write back
        dst = skel_labels[z0:z1, y0:y1, x0:x1]
        dst[sk] = lab
        # métricas mínimas
        skel_vox = int(np.count_nonzero(sk))
        # aproximación de longitud: usar el menor paso (conservador)
        min_step = float(min(spacing_um)) if spacing_um else 1.0
        approx_len_um = float(skel_vox) * float(min_step)
        metrics.append({"label": int(lab), "skeleton_voxels": skel_vox, "approx_length_um": approx_len_um})
    out_path = out_dir / "05_skeleton_labels.tif"
    tifffile.imwrite(out_path, skel_labels)
    # guardar métricas mínimas si hay
    if metrics:
        try:
            import pandas as _pd
            skel_dir = out_dir / "skeletons"
            skel_dir.mkdir(parents=True, exist_ok=True)
            _pd.DataFrame(metrics).to_csv(skel_dir / "summary.csv", index=False)
        except Exception:
            pass
    return out_path


def _merge_params(out_dir: Path, updates: dict):
    """Deprecated: parámetros por preparado deshabilitados para asegurar reproducibilidad."""
    return None

