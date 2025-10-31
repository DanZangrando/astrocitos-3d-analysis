from __future__ import annotations
from pathlib import Path
from typing import Tuple
import json

import numpy as np
import tifffile
from skimage.morphology import binary_closing, binary_dilation, ball
from skimage.measure import regionprops, label as cc_label
from scipy.ndimage import distance_transform_edt
from skimage.transform import resize

# skeletonize fallbacks
try:
    from skimage.morphology import skeletonize_3d as _skeletonize_3d  # type: ignore
except Exception:  # pragma: no cover
    _skeletonize_3d = None
try:
    from skimage.morphology import skeletonize as _skeletonize_2d  # type: ignore
except Exception:  # pragma: no cover
    _skeletonize_2d = None


def _zoom_factors(z_um: float, y_um: float, x_um: float, target_iso_um: float):
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


def _um_to_vox(um_val: float, voxel_um: float) -> int:
    return int(round(max(0.0, float(um_val)) / max(voxel_um, 1e-6)))


def _prune_short_spurs(skel: np.ndarray, min_len_um: float, vox_um: float) -> np.ndarray:
    from scipy.ndimage import convolve
    sk = skel.copy()
    if min_len_um <= 0:
        return sk
    kernel = np.ones((3, 3, 3), dtype=np.int16)
    kernel[1, 1, 1] = 0
    offsets = []
    for dz in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dz == 0 and dy == 0 and dx == 0:
                    continue
                offsets.append((dz, dy, dx))
    offsets = np.array(offsets, dtype=np.int8)

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
                if deg_curr <= 1:
                    break
                if deg_curr >= 3:
                    break
                if len(nbrs) == 0:
                    break
                nxt = nbrs[0]
                dz = abs(nxt[0] - curr[0]); dy = abs(nxt[1] - curr[1]); dx = abs(nxt[2] - curr[2])
                step_um = float(np.sqrt(dz*dz + dy*dy + dx*dx)) * float(vox_um)
                length_um += step_um
                path.append(nxt)
                prev, curr = curr, nxt
            if length_um < float(min_len_um):
                if int(deg[curr]) >= 3 and len(path) > 0:
                    to_remove.extend(path[:-1])
                else:
                    to_remove.extend(path)
        if not to_remove:
            break
        rr = tuple(np.array(to_remove).T)
        sk[rr] = False
    return sk


def run_advanced_skeleton_and_save(
    img_path: Path,
    out_dir: Path,
    gfap_idx: int,
    cal: dict,
    conflict_resolve: bool = True,
) -> Path:
    """Runs advanced skeletonization with isotropic resampling, ROI-based threshold, Voronoi territories (optional),
    pruning (optional), tube-based GFAP intensity, and summary metrics. Writes:
      - 05_skeleton_labels.tif
      - skeletons/summary.csv
    Uses parameters from cal (global or merged), with keys like SKELETON_*.
    """
    out_dir = Path(out_dir)
    final_mask_path = out_dir / "04_final_astrocytes_mask.tif"
    if not final_mask_path.exists():
        raise FileNotFoundError("04_final_astrocytes_mask.tif no encontrado")

    # Load image and channels
    from .pipeline import load_image_any, reorder_to_zcyx
    arr, axes = load_image_any(img_path)
    vol = reorder_to_zcyx(arr, axes)  # Z,C,Y,X
    gfap = vol[:, int(gfap_idx), :, :].astype(np.float32)
    mask = tifffile.imread(final_mask_path)
    if mask.shape != gfap.shape:
        raise RuntimeError("Dimensiones incompatibles entre GFAP y máscara final")

    # Calibration and skeleton params
    z_um = float(cal.get('z', 1.0)); y_um = float(cal.get('y', 1.0)); x_um = float(cal.get('x', 1.0))
    target_iso_um = float(cal.get("SKELETON_TARGET_ISO_UM", min(z_um, y_um, x_um)))
    padding_um = float(cal.get("SKELETON_PADDING_UM", 2.0))
    seed_dilate_um = float(cal.get("SKELETON_SEED_DILATE_UM", 2.0))
    connectivity = int(cal.get("SKELETON_CONNECTIVITY", 26))
    closing_um = float(cal.get("SKELETON_CLOSING_UM", 0.8))
    max_radius_um = float(cal.get("SKELETON_MAX_RADIUS_UM", 10.0 * target_iso_um))
    thresh_mode = str(cal.get("SKELETON_THRESHOLD_MODE", "otsu")).lower()
    manual_thr = float(cal.get("SKELETON_MANUAL_THRESHOLD", 50.0))
    territory_voronoi = bool(cal.get("SKELETON_TERRITORY_VORONOI", False))
    territory_excl_um = float(cal.get("SKELETON_TERRITORY_EXCLUSION_UM", 1.0))
    domain_volume_source = str(cal.get("SKELETON_DOMAIN_VOLUME_SOURCE", "gfap")).lower()
    prune_enable = bool(cal.get("SKELETON_PRUNE_ENABLE", False))
    prune_min_len_um = float(cal.get("SKELETON_PRUNE_MIN_LEN_UM", max(2.0, 2.0 * target_iso_um)))
    tube_radius_um = float(cal.get("SKELETON_TUBE_RADIUS_UM", 1.5))

    skel_dir = out_dir / "skeletons"
    skel_dir.mkdir(parents=True, exist_ok=True)
    combined_labels = np.zeros_like(mask, dtype=np.uint16)
    combined_dist_um = np.full(mask.shape, np.inf, dtype=np.float32) if conflict_resolve else None

    # nucleus centroids for conflict resolution
    label_centroid_um = {}
    if conflict_resolve:
        for p in regionprops(mask.astype(np.int32)):
            lab_id = int(p.label)
            cz, cy, cx = p.centroid
            label_centroid_um[lab_id] = (
                float(cz) * float(z_um),
                float(cy) * float(y_um),
                float(cx) * float(x_um),
            )

    zf, yf, xf = _zoom_factors(z_um, y_um, x_um, target_iso_um)
    closing_r_vox = _um_to_vox(closing_um, target_iso_um) if closing_um > 0 else 0
    seed_r_vox = _um_to_vox(seed_dilate_um, target_iso_um) if seed_dilate_um > 0 else 0
    max_r_vox = _um_to_vox(max_radius_um, target_iso_um) if max_radius_um > 0 else 0
    excl_r_vox = _um_to_vox(territory_excl_um, target_iso_um) if territory_excl_um > 0 else 0
    tube_r_vox = _um_to_vox(tube_radius_um, target_iso_um) if tube_radius_um > 0 else 0

    pad_z = _um_to_vox(padding_um + max_radius_um, z_um)
    pad_y = _um_to_vox(padding_um + max_radius_um, y_um)
    pad_x = _um_to_vox(padding_um + max_radius_um, x_um)

    labels = np.unique(mask)
    labels = labels[labels > 0]

    metrics = []

    for lab in labels:
        terr_lab = None
        reg_mask = (mask == lab)
        if not np.any(reg_mask):
            continue
        props = regionprops(reg_mask.astype(np.uint8))
        if not props:
            continue
        z0, y0, x0, z1, y1, x1 = (*props[0].bbox[0:3], *props[0].bbox[3:6])
        z0 = max(0, z0 - pad_z); y0 = max(0, y0 - pad_y); x0 = max(0, x0 - pad_x)
        z1 = min(mask.shape[0], z1 + pad_z); y1 = min(mask.shape[1], y1 + pad_y); x1 = min(mask.shape[2], x1 + pad_x)

        roi_mask = reg_mask[z0:z1, y0:y1, x0:x1]
        roi_gfap = gfap[z0:z1, y0:y1, x0:x1].astype(np.float32)
        roi_all_labels = mask[z0:z1, y0:y1, x0:x1]

        roi_mask_iso = _resize_to_isotropic(roi_mask.astype(np.uint8), (zf, yf, xf), order=0) > 0.5
        roi_gfap_iso = _resize_to_isotropic(roi_gfap, (zf, yf, xf), order=1)
        roi_all_labels_iso = _resize_to_isotropic(roi_all_labels.astype(np.int32), (zf, yf, xf), order=0).astype(np.int32)

        if thresh_mode.startswith("manual"):
            thr = float(manual_thr)
        else:
            from skimage.filters import threshold_otsu
            thr = threshold_otsu(roi_gfap_iso[roi_mask_iso]) if np.any(roi_mask_iso) else threshold_otsu(roi_gfap_iso)
        bin_gfap_iso = (roi_gfap_iso >= thr)

        if closing_r_vox > 0:
            try:
                bin_gfap_iso = binary_closing(bin_gfap_iso, ball(int(closing_r_vox)))
            except Exception:
                pass

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

        if territory_voronoi:
            labs_here = np.unique(roi_all_labels_iso)
            labs_here = labs_here[labs_here > 0]
            if labs_here.size > 0:
                props_here = regionprops(roi_all_labels_iso)
                centroids = {int(p.label): p.centroid for p in props_here}
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
                    from scipy.ndimage import distance_transform_edt as _edt
                    interior = _edt(terr_lab) >= float(excl_r_vox)
                    terr_lab = terr_lab & interior
                bin_iso &= terr_lab

        if max_r_vox > 0 and np.any(bin_iso):
            dist = distance_transform_edt(~seed_iso)
            bin_iso &= (dist <= float(max_r_vox))

        if _skeletonize_3d is not None:
            skel_iso = _skeletonize_3d(bin_iso.astype(np.uint8)) > 0
        elif _skeletonize_2d is not None:
            skel_iso = np.zeros_like(bin_iso, dtype=bool)
            for zz in range(bin_iso.shape[0]):
                try:
                    skel_iso[zz] = _skeletonize_2d(bin_iso[zz].astype(bool))
                except Exception:
                    skel_iso[zz] = bin_iso[zz]
        else:
            skel_iso = bin_iso

        if prune_enable and float(prune_min_len_um) > 0:
            try:
                skel_iso = _prune_short_spurs(skel_iso, float(prune_min_len_um), float(target_iso_um))
            except Exception:
                pass

        skel_voxels = int(np.count_nonzero(skel_iso))
        approx_length_um = skel_voxels * float(target_iso_um)
        gfap_connected_vox = int(np.count_nonzero(bin_iso))
        gfap_connected_um3 = float(gfap_connected_vox) * float(target_iso_um) ** 3
        voronoi_um3 = np.nan
        if territory_voronoi:
            try:
                vor_vox = int(np.count_nonzero(terr_lab)) if terr_lab is not None else 0
                voronoi_um3 = float(vor_vox) * float(target_iso_um) ** 3
            except Exception:
                voronoi_um3 = np.nan
        domain_um3 = gfap_connected_um3 if domain_volume_source.startswith("gfap") else voronoi_um3

        tube_sum_intensity = np.nan
        tube_mean_intensity = np.nan
        tube_voxels = 0
        tube_volume_um3 = np.nan
        if tube_r_vox > 0 and np.any(skel_iso):
            try:
                dist_to_skel = distance_transform_edt(~skel_iso)
                tube_mask = dist_to_skel <= float(tube_r_vox)
                tube_sel = tube_mask & (roi_gfap_iso >= float(thr))
                tube_voxels = int(np.count_nonzero(tube_sel))
                tube_sum_intensity = float(np.sum(roi_gfap_iso[tube_sel])) if tube_voxels > 0 else 0.0
                tube_mean_intensity = float(np.mean(roi_gfap_iso[tube_sel])) if tube_voxels > 0 else 0.0
                tube_volume_um3 = float(tube_voxels) * float(target_iso_um) ** 3
            except Exception:
                pass

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

        skel_roi = _resize_to_original(skel_iso.astype(np.uint8), roi_mask.shape, order=0) > 0.5

        target_slice = combined_labels[z0:z1, y0:y1, x0:x1]
        if conflict_resolve:
            dist_slice = combined_dist_um[z0:z1, y0:y1, x0:x1]
            cz_um, cy_um, cx_um = label_centroid_um.get(int(lab), (None, None, None))
            if cz_um is None:
                target_slice[skel_roi] = lab
            else:
                zz = (np.arange(z0, z1, dtype=np.float32) * float(z_um))[:, None, None]
                yy = (np.arange(y0, y1, dtype=np.float32) * float(y_um))[None, :, None]
                xx = (np.arange(x0, x1, dtype=np.float32) * float(x_um))[None, None, :]
                dist_um = np.sqrt((zz - cz_um) ** 2 + (yy - cy_um) ** 2 + (xx - cx_um) ** 2)
                cand = skel_roi
                better = cand & ((target_slice == 0) | (dist_um < dist_slice))
                target_slice[better] = lab
                dist_slice[better] = dist_um[better]
        else:
            target_slice[skel_roi] = lab

        tifffile.imwrite(skel_dir / f"astro_{int(lab)}_skeleton_roi.tif", skel_roi.astype(np.uint8))
        metrics.append({
            "label": int(lab),
            "skeleton_voxels_iso": int(skel_voxels),
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

    # write combined labels and metrics
    out_skel_labels = out_dir / "05_skeleton_labels.tif"
    tifffile.imwrite(out_skel_labels, combined_labels.astype(np.uint16))
    if metrics:
        import pandas as pd
        df = pd.DataFrame(metrics)
        # Enriquecimiento opcional con Skan (longitudes por rama, endpoints, junctions)
        try:
            from skan import Skeleton, summarize as sk_summarize
            # construir resumen por etiqueta usando combined_labels
            skel_lab_img = combined_labels
            labels_all = [int(v) for v in np.unique(skel_lab_img) if v > 0]
            rows = []
            for lab in labels_all:
                vox = (skel_lab_img == lab)
                # Afinar por si acaso
                vox_thin = vox
                try:
                    if _skeletonize_3d is not None:
                        vox_thin = (_skeletonize_3d(vox.astype(np.uint8)) > 0)
                except Exception:
                    pass
                try:
                    sk = Skeleton(vox_thin.astype(np.uint8), spacing=(z_um, y_um, x_um))
                    dfb = sk_summarize(sk, separator="_")
                    if dfb is not None and not dfb.empty:
                        # Helper: try multiple possible column names for branch length and branch type
                        def _try_get_series(dframe, candidates):
                            for c in candidates:
                                if c in dframe.columns:
                                    return dframe[c].to_numpy(), c
                            # fallback: try to find a column that looks like a branch-length column
                            for c in dframe.columns:
                                low = str(c).lower()
                                if "branch" in low and ("dist" in low or "length" in low):
                                    return dframe[c].to_numpy(), c
                            return np.array([], dtype=float), None

                        dist, dist_col = _try_get_series(dfb, ["branch-distance", "branch_distance", "branch_dist", "branch_length", "branch_length_um"])
                        keep = (dist >= float(prune_min_len_um)) if dist.size else np.array([], dtype=bool)
                        # If we have distances, sum them; otherwise fallback to voxel-count-based estimate
                        if dist.size:
                            total_len_um = float(np.nansum(dist[keep])) if keep.size else float(np.nansum(dist))
                            n_branches = int(np.count_nonzero(keep)) if keep.size else int(dfb.shape[0])
                            mean_branch_len_um = float(np.nanmean(dist[keep])) if np.any(keep) else (float(np.nanmean(dist)) if dist.size else np.nan)
                        else:
                            # fallback: estimate from voxels in the original binary mask for this label
                            try:
                                vox_count = int(np.count_nonzero(vox))
                            except Exception:
                                vox_count = 0
                            total_len_um = float(vox_count) * float(target_iso_um)
                            n_branches = int(dfb.shape[0]) if dfb is not None else 0
                            mean_branch_len_um = float(total_len_um / n_branches) if n_branches else np.nan

                        endpoints = None; junctions = None
                        bt, bt_col = _try_get_series(dfb, ["branch-type", "branch_type", "branch_type_id"])
                        if bt is not None and bt.size:
                            if keep.size:
                                endpoints = int(np.count_nonzero((bt == 1) & keep))
                                junctions = int(np.count_nonzero((bt >= 2) & keep))
                            else:
                                endpoints = int(np.count_nonzero(bt == 1))
                                junctions = int(np.count_nonzero(bt >= 2))

                        if total_len_um and total_len_um > 0 and n_branches:
                            density = (n_branches / total_len_um * 100.0)
                        else:
                            density = np.nan

                        # If expected columns were missing, log available columns for debugging
                        if dist_col is None or bt_col is None:
                            try:
                                print(f"[skan] label={lab}: available columns: {list(dfb.columns)}")
                            except Exception:
                                pass

                        rows.append({
                            "label": int(lab),
                            "total_length_um": float(total_len_um),
                            "mean_branch_len_um": float(mean_branch_len_um) if mean_branch_len_um == mean_branch_len_um else np.nan,
                            "n_branches": int(n_branches) if n_branches is not None else np.nan,
                            "endpoints": int(endpoints) if endpoints is not None else np.nan,
                            "junctions": int(junctions) if junctions is not None else np.nan,
                            "branch_density_per_100um": float(density) if not (isinstance(density, float) and np.isnan(density)) else np.nan,
                        })
                except Exception:
                    continue
            if rows:
                import pandas as pd
                df_skan = pd.DataFrame(rows)
                try:
                    df = df.merge(df_skan, on="label", how="left")
                except Exception:
                    pass
        except Exception:
            pass
        df.to_csv(skel_dir / "summary.csv", index=False)
    return out_skel_labels
