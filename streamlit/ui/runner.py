from __future__ import annotations
from pathlib import Path
from typing import Iterable, Literal
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import tifffile

# Reuse pipeline utilities
from .pipeline import (
    load_image_any,
    reorder_to_zcyx,
    run_otsu_and_save,
    run_cellpose_and_save,
    run_filter_and_save,
    run_size_filter_and_save,
)
from .skeleton import run_advanced_skeleton_and_save
from .sholl import compute_sholl, ShollParams

Step = Literal["01", "02", "03", "04", "05"]


def read_calibration(calib_path: Path) -> dict:
    try:
        return json.loads(calib_path.read_text()) if calib_path.exists() else {}
    except Exception:
        return {}


def detect_group(p: Path, root: Path) -> str:
    try:
        rel = str(p.relative_to(root)).lower()
    except Exception:
        rel = str(p).lower()
    if "/hip/" in rel:
        return "Hipoxia"
    return "CTL"


def list_raw_images(root: Path) -> list[Path]:
    raw_dir = root / "data" / "raw"
    return sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])


def _ensure_out_dir(root: Path, img_path: Path) -> Path:
    out_dir = root / "data" / "processed" / img_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir


def run_pipeline_for(
    root: Path,
    img_path: Path,
    cal: dict,
    start_step: Step,
    overwrite_from_step: bool = True,
) -> dict:
    """Run the pipeline from a given step onward for one image.
    Returns a status dict of files produced.
    """
    out_dir = _ensure_out_dir(root, img_path)

    # Read image and split channels
    arr, axes = load_image_any(img_path)
    vol = reorder_to_zcyx(arr, axes)  # Z,C,Y,X
    nC = vol.shape[1]
    dapi_idx = int(cal.get("DAPI_CHANNEL_INDEX", 0))
    gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1 if nC > 1 else 0))
    micro_idx = int(cal.get("MICROGLIA_CHANNEL_INDEX", 2 if nC > 2 else max(1, nC-1)))
    dapi = vol[:, dapi_idx, :, :]
    gfap = vol[:, gfap_idx, :, :]
    micro = vol[:, micro_idx, :, :]

    spacing_um = (
        float(cal.get("z", 1.0)),
        float(cal.get("y", 1.0)),
        float(cal.get("x", 1.0)),
    )

    # Parameters
    nucleus_diam = int(cal.get("NUCLEUS_DIAMETER", cal.get("nucleus_diameter_px", 30)))
    use_gpu = bool(cal.get("CELLPOSE_USE_GPU", cal.get("cellpose_gpu", True)))
    gfap_thr = float(cal.get("GFAP_INTENSITY_THRESHOLD", 40))
    micro_thr = float(cal.get("MICROGLIA_INTENSITY_THRESHOLD", 200))
    min_volume_um3 = float(cal.get("MIN_VOLUME_UM3", 75))

    # Output paths
    p_otsu = out_dir / "01_otsu_mask.tif"
    p_cp = out_dir / "02_cellpose_mask.tif"
    p_gfap = out_dir / "03_gfap_microglia_filtered_mask.tif"
    p_final = out_dir / "04_final_astrocytes_mask.tif"
    p_skel = out_dir / "05_skeleton_labels.tif"

    def need(step: Step) -> bool:
        if not overwrite_from_step:
            return False
        order = {"01":1, "02":2, "03":3, "04":4, "05":5}
        return order[step] >= order[start_step]

    # 01 Otsu
    if need("01") or not p_otsu.exists():
        thr = run_otsu_and_save(dapi, out_dir)
        # optionally return thr if needed

    # 02 Cellpose (use Otsu masked DAPI if available)
    if need("02") or not p_cp.exists():
        if p_otsu.exists():
            otsu = tifffile.imread(p_otsu).astype(bool)
            dapi_in = np.where(otsu, dapi, 0)
        else:
            dapi_in = dapi
        run_cellpose_and_save(dapi_in, out_dir, nucleus_diameter=nucleus_diam, use_gpu=use_gpu)

    # 03 GFAP/Microglía filter
    if need("03") or not p_gfap.exists():
        if not p_cp.exists():
            raise FileNotFoundError("02_cellpose_mask.tif no encontrado")
        cp_masks = tifffile.imread(p_cp)
        run_filter_and_save(cp_masks, gfap, micro, gfap_thr, micro_thr, shell_um=0.0, spacing_um=spacing_um, out_dir=out_dir)
        # 03b Size filter
        run_size_filter_and_save(tifffile.imread(p_gfap), min_volume_um3=min_volume_um3, spacing_um=spacing_um, out_dir=out_dir)

    # 04 Skeleton
    if need("04") or not p_skel.exists():
        merged_cal = {**cal}
        run_advanced_skeleton_and_save(img_path, out_dir, gfap_idx, merged_cal, conflict_resolve=True)

    # 05 Sholl
    sh_min = float(cal.get("SHOLL_MIN_RADIUS_UM", 0.0))
    sh_max = float(cal.get("SHOLL_MAX_RADIUS_UM", 50.0))
    sh_step = float(cal.get("SHOLL_STEP_UM", 20.0))
    if sh_step > 0 and (need("05") or not (out_dir/"sholl.csv").exists()):
        params = ShollParams(sh_min, sh_max, sh_step)
        compute_sholl(out_dir, spacing_um, params, restrict_to_final=True, save_rings_json=False)

    return {
        "otsu": p_otsu.exists(),
        "cellpose": p_cp.exists(),
        "gfap": p_gfap.exists(),
        "final": p_final.exists(),
        "skeleton": p_skel.exists(),
        "sholl": (out_dir/"sholl.csv").exists(),
    }


def run_scope(
    root: Path,
    scope: Literal["selected", "group", "all"],
    start_step: Step,
    cal: dict,
    selected: Path | None = None,
    group: str | None = None,
    overwrite_from_step: bool = True,
) -> list[tuple[Path, dict]]:
    files = list_raw_images(root)
    if scope == "selected":
        targets = [selected] if selected is not None else []
    elif scope == "group" and group in {"CTL", "Hipoxia"}:
        targets = [p for p in files if detect_group(p, root) == group]
    else:
        targets = files

    # Prefer parallel execution for CPU-bound steps 04–05 when processing many images
    results: list[tuple[Path, dict]] = []
    do_parallel = (len(targets) > 1) and (start_step in ("04", "05"))
    max_workers = int(cal.get("MAX_WORKERS", max(1, (os.cpu_count() or 4) - 1)))

    if do_parallel:
        out_map: dict[Path, dict] = {}
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futs = {
                ex.submit(
                    run_pipeline_for,
                    root,
                    p,
                    cal,
                    start_step,
                    overwrite_from_step,
                ): p for p in targets
            }
            for fut in as_completed(futs):
                p = futs[fut]
                try:
                    out_map[p] = fut.result()
                except Exception:
                    out_map[p] = {"error": True}
        # preserve input order
        results = [(p, out_map.get(p, {"error": True})) for p in targets]
    else:
        for p in targets:
            try:
                st = run_pipeline_for(root, p, cal, start_step=start_step, overwrite_from_step=overwrite_from_step)
                results.append((p, st))
            except Exception:
                results.append((p, {"error": True}))
    return results
