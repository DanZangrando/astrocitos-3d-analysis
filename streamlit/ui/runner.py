from __future__ import annotations
from pathlib import Path
from typing import Iterable, Literal
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import tifffile
import streamlit as st # Para mostrar progreso en el recalculado

# --- IMPORTAR EL NUEVO PIPELINE UNIFICADO ---
from . import pipeline
# -------------------------------------------

Step = Literal["01", "02", "03", "04"]


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
    """
    Ejecuta el pipeline unificado desde un paso dado para una imagen.
    """
    out_dir = _ensure_out_dir(root, img_path)
    
    def get_path(step_id: str) -> Path:
        name_map = {
            "01": "01_otsu_mask.tif",
            "02": "02_cellpose_mask.tif",
            "03": "03_gfap_filtered_mask.tif",
            "03_metrics": "03_nucleus_metrics.csv",
            "04": "04_final_astrocytes_mask.tif",
            "skeleton_2d": "05_skeleton_labels_2d.tif",
            "skeleton_summary": "skeletons/summary.csv",
            "sholl_curves": "sholl_2d.csv",
            "sholl_summary": "sholl_summary.csv"
        }
        return out_dir / name_map[step_id]

    def need(step: Step) -> bool:
        order = {"01":1, "02":2, "03":3, "04":4}
        
        # Para paso 03: verificar tanto la máscara como las métricas
        if step == "03":
            mask_exists = get_path("03").exists()
            metrics_exists = get_path("03_metrics").exists()
            if not mask_exists or not metrics_exists:
                return order[step] >= order[start_step]
            if overwrite_from_step:
                return order[step] >= order[start_step]
            return False
        
        # Para paso 04: verificar máscara, skeleton summary Y sholl summary
        if step == "04":
            mask_exists = get_path("04").exists()
            skeleton_exists = get_path("skeleton_summary").exists()
            sholl_exists = get_path("sholl_summary").exists()
            
            if not mask_exists or not skeleton_exists or not sholl_exists:
                return order[step] >= order[start_step]
            if overwrite_from_step:
                return order[step] >= order[start_step]
            return False
        
        # Para otros pasos: comportamiento original
        if not get_path(step).exists():
            return order[step] >= order[start_step]
        if overwrite_from_step:
            return order[step] >= order[start_step]
        return False

    # --- Cargar Imagen y Canales ---
    arr, axes = pipeline.load_image_any(img_path)
    vol = pipeline.reorder_to_zcyx(arr, axes)  # Z,C,Y,X
    nC = vol.shape[1]
    
    dapi_idx = int(cal.get("DAPI_CHANNEL_INDEX", 0))
    gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1 if nC > 1 else 0))
    
    dapi = vol[:, dapi_idx, :, :]
    gfap = vol[:, gfap_idx, :, :]
    
    # --- Parámetros ---
    nucleus_diam = int(cal.get("NUCLEUS_DIAMETER", 30))
    use_gpu = bool(cal.get("CELLPOSE_USE_GPU", True))

    # --- 01 Otsu ---
    otsu_mask = None
    if need("01") or (need("03") and not get_path("01").exists()): # Paso 03 necesita Otsu
        otsu_mask, _ = pipeline.run_otsu_and_save(dapi, out_dir)
    
    # --- 02 Cellpose ---
    if need("02"):
        if otsu_mask is None and get_path("01").exists():
            otsu_mask = tifffile.imread(get_path("01")).astype(np.uint8)
        
        dapi_in = dapi
        if otsu_mask is not None:
            dapi_in = np.where(otsu_mask > 0, dapi, 0)
        
        pipeline.run_cellpose_and_save(dapi_in, out_dir, nucleus_diameter=nucleus_diam, use_gpu=use_gpu)

    # --- 03 Filtrado (GFAP/Tamaño) ---
    if need("03"):
        if not get_path("02").exists():
            raise FileNotFoundError(f"Falta 02_cellpose_mask.tif para {img_path.stem}")
        if otsu_mask is None:
            if not get_path("01").exists():
                print(f"Generando Otsu temporal para filtrado en {img_path.stem}")
                otsu_mask, _ = pipeline.run_otsu_and_save(dapi, out_dir)
            else:
                otsu_mask = tifffile.imread(get_path("01")).astype(np.uint8)
        
        cp_masks = tifffile.imread(get_path("02"))
        
        # 03a - Filtro GFAP (relativo)
        gfap_filtered_mask, _ = pipeline.run_filter_and_save(
            cp_masks, gfap, otsu_mask, cal, out_dir
        )
        
        # 03b - Filtro por Tamaño
        pipeline.run_size_filter_and_save(gfap_filtered_mask, cal, out_dir)

    # --- 04 Esqueletización + Sholl (Unificado 2D) ---
    if need("04"):
        if not get_path("04").exists():
             raise FileNotFoundError(f"Falta 04_final_astrocytes_mask.tif para {img_path.stem}")
        
        # Pipeline unificado: esqueletización + Sholl en un solo paso
        pipeline.run_unified_2d_skeleton_and_sholl(
            img_path=img_path,
            mask_path=get_path("04"),
            out_dir=out_dir,
            cal=cal
        )

    return {
        "otsu": get_path("01").exists(),
        "cellpose": get_path("02").exists(),
        "gfap_filtered": get_path("03").exists(),
        "nucleus_metrics": get_path("03_metrics").exists(),
        "final_mask": get_path("04").exists(),
        "skeleton_2d": get_path("skeleton_2d").exists(),
        "skeleton_summary": get_path("skeleton_summary").exists(),
        "sholl_curves": get_path("sholl_curves").exists(),
        "sholl_summary": get_path("sholl_summary").exists(),
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

    if not targets:
        return []

    # --- Ejecución ---
    results: list[tuple[Path, dict]] = []
    
    # Crear una barra de progreso en Streamlit
    progress_bar = st.progress(0.0)
    status_text = st.empty()
    
    for i, p in enumerate(targets):
        status_text.text(f"Procesando {i+1}/{len(targets)}: {p.stem}...")
        try:
            st_res = run_pipeline_for(root, p, cal, start_step=start_step, overwrite_from_step=overwrite_from_step)
            results.append((p, st_res))
        except Exception as e:
            print(f"ERROR en {p.stem}: {e}")
            results.append((p, {"error": str(e)}))
        
        progress_bar.progress((i + 1) / len(targets))
    
    status_text.text(f"¡Procesamiento completado para {len(targets)} preparados!")
    return results