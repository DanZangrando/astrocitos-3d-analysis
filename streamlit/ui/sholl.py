from __future__ import annotations
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi


@dataclass
class ShollParams:
    min_radius_um: float
    max_radius_um: float
    step_um: float

    def radii(self) -> np.ndarray:
        if self.step_um <= 0:
            raise ValueError("step_um debe ser > 0")
        if self.max_radius_um <= 0:
            return np.array([], dtype=float)
        r = np.arange(float(self.min_radius_um), float(self.max_radius_um) + 1e-6, float(self.step_um), dtype=float)
        return r[r >= 0]


def _centroids_from_labels(lbl: np.ndarray, spacing: tuple[float, float, float]) -> dict[int, tuple[float, float, float]]:
    # centroid en unidades físicas (µm)
    z_um, y_um, x_um = spacing
    ids = [int(v) for v in np.unique(lbl) if v > 0]
    cz = {}
    for lab in ids:
        zz, yy, xx = np.where(lbl == lab)
        if zz.size == 0:
            continue
        cz_um = float(np.mean(zz)) * z_um
        cy_um = float(np.mean(yy)) * y_um
        cx_um = float(np.mean(xx)) * x_um
        cz[lab] = (cz_um, cy_um, cx_um)
    return cz


def _dist_um_grid(shape: tuple[int, int, int], center_um: tuple[float, float, float], spacing: tuple[float, float, float]) -> np.ndarray:
    z_um, y_um, x_um = spacing
    cz, cy, cx = center_um
    z = (np.arange(shape[0], dtype=np.float32) * z_um - cz)
    y = (np.arange(shape[1], dtype=np.float32) * y_um - cy)
    x = (np.arange(shape[2], dtype=np.float32) * x_um - cx)
    zz, yy, xx = np.meshgrid(z, y, x, indexing='ij', sparse=True)
    return np.sqrt(zz**2 + yy**2 + xx**2)


def compute_sholl(out_dir: Path, spacing_um: tuple[float, float, float], params: ShollParams, 
                   restrict_to_final: bool = True,
                   save_rings_json: bool = True) -> pd.DataFrame:
    """Compute Sholl intersections per nucleus label using skeleton labels.

    Saves CSV to out_dir / 'sholl.csv' and returns the DataFrame.
    Optionally writes a rings JSON for Napari at out_dir / 'sholl_rings.json'.
    """
    out_dir = Path(out_dir)
    skel_path = out_dir / "05_skeleton_labels.tif"
    nuc_path = out_dir / "02_cellpose_mask.tif"
    final_path = out_dir / "04_final_astrocytes_mask.tif"
    if not skel_path.exists() or not nuc_path.exists():
        raise FileNotFoundError("Faltan archivos requeridos: 05_skeleton_labels.tif y/o 02_cellpose_mask.tif")

    skel = tifffile.imread(skel_path)
    nucs = tifffile.imread(nuc_path)
    finals = tifffile.imread(final_path) if (restrict_to_final and final_path.exists()) else None

    # Limitar a labels presentes en skeleton (si está por célula) y finales
    labs = [int(v) for v in np.unique(nucs) if v > 0]
    if finals is not None:
        labs_final = set(int(v) for v in np.unique(finals) if v > 0)
        labs = [l for l in labs if l in labs_final]
    if len(labs) == 0:
        return pd.DataFrame()

    spacing = tuple(float(v) for v in spacing_um)
    centroids = _centroids_from_labels(nucs, spacing)

    # Pre-cálculos
    radii = params.radii()
    if radii.size == 0:
        return pd.DataFrame()
    half = float(params.step_um) / 2.0

    rows: list[dict] = []
    # Anillos para Napari (en px, sobre el corte Z más cercano al centro)
    rings_for_viewer: list[dict] = []

    for lab in labs:
        if lab not in centroids:
            continue
        center_um = centroids[lab]
        # Distancias en UM por vóxel desde el centro
        d_um = _dist_um_grid(skel.shape, center_um, spacing)
        # Voxeles del esqueleto de esta célula (si labels por célula) o todo skeleton
        # Preferimos voxels de skeleton que correspondan a la célula, si no entonces todo skeleton
        vox = (skel == lab) if (np.max(skel) > 1) else (skel > 0)
        if not np.any(vox):
            continue
        for r in radii:
            shell = (d_um >= (r - half)) & (d_um < (r + half)) & vox
            # Aproximación: # de componentes conectados en el cascarón como # intersecciones
            n_intersections = int(ndi.label(shell.astype(np.uint8), structure=np.ones((3,3,3), dtype=np.uint8))[1]) if np.any(shell) else 0
            rows.append({
                "label": int(lab),
                "radius_um": float(r),
                "intersections": n_intersections,
            })
        # Para Napari: círculos 2D sobre el plano Z más cercano al centro
        z_um, y_um, x_um = spacing
        z_idx = int(round(center_um[0] / z_um))
        cy_px = float(center_um[1] / y_um)
        cx_px = float(center_um[2] / x_um)
        rings_for_viewer.append({
            "z_index": int(np.clip(z_idx, 0, skel.shape[0]-1)),
            "y_px": float(cy_px),
            "x_px": float(cx_px),
            "radii_px": [float(r / y_um) for r in radii.tolist()],
            "label": int(lab),
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(out_dir / "sholl.csv", index=False)
        # Solo guardar junto al preparado; no duplicar en results/tables
        if save_rings_json and rings_for_viewer:
            (out_dir / "sholl_rings.json").write_text(json.dumps({"items": rings_for_viewer}, indent=2))
    return df

