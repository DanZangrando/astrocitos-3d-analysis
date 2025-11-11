import argparse
from pathlib import Path
import numpy as np
import tifffile
import napari
import json as _json

def reorder_to_zcyx(arr, axes: str | None):
    """Reordena un array con ejes reportados por tifffile.series a (Z, C, Y, X)."""
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
        raise ValueError(f"Ejes insuficientes tras seleccionar T=0: {ax_list} (se requieren Z,C,Y,X)")
    src_order = [ax_list.index(a) for a in needed]
    return np.transpose(arr, axes=src_order)


def load_image_any(path: Path):
    """Carga .tif/.tiff con tifffile y devuelve (image, axes str or None)."""
    suffix = path.suffix.lower()
    if suffix in ('.tif', '.tiff'):
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            axes = getattr(series, 'axes', None)
            arr = series.asarray()
        return arr, axes
    else:
        raise ValueError(f"Extensión no soportada: {suffix}")


def open_napari_3d(
    image_path: Path, 
    scale, 
    dapi_idx: int = 0,
    gfap_idx: int = 1,
    otsu_path: Path | None = None, 
    cellpose_path: Path | None = None, 
    gfap_path: Path | None = None, 
    final_path: Path | None = None, 
    skeleton_path: Path | None = None
):
    image, axes = load_image_any(image_path)
    image = reorder_to_zcyx(image, axes)

    # Extraer solo los canales DAPI y GFAP
    dapi_channel = image[:, dapi_idx:dapi_idx+1, :, :]
    gfap_channel = image[:, gfap_idx:gfap_idx+1, :, :]
    image_filtered = np.concatenate([dapi_channel, gfap_channel], axis=1)
    
    ch_names = ["DAPI (3D)", "GFAP (3D)"]

    v = napari.Viewer(title=f"{image_path.name} - Análisis 3D")
    v.add_image(
        image_filtered,
        channel_axis=1,
        name=ch_names,
        colormap=["blue", "green"],
        scale=scale,
    )

    # --- Cargar Máscaras 3D (con escala) ---
    mask_paths = {
        "01 - Máscara Otsu (3D)": otsu_path,
        "02 - Máscara Cellpose (3D)": cellpose_path,
        "03 - Filtro GFAP (3D)": gfap_path,
        "04 - Astrocitos Finales (3D)": final_path,
        "05 - Skeleton Labels (3D)": skeleton_path
    }
    
    for name, path in mask_paths.items():
        if path is not None and Path(path).exists():
            try:
                mask = tifffile.imread(str(path))
                v.add_labels(mask, name=name, scale=scale, visible=(name == "04 - Astrocitos Finales (3D)"))
            except Exception as e:
                print(f"Advertencia: no se pudo cargar {name} {path}: {e}")

    napari.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="Ruta al archivo de imagen")
    parser.add_argument("--z", type=float, required=True, help="Tamaño de vóxel Z en µm")
    parser.add_argument("--y", type=float, required=True, help="Tamaño de vóxel Y en µm")
    parser.add_argument("--x", type=float, required=False, help="Tamao de vóxel X en µm")
    parser.add_argument("--dapi_idx", type=int, default=0, help="Índice del canal DAPI")
    parser.add_argument("--gfap_idx", type=int, default=1, help="Índice del canal GFAP")
    parser.add_argument("--otsu", type=str, required=False, help="Ruta a 01_otsu_mask.tif")
    parser.add_argument("--cellpose", type=str, required=False, help="Ruta a 02_cellpose_mask.tif")
    parser.add_argument("--gfap", type=str, required=False, help="Ruta a 03_gfap_filtered_mask.tif")
    parser.add_argument("--final", type=str, required=False, help="Ruta a 04_final_astrocytes_mask.tif")
    parser.add_argument("--skeleton", type=str, required=False, help="Ruta a 05_skeleton_labels.tif")
    args = parser.parse_args()

    p = Path(args.path)
    x_um = args.x if args.x is not None else args.y
    scale = (args.z, args.y, x_um)
    
    open_napari_3d(
        p,
        scale,
        dapi_idx=args.dapi_idx,
        gfap_idx=args.gfap_idx,
        otsu_path=Path(args.otsu) if args.otsu else None,
        cellpose_path=Path(args.cellpose) if args.cellpose else None,
        gfap_path=Path(args.gfap) if args.gfap else None,
        final_path=Path(args.final) if args.final else None,
        skeleton_path=Path(args.skeleton) if args.skeleton else None,
    )