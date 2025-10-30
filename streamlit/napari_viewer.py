import argparse
from pathlib import Path
import numpy as np
import tifffile
import napari


def reorder_to_zcyx(arr, axes: str | None):
    """Reordena un array con ejes reportados por tifffile.series a (Z, C, Y, X)."""
    if axes is None:
        if arr.ndim != 4:
            raise ValueError(f"Forma inesperada sin ejes: {arr.shape}")
        # Heurística: canal es el eje menor
        chan_axis = int(np.argmin(arr.shape))
        if chan_axis != 1:
            arr = np.moveaxis(arr, chan_axis, 1)
        return arr
    ax_list = list(axes)
    # Seleccionar T=0 si existe
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


def open_napari(image_path: Path, scale, otsu_path: Path | None = None, cellpose_path: Path | None = None, gfap_path: Path | None = None, final_path: Path | None = None, skeleton_path: Path | None = None, rings_json: Path | None = None):
    image, axes = load_image_any(image_path)
    image = reorder_to_zcyx(image, axes)

    # Nombres de canales típicos
    ch_names = None
    if image.shape[1] == 3:
        ch_names = ["DAPI", "GFAP", "Microglia"]
    elif image.shape[1] == 2:
        ch_names = ["Ch0", "Ch1"]

    v = napari.Viewer(title=f"{image_path.name} - Calibración {scale}")
    v.add_image(
        image,
        channel_axis=1,
        name=ch_names if ch_names else None,
        colormap=["blue","green","magenta","yellow","cyan"][:image.shape[1]],
        scale=scale,
    )

    # Agregar máscaras si fueron provistas
    if otsu_path is not None and Path(otsu_path).exists():
        try:
            otsu_mask = tifffile.imread(str(otsu_path))
            v.add_labels(otsu_mask, name="01 - Máscara Otsu", scale=scale, visible=False)
        except Exception as e:
            print(f"Advertencia: no se pudo cargar Otsu {otsu_path}: {e}")

    if cellpose_path is not None and Path(cellpose_path).exists():
        try:
            cellpose_mask = tifffile.imread(str(cellpose_path))
            v.add_labels(cellpose_mask, name="02 - Máscara Cellpose", scale=scale, visible=True)
        except Exception as e:
            print(f"Advertencia: no se pudo cargar Cellpose {cellpose_path}: {e}")

    if gfap_path is not None and Path(gfap_path).exists():
        try:
            gfap_mask = tifffile.imread(str(gfap_path))
            v.add_labels(gfap_mask, name="03 - Filtro GFAP/Microglía", scale=scale, visible=True)
        except Exception as e:
            print(f"Advertencia: no se pudo cargar GFAP filtrado {gfap_path}: {e}")

    if final_path is not None and Path(final_path).exists():
        try:
            final_mask = tifffile.imread(str(final_path))
            v.add_labels(final_mask, name="04 - Astrocitos Finales", scale=scale, visible=True)
        except Exception as e:
            print(f"Advertencia: no se pudo cargar Máscara Final {final_path}: {e}")

    if skeleton_path is not None and Path(skeleton_path).exists():
        try:
            skel = tifffile.imread(str(skeleton_path))
            v.add_labels(skel, name="05 - Skeleton (labels)", scale=scale, visible=True)
        except Exception as e:
            print(f"Advertencia: no se pudo cargar Skeleton {skeleton_path}: {e}")
    # Círculos de Sholl (Shapes 2D sobre planos Z dados)
    if rings_json is not None and Path(rings_json).exists():
        try:
            import json as _json
            data = _json.loads(Path(rings_json).read_text())
            items = data.get("items", [])
            shapes = []
            for it in items:
                z = int(it.get("z_index", 0))
                cy = float(it.get("y_px", 0.0))
                cx = float(it.get("x_px", 0.0))
                radii = it.get("radii_px", [])
                for r in radii:
                    r = float(r)
                    # polígono aproximando el círculo en el plano z
                    theta = np.linspace(0, 2*np.pi, 64, endpoint=True)
                    yy = cy + r * np.sin(theta)
                    xx = cx + r * np.cos(theta)
                    # Napari espera data como (N, D) con D=3 para 3D; fijamos z constante
                    poly = np.stack([np.full_like(yy, z), yy, xx], axis=1)
                    # Para 'path', repetir el primer punto al final para cerrar el anillo
                    poly = np.concatenate([poly, poly[:1]], axis=0)
                    shapes.append(poly)
            if shapes:
                # Usar 'path' en lugar de 'polygon' para evitar triangulación en 3D
                # y errores al activar la vista 3D de napari.
                v.add_shapes(
                    shapes,
                    shape_type='path',
                    name='Sholl rings',
                    edge_color='#FFFF00',  # hex color to avoid property lookup
                    edge_width=1,
                    face_color='transparent',
                    scale=scale,  # alinear con imagen
                    blending='translucent',
                    opacity=0.7,
                )
        except Exception as e:
            print(f"Advertencia: no se pudieron dibujar anillos de Sholl: {e}")
    napari.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="Ruta al archivo de imagen")
    parser.add_argument("--z", type=float, required=True, help="Tamaño de vóxel Z en µm")
    parser.add_argument("--y", type=float, required=True, help="Tamaño de vóxel Y en µm")
    parser.add_argument("--x", type=float, required=True, help="Tamaño de vóxel X en µm")
    parser.add_argument("--otsu", type=str, required=False, help="Ruta a 01_otsu_mask.tif")
    parser.add_argument("--cellpose", type=str, required=False, help="Ruta a 02_cellpose_mask.tif")
    parser.add_argument("--gfap", type=str, required=False, help="Ruta a 03_gfap_microglia_filtered_mask.tif")
    parser.add_argument("--final", type=str, required=False, help="Ruta a 04_final_astrocytes_mask.tif")
    parser.add_argument("--skeleton", type=str, required=False, help="Ruta a 05_skeleton_labels.tif")
    parser.add_argument("--rings", type=str, required=False, help="Ruta a JSON con anillos de Sholl")
    args = parser.parse_args()

    p = Path(args.path)
    scale = (args.z, args.y, args.x)
    open_napari(
        p,
        scale,
        Path(args.otsu) if args.otsu else None,
        Path(args.cellpose) if args.cellpose else None,
        Path(args.gfap) if args.gfap else None,
        Path(args.final) if args.final else None,
        Path(args.skeleton) if args.skeleton else None,
        Path(args.rings) if args.rings else None,
    )
