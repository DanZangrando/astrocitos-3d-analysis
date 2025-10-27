import argparse
from pathlib import Path
import numpy as np
import tifffile
import napari
import readlif


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
    """Carga .tif/.tiff con tifffile o .lif con readlif y devuelve (image Z,C,Y,X, axes str or None)."""
    suffix = path.suffix.lower()
    if suffix in ('.tif', '.tiff'):
        with tifffile.TiffFile(str(path)) as tf:
            series = tf.series[0]
            axes = getattr(series, 'axes', None)
            arr = series.asarray()
        return arr, axes
    elif suffix == '.lif':
        rdr = readlif.Reader(str(path))
        # Tomamos la primera imagen/serie
        try:
            img0 = rdr.get_image(0)
            arr = np.asarray(img0)
        except Exception:
            # Fallback: algunas versiones exponen getSeries
            series_list = rdr.getSeries() if hasattr(rdr, 'getSeries') else []
            if series_list:
                s0 = series_list[0]
                # Intentamos obtener un volumen 4D aproximado
                try:
                    arr = np.asarray(s0)
                except Exception:
                    raise RuntimeError("No se pudo leer la primera serie del LIF")
            else:
                raise RuntimeError("No se encontraron series en el LIF")

        # Intento de deducir ejes para LIF: comúnmente (Z, Y, X) o (Z, Y, X, C)
        if arr.ndim == 3:
            # Z, Y, X
            arr = arr[:, None, :, :]
            axes = 'ZCYX'
        elif arr.ndim == 4:
            # Ubicar canal como el menor eje
            chan_axis = int(np.argmin(arr.shape))
            if chan_axis != 1:
                arr = np.moveaxis(arr, chan_axis, 1)
            axes = None  # ya lo dejamos en (Z, C, Y, X) heurísticamente
        else:
            # Otras combinaciones (p.ej., T presente): tomamos primer T si existe
            if arr.ndim == 5:
                # Asumimos (T, Z, Y, X, C) o similares: tomar T=0
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


def open_napari(image_path: Path, scale, otsu_path: Path | None = None, cellpose_path: Path | None = None, gfap_path: Path | None = None, final_path: Path | None = None):
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
    )
