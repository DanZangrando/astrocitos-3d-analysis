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


def open_napari(
    image_path: Path, 
    scale, 
    otsu_path: Path | None = None, 
    cellpose_path: Path | None = None, 
    gfap_path: Path | None = None, 
    final_path: Path | None = None, 
    skeleton_path: Path | None = None, 
    rings_json: Path | None = None
):
    image, axes = load_image_any(image_path)
    image = reorder_to_zcyx(image, axes)

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

    # --- Cargar Máscaras (con escala) ---
    mask_paths = {
        "01 - Máscara Otsu": otsu_path,
        "02 - Máscara Cellpose": cellpose_path,
        "03 - Filtro GFAP/Microglía": gfap_path,
        "04 - Astrocitos Finales": final_path,
        "05 - Skeleton (labels)": skeleton_path
    }
    
    for name, path in mask_paths.items():
        if path is not None and Path(path).exists():
            try:
                mask = tifffile.imread(str(path))
                v.add_labels(mask, name=name, scale=scale, visible=(name == "04 - Astrocitos Finales"))
            except Exception as e:
                print(f"Advertencia: no se pudo cargar {name} {path}: {e}")

    # --- Cargar Anillos de Sholl (convertir µm -> voxels y dibujar) ---
    if rings_json is not None and Path(rings_json).exists():
        try:
            data = _json.loads(Path(rings_json).read_text())

            all_ellipsoid_boxes = []
            all_ring_polygons = []

            # Asegurarse de que `scale` esté disponible: esperado (z_um, y_um, x_um)
            try:
                scale_z, scale_y, scale_x = scale
            except Exception:
                # Si scale no tiene 3 elementos, asumir 1.0
                scale_z = scale_y = scale_x = 1.0

            # El formato es {"label": {"centroid_um": [z,y,x], "radii_um": [r1, r2...]}}
            for label, info in data.items():
                center_um = info.get("centroid_um")
                radii_um_list = info.get("radii_um", [])
                if center_um is None or not radii_um_list:
                    continue
                center_um = np.asarray(center_um, dtype=float)

                for r in radii_um_list:
                    try:
                        r = float(r)
                    except Exception:
                        continue
                    if r <= 0:
                        continue

                    # Convertir centro y radios de µm a coordenadas en vóxeles (Z,Y,X)
                    # Pero detectamos si `center_um` ya está en voxeles: si dividir por scale
                    # produce coordenadas fuera de la imagen, preferimos tratarlas como voxels.
                    try:
                        cand_vox = center_um / np.array([scale_z, scale_y, scale_x], dtype=float)
                    except Exception:
                        cand_vox = center_um

                    # Imagen dims en voxeles
                    z_dim = image.shape[0]
                    y_dim = image.shape[2]
                    x_dim = image.shape[3]

                    def in_bounds(v):
                        return (v[0] >= 0 and v[0] < z_dim) and (v[1] >= 0 and v[1] < y_dim) and (v[2] >= 0 and v[2] < x_dim)

                    if in_bounds(cand_vox):
                        center_vox = cand_vox
                    else:
                        # If the raw center_um already looks like voxel coords, use it directly
                        if in_bounds(center_um):
                            center_vox = center_um.astype(float)
                        else:
                            # As a last resort, clamp the cand_vox into image bounds
                            center_vox = np.clip(cand_vox, [0,0,0], [z_dim-1, y_dim-1, x_dim-1])

                    # For ring polygons we need separate radii in Y and X voxels
                    ry = r / scale_y if scale_y > 0 else r
                    rx = r / scale_x if scale_x > 0 else r

                    # Build a circular polygon in the Y-X plane at constant Z=center_vox[0]
                    npts = 128
                    thetas = np.linspace(0, 2*np.pi, npts, endpoint=False)
                    ys = center_vox[1] + ry * np.sin(thetas)
                    xs = center_vox[2] + rx * np.cos(thetas)
                    zs = np.full_like(xs, fill_value=center_vox[0])
                    poly = np.vstack([zs, ys, xs]).T  # shape (npts, 3) as (z,y,x)
                    # Only keep polygons that are at least partially inside the image bounds
                    if np.any((ys >= 0) & (ys < y_dim) & (xs >= 0) & (xs < x_dim)):
                        all_ring_polygons.append(poly)

                    # Also keep a simple ellipse-like shape for compatibility with older napari
                    shape_data = np.asarray([center_vox.tolist(), [r/scale_z if scale_z>0 else r, r/scale_y if scale_y>0 else r, r/scale_x if scale_x>0 else r]], dtype=float)
                    if shape_data.shape == (2, 3):
                        all_ellipsoid_boxes.append(shape_data)

            # Prefer drawing ring polygons (concentric rings per astrocyte)
            if all_ring_polygons:
                try:
                    v.add_shapes(all_ring_polygons, shape_type='polygon', name='Sholl rings (polygons)', face_color='transparent', edge_color='#FFFF00', edge_width=1.0, opacity=0.8, scale=scale)
                except Exception as e:
                    print(f"Advertencia: add_shapes falló con polygons ({e}), intentando ellipses...")
                    # Fall back to ellipses if polygons are problematic
                    if all_ellipsoid_boxes:
                        try:
                            v.add_shapes(all_ellipsoid_boxes, shape_type='ellipse', name='Sholl rings (ellipses)', face_color='transparent', edge_color='#FFFF00', edge_width=1.0, opacity=0.6, scale=scale)
                        except Exception as ex:
                            print(f"Advertencia: add_shapes falló con ellipses también ({ex}), usando puntos centrales.")
                            # Final fallback: points at centers
                            try:
                                centers = np.array([s[0] for s in all_ellipsoid_boxes])
                                if centers.size:
                                    v.add_points(centers, name='Sholl centers (fallback)', size=5.0, face_color='#FFFF00', scale=scale)
                            except Exception as exc:
                                print(f"Advertencia: no se pudieron añadir puntos de fallback: {exc}")
                    else:
                        print("Advertencia: no hay datos de elipses para intentar fallback.")
            elif all_ellipsoid_boxes:
                # If no polygons available, try ellipses directly
                try:
                        v.add_shapes(all_ellipsoid_boxes, shape_type='ellipse', name='Sholl rings (ellipses)', face_color='transparent', edge_color='#FFFF00', edge_width=1.0, opacity=0.6, scale=scale)
                except Exception as e:
                    print(f"Advertencia: add_shapes falló con ellipses ({e}), usando puntos centrales.")
                    try:
                        centers = np.array([s[0] for s in all_ellipsoid_boxes])
                        if centers.size:
                            v.add_points(centers, name='Sholl centers (fallback)', size=5.0, face_color='#FFFF00', scale=scale)
                    except Exception as exc:
                        print(f"Advertencia: no se pudieron añadir puntos de fallback: {exc}")

        except Exception as e:
            print(f"Advertencia: no se pudieron dibujar anillos de Sholl 3D: {e}")
            
    napari.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--path", required=True, help="Ruta al archivo de imagen")
    parser.add_argument("--z", type=float, required=True, help="Tamaño de vóxel Z en µm")
    parser.add_argument("--y", type=float, required=True, help="Tamaño de vóxel Y en µm")
    parser.add_argument("--x", type=float, required=False, help="Tamao de vóxel X en µm")
    parser.add_argument("--otsu", type=str, required=False, help="Ruta a 01_otsu_mask.tif")
    parser.add_argument("--cellpose", type=str, required=False, help="Ruta a 02_cellpose_mask.tif")
    parser.add_argument("--gfap", type=str, required=False, help="Ruta a 03_gfap_microglia_filtered_mask.tif")
    parser.add_argument("--final", type=str, required=False, help="Ruta a 04_final_astrocytes_mask.tif")
    parser.add_argument("--skeleton", type=str, required=False, help="Ruta a 05_skeleton_labels.tif")
    parser.add_argument("--rings", type=str, required=False, help="Ruta a JSON con anillos de Sholl")
    args = parser.parse_args()

    p = Path(args.path)
    x_um = args.x if args.x is not None else args.y
    scale = (args.z, args.y, x_um)
    
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