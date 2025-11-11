import argparse
from pathlib import Path
import numpy as np
import tifffile
import napari
import json as _json

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

def open_napari_2d(
    image_path: Path,
    scale_2d,
    dapi_idx: int = 0,
    gfap_idx: int = 1,
    gfap_2d_projection: Path | None = None,
    gfap_2d_mask: Path | None = None,
    final_2d_mask: Path | None = None,
    skeleton_2d_mask: Path | None = None,
    rings_json: Path | None = None
):
    """
    Visualizador especializado para proyecciones 2D.
    Solo carga las capas 2D para evitar confusión.
    """
    
    # Cargar imagen original (para referencia de canales)
    image, axes = load_image_any(image_path)
    image = reorder_to_zcyx(image, axes)
    
    # Crear proyección máxima de la imagen original para contexto
    if image.ndim == 4:  # (Z, C, Y, X)
        # Extraer solo los canales DAPI y GFAP y proyectar en Z
        dapi_channel = image[:, dapi_idx, :, :]  # (Z, Y, X)
        gfap_channel = image[:, gfap_idx, :, :]  # (Z, Y, X)
        
        # Proyección máxima en Z para cada canal
        dapi_proj = np.max(dapi_channel, axis=0)  # (Y, X)
        gfap_proj = np.max(gfap_channel, axis=0)  # (Y, X)
        
        # Combinar en (C, Y, X) manteniendo canales separados
        image_proj = np.stack([dapi_proj, gfap_proj], axis=0)  # (2, Y, X)
    else:
        image_proj = image
    
    ch_names = ["DAPI (Proj)", "GFAP (Proj)"]

    # Determinar si las máscaras incluyen conexiones
    connection_info = ""
    if skeleton_2d_mask and skeleton_2d_mask.exists():
        # Buscar archivo de parámetros para verificar conexión
        parent_dir = skeleton_2d_mask.parent
        calibration_file = parent_dir.parent / "calibration.json"
        if calibration_file.exists():
            try:
                cal_data = _json.loads(calibration_file.read_text())
                if cal_data.get('CONNECT_SKELETON_FRAGMENTS', False):
                    connection_radius = cal_data.get('CONNECTION_RADIUS_UM', 0.0)
                    connection_info = f" [Fragmentos Conectados: {connection_radius:.1f}µm]"
            except:
                pass

    v = napari.Viewer(title=f"{image_path.name} - Análisis 2D{connection_info}")
    
    # Agregar imagen proyectada como referencia
    v.add_image(
        image_proj,
        channel_axis=0,
        name=ch_names,
        colormap=["blue", "green"],
        scale=scale_2d,
        opacity=0.7,
        visible=False  # Inicialmente oculta para no saturar
    )
    
    # --- Cargar Proyección GFAP 2D como imagen principal ---
    if gfap_2d_projection and gfap_2d_projection.exists():
        try:
            gfap_2d = tifffile.imread(str(gfap_2d_projection))
            v.add_image(
                gfap_2d, 
                name="GFAP Proyección 2D", 
                scale=scale_2d, 
                colormap='green',
                visible=True
            )
        except Exception as e:
            print(f"Advertencia: no se pudo cargar proyección GFAP 2D: {e}")

    # --- Cargar Máscaras 2D ---
    mask_2d_paths = {
        "GFAP Filtrado (2D)": gfap_2d_mask,
        "Astrocitos Finales (2D)": final_2d_mask,
        "Skeleton Labels (2D)": skeleton_2d_mask
    }
    
    # Detectar si hay conexión de fragmentos activa
    fragments_connected = False
    connection_radius = 0.0
    parent_dir = None
    if skeleton_2d_mask and skeleton_2d_mask.exists():
        parent_dir = skeleton_2d_mask.parent
        calibration_file = parent_dir.parent / "calibration.json"
        if calibration_file.exists():
            try:
                cal_data = _json.loads(calibration_file.read_text())
                fragments_connected = cal_data.get('CONNECT_SKELETON_FRAGMENTS', False)
                connection_radius = cal_data.get('CONNECTION_RADIUS_UM', 0.0)
            except:
                pass
    
    for name, path in mask_2d_paths.items():
        if path and path.exists():
            try:
                mask_2d = tifffile.imread(str(path))
                # Para máscaras 2D, verificar si son labels o binarias
                is_labels = len(np.unique(mask_2d)) > 2
                
                # Modificar nombre si incluye conexiones
                display_name = name
                if name == "Skeleton Labels (2D)" and fragments_connected:
                    display_name = f"🔗 Skeleton Conectado (r={connection_radius:.1f}µm)"
                
                if is_labels:
                    v.add_labels(
                        mask_2d, 
                        name=display_name, 
                        scale=scale_2d,
                        visible=(name == "Astrocitos Finales (2D)" or name == "Skeleton Labels (2D)")
                    )
                else:
                    # Convertir máscara binaria a imagen
                    v.add_image(
                        mask_2d.astype(np.float32), 
                        name=display_name, 
                        scale=scale_2d,
                        colormap='red',
                        opacity=0.5,
                        visible=False
                    )
            except Exception as e:
                print(f"Advertencia: no se pudo cargar {name}: {e}")

    # --- Cargar máscara original (sin conectar) para comparación ---
    if fragments_connected and parent_dir:
        # Buscar archivo original de skeleton (sin conexión)
        original_skel_candidates = [
            parent_dir / "05_skeleton_labels_2d_original.tif",
            parent_dir / "skeletons" / "skeleton_labels_original.tif"
        ]
        
        for orig_path in original_skel_candidates:
            if orig_path.exists():
                try:
                    orig_skel = tifffile.imread(str(orig_path))
                    v.add_labels(
                        orig_skel,
                        name="Skeleton Original (sin conectar)",
                        scale=scale_2d,
                        visible=False,
                        opacity=0.7
                    )
                    print(f"Cargada máscara original para comparación: {orig_path.name}")
                    break
                except Exception as e:
                    print(f"Advertencia: no se pudo cargar skeleton original: {e}")

    # --- Cargar Anillos de Sholl 2D (si existen) ---
    if rings_json and rings_json.exists():
        try:
            data = _json.loads(rings_json.read_text())
            all_ring_polygons = []
            
            # Para 2D, usamos solo Y,X de la escala
            try:
                scale_y, scale_x = scale_2d
            except Exception:
                scale_y = scale_x = 1.0
            
            # Procesar anillos para 2D (solo en el plano Y-X)
            for label, info in data.items():
                center_um = info.get("centroid_um")
                radii_um_list = info.get("radii_um", [])
                if center_um is None or not radii_um_list:
                    continue
                
                center_um = np.asarray(center_um, dtype=float)
                
                # Para 2D, usar solo Y,X del centroide (orden Z,Y,X -> tomar Y,X)
                if len(center_um) >= 3:
                    center_y_um, center_x_um = center_um[1], center_um[2]  # Y, X
                elif len(center_um) >= 2:
                    center_y_um, center_x_um = center_um[0], center_um[1]  # Y, X
                else:
                    continue
                
                for r in radii_um_list:
                    try:
                        r = float(r)
                    except Exception:
                        continue
                    if r <= 0:
                        continue
                    
                    # Convertir centroide de µm a voxeles
                    center_y_vox = center_y_um / scale_y if scale_y > 0 else center_y_um
                    center_x_vox = center_x_um / scale_x if scale_x > 0 else center_x_um
                    
                    # Convertir radio de µm a voxeles
                    ry = r / scale_y if scale_y > 0 else r
                    rx = r / scale_x if scale_x > 0 else r
                    
                    # Crear polígono circular en 2D
                    npts = 64
                    thetas = np.linspace(0, 2*np.pi, npts, endpoint=False)
                    ys = center_y_vox + ry * np.sin(thetas)
                    xs = center_x_vox + rx * np.cos(thetas)
                    poly = np.vstack([ys, xs]).T  # shape (npts, 2) como (y,x)
                    all_ring_polygons.append(poly)
            
            if all_ring_polygons:
                try:
                    v.add_shapes(
                        all_ring_polygons,
                        shape_type='polygon',
                        name='Anillos Sholl (2D)',
                        face_color='transparent',
                        edge_color='#FFFF00',
                        edge_width=2.0,
                        opacity=0.8,
                        scale=scale_2d
                    )
                except Exception as e:
                    print(f"Advertencia: no se pudieron agregar anillos 2D: {e}")
                    
        except Exception as e:
            print(f"Advertencia: no se pudieron procesar anillos de Sholl 2D: {e}")

    napari.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualizador Napari especializado para análisis 2D")
    parser.add_argument("--path", required=True, help="Ruta al archivo de imagen original")
    parser.add_argument("--y", type=float, required=True, help="Tamaño de vóxel Y en µm")
    parser.add_argument("--x", type=float, required=False, help="Tamaño de vóxel X en µm")
    parser.add_argument("--dapi_idx", type=int, default=0, help="Índice del canal DAPI")
    parser.add_argument("--gfap_idx", type=int, default=1, help="Índice del canal GFAP")
    
    # Archivos específicos 2D
    parser.add_argument("--gfap_proj", type=str, required=False, help="Ruta a gfap_projection_2d.tif")
    parser.add_argument("--gfap_2d", type=str, required=False, help="Ruta a 03_gfap_filtered_mask_2d.tif")
    parser.add_argument("--final_2d", type=str, required=False, help="Ruta a 04_final_astrocytes_mask_2d.tif")
    parser.add_argument("--skeleton_2d", type=str, required=False, help="Ruta a 05_skeleton_labels_2d.tif")
    parser.add_argument("--rings", type=str, required=False, help="Ruta a JSON con anillos de Sholl")
    
    args = parser.parse_args()

    p = Path(args.path)
    x_um = args.x if args.x is not None else args.y
    scale_2d = (args.y, x_um)  # Solo Y,X para 2D
    
    open_napari_2d(
        p,
        scale_2d,
        dapi_idx=args.dapi_idx,
        gfap_idx=args.gfap_idx,
        gfap_2d_projection=Path(args.gfap_proj) if args.gfap_proj else None,
        gfap_2d_mask=Path(args.gfap_2d) if args.gfap_2d else None,
        final_2d_mask=Path(args.final_2d) if args.final_2d else None,
        skeleton_2d_mask=Path(args.skeleton_2d) if args.skeleton_2d else None,
        rings_json=Path(args.rings) if args.rings else None,
    )