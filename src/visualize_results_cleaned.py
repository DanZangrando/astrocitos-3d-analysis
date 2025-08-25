import napari
import tifffile
from pathlib import Path

# --- 1. Define qué archivo quieres visualizar ---
project_root = Path.cwd().parent

# Elige el nombre base del archivo que procesaste
base_filename = "Inmuno 26-07-23.lif - CTL 1-2 a"

# Construimos las rutas a la imagen original y a su máscara LIMPIA
original_image_path = project_root / f"data/raw/CTL/CTL 1-2/{base_filename}.tif"
mask_path = project_root / f"data/processed/{base_filename}_mask_cleaned.tif"

# --- 2. Verificar que los archivos existan ---
if not original_image_path.exists() or not mask_path.exists():
    print("Error: No se encontró la imagen original o el archivo de máscara limpia.")
    print(f"Buscando imagen en: {original_image_path}")
    print(f"Buscando máscara en: {mask_path}")
else:
    # --- 3. Cargar Datos ---
    print("Cargando imagen original y máscara limpia...")
    original_image = tifffile.imread(original_image_path)
    cleaned_mask = tifffile.imread(mask_path)

    # --- 4. Visualizar en Napari ---
    viewer = napari.Viewer()

    # Añadimos la imagen original con sus colores
    channel_names = ["DAPI", "GFAP", "Microglia"]
    channel_colors = ["blue", "green", "red"]
    viewer.add_image(
        original_image,
        name=channel_names,
        channel_axis=1,
        colormap=channel_colors
    )

    # Añadimos la capa de máscaras limpias
    viewer.add_labels(cleaned_mask, name='Segmentación Limpia (por tamaño)')
    
    print("Iniciando Napari. Cierra la ventana para terminar.")
    napari.run()