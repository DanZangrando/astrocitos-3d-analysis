import napari
import tifffile
from pathlib import Path

# --- 1. Define qué archivo RAW quieres visualizar ---
project_root = Path.cwd().parent

# Cambia este nombre para explorar diferentes imágenes
base_filename = "Inmuno 26-07-23.lif - CTL 1-2 a"
# Asegúrate de que la ruta coincida con tu estructura (ej. CTL o hip)
subfolder = "CTL/CTL 1-2" 

original_image_path = project_root / f"data/raw/{subfolder}/{base_filename}.tif"

# --- 2. Verificar que el archivo exista ---
if not original_image_path.exists():
    print(f"Error: No se encontró la imagen en: {original_image_path}")
else:
    # --- 3. Cargar Datos ---
    print(f"Cargando imagen original: {original_image_path.name}...")
    original_image = tifffile.imread(original_image_path)
    print(f"Dimensiones: {original_image.shape}")

    # --- 4. Visualizar en Napari ---
    viewer = napari.Viewer()

    # Añadimos la imagen con sus canales y colores
    channel_names = ["DAPI", "GFAP", "Microglia"]
    channel_colors = ["blue", "green", "red"]
    viewer.add_image(
        original_image,
        name=channel_names,
        channel_axis=1, # Asume que el eje 1 son los canales (Z, C, Y, X)
        colormap=channel_colors
    )
    
    print("Iniciando Napari para exploración. Cierra la ventana para terminar.")
    napari.run()