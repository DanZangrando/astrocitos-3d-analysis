import napari
import tifffile
from pathlib import Path

# 1. Definir la ruta a la imagen
project_root = Path.cwd().parent
image_path = project_root / "data/raw/CTL/CTL 1-2/Inmuno 26-07-23.lif - CTL 1-2 a.tif"

# 2. Cargar la imagen
image_stack = tifffile.imread(image_path)
print(f"Imagen cargada. Dimensiones (Z, Canales, Y, X): {image_stack.shape}")

# Asumimos que el orden de los canales es: DAPI, GFAP, Microglia
# Si tu orden es diferente, ajusta esto.
channel_names = ["DAPI", "GFAP", "Microglia"]
channel_colors = ["blue", "green", "red"]

# 3. Iniciar el visualizador de Napari
viewer = napari.Viewer()

# 4. Añadir la imagen al visualizador, especificando colores y nombres
viewer.add_image(
    image_stack,
    name=channel_names,
    channel_axis=1,  # Asumimos que los canales están en el segundo eje (índice 1)
    colormap=channel_colors,
    contrast_limits=[(0, 1000), (0, 1500), (0, 1200)] # Ajusta estos límites según tus datos
)

# Ajustar brillo y contraste inicial (opcional)
if len(viewer.layers) == 3:
    viewer.layers["DAPI"].contrast_limits = (0, image_stack[:, 0].max() * 0.3) # Ejemplo
    viewer.layers["GFAP"].contrast_limits = (0, image_stack[:, 1].max() * 0.5) # Ejemplo
    viewer.layers["Microglia"].contrast_limits = (0, image_stack[:, 2].max() * 0.4) # Ejemplo

# 5. INSTRUCCIÓN CLAVE: Mantener la ventana abierta
print("Iniciando la aplicación de Napari. Cierra la ventana para continuar.")
napari.run()

print("Ventana de Napari cerrada. El script ha finalizado.")