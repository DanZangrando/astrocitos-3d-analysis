import tifffile
from pathlib import Path
from cellpose import models
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table
from skimage.filters import threshold_otsu # Importamos la función de Otsu

# --- 1. Definir Rutas ---
project_root = Path.cwd().parent
image_path = project_root / "data/raw/CTL/CTL 1-2/Inmuno 26-07-23.lif - CTL 1-2 a.tif"

processed_dir = project_root / "data/processed"
processed_dir.mkdir(exist_ok=True)
output_mask_path = processed_dir / f"{image_path.stem}_mask_cleaned.tif"

print(f"Procesando imagen: {image_path.name}")

# --- 2. Cargar y Preparar Imagen ---
image_stack = tifffile.imread(image_path)
dapi_channel = image_stack[:, 0, :, :] # Asumimos DAPI en el canal 0

# --- 3. Pre-procesamiento: Limpieza de Ruido con Umbral de Intensidad ---
print("Aplicando umbral de intensidad para limpiar la señal DAPI...")

# Calculamos el umbral óptimo usando el método de Otsu
otsu_threshold = threshold_otsu(dapi_channel)

# Creamos una máscara binaria: los vóxeles por debajo del umbral se hacen 0 (fondo).
# Multiplicamos esta máscara por la imagen original para "apagar" el ruido.
dapi_channel_cleaned = np.where(dapi_channel > otsu_threshold, dapi_channel, 0)

print(f"Umbral de Otsu calculado: {otsu_threshold:.2f}. Señal DAPI limpiada.")

# --- 4. Ejecutar Cellpose sobre la Imagen Limpia ---
# Ahora le pasamos la imagen 'dapi_channel_cleaned' a Cellpose.
model = models.CellposeModel(gpu=True)
print("Ejecutando segmentación con Cellpose sobre la imagen limpia...")

masks, flows, styles = model.eval(
    dapi_channel_cleaned, # Usamos la imagen filtrada por intensidad
    diameter=30,
    z_axis=0,
    do_3D=True,
    cellprob_threshold=0.0
)
print(f"Segmentación inicial completada. Se encontraron {masks.max()} objetos.")

# --- 5. Post-procesamiento para Limpiar Máscaras por Tamaño ---
print("Iniciando post-procesamiento para eliminar artefactos por tamaño...")

props = regionprops_table(masks, properties=('label', 'area'))
props_df = pd.DataFrame(props)
min_volume_threshold = 500
labels_to_keep = props_df[props_df['area'] >= min_volume_threshold]['label']
cleaned_masks = np.where(np.isin(masks, labels_to_keep), masks, 0)

print(f"Post-procesamiento completado. Quedan {len(labels_to_keep)} núcleos después del filtrado.")

# --- 6. Guardar las Máscaras LIMPIAS ---
tifffile.imwrite(output_mask_path, cleaned_masks.astype(np.uint16))

print(f"¡Resultados limpios guardados exitosamente en: {output_mask_path}!")