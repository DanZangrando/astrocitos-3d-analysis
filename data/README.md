# Carpeta de Datos 🔬

Esta carpeta almacena todos los datos de imágenes y resultados del proyecto.

## Subcarpetas

- `raw/`: Archivos de imagen originales del microscopio en formato `.tif`. La app actualmente trabaja en modo TIF‑only por estabilidad.

- `processed/`: Resultados generados por el pipeline, organizados por preparado. Incluye, por cada imagen:
  - `01_otsu_mask.tif`, `02_cellpose_mask.tif`, `03_gfap_microglia_filtered_mask.tif`, `04_final_astrocytes_mask.tif`, `05_skeleton_labels.tif`.
  - `skeletons/summary.csv` con métricas por célula.
  - `sholl.csv` y `sholl_rings.json` con resultados del análisis de Sholl y geometría de anillos.
  - `params.json` con los parámetros usados en ese preparado.

Además, se genera `data/processed/master_morphology_results.csv` como agregación global de métricas por preparado.

Nota: Por tamaño, los datos están excluidos de Git vía `.gitignore`. Mantener y respaldar localmente.
