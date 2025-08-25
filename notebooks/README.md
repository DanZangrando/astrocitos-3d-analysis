# Carpeta de Cuadernos de Jupyter 🧠

Esta carpeta contiene los cuadernos de Jupyter que detallan el proceso de análisis, desde la segmentación hasta la cuantificación final de los astrocitos.

## Cuadernos

- **`01_segmentation_with_cellpose.ipynb`**:

  - **Objetivo**: Cargar las imágenes `.lif` y realizar la segmentación 3D de los astrocitos.
  - **Pasos**:
    1.  Carga de archivo `.lif` usando la librería `readlif`.
    2.  Visualización de los canales (GFAP, DAPI, etc.) en Napari.
    3.  Uso interactivo del plugin **`cellpose-napari`** para segmentar los astrocitos en 3D.
    4.  Corrección manual de la segmentación si es necesario.
    5.  Guardado de las máscaras de etiquetas (`label images`) en la carpeta `/data/processed`.

- **`02_quantification_and_sholl_analysis.ipynb`**:
  - **Objetivo**: Realizar el análisis cuantitativo y morfológico a partir de las máscaras de astrocitos segmentados.
  - **Pasos**:
    1.  Carga de las máscaras de segmentación.
    2.  Cálculo de propiedades de los objetos 3D (ej. **volumen**) con `scikit-image`.
    3.  Implementación del **análisis de Sholl** para cada astrocito segmentado.
    4.  Generación de visualizaciones 3D de los resultados en Napari.
    5.  Exportación de los datos cuantitativos a la carpeta `/results/tables`.
