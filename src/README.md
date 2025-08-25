# Carpeta de Cuadernos y Documentación de Procesos 🧠

Esta sección documenta en detalle el flujo de trabajo implementado en los scripts de la carpeta `/src`.

## Flujo de Trabajo de Segmentación

El proceso se divide en tres scripts principales que siguen un orden lógico: **Explorar, Procesar y Validar**.

1.  **`explore_raw_image.py` (Explorar)**:

    - **Propósito**: Este script es el punto de partida. Carga una imagen original (`.tif`) y la muestra en Napari.
    - **Uso**: Se utiliza para la inspección visual inicial de los datos. Es fundamental para estimar parámetros clave antes del procesamiento, como el **`diameter`** de los núcleos, midiendo algunos ejemplos con la herramienta de regla de Napari.

2.  **`segment_and_save.py` (Procesar)**:

    - **Propósito**: Este script realiza el trabajo pesado. Carga una imagen, aplica la cadena de filtros y el modelo Cellpose para segmentar los núcleos, y guarda la máscara de etiquetas resultante en la carpeta `/data/processed`.
    - **Uso**: Se ejecuta después de haber explorado la imagen y tener una idea de los parámetros necesarios (especialmente el `diameter`).

3.  **`visualize_results.py` (Validar)**:
    - **Propósito**: Es un visor que carga una imagen original y su máscara procesada correspondiente.
    - **Uso**: Permite una inspección visual detallada para validar la calidad de la segmentación final, comparando la máscara limpia con los canales originales.

## Estrategia de Filtrado y Parámetros

Para asegurar una segmentación de alta calidad, se utiliza una cadena de tres filtros principales en el script `segment_and_save.py`:

### Filtro 1: Umbral de Intensidad (Pre-procesamiento)

- **Método**: Se utiliza la función `skimage.filters.threshold_otsu`.
- **Propósito**: Este filtro se aplica directamente sobre el canal DAPI antes de la segmentación para eliminar el ruido de fondo de baja señal, proporcionando una imagen más limpia a Cellpose.
- **Parámetros**:
  - `otsu_threshold`: **Este valor es calculado automáticamente por el algoritmo**, por lo que no requiere ajuste manual.

### Filtro 2: Parámetros de Segmentación de Cellpose

- **Método**: Modelo `CellposeModel` de la librería `cellpose`.
- **Propósito**: Realizar la segmentación de los núcleos en la imagen DAPI ya filtrada.
- **Parámetros Clave**:
  - `diameter`: **El parámetro más crítico para Cellpose**. Define el diámetro promedio (en píxeles) de los objetos a segmentar.
    - **Rango Apropiado**: Este valor depende totalmente de la magnificación de la imagen. Un rango típico puede ser de **`20` a `80` píxeles**. Se debe estimar usando el script `explore_raw_image.py`.
  - `cellprob_threshold`: Define el umbral de confianza.
    - **Rango Apropiado**: Dado que usamos un potente filtro de post-procesamiento por tamaño, podemos mantener este valor bajo (ej. **`0.0` a `0.5`**) para ser inclusivos.

### Filtro 3: Umbral de Tamaño/Volumen (Post-procesamiento)

- **Método**: Se utilizan las funciones de `skimage.measure.regionprops_table` y `pandas`.
- **Propósito**: Este es el filtro final contra los artefactos. Elimina los objetos que son demasiado pequeños para ser núcleos reales.
- **Parámetros Clave**:
  - `min_volume_threshold`: El volumen mínimo (en vóxeles) que un objeto debe tener para ser considerado un núcleo válido.
    - **Rango Apropiado**: Este valor es específico de los datos. La mejor estrategia es ejecutar la segmentación, visualizarla, y usar Napari para inspeccionar los volúmenes de los núcleos reales versus los artefactos para elegir un umbral que los separe eficazmente. Un valor inicial de **`500` a `1000` vóxeles** suele ser un buen punto de partida.
