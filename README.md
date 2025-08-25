# Análisis Morfológico 3D de Astrocitos a partir de Imágenes Confocales

**Proyecto:** Reconstrucción, segmentación y análisis cuantitativo de astrocitos (GFAP), otros marcadores y núcleos (DAPI) a partir de imágenes `.lif`.
**Tecnologías:** Python, Napari, Cellpose, Scikit-image, Readlif

## 1. Objetivo del Proyecto

El objetivo principal de este proyecto es desarrollar un flujo de trabajo para segmentar y analizar la morfología tridimensional de **astrocitos** en imágenes de microscopía confocal multi-canal. El análisis incluye la cuantificación de volúmenes celulares y un análisis de Sholl para caracterizar la complejidad de sus procesos.

## 2. Metodología Propuesta

1.  **Carga de Datos `.lif`**: Las imágenes en formato propietario de Leica (`.lif`) se cargarán en Python utilizando la librería `readlif` para acceder a los z-stacks y a los metadatos de la imagen.
2.  **Segmentación con Cellpose**: Se utilizará el modelo pre-entrenado de **Cellpose** para la segmentación de los núcleos (canal DAPI) y los cuerpos celulares de los **astrocitos** (canal GFAP). Este proceso se realizará de forma interactiva en `napari` a través del plugin `napari-cellpose`, permitiendo ajustes y correcciones manuales para garantizar la máxima precisión.
3.  **Cuantificación de Volúmenes**: A partir de las máscaras de segmentación 3D, se calculará el volumen de cada **astrocito**/núcleo individual utilizando `scikit-image` (`regionprops_table`).
4.  **Análisis de Sholl**: Para cuantificar la complejidad morfológica de los **astrocitos**, se implementará un análisis de Sholl. Este método consiste en trazar esferas concéntricas a partir del centroide de la célula y contar el número de intersecciones de los procesos astrocitarios con la superficie de cada esfera.
5.  **Visualización 3D**: Se utilizará `napari` para generar visualizaciones 3D interactivas, superponiendo la imagen original, las máscaras de segmentación y los resultados del análisis de Sholl.

## 3. Estructura del Repositorio

- **/data**: Almacena los datos de imágenes (`.lif` en `raw` y máscaras en `processed`). Excluida del repositorio vía `.gitignore`.
- **/notebooks**: Cuadernos de Jupyter para el análisis interactivo.
- **/results**: Contiene los resultados finales (tablas `.csv` y figuras `.png`).
- **/src**: Código Python reutilizable, como funciones para el análisis de Sholl.
- `requirements.txt`: Lista de dependencias de Python del proyecto.

## 4. Cómo Empezar

1.  **Clonar el repositorio.**
2.  **Crear y activar el entorno virtual (`venv`).**
3.  **Instalar las dependencias:**
    ```bash
    pip install -r requirements.txt
    ```
4.  **Instalar el plugin de Cellpose en Napari:**
    Dentro de Napari, ir a `Plugins > Install/Uninstall Package(s)` y buscar e instalar `cellpose-napari`.
5.  **Ejecutar los cuadernos en VS Code.**
