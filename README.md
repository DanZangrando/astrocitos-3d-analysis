# Análisis Morfológico 3D de Astrocitos a partir de Imágenes Confocales

**Proyecto:** Reconstrucción, segmentación y análisis cuantitativo de astrocitos (GFAP), otros marcadores y núcleos (DAPI) a partir de imágenes `.lif`.
**Tecnologías:** Python, Napari, Cellpose, Scikit-image, Readlif

## 1. Objetivo del Proyecto

El objetivo es desarrollar un flujo de trabajo para segmentar y analizar la morfología tridimensional de astrocitos en imágenes de microscopía confocal multi-canal. El análisis se centra en la identificación de núcleos de astrocitos para su posterior caracterización morfológica.

## 2. Metodología Propuesta

El flujo de trabajo actual se basa en una estrategia de filtrado en múltiples etapas para aislar los núcleos de interés a partir del canal DAPI:

1.  **Carga de Datos**: Las imágenes `.lif` se cargan utilizando la librería `readlif`.
2.  **Pre-procesamiento (Filtrado de Intensidad)**: Se aplica un umbral de intensidad automático (**Otsu's Method**) al canal DAPI. Esto elimina el ruido de fondo de baja señal y proporciona una imagen más limpia para la segmentación.
3.  **Segmentación con Cellpose**: El modelo de deep learning **Cellpose** se ejecuta sobre el canal DAPI pre-procesado para generar un "mapa" de etiquetas (`label mask`), donde cada objeto detectado es un núcleo potencial.
4.  **Post-procesamiento (Filtrado por Tamaño)**: Se analizan las propiedades de cada objeto segmentado. Aquellos objetos cuyo volumen (en vóxeles) es inferior a un **umbral mínimo (`min_volume_threshold`)** son descartados, eliminando eficazmente los artefactos y segmentaciones de ruido.
5.  **Almacenamiento y Visualización**: Las máscaras de etiquetas finales y limpias se guardan como archivos `.tif` para su posterior análisis y visualización en Napari.

## 3. Estructura del Repositorio

- **/data**: Almacena los datos de imágenes. Excluida del repositorio vía `.gitignore`.
- **/notebooks**: Contiene la documentación detallada del flujo de trabajo, parámetros y guías de uso.
- **/results**: Almacena los resultados finales como tablas y figuras.
- **/src**: Contiene los scripts de Python para el procesamiento y la visualización.
- `requirements.txt`: Lista de dependencias de Python del proyecto.

## 4. Cómo Empezar

1.  **Clonar el repositorio.**
2.  **Crear y activar el entorno virtual (`venv`).**
3.  **Instalar las dependencias (`pip install -r requirements.txt`).**
4.  **Ejecutar los scripts** desde la carpeta `/src` para procesar las imágenes y visualizar los resultados.
