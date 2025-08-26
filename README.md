# Análisis Morfológico 3D de Astrocitos a partir de Imágenes Confocales

**Proyecto:** Reconstrucción, segmentación y análisis cuantitativo de astrocitos (GFAP), otros marcadores y núcleos (DAPI) a partir de imágenes `.lif`.
**Tecnologías:** Python, Napari, Cellpose, Scikit-image, Scipy, Skan

## 1. Objetivo del Proyecto

El objetivo es desarrollar un flujo de trabajo para segmentar y analizar la morfología tridimensional de astrocitos en imágenes de microscopía confocal multi-canal. El análisis se centra en la identificación de astrocitos completos para cuantificar su complejidad estructural.

## 2. Metodología Propuesta

El flujo de trabajo se divide en dos fases principales, implementadas en cuadernos de Jupyter:

**Fase 1: Identificación de Astrocitos**

1.  **Segmentación de Núcleos**: Se utiliza el modelo de deep learning **Cellpose** sobre el canal DAPI para obtener una máscara de todos los núcleos potenciales.
2.  **Filtrado por Co-localización**: Para cada núcleo, se mide la intensidad promedio de la señal GFAP en un "anillo" peri-nuclear. Solo se conservan los núcleos que superan un umbral de intensidad, identificándolos como núcleos de astrocitos.
3.  **Limpieza de Máscaras**: Se aplican filtros de pre-procesamiento (umbral de Otsu) y post-procesamiento (tamaño mínimo) para asegurar la robustez de la segmentación.

**Fase 2: Análisis Morfológico Avanzado**

1.  **Identificación del Astrocito Completo**: Para cada núcleo de astrocito identificado, se localiza la estructura GFAP conectada a él para obtener la máscara completa de la célula.
2.  **Skeletonización 3D**: Los procesos de cada astrocito se reducen a un esqueleto de 1 vóxel de grosor para analizar su topología.
3.  **Análisis de Ramificaciones**: Se utiliza la librería `skan` sobre el esqueleto para cuantificar la **longitud total de los procesos** y el **número de ramas**.
4.  **Análisis de Sholl**: Se mide la complejidad radial de las ramificaciones contando las intersecciones con esferas concéntricas.
5.  **Cálculo del Territorio Celular**: Se determina el volumen de la envolvente convexa (`convex hull`) para medir el espacio total de influencia de la célula.
6.  **Exportación de Datos**: Todas las métricas morfológicas se guardan en formato `.csv` para su análisis estadístico.

## 3. Estructura del Repositorio

- **/data**: Almacena los datos de imágenes. Excluida del repositorio vía `.gitignore`.
- **/notebooks**: Contiene los cuadernos de Jupyter con el flujo de trabajo detallado.
- **/results**: Almacena los resultados finales como tablas (`.csv`) y figuras.
- `requirements.txt`: Lista de dependencias de Python del proyecto.

## 4. Cómo Empezar

1.  **Clonar el repositorio.**
2.  **Crear y activar el entorno virtual (`venv`).**
3.  **Instalar las dependencias (`pip install -r requirements.txt`).**
4.  **Ejecutar los scripts** desde la carpeta `/src` para procesar las imágenes y visualizar los resultados.
