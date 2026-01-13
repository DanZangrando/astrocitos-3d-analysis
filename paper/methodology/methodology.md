
# Methodology


### Librerías Principales

-   **`numpy`**: Operaciones matriciales para el cálculo de sumas de intensidad y recorte de arrays 4D.
-   **`tifffile`**: Lectura y escritura de imágenes compatible con ImageJ Hyperstacks.

## 1. Calibración y Visualización de los Preparados

Este paso inicial tiene como objetivo definir la escala física del experimento (resolución en micrones por vóxel) y visualizar los datos crudos para verificar su integridad.

### Descripción del Proceso

El pipeline comienza con la carga de imágenes de microscopía confocal (formato `.tif` o `.tiff`). Debido al pre-procesamiento previo, los metadatos físicos originales no siempre son detectables automáticamente.

Por esta razón, se implementó un sistema de **Calibración Global Manual permitiendo al usuario ingresar las dimensiones físicas del vóxel (ancho, alto y profundidad en µm). Estos valores se guardan en un archivo de configuración (`calibration.json`) y se reutilizan en todos los pasos subsiguientes del análisis para garantizar la consistencia de las mediciones morfométricas.

### Librerías Principales

-   **`streamlit`**: Framework para la interfaz de usuario interactiva.
-   **`tifffile`**: Lectura y manejo de archivos.
-   **`json`**: Persistencia de la configuración.
-   **`napari`**: Visualización multidimensional.

### Funciones Clave

-   `read_calibration_from_tiff(path)`: Intenta recuperar automáticamente la escala física.
-   `save_override_global(z, y, x)`: Guarda la calibración manual.
-   `detect_group(path)`: Clasifica muestras en grupos (CTL/Hipoxia).

### Consideraciones Metodológicas

Calibración estándar utilizada:
-   **Resolución XY**: 0.3788 µm/pixel.
-   **Resolución Z**: 1.0071 µm/slice.

## 2. Umbral de Otsu y Segmentación Nuclear (Cellpose)

Este paso se enfoca en la identificación y segmentación de los núcleos celulares (DAPI).

### Descripción del Proceso

1.  **Umbral de Otsu (Global)**: Cálculo automático de umbral para separar fondo de señal biológica.
2.  **Segmentación con Cellpose (Instancia)**: Algoritmo de deep learning (`cyto` o `nuclei`) para separar núcleos individuales, ideal para clusters densos.

### Librerías Principales

-   **`skimage.filters.threshold_otsu`**: Cálculo de umbral global.
-   **`cellpose`**: Segmentación basada en redes neuronales.
-   **`tifffile`**: Manejo de máscaras.

### Funciones Clave

-   `run_otsu_and_save`: Genera máscara binaria de Otsu.
-   `run_cellpose_and_save`: Ejecuta inferencia de Cellpose.

### Parámetros Críticos

-   **`NUCLEUS_DIAMETER`**: Diámetro característico (e.g., 30 px).
-   **`CELLPOSE_USE_GPU`**: Aceleración por GPU.

## 3. Filtrado de Núcleos de Astrocitos (Colocalización GFAP)

Una vez segmentados todos los núcleos presentes en la muestra, es necesario discriminar cuáles corresponden a astrocitos y cuáles a otras células (neuronas, microglia, etc.).

### Descripción del Proceso

El criterio principal de clasificación es la **positividad a GFAP** en la perinuclea. Se asume que los núcleos rodeados por una "cáscara" o *shell* de señal GFAP intensa corresponden a astrocitos.

1.  **Construcción del Shell**: Para cada núcleo segmentado, se genera una región de interés tridimensional extendiendo sus bordes una distancia definida (`SHELL_RADIUS_UM`, típicamente 2.0 µm).
2.  **Cuantificación de Señal GFAP**: Se mide la intensidad media del canal de GFAP dentro de este *shell*, excluyendo la región nuclear propiamente dicha.
3.  **Filtrado Relativo (Robustez)**: En lugar de usar un umbral de intensidad absoluto (variable entre imágenes), se calcula el fondo local y la desviación estándar de la señal GFAP. Un núcleo se considera "Positivo" si la intensidad de medio de su shell supera al fondo en un número determinado de desviaciones estándar (`GFAP_STD_DEV_THRESHOLD`, e.g., > 3.0 SD).
4.  **Filtrado por Tamaño (Volumen)**: Se aplica un segundo criterio morfológico para descartar artefactos o núcleos fragmentados, exigiendo un volumen mínimo (`MIN_VOLUME_UM3`, e.g., > 75 µm³).

### Librerías Principales

-   **`skimage.measure.regionprops`**: Utilizado para iterar sobre cada núcleo y extraer sus coordenadas y propiedades geométricas.
-   **`numpy`**: Operaciones de enmascaramiento y cálculos estadísticos (media, desviación estándar) sobre los vóxeles de imagen.
-   **`pandas`**: Tabulación de métricas por núcleo, permitiendo el filtrado y análisis posterior.

### Funciones Clave

-   `pipeline.run_filter_and_save(...)`: Orquesta el proceso de medición de intensidad en shells y clasificación binaria (Candidato / No Candidato). Genera `03_gfap_filtered_mask.tif` con solo los núcleos candidatos.
-   `pipeline.run_size_filter_and_save(...)`: Aplica el filtro de volumen físico sobre los candidatos y genera la máscara final `04_final_astrocytes_mask.tif`.

### Parámetros Críticos

-   **`SHELL_RADIUS_UM`**: Grosor de la zona perinuclear a evaluar.
-   **`GFAP_STD_DEV_THRESHOLD`**: Rigurosidad del filtro de intensidad. Valores más altos aumentan la especificidad pero reducen la sensibilidad.
-   **`MIN_VOLUME_UM3`**: Criterio de exclusión para detritos o fragmentos de segmentación.

## 4. Esqueletización y Análisis Sholl 2D

Para caracterizar la complejidad morfológica de los astrocitos, se realiza un análisis topológico sobre una proyección bidimensional (2D) de la señal. Este enfoque maximiza la resolución lateral (XY) y simplifica la interpretación de las ramificaciones.

### Descripción del Proceso

1.  **Proyección 2D**: Se colapsa el stack 3D (10 planos) en una imagen 2D utilizando la proyección de máxima intensidad (`Max Projection`) sobre el eje Z. Esto conserva las estructuras más brillantes y continuas de las ramificaciones.
2.  **Definición de Territorios (Voronoi)**: Para evitar que las ramificaciones de astrocitos adyacentes se mezclen, se particiona el espacio utilizando un diagrama de Voronoi sembrado en los centroides de los núcleos. Se aplica una zona de exclusión (`TERRITORY_EXCLUSION_UM`) en las fronteras para evitar asignaciones ambiguas.
3.  **Esqueletización**: Dentro del territorio de cada astrocito, la señal GFAP segmentada se reduce a un "esqueleto" de 1 píxel de grosor, preservando la topología de la red de fibras.
4.  **Análisis de Sholl**: Se cuantifica la complejidad de la arborización contando el número de intersecciones del esqueleto con círculos concéntricos centrados en el soma celular, en pasos incrementales (`SHOLL_STEP_UM`).

### Librerías Principales

-   **`skimage.morphology.skeletonize`**: Algoritmo de adelgazamiento (thinning) para obtener el eje medial de las ramificaciones.
-   **`skan`**: Librería especializada para el análisis de esqueletos. Permite descomponer la red en ramas (`branches`), uniones (`junctions`) y puntos finales (`endpoints`), calculando longitudes y jerarquías.
-   **`scipy.spatial.Voronoi`**: Para la segmentación geométrica de los dominios celulares.

### Funciones Clave

-   `pipeline.run_unified_2d_skeleton_and_sholl(...)`: Función maestra que orquesta la proyección, partición de territorios, esqueletización y cuantificación de métricas en un solo paso eficiente.

### Parámetros Críticos

-   **`PROJECTION_2D_METHOD`**: Método de colapso Z (usualmente "max").
-   **`TERRITORY_EXCLUSION_UM`**: Separación física forzada entre células vecinas.
-   **`MAX_RADIUS_FROM_NUCLEUS_UM`**: Límite distal de análisis para evitar ruido lejano.
-   **`SHOLL_STEP_UM`**: Resolución espacial del análisis de Sholl.

## 5. Análisis Detallado por Preparado

Este módulo integra y visualiza todos los datos generados para una muestra individual, permitiendo el control de calidad y la exploración profunda de las métricas antes del análisis grupal.

### Descripción del Proceso

El sistema carga los archivos de resultados generados en las etapas previas para construir un perfil morfométrico completo del preparado:

1.  **Integración de Datos**: Se leen y correlacionan las máscaras de segmentación, tablas de métricas nucleares, datos de esqueletos y perfiles de Sholl.
2.  **Validación de Conteos**: Se reportan los números absolutos de núcleos detectados, candidatos validos GFAP+ y astrocitos finales con morfología resuelta, cuantificando la eficiencia del filtrado.
3.  **Visualización Estadística**: Se generan histogramas interactivos y matrices de correlación para evaluar la población celular (e.g., distribución de volúmenes, tortuosidad vs complejidad).
4.  **Inspección Visual (Napari)**: Se provee acceso directo a la visualización 3D y 2D de las capas procesadas (máscaras, esqueletos, anillos de Sholl) superpuestas a la imagen original, permitiendo verificar manualmente la precisión de la segmentación.

### Librerías Principales

-   **`pandas`**: Agregación y manipulación de tablas de datos CSV.
-   **`altair`**: Generación de gráficos interactivos (scatter plots, histogramas) integrados en la interfaz web.
-   **`napari`**: Renderizado avanzado para inspección visual cualitativa.

### Métricas Principales (Biomarcadores)

Se seleccionaron cuatro métricas clave para caracterizar la reactividad astrocitaria:

1.  **Índice de Ramificación (Ramification Index)**: Relación entre el número de ramas y uniones. Refleja la complejidad topológica.
2.  **Longitud Total del Esqueleto (Total Branch Length)**: Suma de la longitud de todos los procesos. Indica el volumen de exploración del astrocito.
3.  **Número de Terminaciones (Endpoints)**: Cantidad total de puntos finales de las ramas. Indica el grado de división terminal.
4.  **Radio Crítico (Sholl)**: Distancia desde el núcleo donde ocurre la máxima arborización. Indica la distribución espacial de los procesos.

Estas métricas constituyen el núcleo del análisis fenotípico presentado en este estudio.

## 6. Comparación entre Grupos y Análisis Multivariado

El objetivo final del pipeline es identificar las alteraciones fenotípicas inducidas por la condición experimental (Hipoxia) frente al control (CTL), asegurando el rigor estadístico.

### Descripción del Proceso

El módulo de comparación realiza una carga masiva de todos los datos procesados disponibles y ejecuta un flujo de análisis estadístico automatizado:

1.  **Agregación Robusta**: Para evitar la *pseudoreplicación* (un error común al tratar múltiples células de un mismo animal como independientes), el sistema calcula la **mediana por preparado** de todas las métricas. Esta mediana se utiliza como la unidad experimental para las pruebas de hipótesis.
    -   *Datos por Célula*: Se conservan para visualizar la variabilidad intrínseca y las distribuciones poblacionales.
    -   *Datos por Preparado*: Se utilizan para el cálculo de p-valores y distancias entre grupos.

2.  **Pruebas de Hipótesis Automatizadas**:
    -   Para cada métrica, se evalúa primero la normalidad de la distribución en ambos grupos mediante el test de **Shapiro-Wilk**.
    -   Si ambos grupos son normales: Se aplica el **Test t de Welch** (no asume varianzas iguales).
    -   Si alguno no es normal: Se aplica el test no paramétrico **U de Mann-Whitney**.
    -   El sistema reporta automáticamente el estadístico, el p-valor y el tamaño muestral (N).

### Librerías Principales

-   **`scipy.stats`**: Implementación de tests estadísticos (Shapiro, T-test, Mann-Whitney).
-   **`altair`**: Visualización de distribuciones.

### Estrategia de Análisis

Este enfoque jerárquico (Célula -> Preparado -> Grupo) garantiza que las conclusiones biológicas sean robustas y no estén sesgadas por el número variable de células detectadas en cada imagen.

## 7. Recursos de Software

Para facilitar la reproducibilidad de este estudio, a continuación se detallan las versiones específicas de las librerías principales utilizadas en el pipeline de análisis.

**Entorno de Desarrollo:**
*   **Lenguaje:** Python 3.10+
*   **Sistema Operativo:** Linux (Ubuntu 22.04 LTS o superior recomendado)

**Librerías de Procesamiento de Imagen:**
*   **`numpy`** (v2.2.6): Computación numérica y manejo de arrays multidimensionales.
*   **`scikit-image`** (v0.25.2): Preprocesamiento, filtrado, umbralización y operaciones morfológicas.
*   **`scipy`** (v1.16.1): Algoritmos de optimización y estructuras espaciales (Voronoi, KDTree).
*   **`tifffile`** (v2023.2.28): Lectura/escritura de alta fidelidad para formatos de microscopía (OME-TIFF, ImageJ Hyperstacks).

**Biología Computacional y Segmentación:**
*   **`cellpose`** (v4.0.6): Segmentación nuclear basada en Deep Learning (modelo `cyto2`/`nuclei`).
*   **`skan`** (v0.13.0): Análisis esqueletal de grafos y topología de redes.

**Análisis de Datos y Estadística:**
*   **`pandas`** (v2.3.2): Estructuración de datos, agregación y manejo de DataFrames.
*   **`scipy.stats`**: Pruebas de hipótesis (Shapiro-Wilk, Mann-Whitney U, t-Test de Welch).

**Visualización e Interfaz:**
*   **`napari`** (v0.6.4): Visualizador multidimensional n-dimensional para validación cualitativa.
*   **`streamlit`** (v1.39.0): Framework para la interfaz gráfica de usuario (GUI) interactiva.
*   **`altair`** (v5.5.0): Visualización estadística declarativa y gráficos interactivos.

El código fuente completo y los scripts de reproducibilidad están disponibles en el repositorio del proyecto.
