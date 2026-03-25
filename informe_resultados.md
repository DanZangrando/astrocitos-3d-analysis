# Informe Técnico: Pipeline Computacional para el Análisis Morfológico de Astrocitos

## 1. Introducción
El presente informe documenta el flujo de trabajo computacional desarrollado para la cuantificación objetiva de la morfología en astrocitos. El conjunto de datos de entrada consiste en imágenes tridimensionales de microscopía (*z-stacks*), en las cuales se adquirieron dos canales fluorescentes:
- **DAPI**: Marcador nuclear utilizado para la identificación de somas celulares.
- **GFAP**: Marcador del citoesqueleto empleado para visualizar la arborización astrocitaria.

Las imágenes corresponden a dos grupos experimentales distintos: un grupo **Control (CTL)** y un grupo bajo condiciones de **Hipoxia**. El objetivo de este documento es describir paso a paso la metodología automatizada empleada para segmentar individualmente cada astrocito, extraer métricas topológicas precisas y evitar sesgos de medición, proporcionando el marco técnico y estadístico de los resultados obtenidos.

## 2. Herramientas Computacionales
El análisis estructurado se ejecutó localmente utilizando el lenguaje **Python**, integrando algoritmos de propósito específico:
- **Cellpose**: Modelo de *Deep Learning* de estado del arte adoptado para la segmentación automatizada e independiente de núcleos a partir de la señal DAPI.
- **SKAN (Skeleton Analysis)**: Librería científica matemática especializada en el análisis rápido de grafos (branching). Fue utilizada como *gold-standard* para la evaluación morfológica, ya que permite cuantificar "ramas" conectadas físicamente en vez de realizar un simple conteo de píxeles, maximizando la pureza del análisis de Sholl y la topología 2D.

## 3. Metodología: Pipeline Computacional
El flujo de procesamiento (*pipeline*) operó de la siguiente forma secuencial para garantizar la individualización precisa de cada red de citoesqueleto astrocitario:

### 3.1. Estandarización Paramétrica (Calibración) y Procesamiento Predictivo en Lote
Previo a la ejecución cuantitativa a gran escala, los parámetros analíticos (umbrales de intensidad GFAP, márgenes de Voronoi y umbrales de tamaño nuclear) fueron calibrados de forma interactiva sobre un subconjunto representativo y aleatorio de imágenes (*z-stacks*) provenientes de ambos grupos. Una vez validados visualmente, estos parámetros fueron fijados globalmente para todo el experimento. A continuación, el volumen total de base de datos fue analizado mediante un sistema de procesamiento desatendido en lote (*batch processing*). Esto asegura de manera estricta que todas las imágenes—independientemente de su grupo de origen—reciban un idéntico criterio de medición bioinformática iterativa, previniendo sesgos de experimentador (*blinded processing*).

### 3.2. Segmentación Nuclear y Filtrado de Población
1. **Segmentación 3D (DAPI)**: Se aislaron las máscaras nucleares individuales. Se aplicaron filtros de exclusión volumétrica para descartar partículas anómalas por tamaño. 
2. **Asignación de GFAP**: Se inspeccionó el vecindario volumétrico inmediato de cada núcleo segmentado; aquellos que validaron colocalización significativa con GFAP fueron catalogados positivamente como astrocitos.
3. **Filtro de Exclusión de Bordes (Border Clearance)**: Con el fin de evitar sesgos analíticos por células truncadas, **se descartó todo astrocito cuyo núcleo estuviera situado a menos de 40 µm** de los márgenes geométricos (XY o Z) de la imagen. Esto garantiza que todos los dominios celulares analizados estén completamente contenidos en el campo de visión.

### 3.2. Proyección 2D y Regionalización (Voronoi)
Para el análisis morfológico detallado, se implementó un enfoque 2D que mantiene la fidelidad de la resolución nativa lateral (0.38 µm/px):
1. **Proyección (MIP)**: Las imágenes volumétricas (núcleos y GFAP) se colapsaron mediante una Proyección de Máxima Intensidad en el plano 2D.
2. **Territorios Funcionales de Voronoi**: Dado el alto nivel de superposición natural de los procesos intercelulares, se calculó una **tesselación de Voronoi** usando los centroides nucleares detectados. Este algoritmo asigna un dominio bidimensional inequívoco a cada célula.
3. **Confinamiento Espacial**: Al diagrama de Voronoi se le añadió un margen geométrico mínimo de exclusión (*gap* de 0.1 µm) entre fronteras. Adicionalmente, el territorio de cualquier célula fue restringido numéricamente a un radio máximo circular de 25 µm desde su suma celular. Toda señal GFAP fuera de esta máscara individual fue ignorada.

### 3.3. Obtención del Esqueleto Individual
Dentro del límite asignado para cada astrocito:
1. **Umbralización Bimodal Local**: La señal fluorescente GFAP de cada región fue segmentada de forma adaptativa.
2. **Reconexión de Fragmentos**: Se aplicó una operación morfológica (radio de conexión cerrado de ~0.5 µm) para compensar interrupciones de señal y unir fragmentos contiguos locales propios del citoesqueleto GFAP.
3. **Esqueletización Topológica**: Se extrajo la línea constitutiva central del dominio GFAP, resultando en ramas de 1 píxel de espesor.
4. **Filtrado de Conectividad Estricta**: Para asegurar medir íntegramente a cada célula, se identificaron todas las componentes desconectadas del esqueleto y se filtraron, conservando **únicamente** aquella cuya ramificación principal posea conectividad física directa con la máscara central nuclear.

### 3.4. Extracción de Métricas (SKAN y Sholl 2D)
El esqueleto limpio y categorizado se analizó matemáticamente.
- **Topología General**:
  - `total_branch_length_um`: Longitud física (µm) agregada de todas las ramificaciones esqueléticas.
  - `n_branches`: Número total de segmentos o cortes entre nodos (conteo de ramas individuales).
  - `n_junctions`: Número de nodos o puntos de convergencia/bifurcación de la arborización.
  - `ramification_index`: Ratio del número de ramas sobre ramas divido número de nodos, indicando la complejidad del tendido.
  - `tortuosity_mean`: Relación de curvatura calculada entre la longitud euclidiana y la longitud del segmento para describir la rectitud intrínseca del citoesqueleto.
- **Análisis de Sholl Clásico 2D**: Se configuraron circunferencias concéntricas empezando desde 0 µm (zona próxima al núcleo) hasta abarcar los 25 µm, a intervalos estrictos de 0.5 µm. El algoritmo nativo de SKAN evaluó el corte de **Ramas** genuinas. Fueron derivadas las siguientes métricas: el número máximo de intersecciones detectado (`peak_intersections`), el radio de ocurrencia de dicho pico (`critical_radius_um`), y la matriz de Área Bajo la Curva completa (`auc`).

## 4. Diseño del Análisis Estadístico
Para evitar pseudoreplicación y asegurar una representación poblacional:
1. **Agrupación Biológica Central**: Todos los resultados fueron condensados y testeados a nivel estadístico de "Preparado Individual" (Animal), utilizando el valor general (**Mediana**) de cada población celular. El "N" efectivo final corresponde a N=23 (Control) y N=21 (Hipoxia), respectivamente.
2. **Análisis Univariado (Decisión Paramétrica/No-Paramétrica)**: Para cada métrica morfológica evaluada, se analizó inicialmente la distribución de las medianas muestrales para determinar su normalidad. Basado en este perfil estructural, se decidió el algoritmo estadístico adecuado caso por caso:
   - Se aplicó el **Test t de Welch** (Prueba Paramétrica para varianzas desiguales) cuando los datos cumplían asunciones de normalidad.
   - Se aplicó la **Prueba U de Mann-Whitney** (Prueba No Paramétrica en bases de rangos) cuando los datos no presentaban una distribución normal clásica.
   - **Umbral de Significancia**: Se definió un $p < 0.05$ como diferencia estadísticamente significativa. Para aquellos resultados con un valor $p$ situado entre $0.05$ y $0.10$, se documentaron y analizaron como una **tendencia biológica a la diferencia**.
3. **Test Multivariado de Fenotipos Globales (PCA y Hotelling T²)**: Para descartar la redundancia de las variables individuales y validar la existencia de un cambio morfológico sistémico, se aplicó un **Análisis de Componentes Principales (PCA)** configurado analíticamente con cuatro métricas de la red citoesquelética directa (`total_branch_length_um`, `peak_intersections`, `auc`, `critical_radius_um`). La significancia estadística de la separación global entre fenotipos en este espacio dimensional se calculó numéricamente con el **Test T² de Hotelling**. Este algoritmo es el análogo multivariado de la prueba T de Student: en lugar de comparar un único promedio geométrico, contrasta simultáneamente los vectores de medias (es decir, la posición de los centroides) de ambos grupos para múltiples variables a la vez. Al trabajar en todo nuestro diseño con valores aglutinados por preparación/animal, el influjo del Teorema del Límite Central permite asumir una tendencia natural hacia la normalidad multivariante en las distribuciones, cumpliendo satisfactoriamente la asunción fundamental requerida por el estadístico T² para medir la probabilidad de que la distancia entre los centroides de Control e Hipoxia no se deba enteramente al azar exploratorio.

## 5. Resultados Estadísticos (Medianas por Preparado)

### 5.1. Impacto Morfológico Global (Retracción del Citoesqueleto)
El análisis univariado del citoesqueleto GFAP reveló un fenotipo celular de dimensiones reducidas en tejido bajo Hipoxia en comparación con los controles (N=23 CTL, N=21 Hipoxia):
- **Longitud Total del Esqueleto**: Se midió la suma total de las ramas del esqueleto GFAP, lo que funciona como un claro estimador de la masa citoesquelética o gasto energético expansivo. Se observó una disminución muy significativa en el grupo hipóxico (**Mediana Hipoxia = 46.62 µm** [rango intercuartílico IQR: 37.74 - 68.69] frente a **Mediana CTL = 72.89 µm** [IQR: 42.24 - 96.54]). Esta caída ($t = 2.497, p = 0.018$, Test t de Welch) ilustra una clara retracción en el campo ramificacional de la célula.
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Longitud Total Esqueleto]***

- **Nodos y Bifurcaciones**: Complementario a la retracción en longitud, se cuantificó un decremento estadísticamente significativo en la cantidad de nodos de ramificación o bifurcaciones en el árbol GFAP (**Mediana Hipoxia = 6.5** [IQR: 5.0 - 9.0] frente a **Mediana CTL = 9.0** [IQR: 5.0 - 13.25]; $t = 2.255, p = 0.031$, Test t de Welch). Esto señala que los astrocitos bajo hipoxia no solo cubren menos territorio de longitud, sino que son funcionalmente menos arborescentes.
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Número de Nodos/Bifurcaciones]***

- **Número Total de Ramas**: El recuento absoluto de los distintos gajos independientes que forman la célula mostró una tendencia fuerte a la reducción en los individuos hipóxicos (**Mediana Hipoxia = 13.0** [IQR: 10.0 - 17.0] frente a **Mediana CTL = 19.0** [IQR: 10.5 - 24.0]; $U = 314.50, p = 0.088$, Prueba U de Mann-Whitney). Esta tendencia consolida el hallazgo previo de pérdida de complejidad topológica.
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Número Total de Ramas]***

- **Índice de Ramificación (Ratio Ramas/Nodos)**: Resultó de particular interés que, a pesar de tener menos masa total, los astrocitos hipóxicos presentaron un incremento estadísticamente significativo en su Índice de Ramificación (**Mediana Hipoxia = 2.06** [IQR: 2.00 - 2.17] frente a **Mediana CTL = 2.00** [IQR: 1.90 - 2.00]; $U = 156.000, p = 0.044$, Prueba U de Mann-Whitney). Esta relación (más ramas remanentes por cada nodo sobreviviente) sugiere una reconfiguración geométrica que tendería a perder los "hubs" arborescentes más complejos, propiciando tramos o cadenas más simples y lineales.
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Índice de Ramificación]***

- **Tortuosidad Media**: El parámetro de tortuosidad a nivel de ramas no presentó diferencias estadísticamente constatables entre grupos (**Mediana Hipoxia = 1.09** [IQR: 1.09 - 1.10] y **Mediana CTL = 1.10** [IQR: 1.09 - 1.10]; $t = 0.627, p = 0.534$, Test t de Welch). Concluimos que, a pesar de retraerse, lo que resta del esqueleto preserva sus propiedades geométricas fundamentales (sinuosidad intrínseca fisiológica) y descarta una deformación o arrugamiento pasivo de los filamentos.
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Tortuosidad Media]***

### 5.2. Reorganización Espacial del Árbol (Análisis Sholl)
El decaimiento radial de la arborización somatofugal demostró que la pérdida de masa no es un efecto uniformemente esparcido o un simple error de escalado:
- **Área Bajo la Curva (AUC)**: Empleado como resumen holístico integral de la red, la función total del test de Sholl evidenció una pérdida altamente significativa del volumen de distribución de procesos a favor de un astrocito desprovisto en hipoxia (**Mediana Hipoxia = 40.43** [IQR: 28.27 - 53.21] vs **Mediana CTL = 53.50** [IQR: 36.14 - 81.05]; $t = 2.576, p = 0.014$, Test t de Welch). 
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Área Bajo la Curva Sholl]***

- **Perfil de Curva Sholl (Intersecciones Espaciales)**: Visualizando directamente los promedios espaciales (con error estándar de la media poblacional), la traza hipóxica permaneció invariablemente deprimida ubicándose debajo de la curva del Control a lo largo de prácticamente todo el dominio territorial investigado (desde 4 µm hasta el fin del territorio de 25 µm), corroborando fenotípicamente la pérdida progresiva y severa de intersecciones GFAP.
  > ***[INSERTAR GRÁFICA AQUÍ: Gráfico de Líneas - Análisis de Curvas de Sholl]***

- **Radio Crítico (Distancia de Máxima Complejidad)**: Esta distancia representa el perímetro alrededor del astrocito donde se intercepta la mayor densidad simultánea de ramas. El grupo con Hipoxia presentó una **tendencia a la contracción proximal** de su radio crítico respecto al Control (**Mediana Hipoxia = 4.50 µm** [IQR: 4.0 - 5.0] versus **Mediana CTL = 5.75 µm** [IQR: 4.0 - 8.0]; $U = 320.00, p = 0.066$, Prueba U de Mann-Whitney). Esto significa que bajo hipoxia, el astrocito tiende a compactar su zona central de máxima actividad arquitectónica, aproximándola forzosamente hacia su soma.
  > ***[INSERTAR GRÁFICA AQUÍ: Boxplot Vertical - Radio Crítico Sholl]***

### 5.3. Segregación Multivariante (PCA Global del Citoesqueleto)
Habiendo descartado la influencia de variables somato-nucleares para enfocarnos estrictamente en el "hardware" de la arborización astrocitaria, se ejecutó un Análisis de Componentes Principales (PCA) ingresando exclusivamente cuatro macro-variables morfológicas:
1. `total_branch_length_um`
2. `peak_intersections` (Densidad máxima de Sholl)
3. `auc` (Área bajo la curva Sholl)
4. `critical_radius_um` (Distancia al pico de ramificación)

Bajo este modelo focalizado, la reducción de dimensionalidad indicó que la primera Componente Principal (**PC1**) logra explicar de manera abrumadora el **92.1%** de toda la varianza morfológica de la población, relegando un ínfimo 4.7% al eje ortogonal secundario (PC2). Este acaparamiento rotundo de la varianza en una única dimensión unidireccional (PC1) revela que las cuatro métricas individuales de atrofia (longitudes, intersecciones, radio de alcance y área) decaen de manera conjunta, simultánea y altamente correlacionada.

Al proyectar los centroides de ambos grupos experimentales en el hiperplano PC1 vs PC2, se evidenció la separación formal de los fenotipos. Esta divergencia multivariada fue validada probabilísticamente mediante el **Test T² de Hotelling** sobre las distribuciones combinadas, estableciendo una definida significancia en la distancia multivariada biológica ($T^2 = 7.58, \text{Distancia} = 1.289, p = 0.033$). 

De esta manera se confirma formalmente la premisa central del análisis tridimensional propuesto: la reacción astrocitaria hipóxica no se ciñe a una simple disminución de una métrica aislada, sino a un colapso sistémico y orquestado de toda su red arborescente (cuantificada íntegramente de manera indirecta sobre PC1).
  > ***[INSERTAR GRÁFICA AQUÍ: Proyección PCA 2D (Elipses de Confianza) con Estadísticas Multivariadas]***
