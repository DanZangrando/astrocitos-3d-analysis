# Pipeline de Análisis Morfológico de Astrocitos — Documentación Completa

**Última actualización:** 2025-11-11  
**Versión:** 2.0 (Pipeline 2D Unificado + Métricas Avanzadas)

---

## 📋 Tabla de Contenidos

1. [Arquitectura General](#arquitectura-general)
2. [Paso 01: Calibración](#paso-01-calibración)
3. [Paso 02: Segmentación Nuclear (Cellpose)](#paso-02-segmentación-nuclear)
4. [Paso 03: Filtrado de Astrocitos](#paso-03-filtrado-de-astrocitos)
5. [Paso 04: Pipeline 2D Unificado](#paso-04-pipeline-2d-unificado)
6. [Paso 05: Análisis por Preparado](#paso-05-análisis-por-preparado)
7. [Paso 06: Comparación entre Grupos](#paso-06-comparación-entre-grupos)
8. [Métricas Calculadas: Resumen Completo](#métricas-calculadas)
9. [Tests Estadísticos Aplicados](#tests-estadísticos)
10. [Troubleshooting](#troubleshooting)

---

## Arquitectura General

El pipeline consta de **4 pasos de procesamiento** (generan datos) y **2 páginas de análisis** (visualizan y comparan):

### Pasos de Procesamiento (01-04)
1. **Calibración** → Escala física (µm)
2. **Segmentación Nuclear** → Detectar núcleos (Cellpose)
3. **Filtrado de Astrocitos** → Selección por GFAP y tamaño
4. **Pipeline 2D Unificado** → Esqueletización + Sholl + Métricas topológicas

### Páginas de Análisis (05-06)
5. **Análisis por Preparado** → Revisión individual
6. **Comparación entre Grupos** → Estadística CTL vs Hipoxia

**Flujo completo:**
```
data/raw/*.tif
    ↓
[01] Calibración → calibration.json
    ↓
[02] Cellpose → 02_cellpose_mask.tif
    ↓
[03] Filtrado → 03_nucleus_metrics.csv + 04_final_astrocytes_mask.tif
    ↓
[04] Pipeline 2D → skeletons/summary.csv (22 columnas)
                → sholl_summary.csv + sholl_2d_native.csv
    ↓
[05] Análisis Individual
    ↓
[06] Comparación Grupos (estadística rigurosa)
```

---

## Paso 01: Calibración

**Objetivo:** Establecer escala física global (µm/píxel) para todo el experimento.

### Metodología

#### Detección Automática de Metadatos
- **TIFF:** Lee OME-XML, metadatos ImageJ, tags de resolución
- **LIF:** Usa `readlif` para extraer `scale` en (Z,Y,X) µm

#### Conversión de Unidades
Normaliza automáticamente: nm → µm, mm → µm, m → µm

#### Persistencia
Guarda en `streamlit/calibration.json`:
```json
{
  "z": 0.3788,
  "y": 0.3788,
  "x": 0.3788,
  "DAPI_CHANNEL_INDEX": 0,
  "GFAP_CHANNEL_INDEX": 1,
  ...
}
```

### Parámetros
- `z`, `y`, `x` — Resolución espacial (µm/píxel)
- `DAPI_CHANNEL_INDEX` — Índice del canal DAPI (típicamente 0)
- `GFAP_CHANNEL_INDEX` — Índice del canal GFAP (típicamente 1)

### Entradas/Salidas
- **Entrada:** Imágenes en `data/raw/CTL/` y `data/raw/hip/`
- **Salida:** `streamlit/calibration.json`

---

## Paso 02: Segmentación Nuclear

**Objetivo:** Generar máscaras de núcleos individuales con Cellpose 3D.

### Metodología

#### 1. Normalización de Ejes
Convierte todas las imágenes a formato estándar `(Z, C, Y, X)`

#### 2. Máscara de Otsu (Opcional)
- Umbral global sobre volumen DAPI
- Genera `01_otsu_mask.tif` como pre-filtro
- Elimina fondo antes de Cellpose

#### 3. Segmentación con Cellpose
- Modelo: `cyto2` o `nuclei`
- Modo 3D: `do_3D=True`, `z_axis=0`
- Genera labels únicos: `02_cellpose_mask.tif`
- Limpieza opcional: enmascara DAPI por Otsu previo

### Parámetros
- `NUCLEUS_DIAMETER` — Diámetro esperado (px, típicamente 30)
- `CELLPOSE_USE_GPU` — Usar GPU (true/false)

### Entradas/Salidas
- **Entrada:** Canal DAPI
- **Salidas:**
  - `01_otsu_mask.tif` — Máscara binaria 3D (opcional)
  - `02_cellpose_mask.tif` — Labels de núcleos
  - `params.json` — Parámetros usados

---

## Paso 03: Filtrado de Astrocitos

**Objetivo:** Seleccionar candidatos astrocitarios por señal GFAP perinuclear y tamaño físico.

### ⚠️ Nota Importante
El canal **Iba-1 ha sido eliminado** debido a contaminación detectada con GFAP. Filtrado basado **exclusivamente en GFAP**.

### Metodología

#### Paso 3a: Filtrado Biológico (GFAP Perinuclear)

**Estrategia:**
1. Expandir cada núcleo como "shell" mediante dilatación morfológica
2. Medir intensidad media GFAP en el shell (región expandida - núcleo)
3. Calcular fondo GFAP fuera de máscara Otsu: `media_fondo ± std_fondo`
4. **Regla:** Si `GFAP_shell > (media_fondo + N × std_fondo)` → candidato astrocito

**Ventajas del umbral relativo:**
- Robusto ante variaciones de intensidad entre preparados
- Se adapta automáticamente al nivel de señal
- Justificación estadística (desviaciones estándar)

#### Paso 3b: Filtrado por Tamaño Físico
1. Calcular volumen mínimo en voxels: `MIN_VOLUME_UM3 / (z × y × x)`
2. Retener solo núcleos con área ≥ mínimo
3. Elimina artefactos (muy pequeños) y fusiones (muy grandes)

### Métricas Generadas (`03_nucleus_metrics.csv`)

| Columna | Descripción | Unidades |
|---------|-------------|----------|
| `label` | ID del núcleo | - |
| `nucleus_volume_um3` | Volumen físico | µm³ |
| `nucleus_sphericity` | **Esfericidad 2D** (circularidad) | 0-1 |
| `shell_gfap_mean` | Intensidad GFAP en shell | a.u. |
| `is_astrocyte_candidate` | Pasó filtro GFAP | bool |

**Esfericidad 2D (Circularidad):**
- Fórmula: `4π × area / perímetro²`
- Calculada en proyección MIP (preparados planos)
- **Interpretación:**
  - 1.0 = círculo perfecto
  - 0.8-0.9 = núcleos redondeados típicos
  - 0.6-0.7 = alargados/irregulares
  - <0.6 = muy irregulares/fragmentados

### Parámetros

#### Filtro GFAP
- `SHELL_RADIUS_UM` — Radio del shell perinuclear (típicamente 2.0 µm)
- `GFAP_STD_DEV_THRESHOLD` — N° desviaciones sobre fondo (típicamente 3.0)
- `MAX_DILATION_ITERATIONS` — Máximo iteraciones (seguridad, 20)

#### Filtro por Tamaño
- `MIN_VOLUME_UM3` — Volumen mínimo núcleo válido (típicamente 75 µm³)

### Entradas/Salidas
- **Entradas:**
  - `02_cellpose_mask.tif`
  - Canal GFAP
  - `01_otsu_mask.tif`
- **Salidas:**
  - `03_gfap_filtered_mask.tif` — Pasaron filtro GFAP
  - `03_nucleus_metrics.csv` — Tabla de métricas
  - `04_final_astrocytes_mask.tif` — Máscara final tras filtro tamaño

---

## Paso 04: Pipeline 2D Unificado

**Objetivo:** Generar esqueletos 2D, calcular métricas topológicas avanzadas y análisis de Sholl en una ejecución integrada.

### 🎯 Filosofía 2D vs 3D

**¿Por qué 2D?**

Preparados de hipocampo son planos (~15-20 µm grosor total). Territorios astrocitarios se extienden principalmente en plano XY.

**Ventajas del enfoque 2D:**
- ✅ Mantiene resolución XY completa (0.38 µm, sin remuestreo)
- ✅ Sholl 2D nativo más preciso
- ✅ Voronoi 2D robusto y simple
- ✅ ~10x más rápido que 3D isotrópico
- ✅ Científicamente justificado para morfología territorial

### Metodología: 5 Sub-pasos Integrados

#### [1] Proyección 2D

**Máscaras nucleares:**
- MIP (Max Intensity Projection) preservando labels
- Para cada (y,x), toma label del slice Z con mayor área

**Señal GFAP:**
- MIP estándar (máximo por píxel en Z)
- Mantiene intensidad máxima por posición XY

#### [2] Territorios Voronoi

**Objetivo:** Definir territorios no solapados por célula

**Método:**
1. Centroides 2D de núcleos proyectados
2. Diagrama de Voronoi en plano XY
3. Zona de exclusión (gap) entre territorios

**Parámetro:**
- `TERRITORY_EXCLUSION_UM` — Gap entre territorios (típicamente 1-2 µm)

**Ventaja:** Evita entrelazamiento entre células vecinas

#### [3] Esqueletización + Filtrado de Conectividad

**Por cada territorio:**

1. **Restricción radial desde núcleo:**
   - Máscara circular: radio = `MAX_RADIUS_FROM_NUCLEUS_UM` (típicamente 100 µm)
   - Dominio = Voronoi ∩ círculo
   - **Justificación:** Astrocitos no extienden procesos indefinidamente

2. **Umbralización GFAP:** Otsu local en dominio restringido

3. **Conexión de fragmentos** (opcional):
   - Conecta fragmentos cercanos si distancia < `CONNECTION_RADIUS_UM`
   - No cruza entre territorios

4. **Esqueletización 2D:** Algoritmo de adelgazamiento morfológico

5. **🔥 FILTRADO DE CONECTIVIDAD (CRÍTICO):**
   - Detecta componentes conectadas del esqueleto
   - **RETIENE SOLO** la componente que toca el núcleo
   - **ELIMINA** fragmentos aislados dentro del territorio
   - **Justificación:** Fragmentos aislados = señal de fondo, procesos de células vecinas o artefactos; NO pertenecen al astrocito

**Salida:** `05_skeleton_labels_2d.tif` — Solo componentes conectadas al soma

#### [4] Análisis Topológico con SKAN

**SKAN (Skeleton Analysis)** construye grafo topológico y extrae métricas por rama.

**Base de datos interna (`summarize()`):**
- `skeleton_id` — ID componente conectada
- `branch_distance` — Longitud real rama (µm)
- `euclidean_distance` — Distancia directa extremos (µm)
- `branch_type` — 0=path, 1=endpoint, 2=junction, 3=cycle

**Métricas Calculadas (22 columnas en `skeletons/summary.csv`):**

##### A. Métricas Básicas
| Columna | Descripción | Unidades |
|---------|-------------|----------|
| `label` | ID astrocito | - |
| `n_branches` | Número total de ramas | count |
| `total_branch_length_um` | Suma longitudes | µm |
| `mean_branch_length_um` | Longitud promedio | µm |
| `n_endpoints` | Terminaciones libres | count |
| `n_junctions` | Bifurcaciones | count |

##### B. Tortuosidad (NUEVO - Relevante para Astrocitos en Hipoxia)

**Definición:** `tortuosidad = branch_distance / euclidean_distance`
- 1.0 = rama perfectamente recta
- >1.0 = rama sinuosa/curvada

| Columna | Descripción | Interpretación | Hipótesis Hipoxia |
|---------|-------------|----------------|-------------------|
| `tortuosity_mean` | Tortuosidad promedio | 1.05-1.15 normal, >1.2 desorganización | ↑ por estrés citoesquelético |
| `tortuosity_max` | Máxima tortuosidad | Detecta procesos extremos | ↑ heterogeneidad |
| `tortuosity_std` | Desv. std tortuosidad | Heterogeneidad morfológica | ↑ respuesta no uniforme |

**Relevancia biológica en hipoxia:**
- ↑ Tortuosidad → Desorganización citoesquelética, edema celular
- ↑ Desv. std → Respuesta morfológica heterogénea

##### C. Índices de Complejidad (NUEVO - Específicos para Astrocitos)

| Columna | Fórmula | Interpretación | Hipótesis Hipoxia |
|---------|---------|----------------|-------------------|
| `ramification_index` | `n_branches / max(n_junctions, 1)` | Ramas por bifurcación, 2-4 típico | ↓ simplificación/atrofia |
| `termination_index` | `n_endpoints / n_junctions` | Balance extensión/ramificación | ↑ búsqueda O₂ (vasos) |

**Interpretación `termination_index`:**
- Alto (>2): Muchas terminaciones → extensión hacia targets
- Bajo (<2): Ramificación central predominante

##### D. Distribución de Longitudes (NUEVO - Heterogeneidad Morfológica)

| Columna | Descripción | Interpretación |
|---------|-------------|----------------|
| `branch_length_median_um` | Mediana longitudes | Más robusta que media |
| `branch_length_p25_um` | Percentil 25 | Cuartil inferior |
| `branch_length_p75_um` | Percentil 75 | Cuartil superior |
| `branch_length_std_um` | Desviación estándar | Variabilidad absoluta |
| `branch_length_cv` | Coeficiente variación (std/mean) | **Clave:** 0.3-0.5 normal, >0.6 desorganización |

**Relevancia biológica:**
- ↑ CV en hipoxia → Morfología desorganizada, respuesta heterogénea
- Percentiles detectan si afecta procesos cortos vs largos

##### E. Fragmentación

| Columna | Descripción | Valor Esperado |
|---------|-------------|----------------|
| `n_connected_components` | Componentes conectadas | 1 (esqueleto íntegro) |

**Nota:** El filtro de conectividad debería prevenir valores >1, pero se reporta para QC

#### [5] Análisis de Sholl 2D Nativo

**Principio:** Cuantificar complejidad territorial mediante intersecciones con anillos concéntricos

**Método:**
1. Centro = centroide nuclear (coordenadas físicas µm)
2. Anillos: de `SHOLL_MIN_UM` a `SHOLL_MAX_UM` con paso `SHOLL_STEP_UM`
3. SKAN `sholl_analysis()` cuenta **cruces de branches** (método estándar desde Sholl, 1953)
4. Genera perfil: intersecciones vs radio

**Importante:** SKAN cuenta cruces de RAMAS, no píxeles. Un branch puede cruzar múltiples anillos.

**Métricas Agregadas (`sholl_summary.csv`):**

| Columna | Descripción | Unidades | Interpretación | Hipótesis Hipoxia |
|---------|-------------|----------|----------------|-------------------|
| `auc` | Área bajo curva Sholl | µm·inter. | Complejidad territorial total | ↓ pérdida arborización |
| `peak_intersections` | Máximo intersecciones | count | Máxima densidad ramificación | ↓ menor ramificación |
| `critical_radius_um` | Radio del pico | µm | Distancia máxima arborización | Puede cambiar con reactividad |

**Perfiles Completos (`sholl_2d_native.csv`):**
- Columnas: `label`, `radius_um`, `intersections`
- Una fila por (astrocito, radio)
- Para graficar perfiles individuales

**Visualización (`sholl_rings_2d_native.json`):**
```json
{
  "1": {
    "centroid_um": [50.5, 75.2],
    "radii_um": [5.0, 10.0, 15.0, ...]
  },
  ...
}
```

### Parámetros Configurables

| Parámetro | Descripción | Típico |
|-----------|-------------|--------|
| `PROJECTION_2D_METHOD` | Método proyección | 'max' |
| `TERRITORY_EXCLUSION_UM` | Gap Voronoi | 1-2 µm |
| `MAX_RADIUS_FROM_NUCLEUS_UM` | Radio máximo análisis | 100 µm |
| `CONNECT_SKELETON_FRAGMENTS` | Habilitar conexión | true/false |
| `CONNECTION_RADIUS_UM` | Radio máximo conexión | 0.5-1.0 µm |
| `SHOLL_MIN_UM` | Radio mínimo Sholl | 5.0 µm |
| `SHOLL_MAX_UM` | Radio máximo Sholl | 20-30 µm |
| `SHOLL_STEP_UM` | Separación anillos | 5.0 µm |

### Entradas/Salidas
- **Entradas:**
  - `04_final_astrocytes_mask.tif`
  - Canal GFAP 3D
  - `calibration.json`
- **Salidas:**
  - `gfap_projection_2d.tif`
  - `05_skeleton_labels_2d.tif`
  - `skeletons/summary.csv` (**22 columnas**)
  - `sholl_2d_native.csv`
  - `sholl_summary.csv`
  - `sholl_rings_2d_native.json`

---

## Paso 05: Análisis por Preparado

**Objetivo:** Revisión individual de resultados por preparado con visualización completa.

### Contenido

#### Dashboard de Conteo
- #Núcleos detectados (Cellpose)
- #Candidatos GFAP
- #Astrocitos finales
- #Células con Sholl completo

#### Métricas Nucleares
**Fuente:** `03_nucleus_metrics.csv`

**Visualizaciones:**
- Tabla detallada por núcleo
- Histograma volumen (coloreado por retención)
- Histograma esfericidad 2D
- Scatter: GFAP vs Volumen (umbral de decisión)

#### Métricas Topológicas
**Fuente:** `skeletons/summary.csv`

**Organizado en 3 tabs:**
1. **📏 Básicas:** Longitud, endpoints, junctions, ramas
2. **🌀 Tortuosidad:** Media, máx, std + histograma
3. **🔀 Complejidad:** Índices + scatter (ramificación vs tortuosidad)

#### Métricas de Sholl
**Fuente:** `sholl_summary.csv` y `sholl_2d_native.csv`

**Visualizaciones:**
- Métricas resumen: AUC mediana, pico mediano, radio crítico
- Tabla detallada por célula
- Perfiles de Sholl: líneas intersecciones vs radio (interactivo)

#### Visualización Napari
**Modo 2D (recomendado):**
- GFAP proyección (fondo)
- Esqueletos 2D coloreados
- Anillos Sholl superpuestos
- Escala física aplicada

### Archivos Leídos
- `02_cellpose_mask.tif`
- `03_nucleus_metrics.csv`
- `04_final_astrocytes_mask.tif`
- `05_skeleton_labels_2d.tif`
- `skeletons/summary.csv`
- `sholl_2d_native.csv`
- `sholl_summary.csv`
- `sholl_rings_2d_native.json`

---

## Paso 06: Comparación entre Grupos

**Objetivo:** Comparar métricas entre grupos experimentales con tests estadísticos rigurosos que evitan pseudoreplicación.

### Problema: Pseudoreplicación

**Definición:** Tratar células individuales como unidades independientes cuando provienen del mismo preparado.

**Consecuencia:** Infla N artificialmente → falsos positivos

**Solución:**
1. **DataFrame para gráficos (`df_plot`):** Todas las células (visualización)
2. **DataFrame para estadística (`df_stats`):** Mediana por preparado
   - N = número de preparados (unidad experimental verdadera)

### Metodología Estadística

#### 1. Consolidación
Agrega datos de todos los preparados:
- `sholl_summary.csv`
- `skeletons/summary.csv`
- `03_nucleus_metrics.csv`

#### 2. Agregación por Preparado
Para cada métrica: **mediana por preparado**

```
Preparado_A_CTL → [cel1, cel2, ..., celM] → mediana_A
Preparado_B_CTL → [...] → mediana_B
...
Preparado_X_Hip → [...] → mediana_X
```

Resultado: 1 valor por preparado (evita pseudoreplicación)

#### 3. Test de Normalidad
**Shapiro-Wilk** por grupo (α = 0.05):
- H₀: Distribución normal
- p > 0.05 → normal

#### 4. Selección de Test

**Caso 1: Ambos grupos normales**
→ **Welch's t-test**
- Compara medias
- No asume varianzas iguales
- Paramétrico

**Caso 2: Al menos uno no normal**
→ **Mann-Whitney U**
- Compara medianas
- No paramétrico
- Robusto a outliers

#### 5. Reporte
- Test usado
- Estadístico
- p-valor
- Interpretación (p < 0.05 = significativo)
- N por grupo (preparados)
- Medias/medianas por grupo

### Métricas Disponibles

#### Métricas Topológicas (SKAN)
- `total_branch_length_um` — Longitud total esqueleto
- `tortuosity_mean` — Tortuosidad promedio
- `tortuosity_max` — Tortuosidad máxima
- `tortuosity_std` — Desv. std tortuosidad
- `ramification_index` — Ramas / Junctions
- `termination_index` — Endpoints / Junctions
- `branch_length_cv` — Coef. variación longitudes
- `n_endpoints` — Terminaciones
- `n_junctions` — Bifurcaciones

#### Métricas de Sholl
- `auc` — Área bajo curva
- `peak_intersections` — Pico intersecciones
- `critical_radius_um` — Radio crítico

#### Métricas Nucleares
- `nucleus_volume_um3` — Volumen soma
- `nucleus_sphericity` — Esfericidad 2D

### Visualizaciones

**Boxplots + Puntos Individuales:**
- Box = cuartiles agregados
- Puntos = células individuales (coloreados por grupo)
- Muestra variabilidad intra-preparado

**Tablas Exportables:**
- Por célula (`df_plot.csv`)
- Por preparado (`df_stats.csv`)

### Interpretación

**p < 0.05:** Diferencia estadísticamente significativa
- Rechazamos H₀
- Evidencia de diferencia real

**p ≥ 0.05:** No hay evidencia suficiente
- No rechazamos H₀
- Puede haber diferencia pero sin poder estadístico

**Nota:** Significancia ≠ Relevancia biológica. Evaluar magnitud del efecto.

---

## Métricas Calculadas: Resumen Completo

### Por Paso

| Paso | Archivo | Métricas |
|------|---------|----------|
| 03 | `03_nucleus_metrics.csv` | 5 columnas: label, volume, sphericity, GFAP, flag |
| 04 | `skeletons/summary.csv` | **22 columnas** topológicas + tortuosidad + complejidad |
| 04 | `sholl_summary.csv` | 4 columnas: label, AUC, peak, critical_radius |
| 04 | `sholl_2d_native.csv` | 3 columnas: label, radius_um, intersections (perfiles) |

### Métricas Topológicas (22 columnas)

**Básicas (6):**
- `label`, `n_branches`, `total_branch_length_um`, `mean_branch_length_um`, `n_endpoints`, `n_junctions`

**Tortuosidad (3):**
- `tortuosity_mean`, `tortuosity_max`, `tortuosity_std`

**Complejidad (2):**
- `ramification_index`, `termination_index`

**Distribución Longitudes (5):**
- `branch_length_median_um`, `branch_length_p25_um`, `branch_length_p75_um`, `branch_length_std_um`, `branch_length_cv`

**Fragmentación (1):**
- `n_connected_components`

**Espaciales (5):**
- `skeleton_pixels`, `nuclear_area_um2`, `centroid_y_um`, `centroid_x_um`, (más columnas según implementación específica)

---

## Tests Estadísticos Aplicados

### Test de Normalidad
**Shapiro-Wilk** (α = 0.05):
- Aplicado a datos agregados por preparado (no por célula)
- Por grupo (CTL, Hipoxia)
- Determina si usar test paramétrico o no paramétrico

### Tests de Comparación

**Welch's t-test** (paramétrico):
- **Cuándo:** Ambos grupos normales
- **Hipótesis nula:** μ_CTL = μ_Hip
- **Compara:** Medias
- **Ventaja:** No asume varianzas iguales

**Mann-Whitney U** (no paramétrico):
- **Cuándo:** Al menos un grupo no normal
- **Hipótesis nula:** Distribuciones idénticas
- **Compara:** Medianas/rangos
- **Ventaja:** Robusto a outliers, no asume distribución

### Nivel de Significancia
- **α = 0.05** (estándar)
- **p < 0.05:** Rechazo H₀ → diferencia significativa
- **p ≥ 0.05:** No rechazo H₀ → sin evidencia suficiente

### Corrección por Tests Múltiples
**Actualmente NO implementada**

**Recomendación para publicación:**
- Aplicar Bonferroni si testear >10 métricas
- α_ajustado = 0.05 / n_tests
- O usar FDR (False Discovery Rate)

---

## Troubleshooting

### Pocos astrocitos detectados
**Causas:**
- Umbral GFAP muy estricto
- Filtro de tamaño demasiado grande

**Soluciones:**
- ↓ `GFAP_STD_DEV_THRESHOLD` (de 3.0 a 2.0)
- ↓ `MIN_VOLUME_UM3` (de 75 a 50)

### Territorios solapados
**Causa:** Gap de exclusión insuficiente

**Solución:**
- ↑ `TERRITORY_EXCLUSION_UM` (de 1.0 a 2-3 µm)

### Esqueletos fragmentados
**Causa:** Señal GFAP discontinua

**Soluciones:**
1. Activar `CONNECT_SKELETON_FRAGMENTS=true`
2. Ajustar `CONNECTION_RADIUS_UM` (0.5-1.5 µm)
3. Verificar calidad imaging (Z-sampling adecuado)

### Sholl sin pico claro
**Causas:**
- Astrocitos sin procesos extendidos (solo soma)
- Rango de radios inadecuado
- MAX_RADIUS muy pequeño

**Soluciones:**
- Verificar que `MAX_RADIUS_FROM_NUCLEUS_UM` ≥ `SHOLL_MAX_UM`
- Ajustar rango: `SHOLL_MIN_UM=5`, `SHOLL_MAX_UM=25`, `SHOLL_STEP_UM=5`
- Visualizar en Napari: ¿hay procesos GFAP extendidos?

### Fragmentos desconectados persisten
**Causa:** Filtro de conectividad no aplicado

**Verificación:**
- Revisar código: `filter_skeleton_by_nucleus_connectivity()` debe estar activo
- Revisar output en terminal: debe reportar componentes removidas

**Solución:**
- Regenerar paso 04 con código actualizado

### Métricas de tortuosidad = NaN
**Causa:** Esqueleto vacío o sin ramas

**Solución:**
- Verificar que astrocitos tienen señal GFAP
- Revisar umbrales de GFAP (paso 03-04)

---

## Resumen de Archivos Generados

```
data/processed/<preparado>/
├── 01_otsu_mask.tif                    [Opcional]
├── 02_cellpose_mask.tif                [3D labels nucleares]
├── 03_gfap_filtered_mask.tif           [Pasaron filtro GFAP]
├── 03_nucleus_metrics.csv              [5 columnas: métricas nucleares]
├── 04_final_astrocytes_mask.tif        [Pasaron filtro tamaño]
├── 05_skeleton_labels_2d.tif           [Esqueletos 2D, solo conectados]
├── gfap_projection_2d.tif              [MIP de GFAP]
├── params.json                         [Parámetros usados]
├── sholl_2d_native.csv                 [Perfiles Sholl completos]
├── sholl_summary.csv                   [AUC, pico, crítico]
├── sholl_rings_2d_native.json          [Coords anillos Napari]
└── skeletons/
    └── summary.csv                     [22 columnas topológicas]
```

---

## Referencias Científicas

### Metodología Sholl
- **Sholl, D.A. (1953).** Dendritic organization in the neurons of the visual and motor cortices of the cat. *J Anat*, 87, 387-406.

### Astrocitos: Morfología y Territorios
- **Bushong et al. (2002).** Protoplasmic astrocytes in CA1 stratum radiatum occupy separate anatomical domains. *Neuron*, 34, 127-138.
- **Oberheim et al. (2012).** Uniquely hominid features of adult human astrocytes. *J Neurosci*, 32, 3176-3187.

### Reactividad Astrocitaria
- **Sofroniew, M.V. (2009).** Molecular dissection of reactive astrogliosis and glial scar formation. *Trends Neurosci*, 32, 638-647.
- **Hauglund et al. (2020).** Cleaning the sleeping brain. *Nat Commun*, 11, 5285.

### Análisis de Esqueletos
- **Nunez-Iglesias et al. (2018).** Skeleton Analysis (SKAN). *Journal of Open Source Software*, 3(24), 1382.

---

**Última actualización:** 2025-11-11  
**Autor:** Pipeline desarrollado para análisis de astrocitos en hipoxia vs normoxia  
**Contacto:** Ver repositorio GitHub para issues y contribuciones
