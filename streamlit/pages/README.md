# Páginas del Pipeline — Metodología Detallada

Este documento describe en detalle la metodología aplicada en cada página de la aplicación Streamlit para el análisis morfológico de astrocitos.

---

## 🎯 Arquitectura del Pipeline (4 Pasos + 2 Análisis)

El pipeline se divide en **4 pasos de procesamiento** que generan datos, seguidos de **2 páginas de análisis** que visualizan y comparan resultados:

### Pasos de Procesamiento (01-04)
1. **Calibración** — Definir escala física (µm)
2. **Segmentación Nuclear** — Detectar núcleos con Cellpose
3. **Filtrado de Astrocitos** — Seleccionar candidatos por GFAP y tamaño
4. **Pipeline 2D Unificado** — Esqueletización + Sholl en una sola ejecución

### Páginas de Análisis (05-06)
5. **Análisis por Preparado** — Revisión individual con métricas y visualización
6. **Comparación entre Grupos** — Estadística CTL vs Hipoxia

---

## 01 — Calibración y Visualización de los Preparados

**Objetivo:** Establecer la escala física global (µm) y verificarla sobre los archivos `.tif`/`.lif`.

### Metodología

#### Detección de Metadatos
- **TIFF:** Lee OME-XML, metadatos ImageJ (`spacing`, `unit`) y tags de resolución
- **LIF:** Usa `readlif` y toma `scale` como (Z,Y,X) en µm si está disponible

#### Conversión de Unidades
Convierte automáticamente: nm → µm / mm → µm / m → µm según corresponda

#### Persistencia
La calibración global (z, y, x) se guarda en `streamlit/calibration.json`, manteniendo otros parámetros globales

#### Visualización
Abre Napari con la imagen y la escala física aplicada para verificación visual

#### UI
- Contadores de archivos encontrados por grupo (CTL/Hipoxia)
- Métricas rápidas de calibración detectada por eje
- Vista previa de dimensiones y metadatos

### Entradas/Salidas
- **Entrada:** Imágenes en `data/raw/` (organizadas por grupo)
- **Salida:** `streamlit/calibration.json` (global para todo el experimento)

### Parámetros Configurables
- `z`, `y`, `x` — Resolución espacial en µm
- `DAPI_CHANNEL_INDEX` — Índice del canal DAPI (típicamente 0)
- `GFAP_CHANNEL_INDEX` — Índice del canal GFAP (típicamente 1)

---

## 02 — Umbral de Otsu y Segmentación con Cellpose (DAPI)

**Objetivo:** Generar máscaras de núcleos mediante segmentación 3D con Cellpose.

### Metodología

#### Re-ordenamiento de Ejes
Todas las imágenes se llevan a formato estándar (Z,C,Y,X) con heurísticas para archivos `.lif`

#### Paso 01: Máscara de Otsu (Opcional)
- Calcula umbral global sobre el volumen DAPI
- Genera `01_otsu_mask.tif` como pre-filtro opcional
- Útil para eliminar fondo antes de Cellpose

#### Paso 02: Segmentación con Cellpose
- Modelo: `CellposeModel` (cyto2 o nuclei según configuración)
- Modo: Segmentación 3D (`do_3D=True`, `z_axis=0`)
- Limpieza opcional: enmascara DAPI por Otsu antes de segmentar
- Genera labels únicos por núcleo en `02_cellpose_mask.tif`

#### Persistencia
Guarda parámetros en `params.json` por preparado:
- Umbral Otsu calculado
- Diámetro de núcleo usado
- Uso de GPU (true/false)

#### Visualización
- Métricas: umbral Otsu, fracción enmascarada, número de núcleos detectados
- Napari: carga DAPI + máscara Otsu + labels Cellpose con escala física

### Entradas/Salidas
- **Entrada:** Canal DAPI desde la imagen raw
- **Salidas:**
  - `01_otsu_mask.tif` — Máscara binaria 3D (opcional)
  - `02_cellpose_mask.tif` — Labels de núcleos 3D
  - `params.json` — Parámetros usados

### Parámetros Configurables
- `NUCLEUS_DIAMETER` — Diámetro esperado del núcleo en píxeles (típicamente 30)
- `CELLPOSE_USE_GPU` — Usar GPU para Cellpose (true/false)

---

## 03 — Filtrado de Núcleos de Astrocitos

**Objetivo:** Seleccionar candidatos astrocitarios basándose en señal de GFAP perinuclear y limpiar por tamaño físico.

### ⚠️ Nota Importante
El canal Iba-1 ha sido **eliminado completamente** del análisis debido a contaminación detectada con la señal de GFAP. El filtrado se basa **exclusivamente en GFAP**.

### Metodología

#### Paso 03a: Filtrado Biológico (GFAP Relativo)

**Estrategia Perinuclear:**
1. Para cada núcleo de Cellpose, se expande un "shell" (anillo) mediante dilatación morfológica
2. Se mide la intensidad media de GFAP en el shell (dilatado - núcleo original)
3. Se calcula el fondo de GFAP fuera de la máscara de Otsu (media ± desviación estándar)
4. **Regla de decisión:** Si `GFAP_shell > (media_fondo + N × std_fondo)` → aceptar como candidato

**Ventajas del umbral relativo:**
- Robusto ante variaciones en intensidad de adquisición entre preparados
- Se adapta automáticamente al nivel de señal de cada imagen
- Umbral científicamente justificable (desviaciones estándar sobre el fondo)

#### Paso 03b: Filtrado por Tamaño Físico

1. Calcula volumen mínimo en voxels según `MIN_VOLUME_UM3` y calibración (z,y,x)
2. Conserva solo objetos con área en voxels ≥ mínimo
3. Elimina núcleos muy pequeños (probablemente artefactos) o muy grandes (probablemente células fusionadas)

#### Métricas Generadas (`03_nucleus_metrics.csv`)

Por cada núcleo detectado:
- `label` — ID del núcleo
- `nucleus_volume_um3` — Volumen físico del núcleo (µm³)
- `nucleus_sphericity` — **Esfericidad 2D** (circularidad): 4π × area / perímetro² (rango 0-1)
  - 1.0 = círculo perfecto
  - 0.8-0.9 = núcleos redondeados típicos
  - 0.6-0.7 = núcleos alargados/irregulares
  - <0.6 = muy irregulares o fragmentados
- `shell_gfap_mean` — Intensidad media de GFAP en el shell perinuclear
- `is_astrocyte_candidate` — Boolean indicando si pasó el filtro GFAP

**Nota:** La esfericidad se calcula en 2D (proyección MIP del núcleo) debido a que los preparados tienen pocas slices Z, lo que hacía inviable el cálculo 3D tradicional.

#### Visualización
- **Contadores por etapa:**
  - #Núcleos Cellpose (total detectado)
  - #Candidatos GFAP (pasaron filtro biológico)
  - #Astrocitos Finales (pasaron filtro de tamaño)
- **Tasa de retención:** Porcentaje retenido en cada etapa
- **Histogramas:** Volumen y esfericidad coloreados por flag de retención
- **Scatter plot:** Espacio GFAP vs Volumen para visualizar umbral de decisión
- **Tabla detallada:** Todas las métricas por núcleo

### Entradas/Salidas
- **Entradas:**
  - `02_cellpose_mask.tif` — Labels de núcleos
  - Canal GFAP de la imagen raw
  - `01_otsu_mask.tif` — Para definir fondo
- **Salidas:**
  - `03_gfap_filtered_mask.tif` — Núcleos que pasaron filtro GFAP
  - `03_nucleus_metrics.csv` — Tabla de métricas por núcleo
  - `04_final_astrocytes_mask.tif` — Máscara final tras filtro de tamaño

### Parámetros Configurables

#### Filtro GFAP
- `SHELL_RADIUS_UM` — Radio del shell perinuclear en µm (típicamente 2.0)
- `GFAP_STD_DEV_THRESHOLD` — Número de desviaciones estándar sobre el fondo (típicamente 3.0)
- `MAX_DILATION_ITERATIONS` — Máximo de iteraciones de dilatación (seguridad, típicamente 20)

#### Filtro por Tamaño
- `MIN_VOLUME_UM3` — Volumen mínimo para considerar un núcleo válido (típicamente 75 µm³)

#### Fallback (Opcional)
- `GFAP_INTENSITY_THRESHOLD` — Umbral absoluto de intensidad GFAP (solo si falla el método relativo)

---

## 04 — Pipeline 2D Unificado: Esqueletización + Sholl

**Objetivo:** Generar esqueletos 2D por célula y calcular análisis de Sholl integrado en una sola ejecución optimizada.

### 🎯 Filosofía del Pipeline 2D

**¿Por qué 2D en vez de 3D?**

Los preparados de astrocitos en hipocampo son relativamente planos (~40-50 slices Z × 0.38 µm = ~15-20 µm de grosor). Los dominios territoriales astrocitarios se extienden principalmente en el plano XY. El enfoque 2D:

- ✅ Mantiene resolución XY completa (0.38 µm) sin degradación por remuestreo
- ✅ Sholl 2D nativo más preciso y eficiente
- ✅ Voronoi 2D simple y robusto para territorios
- ✅ ~10x más rápido que pipeline 3D isotrópico
- ✅ Científicamente justificado para morfología territorial en preparados planos

### Metodología: 5 Sub-pasos Integrados

#### [1] Proyección 2D de Volúmenes

**Máscaras de núcleos:**
- Proyección MIP (Maximum Intensity Projection) preservando labels
- Para cada posición (y,x), toma el label del slice Z con mayor área de ese label

**Señal GFAP:**
- Proyección MIP estándar (máximo por píxel en Z)
- Mantiene la intensidad máxima de señal por posición XY

**Resultado:** Imágenes 2D donde cada píxel representa el valor más intenso/relevante del stack 3D

#### [2] Partición Territorial con Voronoi

**Objetivo:** Definir territorios no solapados para cada astrocito

**Método:**
1. Centroides 2D de cada núcleo proyectado
2. Diagrama de Voronoi en el plano XY
3. Zona de exclusión (gap) entre territorios para evitar asignar regiones ambiguas

**Parámetros:**
- `TERRITORY_EXCLUSION_UM` — Gap de exclusión en µm (típicamente 1.0-2.0)
- Píxeles a distancia < gap de cualquier frontera Voronoi son excluidos (no asignados)

**Ventaja:** Evita entrelazamiento de procesos entre células vecinas, problema común en análisis 3D

#### [3] Esqueletización 2D por Territorio

**Por cada territorio:**
1. Máscara GFAP dentro del territorio Voronoi
2. Esqueletización 2D (algoritmo de adelgazamiento morfológico)
3. Labels preservados: cada esqueleto mantiene el ID de su célula

**Salida:** `05_skeleton_labels_2d.tif` — Imagen 2D con esqueletos etiquetados

#### [4] Conexión de Fragmentos (Opcional)

**Problema:** La esqueletización puede generar fragmentos desconectados por ruido o gaps en GFAP

**Solución:**
- Detecta fragmentos cercanos dentro del mismo territorio
- Conecta mediante línea recta si distancia < `CONNECTION_RADIUS_UM`
- Preserva topología (no conecta fragmentos de distintos territorios)

**Parámetros:**
- `CONNECT_SKELETON_FRAGMENTS` — Habilitar conexión (true/false)
- `CONNECTION_RADIUS_UM` — Radio máximo para conectar (típicamente 0.5 µm)

#### [5] Análisis con SKAN + Sholl 2D Nativo

**SKAN (Skeleton Analysis):**
- Construye grafo topológico del esqueleto
- Métricas: número de endpoints, junctions, ramas, longitud total
- Guarda en `skeletons/summary.csv`

**Sholl 2D Nativo Integrado:**
1. Anillos concéntricos desde el centroide nuclear en el plano 2D
2. Cuenta intersecciones del esqueleto con cada anillo (detección pixel-perfect)
3. Genera perfil completo: intersecciones vs radio

**Resultados:**
- `sholl_2d_native.csv` — Perfiles completos (radio × intersecciones por célula)
- `sholl_summary.csv` — Métricas agregadas:
  - `auc` — Área bajo la curva de Sholl (integración trapezoidal)
  - `peak_intersections` — Máximo de intersecciones
  - `critical_radius_um` — Radio donde ocurre el pico
- `sholl_rings_2d_native.json` — Coordenadas de anillos para visualización Napari

### Visualización en Napari

**Modo 2D (recomendado):**
- Proyección GFAP (fondo)
- `05_skeleton_labels_2d.tif` (esqueletos coloreados por célula)
- Anillos de Sholl (shapes) superpuestos desde JSON
- Escala física correcta aplicada

### Entradas/Salidas
- **Entradas:**
  - `04_final_astrocytes_mask.tif` — Máscara 3D de astrocitos finales
  - Canal GFAP 3D
  - Parámetros de calibración (z,y,x)
- **Salidas:**
  - `05_skeleton_labels_2d.tif` — Esqueletos 2D etiquetados
  - `skeletons/summary.csv` — Métricas topológicas SKAN
  - `sholl_2d_native.csv` — Perfiles de Sholl por célula
  - `sholl_summary.csv` — Métricas agregadas (AUC, pico, radio crítico)
  - `sholl_rings_2d_native.json` — Anillos para visualización

### Parámetros Configurables

#### Proyección
- `PROJECTION_2D_METHOD` — Método de proyección ('max', 'mean', 'sum')

#### Territorio Voronoi
- `TERRITORY_EXCLUSION_UM` — Gap de exclusión entre territorios (µm)

#### Conexión de Fragmentos
- `CONNECT_SKELETON_FRAGMENTS` — Habilitar (true/false)
- `CONNECTION_RADIUS_UM` — Radio máximo de conexión (µm)

#### Análisis de Sholl
- `SHOLL_MIN_UM` — Radio mínimo (típicamente 5.0 µm)
- `SHOLL_MAX_UM` — Radio máximo (típicamente 100.0 µm)
- `SHOLL_STEP_UM` — Separación entre anillos (típicamente 2.0 µm)

---

## 05 — Análisis por Preparado

**Objetivo:** Revisar métricas individuales de un preparado específico y visualizar resultados completos.

### Contenido de la Página

#### Dashboard de Conteo
Resumen rápido de cuántos objetos pasaron cada etapa:
- 📊 **Paso 02:** N° de núcleos detectados (Cellpose)
- 🧪 **Paso 03:** N° de candidatos GFAP
- ✅ **Paso 04:** N° de astrocitos finales
- 📈 **Paso 04 (Sholl):** N° de células con análisis Sholl completo

#### Métricas de Núcleo

**Fuente:** `03_nucleus_metrics.csv`

**Visualizaciones:**
- Tabla detallada con todas las métricas por núcleo
- Histograma de volumen nuclear (coloreado por flag de retención)
- Histograma de esfericidad 2D (coloreado por flag de retención)

**Métricas incluidas:**
- Volumen del núcleo (µm³)
- Esfericidad 2D (circularidad 0-1)
- Intensidad GFAP en shell perinuclear
- Flag de candidato a astrocito

#### Métricas de Sholl 2D

**Fuente:** `sholl_summary.csv` y `sholl_2d_native.csv`

**Visualizaciones:**
- Métricas clave: AUC mediana, pico mediano, radio crítico mediano
- Tabla detallada por célula
- **Perfiles de Sholl:** Gráfico de líneas con intersecciones vs radio por célula (interactivo)

#### Visualización Napari

**Modo 3D:**
- Carga máscaras 3D originales con escala física
- Ideal para verificar segmentación nuclear y filtrado

**Modo 2D (recomendado):**
- Proyección GFAP + esqueletos 2D + anillos de Sholl
- Fuentes: `05_skeleton_labels_2d.tif` + `sholl_rings_2d_native.json`
- Permite verificar territorios, esqueletos y análisis de Sholl superpuestos

### Archivos Leídos
- `02_cellpose_mask.tif`
- `03_nucleus_metrics.csv`
- `04_final_astrocytes_mask.tif`
- `05_skeleton_labels_2d.tif`
- `sholl_2d_native.csv`
- `sholl_summary.csv`
- `sholl_rings_2d_native.json`

### Uso
Esta página es ideal para:
- ✅ Verificar calidad de segmentación y filtrado
- ✅ Identificar preparados con problemas (ej: pocos astrocitos detectados)
- ✅ Explorar morfología individual de células
- ✅ Validar que Sholl captura la complejidad esperada

---

## 06 — Comparación entre Grupos

**Objetivo:** Comparar métricas agregadas entre grupos experimentales (CTL vs Hipoxia) con tests estadísticos rigurosos que evitan pseudoreplicación.

### Problema: Pseudoreplicación

**Definición:** Tratar células individuales como unidades independientes cuando provienen del mismo preparado.

**Consecuencia:** Infla artificialmente el tamaño muestral (N) y aumenta falsos positivos en tests estadísticos.

**Solución implementada:**
1. **DataFrame para gráficos (`df_plot`):** Incluye todas las células para visualizar distribuciones completas
2. **DataFrame para estadística (`df_stats`):** Calcula mediana por preparado, luego compara entre grupos
   - N = número de preparados (no número de células)
   - Preparado es la unidad experimental verdadera

### Metodología Estadística

#### Paso 1: Consolidación de Datos

Agrega datos de todos los preparados procesados:
- `sholl_summary.csv` → Métricas de Sholl por célula
- `03_nucleus_metrics.csv` → Métricas nucleares por célula
- Merge por `label`, `prepared`, `group`

#### Paso 2: Agregación por Preparado

Para cada métrica, calcula **mediana por preparado**:
```
Preparado_A_CTL → [célula1, célula2, ..., célulaM] → mediana_A
Preparado_B_CTL → [célula1, célula2, ..., célulaN] → mediana_B
...
Preparado_X_Hip → [célula1, célula2, ..., célulaK] → mediana_X
```

Resultado: 1 valor por preparado (evita pseudoreplicación)

#### Paso 3: Test de Normalidad

**Shapiro-Wilk** por grupo (CTL y Hipoxia):
- H0: Los datos provienen de una distribución normal
- α = 0.05
- Si p > 0.05 → distribución normal

#### Paso 4: Selección de Test

**Caso 1: Ambos grupos normales**
→ **Welch's t-test** (compara medias, no asume varianzas iguales)

**Caso 2: Al menos un grupo no normal**
→ **Mann-Whitney U** (compara medianas, no paramétrico)

#### Paso 5: Reporte de Resultados

**Salida incluye:**
- Test usado (Welch's t-test o Mann-Whitney U)
- Estadístico del test
- p-valor
- Interpretación (significativo si p < 0.05)
- Tamaños muestrales por grupo (N = preparados)
- Medias/medianas por grupo

### Métricas Disponibles para Comparación

#### Métricas de Sholl
- **AUC (Area Under Curve)** — Integración trapezoidal del perfil de Sholl
  - Unidades: µm·intersecciones
  - Interpreta: Complejidad dendrítica total
- **Pico de Intersecciones** — Máximo número de intersecciones alcanzado
  - Unidades: intersecciones
  - Interpreta: Máxima ramificación en un radio específico
- **Radio Crítico** — Radio donde ocurre el pico
  - Unidades: µm
  - Interpreta: Distancia de máxima arborización desde el soma

#### Métricas Nucleares
- **Volumen del Núcleo** — Volumen físico del núcleo
  - Unidades: µm³
  - Interpreta: Tamaño somático del astrocito
- **Esfericidad del Núcleo** — Circularidad 2D (4π·area/perímetro²)
  - Unidades: 0-1 (adimensional)
  - Interpreta: Regularidad morfológica nuclear (1.0 = círculo perfecto)

### Visualizaciones

#### Gráficos de Distribución
**Boxplots con ticks individuales:**
- Muestra distribución completa por célula (`df_plot`)
- Box = cuartiles, línea = mediana
- Ticks = células individuales coloreadas por grupo
- Permite ver variabilidad intra-preparado

#### Tabla de Datos

**Por Célula (`df_plot`):**
- Todas las células con todas las métricas
- Útil para análisis exploratorio detallado
- Botón de descarga CSV

**Por Preparado (`df_stats`):**
- Medianas por preparado (datos usados en test estadístico)
- N pequeño (típicamente 3-5 preparados por grupo)
- Botón de descarga CSV

### Interpretación de Resultados

**p < 0.05:** Diferencia estadísticamente significativa entre grupos
- Rechazamos H0 (grupos son iguales)
- Existe evidencia de diferencia real

**p ≥ 0.05:** No hay evidencia suficiente de diferencia
- No rechazamos H0
- Puede haber diferencia real pero sin poder estadístico suficiente

**Nota importante:** La significancia estadística no implica relevancia biológica. Evaluar también la magnitud del efecto (diferencia de medias/medianas).

### Archivos Leídos
Recorre todos los subdirectorios en `data/processed/`:
- `sholl_summary.csv` (por preparado)
- `03_nucleus_metrics.csv` (por preparado)
- Detecta automáticamente grupo (CTL/Hipoxia) por estructura de carpetas

### Uso
Esta página es ideal para:
- ✅ Comparar efectos de tratamiento/condición entre grupos
- ✅ Reportar resultados con rigor estadístico
- ✅ Evitar pseudoreplicación en análisis
- ✅ Exportar datos para análisis posteriores o publicación

---

## Idempotencia en Todas las Páginas

**Principio:** Cada paso verifica si su salida ya existe antes de recomputar.

**Implementación:**
- Si `02_cellpose_mask.tif` existe → reutiliza sin re-ejecutar Cellpose
- Si `sholl_summary.csv` existe → muestra resultados sin recalcular
- Botones específicos para **sobrescribir** cuando se necesita regenerar

**Ventajas:**
- ⚡ Navegación rápida entre páginas
- 💾 Ahorro de cómputo (Cellpose es costoso)
- 🔄 Permite ajustar parámetros downstream sin rehacer todo

**Excepciones:**
- Cambios en calibración (z,y,x) → suele requerir regenerar desde Paso 03
- Cambios en parámetros de filtrado → regenerar desde Paso 03
- Cambios en parámetros de Sholl → regenerar solo Paso 04

---

## Orquestación desde Home

La página `Home.py` actúa como dashboard y orquestador:

**Dashboard:**
- Contadores de preparados por grupo
- Estado de procesamiento (cuántos tienen cada paso completo)
- Archivos de salida por paso

**Ejecución Batch:**
- Permite ejecutar pipeline completo (pasos 01-04) en lote
- Opciones de ámbito:
  - **Preparado individual:** Solo la imagen seleccionada
  - **Grupo:** Todos los preparados de CTL o Hipoxia
  - **Todos:** Todo el dataset
- Control de sobrescritura por paso inicial
- Barra de progreso en tiempo real

**Lógica inteligente:**
- Detecta qué pasos faltan por preparado
- Ejecuta solo lo necesario (respeta idempotencia)
- Maneja errores por preparado sin detener el batch completo

---

## Notas Generales

### Persistencia de Configuración
- **Global:** `streamlit/calibration.json` — Parámetros compartidos por todo el experimento
- **Por preparado:** `data/processed/<preparado>/params.json` — Parámetros específicos usados

### Escalas Físicas
- Todas las métricas están en unidades físicas (µm, µm², µm³)
- Napari siempre se abre con `scale=(z,y,x)` aplicada
- Anillos de Sholl alineados correctamente con coordenadas físicas

### Nomenclatura de Archivos
Estandarizada para facilitar interoperabilidad:
- `01_` → Paso 01 (Otsu)
- `02_` → Paso 02 (Cellpose)
- `03_` → Paso 03 (Filtrado GFAP)
- `04_` → Paso 04 (Filtrado tamaño)
- `05_` → Paso 04 (Esqueletos 2D)
- `sholl_*` → Análisis de Sholl

### Troubleshooting Común

**Problema:** Pocos astrocitos detectados
- **Solución:** Reducir `GFAP_STD_DEV_THRESHOLD` o `MIN_VOLUME_UM3`

**Problema:** Territorios solapados
- **Solución:** Aumentar `TERRITORY_EXCLUSION_UM` (típico: 1-3 µm)

**Problema:** Esqueletos fragmentados
- **Solución:** Activar `CONNECT_SKELETON_FRAGMENTS=true` y ajustar `CONNECTION_RADIUS_UM`

**Problema:** Sholl sin pico claro
- **Solución:** Verificar que astrocitos tienen procesos extendidos (no solo soma), ajustar rango de radios

---

## Resumen de Flujo Completo

```
data/raw/*.tif
    ↓
[01] Calibración → calibration.json (z,y,x)
    ↓
[02] Cellpose → 02_cellpose_mask.tif (núcleos 3D)
    ↓
[03] Filtrado GFAP → 03_gfap_filtered_mask.tif + 03_nucleus_metrics.csv
                   → 04_final_astrocytes_mask.tif
    ↓
[04] Pipeline 2D → 05_skeleton_labels_2d.tif
                → sholl_2d_native.csv
                → sholl_summary.csv
    ↓
[05] Análisis Individual → Visualización por preparado
    ↓
[06] Comparación Grupos → CTL vs Hipoxia (estadística)
```

**Salida final:** Conjunto de métricas morfológicas por célula, comparables entre grupos experimentales, con validación estadística rigurosa.

## 01 — Calibración y Visualización de los Preparados

Objetivo: establecer la escala física global (µm) y verificarla sobre los archivos `.tif`/`.lif`.

- Detección de metadatos:
  - TIFF: se intenta leer OME-XML, metadatos ImageJ (`spacing`, `unit`) y tags de resolución.
  - LIF: se usa `readlif` y se toma `scale` como (Z,Y,X) en µm si está disponible.
- Conversión de unidades: nm → µm / mm → µm / m → µm según corresponda.
- Guardado: la calibración global (z,y,x) se persiste en `streamlit/calibration.json`, manteniendo otros parámetros globales.
- Visualización: se puede abrir Napari con la imagen y la escala aplicada.
- UI: contadores de archivos encontrados, métricas rápidas de calibración detectada por eje.

Entradas/Salidas:
- Entrada: imágenes en `data/raw/`.
- Salida: `streamlit/calibration.json` (global).

## 02 — Umbral de Otsu y Segmentación con Cellpose (DAPI)

Objetivo: generar máscaras de núcleos con segmentación 3D.

- Re-ordenamiento de ejes: todas las imágenes se llevan a (Z,C,Y,X) con heurísticas para `.lif`.
- Otsu: se calcula un umbral global sobre el volumen DAPI (opcional) y se guarda `01_otsu_mask.tif`.
- Cellpose: segmentación 3D (`do_3D=True`, `z_axis=0`) con `CellposeModel`.
- Limpieza opcional: el volumen DAPI puede enmascararse por Otsu antes de Cellpose.
- Persistencia de parámetros: umbral Otsu, diámetro, uso de GPU.
- Resumen/Visualización: métricas rápidas (umbral Otsu, fracción en máscara, #núcleos) y apertura en Napari con capas disponibles.

Entradas/Salidas:
- Entrada: DAPI (canal elegido) desde la imagen.
- Salidas: `01_otsu_mask.tif`, `02_cellpose_mask.tif`, `params.json`.

## 03 — Filtrado de Núcleos de Astrocitos (GFAP)

**⚠️ Nota:** El canal Iba-1 ha sido eliminado del análisis debido a contaminación detectada con la señal de GFAP que invalidaba la discriminación astrocito/microglia.

Objetivo: seleccionar candidatos astrocitarios y limpiar por tamaño físico.

- Estrategia perinuclear:
  - Para cada núcleo de Cellpose, se expande un anillo iterativo por dilatación.
  - Se mide la intensidad media de GFAP en el anillo (shell = dilatado - actual).
  - Regla: si GFAP > umbral → aceptar como candidato.
- Filtro de tamaño físico:
  - Se calcula el volumen mínimo en voxeles según `min_volume_um3` y la calibración (z,y,x).
  - Se conservan solo objetos con área en voxeles ≥ mínimo.
- Persistencia y resultados: `03_gfap_filtered_mask.tif` y `04_final_astrocytes_mask.tif`.
- Métricas y visualizaciones:
  - Contadores por etapa: #Cellpose, #Candidatos GFAP, #Final.
  - Tasa de retención (%) desde Cellpose a cada etapa.
  - Tabla por etiqueta con volúmenes físicos y flags de retención.

Entradas/Salidas:
- Entradas: `02_cellpose_mask.tif`, canal GFAP.
- Salidas: `03_gfap_filtered_mask.tif`, `04_final_astrocytes_mask.tif`, `params.json`.

## 04 — Esqueletización GFAP y Calibración para Sholl

Objetivo: generar esqueletos 3D por célula y calcular métricas centradas en el esqueleto.

- Remuestreo isotrópico: se fija `SKELETON_TARGET_ISO_UM` y se re-muestrea a vóxel cúbico.
- Umbral y conectividad:
  - Umbral GFAP por ROI (Otsu o Manual) y cierre morfológico opcional.
  - Selección del GFAP conectado a la semilla (núcleo dilatado) con conectividad 3D (6/26 equivalente a 1/3 en `label`).
- Restricciones espaciales:
  - Radio máximo físico desde el núcleo y padding que incluye dicho radio (evita truncar procesos largos).
  - Territorios Voronoi entre células con zona de exclusión (gap) para evitar entrelazados ambiguos.
  - Resolución de solapamientos de esqueletos por cercanía al núcleo (distancia en µm).
- Esqueletización y limpieza:
  - `skeletonize_3d` (o fallback 2D por slice) para un eje 1-vóxel de grosor.
  - Pruning topológico opcional: elimina ramas terminales más cortas que un umbral (µm) hasta convergencia.
- Métricas centradas en el esqueleto:
  - Longitud aproximada (vóxeles en esqueleto × tamaño de vóxel isotrópico).
  - Volúmenes: GFAP conectado, territorio Voronoi y volumen de dominio elegido.
  - Señal en “tubo” alrededor del esqueleto (radio en µm): intensidad total, media e intensidad por µm.
  - Grosor local por distancia al borde del GFAP conectado (media/mediana/p95).
  - Métricas de Skan (con refinamientos): re-afinamiento a 1 vóxel y filtro de ramas cortas, conteo de endpoints/junctions.
- Salidas y persistencia: `05_skeleton_labels.tif`, `skeletons/summary.csv`, actualización de `params.json` y `calibration.json` (SKELETON_*).
- Visualizaciones: tabla de métricas y gráficos (histogramas, correlaciones), incluyendo volumen de dominio e intensidad/µm.

Entradas/Salidas:
- Entradas: máscara final `04_final_astrocytes_mask.tif`, canal GFAP.
- Salidas: `05_skeleton_labels.tif`, `skeletons/summary.csv`, `params.json`.

### Parámetros SKELETON_* (detalle)

- `SKELETON_TARGET_ISO_UM` — tamaño de vóxel cúbico objetivo para remuestreo.
- `SKELETON_PADDING_UM` — margen físico agregado alrededor del ROI para evitar cortes.
- `SKELETON_SEED_DILATE_UM` — dilatación de la semilla (núcleo) para seleccionar GFAP conectado.
- `SKELETON_CONNECTIVITY` — conectividad 3D (6/18/26) para labeling/selección.
- `SKELETON_CLOSING_UM` — cierre morfológico previo al skeleton.
- `SKELETON_MAX_RADIUS_UM` — radio máximo desde el núcleo para delimitar el dominio por célula.
- `SKELETON_THRESHOLD_MODE` — `"otsu"` o `"manual"`.
- `SKELETON_MANUAL_THRESHOLD` — valor de umbral (si modo manual).
- `SKELETON_TERRITORY_VORONOI` — habilita partición Voronoi por cercanía al núcleo.
- `SKELETON_TERRITORY_EXCLUSION_UM` — gap excluyente entre celdas de Voronoi para reducir solapamientos.
- `SKELETON_DOMAIN_VOLUME_SOURCE` — fuente para volumen de dominio: `gfap`/`voronoi`/`final`.
- `SKELETON_PRUNE_ENABLE` — activa pruning topológico.
- `SKELETON_PRUNE_MIN_LEN_UM` — longitud mínima (µm) de ramas terminales a eliminar.
- `SKELETON_TUBE_RADIUS_UM` — radio (µm) del tubo para integrar señal GFAP alrededor del esqueleto.

## 05 — Test de Sholl y Parámetros

Objetivo: cuantificar la complejidad dendrítica por célula mediante intersecciones del esqueleto con anillos concéntricos.

- Centros por célula: centroides de `02_cellpose_mask.tif`, restringidos a `04_final_astrocytes_mask.tif` cuando corresponde.
- Cálculo 3D: distancias físicas (µm) en la grilla (Z,Y,X); cascarón de ancho `step/2` alrededor de cada radio.
- Intersecciones: número de componentes conectados del esqueleto dentro del cascarón (aproximación robusta en 3D).
- Resultados: `sholl.csv` por preparado y `sholl_rings.json` para visualización.
- Visualización en Napari: anillos por célula en el plano Z más cercano al centro. La capa Shapes usa la misma `scale` de las imágenes; anillos alineados y con color amarillo.

Parámetros SHOLL_*:
- `SHOLL_MIN_RADIUS_UM` — radio mínimo (µm).
- `SHOLL_MAX_RADIUS_UM` — radio máximo (µm).
- `SHOLL_STEP_UM` — separación entre anillos (µm) y semiancho del cascarón.

## 05 — Análisis por Preparado

Objetivo: Revisar métricas individuales de un preparado específico y visualizar resultados completos.

- **Resumen de conteo:** Muestra N° de núcleos detectados (02), candidatos GFAP (03), astrocitos finales (04), y células con análisis Sholl completo (04).
- **Métricas de núcleo:** Tabla detallada de `03_nucleus_metrics.csv` con volumen, esfericidad, flags de retención.
- **Métricas de Sholl 2D:** Tabla de `sholl_summary.csv` con AUC, pico de intersecciones, radio crítico por célula. Gráficos de perfiles de Sholl desde `sholl_2d_native.csv`.
- **Visualización Napari:**
  - **3D:** Carga máscaras 3D originales con escala física.
  - **2D:** Carga proyección + esqueletos 2D + anillos de Sholl superpuestos desde `05_skeleton_labels_2d.tif` y `sholl_rings_2d_native.json`.

Entradas:
- `02_cellpose_mask.tif`, `03_nucleus_metrics.csv`, `04_final_astrocytes_mask.tif`
- `05_skeleton_labels_2d.tif`, `sholl_2d_native.csv`, `sholl_summary.csv`, `sholl_rings_2d_native.json`

## 06 — Comparación entre Grupos

Objetivo: Comparar métricas agregadas entre grupos experimentales (CTL vs Hipoxia) con tests estadísticos rigurosos.

- **Consolidación de datos:** Agrega `sholl_summary.csv` y `03_nucleus_metrics.csv` de todos los preparados.
- **Solución a pseudoreplicación:** 
  - DataFrame para gráficos: incluye todas las células individuales para visualización de distribuciones.
  - DataFrame para estadística: calcula mediana por preparado, luego compara entre grupos (N = preparados, no células).
- **Test estadístico automático:**
  - Shapiro-Wilk para normalidad.
  - Si ambos grupos normales → Welch's t-test (compara medias).
  - Si alguno no normal → Mann-Whitney U (compara medianas).
- **Visualizaciones:** Boxplots + ticks individuales con código de colores por grupo.
- **Métricas disponibles:**
  - Sholl: AUC, pico de intersecciones, radio crítico
  - Núcleo: volumen, esfericidad
- **Exportación:** Botones de descarga CSV para datos por célula y por preparado.

Metodología estadística:
- Evita pseudoreplicación tratando células como unidad experimental.
- Usa preparado como unidad de muestreo (N por grupo = N° de preparados).
- Test de normalidad guía selección de prueba paramétrica/no paramétrica.

## Idempotencia en todas las páginas

Cada paso verifica si su salida ya existe y la reutiliza (evita recomputar). Si falta alguna dependencia, la página genera automáticamente lo necesario en el orden correcto. Este comportamiento también lo implementa el orquestador de `Home` (ejecuta 01–04 según lo faltante por preparado).

## Notas generales

- Todas las páginas respetan `streamlit/calibration.json` y utilizan los parámetros guardados por preparado.
- La apertura de Napari está disponible en cada página con capas y escala física coherentes.
- Los nombres de salida están estandarizados para facilitar el traspaso entre etapas.
