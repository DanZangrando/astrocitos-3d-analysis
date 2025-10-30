# Páginas del Pipeline — Metodología Detallada

Este documento describe en detalle la metodología aplicada en cada página de la app.

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

## 03 — Filtrado de Núcleos de Astrocitos (GFAP/Microglía)

Objetivo: seleccionar candidatos astrocitarios y limpiar por tamaño físico.

- Estrategia perinuclear:
  - Para cada núcleo de Cellpose, se expande un anillo iterativo por dilatación.
  - Se mide la intensidad media de GFAP y de microglía en el anillo (shell = dilatado - actual).
  - Regla: si microglía > umbral → descartar; si GFAP > umbral → aceptar como candidato.
- Filtro de tamaño físico:
  - Se calcula el volumen mínimo en voxeles según `min_volume_um3` y la calibración (z,y,x).
  - Se conservan solo objetos con área en voxeles ≥ mínimo.
- Persistencia y resultados: `03_gfap_microglia_filtered_mask.tif` y `04_final_astrocytes_mask.tif`.
- Métricas y visualizaciones:
  - Contadores por etapa: #Cellpose, #Candidatos GFAP/Microglía, #Final.
  - Tasa de retención (%) desde Cellpose a cada etapa.
  - Tabla por etiqueta con volúmenes físicos y flags de retención.

Entradas/Salidas:
- Entradas: `02_cellpose_mask.tif`, canales GFAP y Microglía.
- Salidas: `03_gfap_microglia_filtered_mask.tif`, `04_final_astrocytes_mask.tif`, `params.json`.

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

## 06 — Resumen por imagen y Comparaciones

Objetivo: consolidar métricas por imagen y comparar grupos (p.ej., CTL vs Hipoxia).

- Consolidación: agrega `skeletons/summary.csv` y picos de Sholl a `data/processed/master_morphology_results.csv`.
- Comparaciones: pruebas no paramétricas de Mann–Whitney U sobre métricas clave (longitud total, #endpoints, pico de Sholl, etc.).
- Gráficos: distribuciones por grupo y anotaciones de p‑valor.

## Idempotencia en todas las páginas

Cada paso verifica si su salida ya existe y la reutiliza (evita recomputar). Si falta alguna dependencia, la página genera automáticamente lo necesario en el orden correcto. Este comportamiento también lo implementa el orquestador de `Home` (ejecuta 01–05 según lo faltante por preparado).

## Notas generales

- Todas las páginas respetan `streamlit/calibration.json` y utilizan los parámetros guardados por preparado.
- La apertura de Napari está disponible en cada página con capas y escala física coherentes.
- Los nombres de salida están estandarizados para facilitar el traspaso entre etapas.
