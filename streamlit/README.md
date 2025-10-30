# App de Streamlit — Pipeline de Astrocitos 3D

La app implementa el pipeline completo: calibración, segmentación, filtrado, esqueletización, Sholl por célula y resumen estadístico. Cada página abre Napari con capas y escala física coherentes.

## Estructura

- `Home.py` — Portada, orquestador del pipeline 01–05 para preparado/grupo/todo; inspector/edición de `calibration.json`; dashboards y resúmenes.
- `pages/` — Páginas del pipeline (01–06). Ver metodología detallada en `pages/README.md`.
- `ui/` — Utilidades (pipeline básico, esqueletización avanzada, Sholl, sidebar, etc.).
- `napari_viewer.py` — Visor estable con plugins deshabilitados; todas las capas usan `scale=(z,y,x)` y los anillos de Sholl se dibujan alineados.
- `calibration.json` — Configuración global unificada (z,y,x, `SKELETON_*`, `SHOLL_*`).

## Flujo dentro de la app

1) Calibración global (página 01) — fija `z,y,x` en µm.
2) Otsu/Cellpose (página 02) — segmentación 3D de núcleos.
3) Filtrado GFAP/Microglía (página 03) — candidatos astrocitarios y filtro de tamaño físico.
4) Esqueletización GFAP (página 04) — esqueletos 3D por célula, métricas centradas en el esqueleto y Skan.
5) Test de Sholl (página 05) — anillos concéntricos por célula y conteo de intersecciones; exporta `sholl.csv` y `sholl_rings.json`.
6) Resumen/Comparaciones (página 06) — tablas por imagen y pruebas Mann–Whitney U por grupo.

Idempotencia: cada paso verifica si la salida esperada existe y la reutiliza; si falta, la genera. Esto aplica tanto en páginas como en el orquestador de `Home`.

## Parámetros y configuración

- `calibration.json` incluye:
  - `z`, `y`, `x` (µm) — escala física global.
  - `MAX_DILATION_ITERATIONS`, `GFAP_INTENSITY_THRESHOLD`, `MICROGLIA_INTENSITY_THRESHOLD`, `MIN_VOLUME_UM3`.
  - `SKELETON_*`: `TARGET_ISO_UM`, `PADDING_UM`, `SEED_DILATE_UM`, `CONNECTIVITY`, `CLOSING_UM`, `MAX_RADIUS_UM`, `THRESHOLD_MODE`, `MANUAL_THRESHOLD`, `TERRITORY_VORONOI`, `TERRITORY_EXCLUSION_UM`, `DOMAIN_VOLUME_SOURCE`, `PRUNE_ENABLE`, `PRUNE_MIN_LEN_UM`, `TUBE_RADIUS_UM`.
  - `SHOLL_*`: `MIN_RADIUS_UM`, `MAX_RADIUS_UM`, `STEP_UM`.
- `data/processed/<preparado>/params.json` guarda los parámetros usados por ese preparado cuando se elige origen "Por preparado".

## Salidas por preparado

- 01_otsu_mask.tif, 02_cellpose_mask.tif, 03_gfap_microglia_filtered_mask.tif, 04_final_astrocytes_mask.tif
- 05_skeleton_labels.tif
- skeletons/summary.csv — métricas por célula
- sholl.csv, sholl_rings.json
- params.json

## Tips y notas

- Anillos de Sholl desplazados: ya corregido agregando `scale=scale` a la capa Shapes; actualizar si fuera necesario.
- Ante entrelazados, activar `SKELETON_TERRITORY_VORONOI` y ajustar `SKELETON_TERRITORY_EXCLUSION_UM`.
- El pruning (`SKELETON_PRUNE_ENABLE` y `SKELETON_PRUNE_MIN_LEN_UM`) ayuda a estabilizar métricas.
- Skan `summarize` se invoca con `separator="_"` para evitar warnings y nombres consistentes de columnas.
