# App de Streamlit — Pipeline de Astrocitos 3D

La app implementa el pipeline completo: calibración, segmentación, filtrado, esqueletización, Sholl por célula y resumen estadístico. Cada página abre Napari con capas y escala física coherentes.

## ⚠️ Cambio Importante: Eliminación del Canal Iba-1

El canal Iba-1 (originalmente usado para identificar microglía) ha sido **completamente eliminado del análisis** debido a contaminación detectada con la señal del canal GFAP. El filtrado de candidatos a astrocitos ahora se basa **exclusivamente en la señal de GFAP** alrededor de los núcleos.

Este cambio afecta principalmente a la **Página 03 (Filtrado)**, donde ya no se configura ni se utiliza el umbral de Iba-1. Para más detalles, ver la sección correspondiente en el README.md principal del proyecto.

## Estructura

- `Home.py` — Portada, orquestador del pipeline 01–05 para preparado/grupo/todo; inspector/edición de `calibration.json`; dashboards y resúmenes.
- `pages/` — Páginas del pipeline (01–06). Ver metodología detallada en `pages/README.md`.
- `ui/` — Utilidades (pipeline básico, esqueletización avanzada, Sholl, sidebar, etc.).
- `napari_viewer.py` — Visor estable con plugins deshabilitados; todas las capas usan `scale=(z,y,x)` y los anillos de Sholl se dibujan alineados.
- `calibration.json` — Configuración global unificada (z,y,x, `SKELETON_*`, `SHOLL_*`).

## Flujo dentro de la app

1. **Página Home**: Ejecuta todos los pasos en lote para un preparado (4 pasos unificados)
2. **Páginas 01-04**: Ejecución manual paso a paso con visualización interactiva
3. **Página 06**: Resumen por imagen y comparaciones
4. **Página 07**: Comparaciones globales entre grupos CTL/hip

**Idempotencia**: Cada paso verifica si la salida esperada existe y la reutiliza; si falta, la genera. Esto aplica tanto en páginas como en el orquestador de `Home`.

## Parámetros y configuración

- `calibration.json` incluye:
  - `z`, `y`, `x` (µm) — escala física global.
  - `DAPI_CHANNEL_INDEX`, `GFAP_CHANNEL_INDEX` — índices de canales en LIF (0 y 1 respectivamente).
  - `NUCLEUS_DIAMETER` — diámetro esperado de núcleos para Cellpose.
  - `MAX_DILATION_ITERATIONS`, `GFAP_INTENSITY_THRESHOLD`, `MIN_VOLUME_UM3` — filtrado de núcleos.
  - `GFAP_STD_DEV_THRESHOLD`, `SHELL_RADIUS_UM` — análisis de intensidad GFAP.
  - **Pipeline 2D nativo**:
    - `PROJECTION_2D_METHOD` — método de proyección ("max", "mean", "sum").
    - `TERRITORY_EXCLUSION_UM` — gap de exclusión en territorios Voronoi (µm).
    - `CONNECT_SKELETON_FRAGMENTS` — activar conexión de fragmentos GFAP discontinuos.
    - `CONNECTION_RADIUS_UM` — radio de conexión para fragmentos (µm).
  - `SHOLL_*`: `MIN_RADIUS_UM`, `MAX_RADIUS_UM`, `STEP_UM` — parámetros de Sholl.
- `data/processed/<preparado>/params.json` guarda los parámetros usados por ese preparado cuando se elige origen "Por preparado".

## Salidas por preparado

- `01_otsu_mask.tif`, `02_cellpose_mask.tif`, `03_gfap_filtered_mask.tif`, `04_final_astrocytes_mask.tif`
- `03_nucleus_metrics.csv` — Métricas detalladas por núcleo (volumen, intensidad GFAP, etc.)
- `05_skeleton_labels_2d.tif` — Labels de esqueletos 2D proyectados
- `skeletons/summary.csv` — Métricas por célula (obsoleto, reemplazado por salida inline)
- `sholl_2d_native.csv` — Perfiles de Sholl por célula (2D nativo con skan)
- `sholl_2d_native_detailed.csv` — Intersecciones completas por radio
- `sholl_rings_2d_native.json` — Coordenadas de anillos para visualización
- `sholl_summary.csv` — Radio crítico, pico, AUC por célula
- `params.json` — Parámetros usados para este preparado

## Tips y notas

- **Pipeline 2D nativo**: Los astrocitos en preparaciones de hipocampo son relativamente planos. El flujo 2D nativo mantiene la resolución XY completa (0.38 µm), es ~10x más rápido, y calcula Sholl directamente sobre el esqueleto 2D.
- **Territorios Voronoi**: `TERRITORY_EXCLUSION_UM` controla el gap de exclusión alrededor de cada núcleo antes de calcular territorios. Valores típicos: 2-5 µm.
- **Conexión de fragmentos**: Si `CONNECT_SKELETON_FRAGMENTS=true`, el pipeline conecta fragmentos cercanos de GFAP usando dilatación morfológica con radio `CONNECTION_RADIUS_UM`.
- Skan `summarize` se invoca con `separator="_"` para evitar warnings y nombres consistentes de columnas.
