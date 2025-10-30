# Astrocitos 3D — Pipeline de Análisis Morfológico

Análisis de astrocitos (GFAP) y núcleos (DAPI) en 3D a partir de imágenes `.tif`, con una app de Streamlit como flujo principal. Incluye integración con Napari, segmentación con Cellpose, filtrado biológico y esqueletización 3D con métricas avanzadas centradas en el esqueleto.

Tecnologías: Python, Streamlit, Napari, Cellpose, scikit-image, SciPy, Skan, Altair.

## Estado actual

- La app de Streamlit es el flujo principal y está actualizada con contadores y visualizaciones por etapa.
- Configuración unificada en `streamlit/calibration.json` (z,y,x y parámetros SKELETON_* y SHOLL_*).
- Resultados por preparado en `data/processed/<preparado>/` y un maestro global en `data/processed/master_morphology_results.csv`.
- `results/` dejó de usarse; toda la salida vive bajo `data/processed/`.

## Estructura

- `streamlit/` — App multipágina (Home, calibración, Otsu/Cellpose, filtrado, esqueletización, Sholl y resumen). Ver `streamlit/README.md`.
- `streamlit/pages/` — Páginas del pipeline. Ver metodología detallada en `streamlit/pages/README.md`.
- `data/raw/` — Dataset crudo (.tif) organizado por carpetas de grupo/preparado.
- `data/processed/<preparado>/` — Salidas por imagen preparada.
- `requirements.txt` — Dependencias unificadas del proyecto.

## Instalación rápida

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Ejecutar la app

```bash
cd streamlit
streamlit run Home.py
```

## Flujo de trabajo

1) Calibración física global (µm): define `z`, `y`, `x` en `streamlit/calibration.json` y/o abre Napari con escala.

2) Otsu y Cellpose (DAPI): umbral Otsu opcional y segmentación 3D con Cellpose.

3) Filtrado GFAP/Microglía: evalúa un anillo perinuclear iterativo (máx. iteraciones) y descarta candidatos por microglía alta. Aplica filtro de tamaño físico (µm³).

4) Esqueletización 3D por célula:
   - Re-muestreo a vóxel isotrópico objetivo; umbral GFAP (Otsu/Manual) y cierre morfológico.
   - Selección del GFAP conectado a la semilla (núcleo dilatado) con conectividad 3D.
   - Radio máximo desde el núcleo y padding que incluye dicho radio para evitar truncamientos.
   - Opcionales: territorios Voronoi por cercanía al núcleo con “gap” en fronteras; resolución de solapamientos por distancia física al núcleo.
   - Pruning topológico (opcional): elimina espículas cortas en el esqueleto.
   - Métrica centrada en el esqueleto: integra la señal GFAP dentro de un “tubo” de radio fijo (µm) alrededor del esqueleto y reporta intensidad por unidad de longitud.

5) Sholl por célula: anillos concéntricos alrededor de cada núcleo/astrocito y conteo de intersecciones por radio. Los anillos se guardan en JSON y se visualizan en Napari.

6) Resumen y estadísticas: tablas por imagen y comparaciones por grupo (p.ej., CTL vs Hipoxia) usando pruebas de Mann–Whitney U para métricas clave (longitud, endpoints, pico de Sholl, etc.).

Idempotencia: cada paso verifica si su salida existe y la reutiliza para evitar recomputar innecesariamente.

## Salidas por preparado (nombres estándar)

- `01_otsu_mask.tif` — Máscara Otsu (DAPI) opcional.
- `02_cellpose_mask.tif` — Etiquetas de núcleos de Cellpose.
- `03_gfap_microglia_filtered_mask.tif` — Candidatos GFAP tras filtrado biológico.
- `04_final_astrocytes_mask.tif` — Máscara final tras filtro de tamaño.
- `05_skeleton_labels.tif` — Esqueletos por etiqueta.
- `skeletons/summary.csv` — Métricas por célula (longitudes, volúmenes, señal en tubo, radios locales, etc.).
- `sholl.csv` y `sholl_rings.json` — Intersecciones por radio y geometría de anillos para Napari.
- `params.json` — Parámetros usados por ese preparado.

## Troubleshooting

- Anillos de Sholl desplazados: corregido añadiendo `scale=scale` a la capa de Shapes en Napari; asegúrate de actualizar el código.
- El esqueleto no se extiende pese a subir el radio: confirma la calibración (µm) y que el padding incluye el radio máximo (ya implementado).
- Procesos entrelazados: activa los territorios Voronoi y ajustá la zona de exclusión; en solapamientos directos, la asignación al núcleo más cercano resuelve conflictos.
- Conteo de ramas inflado: el pipeline re-afinó el esqueleto a 1 vóxel y filtra ramas cortas antes de Skan.
- Performance: usa isotrópico objetivo razonable y radios moderados; Napari puede abrirse desde cada página con escala correcta.

## Licencia

MIT (o actualizar según corresponda).
