# Astrocitos 3D — Pipeline de Análisis Morfológico

Pipeline automatizado para análisis morfológico cuantitativo de astrocitos (GFAP) y núcleos (DAPI) a partir de imágenes de microscopía 3D. Implementa un **flujo de trabajo 2D nativo optimizado** para preparados planos, con interfaz Streamlit interactiva, segmentación con Cellpose, esqueletización territorial y análisis de Sholl integrado.

**Tecnologías:** Python 3.12, Streamlit, Napari, Cellpose, scikit-image, SciPy, Skan, Pandas, Altair

---

## 🎯 Características Principales

### Pipeline 2D Nativo Optimizado
- ✅ **Resolución XY completa** (0.38 µm) sin degradación por remuestreo isotrópico
- ✅ **Sholl 2D nativo** más preciso y eficiente que aproximaciones 3D
- ✅ **Territorios Voronoi 2D** robustos para partición astrocitaria
- ✅ **~10x más rápido** que flujos 3D con remuestreo isotrópico
- ✅ **Científicamente justificado** para morfología territorial en preparados relativamente planos

### Flujo Completo (4 Pasos + 2 Análisis)

**Pasos de Procesamiento:**
1. **Calibración Física** — Detección automática de metadatos y definición de escala (µm)
2. **Segmentación Nuclear** — Cellpose 3D para detectar núcleos (DAPI)
3. **Filtrado de Astrocitos** — Selección por señal GFAP perinuclear (umbrales relativos) + tamaño físico
4. **Pipeline 2D Unificado** — Proyección MIP → Voronoi → Esqueletización → Sholl (integrado)

**Páginas de Análisis:**
5. **Análisis por Preparado** — Visualización individual con métricas completas y Napari
6. **Comparación entre Grupos** — Estadística CTL vs Hipoxia (Mann-Whitney U / Welch's t-test)

### Métricas Cuantitativas

**Métricas Nucleares:**
- Volumen del núcleo (µm³)
- Esfericidad 2D (circularidad: 4π·area/perímetro²)
- Intensidad GFAP en shell perinuclear

**Métricas de Sholl:**
- AUC (área bajo la curva) — Complejidad dendrítica total
- Pico de intersecciones — Máxima ramificación
- Radio crítico (µm) — Distancia de máxima arborización

**Métricas Topológicas (SKAN):**
- Longitud total del esqueleto (µm)
- Número de endpoints y junctions
- Distribución de longitudes de ramas

---

## ⚠️ Nota Importante: Canal Iba-1 Eliminado

**Problema detectado:** Durante el procesamiento se identificó contaminación del canal Iba-1 (marcador de microglía) con señal del canal GFAP. Las estructuras brillantes de astrocitos aparecían replicadas en el canal Iba-1, invalidando su uso para discriminación celular.

**Decisión:** El canal Iba-1 ha sido **completamente eliminado** del pipeline. El filtrado de candidatos a astrocitos se basa **exclusivamente en la señal de GFAP** en el shell perinuclear, usando umbrales relativos (desviaciones estándar sobre el fondo).

**Impacto:** Esto puede incluir algunos núcleos de microglía con señal GFAP basal, pero es preferible a descartar astrocitos verdaderos. El enfoque es ahora 100% basado en morfología astrocitaria (GFAP + tamaño nuclear).

---

## 📁 Estructura del Proyecto

```
astrocitos-3d-analysis/
├── README.md                          # Este archivo
├── requirements.txt                   # Dependencias Python
├── CAMBIOS_ESFERICIDAD_IBA1.md       # Documentación de cambios recientes
│
├── data/
│   ├── raw/                          # Imágenes originales .tif/.lif
│   │   ├── CTL/                      # Grupo control
│   │   └── hip/                      # Grupo hipoxia
│   │
│   └── processed/                    # Resultados por preparado
│       ├── master_morphology_results.csv
│       └── <preparado>/
│           ├── 01_otsu_mask.tif
│           ├── 02_cellpose_mask.tif
│           ├── 03_gfap_filtered_mask.tif
│           ├── 03_nucleus_metrics.csv
│           ├── 04_final_astrocytes_mask.tif
│           ├── 05_skeleton_labels_2d.tif
│           ├── sholl_2d_native.csv
│           ├── sholl_summary.csv
│           ├── sholl_rings_2d_native.json
│           ├── params.json
│           └── skeletons/
│               └── summary.csv
│
├── streamlit/                        # Aplicación Streamlit
│   ├── Home.py                       # Dashboard principal + orquestador batch
│   ├── calibration.json              # Configuración global
│   ├── napari_viewer.py              # Launcher Napari 3D
│   ├── napari_viewer_2d.py           # Launcher Napari 2D
│   │
│   ├── pages/                        # Páginas del pipeline
│   │   ├── README.md                 # Metodología detallada por página
│   │   ├── 01_ Calibración y Visualización de los Preparados.py
│   │   ├── 02_ Umbral de Otsu y Segmentación con Cellpose.py
│   │   ├── 03_ Filtrado de Nucleos de Astrocitos.py
│   │   ├── 04_ Esqueletización GFAP y Calibración para Sholl.py
│   │   ├── 05_ Análisis por Preparado.py
│   │   └── 06_ Comparación entre Grupos.py
│   │
│   └── ui/                           # Módulos de lógica
│       ├── __init__.py
│       ├── pipeline.py               # Pasos 01-03 + redirección Paso 04
│       ├── pipeline_2d_unified.py    # Implementación Pipeline 2D (Paso 04)
│       ├── runner.py                 # Orquestación y ejecución batch
│       ├── sidebar.py                # Configuración global (UI)
│       ├── plots.py                  # Utilidades de visualización
│       └── utils.py                  # Funciones auxiliares
│
└── test_*.py                         # Scripts de prueba y diagnóstico
```

---

## 🚀 Instalación

### Requisitos Previos
- Python 3.12 (recomendado) o Python 3.10+
- CUDA Toolkit (opcional, para GPU en Cellpose)
- 16+ GB RAM (recomendado para volúmenes grandes)

### Instalación con virtualenv

```bash
# Clonar repositorio
git clone https://github.com/DanZangrando/astrocitos-3d-analysis.git
cd astrocitos-3d-analysis

# Crear entorno virtual
python3.12 -m venv venv
source venv/bin/activate  # Linux/Mac
# o: venv\Scripts\activate  # Windows

# Instalar dependencias
pip install --upgrade pip
pip install -r requirements.txt
```

### Verificar Instalación

```bash
python -c "import streamlit, napari, cellpose, skan; print('✅ Todas las dependencias instaladas')"
```

---

## 🎮 Uso Básico

### 1. Ejecutar la Aplicación Streamlit

```bash
cd streamlit
streamlit run Home.py
```

La aplicación se abrirá automáticamente en `http://localhost:8501`

### 2. Flujo de Trabajo Recomendado

#### A. Configuración Inicial (Solo una vez)

1. **Página 01 — Calibración:**
   - Cargar un preparado representativo
   - Verificar que la calibración física (z,y,x en µm) se detecte correctamente
   - Guardar calibración global
   - Configurar índices de canales (DAPI=0, GFAP=1 típicamente)

2. **Configurar Parámetros en Sidebar:**
   - Expander "🔬 Segmentación Nuclear": Ajustar `NUCLEUS_DIAMETER` (típicamente 30 píxeles)
   - Expander "🧪 Filtrado de Astrocitos": 
     - `GFAP_STD_DEV_THRESHOLD` = 3.0 (N° desviaciones estándar sobre fondo)
     - `MIN_VOLUME_UM3` = 75 (volumen mínimo nuclear)
   - Expander "🗺️ Pipeline 2D":
     - `TERRITORY_EXCLUSION_UM` = 1.0-2.0 (gap Voronoi)
     - `SHOLL_MIN/MAX/STEP_UM` = 5.0, 100.0, 2.0 (rango de radios)
   - Guardar configuraciones

#### B. Procesamiento Batch

**Opción 1: Desde Home.py (Recomendado)**
1. Sidebar → Expander "▶️ Ejecución Batch"
2. Seleccionar ámbito: Todos / CTL / Hipoxia
3. Paso inicial: 01 (todo desde cero) o 03 (con núcleos ya segmentados)
4. Marcar "Sobrescribir" si necesario
5. Click "Ejecutar pipeline"
6. Esperar barra de progreso (puede tomar varios minutos por preparado)

**Opción 2: Página por Página (Para ajustar parámetros)**
1. Ir a cada página (01-04) secuencialmente
2. Cargar un preparado de prueba
3. Ejecutar paso individual
4. Verificar resultados
5. Ajustar parámetros según sea necesario
6. Repetir para otros preparados o usar batch

#### C. Análisis de Resultados

**Página 05 — Análisis por Preparado:**
- Revisar métricas nucleares (volumen, esfericidad)
- Visualizar perfiles de Sholl por célula
- Abrir en Napari 2D para verificar esqueletos y territorios
- Identificar preparados con problemas

**Página 06 — Comparación entre Grupos:**
- Seleccionar métrica a comparar (AUC, pico, radio crítico, volumen nuclear, esfericidad)
- Ver distribución por célula (boxplots)
- Ver tabla por preparado (medianas agregadas)
- Leer resultados del test estadístico (p-valor, significancia)
- Exportar datos CSV para análisis posteriores

---

## 📊 Métricas Disponibles

### Métricas Nucleares (`03_nucleus_metrics.csv`)

| Métrica | Unidad | Descripción | Interpretación |
|---------|--------|-------------|----------------|
| `nucleus_volume_um3` | µm³ | Volumen físico del núcleo | Tamaño somático del astrocito |
| `nucleus_sphericity` | 0-1 | Circularidad 2D: 4π·area/perímetro² | 1.0 = círculo perfecto<br>0.8-0.9 = redondeados típicos<br><0.6 = irregulares |
| `shell_gfap_mean` | a.u. | Intensidad media GFAP en shell perinuclear | Nivel de expresión GFAP |
| `is_astrocyte_candidate` | boolean | Pasó filtro GFAP relativo | true = candidato a astrocito |

**Nota sobre esfericidad:** Se calcula en 2D (proyección MIP) debido a limitaciones de resolución Z en los preparados.

### Métricas de Sholl (`sholl_summary.csv`)

| Métrica | Unidad | Descripción | Interpretación |
|---------|--------|-------------|----------------|
| `auc` | µm·intersec | Área bajo la curva de Sholl (integración trapezoidal) | Complejidad dendrítica total |
| `peak_intersections` | intersec | Máximo número de intersecciones | Máxima ramificación alcanzada |
| `critical_radius_um` | µm | Radio donde ocurre el pico | Distancia de máxima arborización desde soma |

**Método:** Sholl 2D nativo con anillos concéntricos desde el centroide nuclear, conteo pixel-perfect de intersecciones con el esqueleto.

### Métricas Topológicas (`skeletons/summary.csv`)

| Métrica | Unidad | Descripción |
|---------|--------|-------------|
| `skeleton_length_um` | µm | Longitud total del esqueleto |
| `n_endpoints` | count | Número de terminaciones |
| `n_junctions` | count | Número de bifurcaciones |
| `mean_branch_length_um` | µm | Longitud media de ramas |

**Método:** Análisis con SKAN (Skeleton Analysis) sobre esqueletos 2D.

---

## 🔧 Configuración Avanzada

### Parámetros Globales (`streamlit/calibration.json`)

#### Calibración Espacial
```json
{
  "z": 0.38,              // Resolución Z en µm
  "y": 0.379,             // Resolución Y en µm
  "x": 0.379,             // Resolución X en µm
  "DAPI_CHANNEL_INDEX": 0,
  "GFAP_CHANNEL_INDEX": 1
}
```

#### Segmentación Nuclear (Paso 02)
```json
{
  "NUCLEUS_DIAMETER": 30,      // Diámetro de núcleo en píxeles
  "CELLPOSE_USE_GPU": true     // Usar GPU en Cellpose
}
```

#### Filtrado de Astrocitos (Paso 03)
```json
{
  "SHELL_RADIUS_UM": 2.0,              // Radio del shell perinuclear
  "GFAP_STD_DEV_THRESHOLD": 3.0,       // N° std sobre fondo
  "MAX_DILATION_ITERATIONS": 20,       // Máx iteraciones dilatación
  "MIN_VOLUME_UM3": 75,                // Volumen mínimo nuclear
  "GFAP_INTENSITY_THRESHOLD": 50       // Umbral absoluto (fallback)
}
```

**Recomendaciones:**
- `GFAP_STD_DEV_THRESHOLD` entre 2.0-4.0 (más bajo = más permisivo)
- `MIN_VOLUME_UM3` ajustar según tamaño esperado de núcleos astrocitarios

#### Pipeline 2D (Paso 04)
```json
{
  "PROJECTION_2D_METHOD": "max",         // 'max', 'mean', 'sum'
  "TERRITORY_EXCLUSION_UM": 1.0,         // Gap Voronoi (evita solapamientos)
  "CONNECT_SKELETON_FRAGMENTS": true,    // Conectar fragmentos cercanos
  "CONNECTION_RADIUS_UM": 0.5,           // Radio máximo de conexión
  "SHOLL_MIN_UM": 5.0,                   // Radio mínimo Sholl
  "SHOLL_MAX_UM": 100.0,                 // Radio máximo Sholl
  "SHOLL_STEP_UM": 2.0                   // Separación entre anillos
}
```

**Recomendaciones:**
- `TERRITORY_EXCLUSION_UM`: 1-3 µm (ajustar según densidad celular)
- `CONNECTION_RADIUS_UM`: 0.3-1.0 µm (más alto = más conexiones)
- Rango Sholl: ajustar según extensión esperada de procesos

### Modificar Parámetros

**Desde la UI (Recomendado):**
1. Abrir Streamlit → Sidebar
2. Expandir sección correspondiente
3. Ajustar valores en controles numéricos
4. Click "Guardar [Sección]"

**Manualmente:**
- Editar `streamlit/calibration.json`
- Reiniciar Streamlit para aplicar cambios
- ⚠️ Regenerar datos si parámetros afectan procesamiento

---

## 🧪 Ejecución en Terminal (Sin UI)

Para procesamiento automatizado o HPC:

```bash
# Activar entorno
source venv/bin/activate

# Ejecutar pipeline completo para un preparado
python -c "
from pathlib import Path
from streamlit.ui import runner
import json

root = Path('.')
img_path = Path('data/raw/CTL/preparado_ejemplo.tif')
cal = json.loads((root / 'streamlit/calibration.json').read_text())

result = runner.run_pipeline_for(
    root=root,
    img_path=img_path,
    cal=cal,
    start_step='01',
    overwrite_from_step=True
)
print(f'Resultado: {result}')
"
```

### Script de Batch Processing

```python
from pathlib import Path
from streamlit.ui import runner
import json

root = Path('.')
cal_path = root / 'streamlit/calibration.json'
cal = json.loads(cal_path.read_text())

# Listar todos los preparados
files = runner.list_raw_images(root)

# Ejecutar batch para grupo CTL
results = runner.run_scope(
    root=root,
    scope='group',
    start_step='01',
    cal=cal,
    group='CTL',
    overwrite_from_step=False  # No sobrescribir existentes
)

# Resumen
ok = sum(1 for _, r in results if not r.get('error'))
print(f'Completados: {ok}/{len(results)}')
```

---

## 📈 Interpretación de Resultados

### Criterios de Calidad

#### Segmentación Nuclear (Paso 02)
- ✅ **Buena:** ~200-500 núcleos/preparado, pocos solapamientos
- ⚠️ **Revisar:** <100 núcleos (sub-segmentación) o >1000 (sobre-segmentación)
- 🔧 **Ajustar:** `NUCLEUS_DIAMETER` si hay problemas sistemáticos

#### Filtrado GFAP (Paso 03)
- ✅ **Buena:** Retención 20-40% de núcleos como candidatos
- ⚠️ **Revisar:** <10% (muy estricto) o >70% (muy permisivo)
- 🔧 **Ajustar:** `GFAP_STD_DEV_THRESHOLD` para calibrar

#### Análisis de Sholl (Paso 04)
- ✅ **Buena:** Perfiles con pico claro, AUC > 100
- ⚠️ **Revisar:** Perfiles planos (sin ramificación) o pico en radio mínimo (solo soma)
- 🔧 **Verificar:** Esqueletos en Napari 2D, ajustar parámetros de conexión si fragmentado

### Comparación Estadística (Página 06)

**p < 0.05:** Diferencia estadísticamente significativa
- Hay evidencia de que los grupos difieren en la métrica evaluada
- **Importante:** Verificar también la magnitud del efecto (diferencia de medias/medianas)

**p ≥ 0.05:** No se detecta diferencia significativa
- No hay evidencia suficiente de diferencia
- **No significa** que los grupos sean iguales (puede haber falta de poder estadístico)

**Consideraciones:**
- N = número de preparados (no células) — Ver tabla por preparado
- Típicamente N=3-5 por grupo (poder estadístico limitado)
- Complementar con análisis de magnitud del efecto y visualización

---

## 🐛 Troubleshooting

### Problemas Comunes

#### "No se detectan núcleos en Cellpose"
**Causas:**
- Calibración incorrecta (z,y,x)
- `NUCLEUS_DIAMETER` muy grande o pequeño
- Canal DAPI incorrecto

**Soluciones:**
1. Verificar calibración en Página 01
2. Ajustar `NUCLEUS_DIAMETER` (típico: 20-40 píxeles)
3. Verificar índice de canal DAPI (típicamente 0)
4. Probar con/sin GPU

#### "Muy pocos astrocitos detectados"
**Causas:**
- `GFAP_STD_DEV_THRESHOLD` muy alto (muy estricto)
- `MIN_VOLUME_UM3` muy alto
- Señal GFAP débil en el preparado

**Soluciones:**
1. Reducir `GFAP_STD_DEV_THRESHOLD` a 2.0-2.5
2. Reducir `MIN_VOLUME_UM3` a 50
3. Verificar histograma de intensidad GFAP en Página 03

#### "Territorios solapados / esqueletos entrelazados"
**Causas:**
- `TERRITORY_EXCLUSION_UM` muy bajo
- Astrocitos muy densos

**Soluciones:**
1. Aumentar `TERRITORY_EXCLUSION_UM` a 2-3 µm
2. Verificar territorios Voronoi en Napari 2D

#### "Esqueletos fragmentados"
**Causas:**
- Señal GFAP discontinua
- `CONNECTION_RADIUS_UM` muy bajo

**Soluciones:**
1. Activar `CONNECT_SKELETON_FRAGMENTS=true`
2. Aumentar `CONNECTION_RADIUS_UM` a 0.8-1.0 µm
3. Verificar calidad de señal GFAP en preparado original

#### "Perfiles de Sholl sin pico claro"
**Causas:**
- Astrocitos con procesos muy cortos (solo soma)
- Rango de radios inapropiado

**Soluciones:**
1. Ajustar `SHOLL_MIN_UM` (más bajo si soma domina)
2. Ajustar `SHOLL_MAX_UM` (más bajo si procesos cortos)
3. Verificar que esqueletos capturen procesos (Napari 2D)

### Errores de Ejecución

#### "FileNotFoundError: 02_cellpose_mask.tif"
**Causa:** Paso 02 no ejecutado o falló

**Solución:** Ir a Página 02 y ejecutar segmentación nuclear primero

#### "MemoryError" durante Cellpose
**Causa:** Volumen demasiado grande para RAM disponible

**Soluciones:**
1. Cerrar otras aplicaciones
2. Procesar preparados de a uno
3. Reducir tamaño de imagen (downsampling previo)

#### "RuntimeError: CUDA out of memory"
**Causa:** Volumen muy grande para GPU

**Soluciones:**
1. Desactivar `CELLPOSE_USE_GPU` (usar CPU)
2. Reducir tamaño de batch interno de Cellpose
3. Procesar en máquina con más VRAM

---

## 📚 Documentación Adicional

### Archivos de Documentación

- **`streamlit/pages/README.md`** — Metodología detallada por página (este documento es más completo)
- **`CAMBIOS_ESFERICIDAD_IBA1.md`** — Cambios recientes: esfericidad 2D y eliminación Iba-1
- **`requirements.txt`** — Lista completa de dependencias con versiones

### Referencias Científicas

**Cellpose:**
- Stringer et al. (2021). "Cellpose: a generalist algorithm for cellular segmentation." *Nature Methods*.

**SKAN (Skeleton Analysis):**
- Nunez-Iglesias et al. (2018). "A multi-level topological analysis of cytoskeletal networks." *PLoS ONE*.

**Sholl Analysis:**
- Sholl, D.A. (1953). "Dendritic organization in the neurons of the visual and motor cortices of the cat." *Journal of Anatomy*.

**Voronoi Tessellation:**
- Aplicación estándar de partición espacial para definir territorios celulares no solapados.

---

## 🤝 Contribuciones

### Estructura de Commits

```bash
# Formato recomendado
git commit -m "[TIPO] Descripción breve

Descripción detallada de los cambios realizados.
Razón del cambio y contexto adicional.
"
```

**Tipos:**
- `[FIX]` — Corrección de bugs
- `[FEAT]` — Nueva funcionalidad
- `[REFACTOR]` — Reestructuración sin cambio funcional
- `[DOCS]` — Actualización de documentación
- `[TEST]` — Adición de tests

### Pruebas Antes de Commit

```bash
# Verificar que no hay errores de sintaxis
python -m py_compile streamlit/ui/pipeline.py

# Ejecutar tests (si existen)
pytest tests/

# Verificar que Streamlit carga sin errores
cd streamlit && timeout 10 streamlit run Home.py --server.headless=true
```

---

## 📄 Licencia

MIT License

Copyright (c) 2024-2025

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## 📧 Contacto y Soporte

**Desarrollador:** Daniel Zangrando  
**Repositorio:** https://github.com/DanZangrando/astrocitos-3d-analysis

**Para reportar problemas:**
- Abrir un Issue en GitHub con:
  - Descripción del problema
  - Pasos para reproducirlo
  - Logs de error completos
  - Versión de Python y sistema operativo

**Para consultas:**
- Email del proyecto o
- Discusiones en GitHub

---

## 🔄 Actualizaciones Recientes

### v2.0 — Pipeline 2D Nativo (Noviembre 2024)
- ✅ Implementación completa de pipeline 2D unificado (Paso 04)
- ✅ Eliminación del canal Iba-1 por contaminación con GFAP
- ✅ Esfericidad 2D (circularidad) calculada desde proyección MIP
- ✅ Territorios Voronoi 2D con zona de exclusión
- ✅ Sholl 2D nativo integrado con SKAN
- ✅ Conexión inteligente de fragmentos de esqueleto
- ✅ ~10x mejora en velocidad vs pipeline 3D isotrópico
- ✅ Resolución XY completa (0.38 µm) sin degradación

### v1.0 — Pipeline Base (Julio 2024)
- Segmentación nuclear con Cellpose 3D
- Filtrado biológico por señal GFAP + Iba-1
- Esqueletización 3D con remuestreo isotrópico
- Análisis de Sholl 3D
- Interfaz Streamlit básica

---

## 🔮 Roadmap

### Corto Plazo
- [ ] Tests automatizados (pytest)
- [ ] Paralelización de procesamiento batch
- [ ] Exportación automática de figuras para publicación

### Mediano Plazo
- [ ] Análisis de expansión territorial (convex hull, área efectiva)
- [ ] Métricas de orientación y anisotropía
- [ ] Clustering de células por fenotipos morfológicos

### Largo Plazo
- [ ] Integración con bases de datos (SQLite/PostgreSQL)
- [ ] API REST para procesamiento remoto
- [ ] Interfaz web (migración de Streamlit a framework más robusto)

---

**¡Gracias por usar este pipeline! 🧠✨**

## 🎯 Enfoque 2D Nativo

El pipeline utiliza un **flujo de trabajo 2D nativo optimizado** para astrocitos en preparados relativamente planos (típico de cortes de hipocampo):

**Ventajas:**
- ✅ **Resolución XY completa** (0.38 µm) sin degradación por remuestreo
- ✅ **Sholl 2D nativo** más preciso y eficiente que aproximaciones 3D
- ✅ **Voronoi 2D simple y robusto** para definir territorios astrocitarios
- ✅ **~10x más rápido** que esqueletización 3D con remuestreo isotrópico
- ✅ **Territorios bien definidos** para células con dominios principalmente planares

**Rationale científico:** Los astrocitos en preparados de hipocampo tienen dominios territoriales principalmente en el plano XY. La componente Z (0.38 µm/slice × ~40-50 slices) es útil para segmentación nuclear, pero la morfología territorial se captura eficientemente mediante proyección máxima.

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

1) **Calibración física global (µm)**: Define `z`, `y`, `x` en `streamlit/calibration.json` y/o abre Napari con escala.

2) **Otsu y Cellpose (DAPI)**: Umbral Otsu opcional y segmentación 3D con Cellpose.

3) **Filtrado GFAP**: Evalúa la señal de GFAP en un anillo perinuclear alrededor de cada núcleo usando umbrales relativos (desviaciones estándar sobre el fondo). Los núcleos con suficiente señal de GFAP son marcados como candidatos a astrocitos. Aplica filtro de tamaño físico (µm³).
   - **Nota**: El filtrado por canal Iba-1 (microglía) ha sido eliminado debido a contaminación detectada con la señal de GFAP.

4) **Esqueletización 3D por célula**:
   - Re-muestreo a vóxel isotrópico objetivo; umbral GFAP (Otsu/Manual) y cierre morfológico.
   - Selección del GFAP conectado a la semilla (núcleo dilatado) con conectividad 3D.
   - Radio máximo desde el núcleo y padding que incluye dicho radio para evitar truncamientos.
   - Opcionales: territorios Voronoi por cercanía al núcleo con “gap” en fronteras; resolución de solapamientos por distancia física al núcleo.
   - Pruning topológico (opcional): elimina espículas cortas en el esqueleto.
   - Métrica centrada en el esqueleto: integra la señal GFAP dentro de un “tubo” de radio fijo (µm) alrededor del esqueleto y reporta intensidad por unidad de longitud.

5) **Sholl por célula**: Anillos concéntricos alrededor de cada núcleo/astrocito y conteo de intersecciones por radio. Los anillos se guardan en JSON y se visualizan en Napari.

6) **Resumen y estadísticas**: Tablas por imagen y comparaciones por grupo (p.ej., CTL vs Hipoxia) usando pruebas de Mann–Whitney U para métricas clave (longitud, endpoints, pico de Sholl, etc.).

**Idempotencia**: Cada paso verifica si su salida existe y la reutiliza para evitar recomputar innecesariamente.

## Salidas por preparado (archivos estándar)

- `01_otsu_mask.tif` — Máscara Otsu (DAPI) opcional
- `02_cellpose_mask.tif` — Etiquetas de núcleos de Cellpose (3D)
- `03_gfap_filtered_mask.tif` — Candidatos GFAP tras filtrado biológico
- `03_nucleus_metrics.csv` — Métricas detalladas por núcleo (volumen, intensidad GFAP, flags de retención)
- `04_final_astrocytes_mask.tif` — Máscara final tras filtro de tamaño (3D)
- `05_skeleton_labels_2d.tif` — Esqueletos por célula proyectados en 2D (1 píxel de grosor)
- `sholl_2d_native.csv` — Perfiles de Sholl por célula (intersecciones por radio)
- `sholl_2d_native_detailed.csv` — Intersecciones completas por radio y célula
- `sholl_rings_2d_native.json` — Coordenadas de anillos para visualización en Napari
- `sholl_summary.csv` — Radio crítico, pico de intersecciones, AUC por célula
- `params.json` — Parámetros usados para ese preparado (persistencia por imagen)

## Troubleshooting

- **Anillos de Sholl desplazados**: Corregido añadiendo `scale=scale` a la capa de Shapes en Napari (ya implementado).
- **Territorios solapados**: Ajustar `TERRITORY_EXCLUSION_UM` (gap de exclusión Voronoi, típico: 2-5 µm).
- **GFAP fragmentado**: Activar `CONNECT_SKELETON_FRAGMENTS=true` y ajustar `CONNECTION_RADIUS_UM` para conectar fragmentos cercanos.
- **Performance lenta**: El pipeline 2D nativo es ~10x más rápido que el anterior flujo 3D isotrópico. Napari puede abrirse desde cada página con escala correcta.
- **Resolución espacial**: El flujo 2D mantiene la resolución XY completa (0.38 µm) sin degradación por remuestreo isotrópico.

## Licencia

MIT (o actualizar según corresponda).
