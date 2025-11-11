# Resumen de Actualizaciones - Pipeline de Análisis de Astrocitos

**Fecha:** 2025-11-11  
**Objetivo:** Documentar métricas, tests estadísticos y flujo completo del pipeline

---

## 🎯 Cambios Implementados

### 1. Nuevas Métricas Topológicas (11 columnas adicionales)

**Archivo:** `streamlit/ui/pipeline_2d_unified.py`

#### A. Tortuosidad (3 métricas)
- `tortuosity_mean` — Tortuosidad promedio (1.0 = recto)
- `tortuosity_max` — Tortuosidad máxima
- `tortuosity_std` — Desv. std tortuosidad

**Relevancia para hipoxia:** ↑ indica desorganización citoesquelética, edema

#### B. Complejidad Morfológica (2 índices)
- `ramification_index` — Ramas / Junctions
- `termination_index` — Endpoints / Junctions

**Relevancia para hipoxia:** ↓ ramificación = simplificación; ↑ terminación = búsqueda O₂

#### C. Distribución de Longitudes (5 métricas)
- `branch_length_median_um` — Mediana
- `branch_length_p25_um` — Percentil 25
- `branch_length_p75_um` — Percentil 75
- `branch_length_std_um` — Desv. estándar
- `branch_length_cv` — Coeficiente de variación

**Relevancia para hipoxia:** ↑ CV indica heterogeneidad/desorganización morfológica

#### D. Fragmentación (1 métrica)
- `n_connected_components` — Componentes conectadas (esperado: 1)

**Total:** De 11 columnas → **22 columnas** en `skeletons/summary.csv`

### 2. Filtro de Componentes Desconectadas

**Archivo:** `streamlit/ui/pipeline_2d_unified.py`

**Función:** `filter_skeleton_by_nucleus_connectivity()`

**Qué hace:**
1. Detecta componentes conectadas del esqueleto
2. **RETIENE SOLO** la que toca el núcleo
3. **ELIMINA** fragmentos aislados dentro del territorio

**Justificación:**
- Fragmentos aislados = señal de fondo, procesos de células vecinas, artefactos
- NO pertenecen a la morfología del astrocito en cuestión

**Impacto:** Esqueletos más limpios y precisos

### 3. Corrección de Visualización Napari

**Archivo:** `streamlit/napari_viewer_2d.py`

**Problema:** Canal DAPI mostraba mezcla con GFAP

**Solución:** Proyección MIP correcta por canal independiente
- DAPI proyectado → canal 0
- GFAP proyectado → canal 1
- Sin mezcla entre ellos

### 4. Actualización de Página 06 (Comparaciones)

**Archivo:** `streamlit/pages/06_ Comparación entre Grupos.py`

**Cambios:**
1. Carga de `skeletons/summary.csv` agregada
2. Merge correcto: esqueleto → Sholl → núcleo
3. **19 métricas nuevas** disponibles para comparación:
   - Todas las topológicas (tortuosidad, complejidad, distribución)
   - Mantenidas: Sholl (AUC, pico, radio crítico)
   - Mantenidas: Nucleares (volumen, esfericidad)

### 5. Documentación Completa

**Archivo:** `PIPELINE_DOCUMENTATION_UPDATED.md`

**Contenido:**
- Descripción completa de cada paso (01-06)
- Todas las 22 métricas documentadas con interpretación biológica
- Hipótesis esperadas para hipoxia vs control
- Tests estadísticos explicados
- Troubleshooting

---

## 📊 Métricas Completas por Archivo

### `03_nucleus_metrics.csv` (5 columnas)
| Columna | Descripción |
|---------|-------------|
| `label` | ID núcleo |
| `nucleus_volume_um3` | Volumen físico (µm³) |
| `nucleus_sphericity` | Esfericidad 2D (0-1) |
| `shell_gfap_mean` | Intensidad GFAP shell |
| `is_astrocyte_candidate` | Pasó filtro GFAP |

### `skeletons/summary.csv` (22 columnas)
| Categoría | Columnas | Total |
|-----------|----------|-------|
| Básicas | label, n_branches, total_branch_length_um, mean_branch_length_um, n_endpoints, n_junctions, skeleton_pixels, nuclear_area_um2, centroid_y/x_um | 10 |
| Tortuosidad | tortuosity_mean, tortuosity_max, tortuosity_std | 3 |
| Complejidad | ramification_index, termination_index | 2 |
| Distribución | branch_length_median/p25/p75/std/cv_um | 5 |
| Fragmentación | n_connected_components | 1 |
| **TOTAL** | | **22** |

### `sholl_summary.csv` (4 columnas)
| Columna | Descripción |
|---------|-------------|
| `label` | ID astrocito |
| `auc` | Área bajo curva (µm·intersec) |
| `peak_intersections` | Máximo intersecciones |
| `critical_radius_um` | Radio del pico (µm) |

### `sholl_2d_native.csv` (3 columnas, N filas)
| Columna | Descripción |
|---------|-------------|
| `label` | ID astrocito |
| `radius_um` | Radio del anillo |
| `intersections` | N° cruces en ese radio |

---

## 🧪 Tests Estadísticos Aplicados

### Flujo de Decisión

```
Para cada métrica:
  ↓
1. Agregación por preparado (mediana)
  ↓
2. Shapiro-Wilk en cada grupo
  ↓
  ├─ Ambos normales → Welch's t-test (compara medias)
  └─ Alguno no normal → Mann-Whitney U (compara medianas)
  ↓
3. Reporte: test usado, estadístico, p-valor, interpretación
```

### Test de Normalidad: Shapiro-Wilk
- **H₀:** Distribución normal
- **α:** 0.05
- **Decisión:** p > 0.05 → normal

### Test Paramétrico: Welch's t-test
- **Cuándo:** Ambos grupos normales
- **H₀:** μ_CTL = μ_Hipoxia
- **Compara:** Medias
- **Ventaja:** No asume varianzas iguales

### Test No Paramétrico: Mann-Whitney U
- **Cuándo:** Al menos un grupo no normal
- **H₀:** Distribuciones idénticas
- **Compara:** Medianas/rangos
- **Ventaja:** Robusto a outliers

### Nivel de Significancia
- **α = 0.05**
- **p < 0.05:** Diferencia significativa (rechazar H₀)
- **p ≥ 0.05:** Sin evidencia suficiente (no rechazar H₀)

### ⚠️ Pseudoreplicación: RESUELTA

**Problema:** Células del mismo preparado no son independientes

**Solución:**
1. **df_plot** (para gráficos): Todas las células individuales
2. **df_stats** (para tests): Mediana por preparado
   - N = número de preparados (unidad experimental verdadera)
   - CTL: típicamente 3-5 preparados
   - Hipoxia: típicamente 3-5 preparados

---

## 📈 Hipótesis Biológicas para Hipoxia vs Control

### Esperado en Hipoxia AGUDA
| Métrica | Cambio | Interpretación |
|---------|--------|----------------|
| `tortuosity_mean` | ↑ | Desorganización citoesquelética |
| `total_branch_length_um` | ↑ | Hipertrofia reactiva inicial |
| `termination_index` | ↑ | Extensión hacia vasos (búsqueda O₂) |
| `branch_length_cv` | ↑ | Heterogeneidad morfológica |
| `auc` (Sholl) | = o ↓ | Complejidad territorial conservada/reducida |

### Esperado en Hipoxia CRÓNICA
| Métrica | Cambio | Interpretación |
|---------|--------|----------------|
| `total_branch_length_um` | ↓ | Atrofia, pérdida de procesos |
| `ramification_index` | ↓ | Simplificación morfológica |
| `tortuosity_mean` | ↑↑ | Desorganización progresiva |
| `branch_length_cv` | ↑ | Gran heterogeneidad (disfunción) |
| `auc` (Sholl) | ↓↓ | Pérdida de complejidad territorial |
| `n_connected_components` | ↑ | Fragmentación (daño celular) |

---

## ✅ Checklist de Validación

### Antes de Análisis Estadístico
- [x] Pipeline 2D unificado implementado
- [x] Filtro de componentes desconectadas activo
- [x] 22 métricas topológicas calculadas
- [x] Visualización Napari corregida
- [x] Página 06 actualizada con todas las métricas
- [x] Documentación completa creada
- [ ] **Regenerar todos los preparados** con pipeline actualizado
- [ ] Validar con 2-4 preparados de cada grupo

### Para Análisis Estadístico
- [ ] Verificar N por grupo (mínimo 3 preparados/grupo)
- [ ] Revisar normalidad por métrica (Shapiro-Wilk)
- [ ] Aplicar tests apropiados (Welch/Mann-Whitney)
- [ ] Considerar corrección por tests múltiples (Bonferroni/FDR)
- [ ] Reportar tamaños de efecto (Cohen's d o similar)

### Para Publicación
- [ ] Tabla de métricas con estadísticas
- [ ] Violin plots o box plots por métrica clave
- [ ] Gráficos de correlación (tortuosidad vs Sholl AUC)
- [ ] PCA para identificar patrones discriminantes
- [ ] Reportar poder estadístico alcanzado

---

## 🚀 Próximos Pasos Recomendados

### 1. Validación con Datos Reales (INMEDIATO)
```
1. Streamlit → Página 04
2. Seleccionar 1-2 preparados de cada grupo
3. Ejecutar pipeline
4. Verificar:
   - Métricas topológicas en skeletons/summary.csv (22 columnas)
   - Valores biológicamente plausibles
   - Componentes desconectadas eliminadas correctamente
   - Visualización Napari muestra esqueletos limpios
```

### 2. Regeneración Batch (SI VALIDACIÓN OK)
```
1. Home.py → Sidebar → "Ejecución Batch"
2. Paso inicial: 04 (regenerar esqueletos y Sholl)
3. Sobrescribir: ✓
4. Ejecutar todos los preparados
```

### 3. Análisis Estadístico Completo (PÁGINA 06)
```
1. Comparar métricas clave:
   - Tortuosidad media
   - Índice de ramificación
   - CV longitudes
   - Sholl AUC
2. Exportar datos por célula y por preparado (CSV)
3. Análisis adicional en Python/R si es necesario:
   - PCA
   - Clustering
   - Correlaciones múltiples
```

### 4. Visualizaciones para Paper
```
- Violin plots: Tortuosidad CTL vs Hip
- Scatter: Ramification index vs Tortuosidad
- Heatmap: Correlaciones entre métricas
- Ejemplo Napari: Astrocito CTL vs Hip (con anillos Sholl)
```

---

## 📝 Archivos de Documentación Creados

1. **PIPELINE_DOCUMENTATION_UPDATED.md** — Documentación técnica completa
2. **METRICAS_ASTROCITOS_RECOMENDADAS.md** — Guía de métricas con contexto biológico
3. **IMPLEMENTACION_METRICAS_ASTROCITOS.md** — Resumen de implementación
4. **test_skeleton_connectivity_filter.py** — Test de validación del filtro
5. **test_new_metrics.py** — Test de nuevas métricas topológicas

---

## 🔗 Referencias Clave

### Sholl Analysis
- Sholl, D.A. (1953). *J Anat*, 87, 387-406.

### Astrocitos
- Bushong et al. (2002). *Neuron*, 34, 127-138.
- Oberheim et al. (2012). *J Neurosci*, 32, 3176-3187.
- Sofroniew, M.V. (2009). *Trends Neurosci*, 32, 638-647.

### SKAN
- Nunez-Iglesias et al. (2018). *Journal of Open Source Software*, 3(24), 1382.

---

**Estado:** ✅ Implementación completa, listo para validación con datos reales

**Próximo hito:** Regenerar todos los preparados → Análisis estadístico → Publicación

