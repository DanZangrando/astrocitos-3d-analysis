# Métricas SKAN Recomendadas para Análisis de Astrocitos

## 🔬 Contexto Biológico: Astrocitos vs Neuronas

### Diferencias Morfológicas Clave:

**NEURONAS:**
- Estructura polarizada: dendritas → soma → axón
- Procesos dendríticos jerárquicos (árbol dendrítico)
- Espinas dendríticas para sinapsis
- Métricas típicas: complejidad dendrítica, longitud axonal

**ASTROCITOS:**
- Estructura **no polarizada**: procesos radiales desde el soma
- **Dominio territorial**: cada astrocito ocupa un territorio definido sin superposición
- Procesos finos que contactan vasos sanguíneos y sinapsis
- **Funciones**: soporte metabólico, homeostasis iónica, barrera hematoencefálica
- En hipoxia: **reactividad astrocitaria** → cambios morfológicos y funcionales

---

## 📊 Métricas SKAN Disponibles (de `summarize()`)

| Columna | Descripción | Unidad |
|---------|-------------|--------|
| `skeleton_id` | ID de componente conectada | - |
| `branch_distance` | Longitud real de la rama | µm |
| `euclidean_distance` | Distancia directa entre extremos | µm |
| `branch_type` | 0=path, 1=endpoint, 2=junction, 3=cycle | - |
| `coord_src/dst` | Coordenadas físicas de inicio/fin | µm |

---

## ✅ Métricas YA Implementadas

1. **Longitud total del esqueleto** (`total_branch_length_um`)
2. **Número de ramas** (`n_branches`)
3. **Endpoints** (`n_endpoints`)
4. **Junctions** (`n_junctions`)
5. **Longitud promedio de rama** (`mean_branch_length_um`)

---

## 🆕 Métricas NUEVAS Agregadas (Relevantes para Astrocitos)

### 1️⃣ **TORTUOSIDAD** (Sinuosidad de Procesos)

**¿Qué mide?** Cuán "curvos" o "rectos" son los procesos astrocitarios.

**Cálculo:** `tortuosidad = branch_distance / euclidean_distance`
- Valor 1.0 = proceso perfectamente recto
- Valor > 1.0 = proceso sinuoso/tortuoso

**Relevancia en Hipoxia:**
- **Hipótesis**: Hipoxia aumenta tortuosidad por:
  - Desorganización citoesquelética
  - Edema celular
  - Retracción/reorientación de procesos

**Métricas calculadas:**
- `tortuosity_mean`: Tortuosidad promedio de todos los procesos
- `tortuosity_max`: Máxima tortuosidad observada
- `tortuosity_std`: Variabilidad en tortuosidad

**Interpretación:**
- ↑ Tortuosidad media → Procesos menos eficientes, posible patología
- ↑ Desv. estándar → Mayor heterogeneidad morfológica

---

### 2️⃣ **ÍNDICE DE RAMIFICACIÓN** (Ramification Index)

**¿Qué mide?** Cuántas ramas tiene el astrocito en relación a sus bifurcaciones.

**Cálculo:** `ramification_index = n_branches / max(n_junctions, 1)`

**Relevancia en Hipoxia:**
- **Astrocitos reactivos** pueden mostrar:
  - Hipertrofia: más procesos → ↑ índice
  - Atrofia: menos procesos → ↓ índice

**Interpretación:**
- Índice alto → Astrocito con muchas ramas (más complejo)
- Índice bajo → Astrocito simple, menos ramificado

---

### 3️⃣ **ÍNDICE DE TERMINALIZACIÓN** (Termination Index)

**¿Qué mide?** Balance entre terminaciones libres y bifurcaciones.

**Cálculo:** `termination_index = n_endpoints / n_junctions`

**Relevancia para Astrocitos:**
- **Alto (>2)**: Muchas terminaciones → procesos que se extienden hacia targets (vasos, sinapsis)
- **Bajo (<2)**: Pocas terminaciones → estructura más ramificada centralmente

**Interpretación:**
- ↑ en hipoxia → Extensión de procesos hacia vasos (buscar oxígeno)
- ↓ en hipoxia → Retracción de procesos (daño celular)

---

### 4️⃣ **DISTRIBUCIÓN DE LONGITUDES DE RAMAS**

**¿Qué mide?** Heterogeneidad en la longitud de los procesos.

**Métricas calculadas:**
- `branch_length_median_um`: Mediana (más robusta que media)
- `branch_length_p25_um`: Percentil 25
- `branch_length_p75_um`: Percentil 75
- `branch_length_std_um`: Desviación estándar
- `branch_length_cv`: Coeficiente de variación (CV = std/mean)

**Relevancia en Hipoxia:**
- **CV alto**: Gran heterogeneidad → algunos procesos largos, otros cortos (morfología desorganizada)
- **CV bajo**: Homogeneidad → procesos uniformes (morfología conservada)

**Interpretación:**
- ↑ CV en hipoxia → Mayor variabilidad morfológica (posible disfunción)
- Cambios en percentiles → Detecta si se afectan procesos cortos vs largos

---

### 5️⃣ **COMPONENTES CONECTADAS** (Fragmentación)

**¿Qué mide?** Cuántos fragmentos independientes tiene el esqueleto.

**Cálculo:** `n_connected_components = skeleton_id.nunique()`

**Valor esperado:** 1 (astrocito íntegro)

**Relevancia en Hipoxia:**
- **>1 componentes**: Esqueleto fragmentado → posibles causas:
  - Daño celular real (necrosis, apoptosis)
  - Artefacto de procesamiento (territorios Voronoi cortando procesos)
  - Separación de endfeet (pies vasculares)

**Interpretación:**
- Fragmentación alta → Revisar calidad de datos o daño real

---

## 📈 Análisis Estadístico Sugerido

### Comparación Hipoxia vs Normoxia:

1. **Test de Hipótesis** (por métrica):
   - Mann-Whitney U test (datos no paramétricos)
   - T-test (si distribución normal)
   - Corrección Bonferroni para tests múltiples

2. **Visualizaciones Recomendadas**:
   - **Violin plots**: Comparar distribuciones completas
   - **Box plots**: Detectar outliers
   - **Scatter plots**: Tortuosidad vs Sholl AUC
   - **Heatmaps**: Correlación entre métricas

3. **Análisis Multivariado**:
   - **PCA**: Identificar patrones morfológicos discriminantes
   - **Clustering**: ¿Existen subpoblaciones astrocitarias?
   - **Random Forest**: ¿Qué métricas mejor predicen grupo?

---

## 🎯 Hipótesis Específicas para Astrocitos en Hipoxia

### Esperado en **HIPOXIA AGUDA**:
- ↑ **Tortuosidad** (desorganización citoesquelética)
- ↑ **Longitud total** (hipertrofia reactiva inicial)
- ↑ **N° endpoints** (extensión hacia vasos)
- ↓ **Sholl AUC** en radios grandes (retracción periférica)
- ↑ **CV longitudes** (heterogeneidad morfológica)

### Esperado en **HIPOXIA CRÓNICA**:
- ↓ **Longitud total** (atrofia, pérdida de procesos)
- ↓ **Ramification index** (simplificación morfológica)
- ↑ **Fragmentación** (daño celular progresivo)
- ↓ **Termination index** (retracción de endfeet)

---

## 💡 Métricas Complementarias (Ya Disponibles)

Estas métricas **YA las tenemos** del pipeline existente y se complementan bien:

1. **Sholl Analysis**:
   - `sholl_auc`: Complejidad territorial total
   - `sholl_peak_intersections`: Máxima densidad de procesos
   - `sholl_critical_radius`: Distancia de mayor ramificación

2. **Métricas Nucleares**:
   - `nuclear_area_um2`: Tamaño del soma (hipertrofia en reactividad)
   - `nuclear_circularity`: Morfología nuclear (daño → irregularidad)

3. **Territoriales**:
   - `territorial_area_um2`: Tamaño del dominio
   - `gfap_positive_area_um2`: Expresión de GFAP (↑ en reactividad)

---

## 🔗 Correlaciones Esperadas

| Métrica 1 | Métrica 2 | Correlación Esperada | Interpretación |
|-----------|-----------|----------------------|----------------|
| Tortuosidad ↑ | Sholl AUC ↓ | Negativa | Procesos tortuosos → menos organización |
| Ramification index ↑ | N° endpoints ↑ | Positiva | Más ramas → más terminaciones |
| Longitud total ↑ | Territorial area ↑ | Positiva | Astrocito grande → más procesos |
| CV longitudes ↑ | Fragmentación ↑ | Positiva | Heterogeneidad → posible daño |

---

## 📝 Reporte Sugerido para Paper

### Tabla de Métricas Morfológicas:

| Métrica | Normoxia (n=X) | Hipoxia (n=Y) | p-value | Interpretación |
|---------|----------------|---------------|---------|----------------|
| Longitud total (µm) | 450 ± 120 | 380 ± 150 | 0.023* | ↓ Atrofia en hipoxia |
| Tortuosidad | 1.08 ± 0.12 | 1.24 ± 0.18 | <0.001*** | ↑ Desorganización |
| Ramification index | 3.2 ± 0.8 | 2.5 ± 0.9 | 0.008** | ↓ Simplificación |
| Sholl AUC | 85 ± 22 | 62 ± 28 | 0.002** | ↓ Complejidad territorial |

*(Valores ilustrativos)*

---

## 🚀 Próximos Pasos

1. ✅ **Implementar métricas** en `pipeline_2d_unified.py`
2. ⏳ **Regenerar todos los datos** con nuevas métricas
3. ⏳ **Validar con datos reales** (1-2 preparados de cada grupo)
4. ⏳ **Análisis estadístico** en página 07 (comparaciones globales)
5. ⏳ **Visualizaciones** adicionales para nuevas métricas

---

## 📚 Referencias Sugeridas

- Bushong et al. (2002) *Neuron* - Territorios astrocitarios
- Oberheim et al. (2012) *J Neurosci* - Morfología astrocitaria en primates
- Sofroniew (2009) *Trends Neurosci* - Reactividad astrocitaria
- Hauglund et al. (2020) *Nat Commun* - Astrocitos en hipoxia

