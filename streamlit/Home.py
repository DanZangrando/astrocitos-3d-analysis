import streamlit as st
from pathlib import Path
from ui.sidebar import render_sidebar
from ui.utils import detect_group
import json
from datetime import datetime
import pandas as pd
import subprocess
import os

# Configuración de página
st.set_page_config(
    page_title="Astrocitos 3D - Análisis", 
    page_icon="🧠", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Título y Sidebar
st.title("🧠 Análisis Morfológico de Astrocitos")
render_sidebar(show_calibration=True)

# Descripción del Proyecto
st.markdown("""
### 📌 Descripción del Proyecto
Plataforma interactiva para la **reconstrucción, segmentación y análisis morfológico** de astrocitos a partir de microscopía confocal. 
El sistema está diseñado para caracterizar la reactividad astrocitaria mediante biomarcadores topológicos robustos.

---

### 🚀 Flujo de Trabajo (Pipeline)

El análisis se estructura en 5 pasos secuenciales:

1.  **Calibración y Visualización (01)**: Definición de la escala física (µm/pixel) y control de calidad de la imagen cruda.
2.  **Segmentación Nuclear (02)**: Identificación de nucleos (DAPI) mediante Deep Learning (**Cellpose**).
3.  **Filtrado de Astrocitos (03)**: Selección de núcleos correspondientes a astrocitos mediante colocalización de señal **GFAP** (filtro de intensidad relativo) y tamaño físico.
4.  **Esqueletización y Sholl** (04): 
    *   Proyección inteligente 3D → 2D.
    *   Definición de territorios celulares (Voronoi con exclusión).
    *   Esqueletización topológica y análisis de **Sholl** nativo.
5.  **Análisis por Preparado (05)**: Validación y exploración de resultados individuales.
6.  **Comparación de Grupos (06)**: Análisis estadístico robusto (**CTL vs Hipoxia**) evitando pseudoreplicación.

---

### 📊 Biomarcadores Clave

Para caracterizar la morfología, el sistema se enfoca en 4 métricas validadas:

| Métrica | Significado Biológico | Tipo |
| :--- | :--- | :--- |
| **Radio Crítico (Sholl)** | Distancia de máxima arborización. Indica la *expansión espacial*. | Sholl |
| **Índice de Ramificación** | Relación Ramas/Uniones. Refleja la *complejidad topológica*. | Topología |
| **Longitud Total del Esqueleto** | Suma de todas las ramas. Indica el *volumen de exploración*. | Topología |
| **Número de Terminaciones** | Puntos finales de las ramas. Indica la *división terminal*. | Topología |

---
""")

st.info("💡 **Tip:** Usá el menú lateral para navegar paso a paso. Los resultados se guardan automáticamente en `data/processed/`.")

# --- Inspector Configuración ---
root = Path(__file__).resolve().parents[1]
calib_path = root / "streamlit" / "calibration.json"

with st.expander("🛠️ Configuración Global (calibration.json)", expanded=False):
    def _load_calib() -> dict:
        if calib_path.exists():
            try: return json.loads(calib_path.read_text())
            except: return {}
        return {}

    cfg = _load_calib()
    st.json(cfg)

# --- Dashboard de Estado ---
st.markdown("### 📋 Estado de Procesamiento por Preparado")

raw_dir = root / "data" / "raw"
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])

if files:
    stem_to_group_all = {p.stem: detect_group(p, root) for p in files}
    
    rows = []
    for p in files:
        od = root / "data" / "processed" / p.stem
        rows.append({
            "Preparado": p.stem,
            "Grupo": stem_to_group_all.get(p.stem, "CTL"),
            "1. Calibrado": (root / "streamlit" / "calibration.json").exists(), # Global por ahora
            "2. Núcleos": (od/"02_cellpose_mask.tif").exists(),
            "3. Filtrado": (od/"03_gfap_filtered_mask.tif").exists(),
            "4. Esqueletos/Sholl": (od/"skeletons"/"summary.csv").exists() and (od/"sholl_summary.csv").exists(),
            "5. Analizado": (od/"skeletons"/"summary.csv").exists()
        })
    
    df_state = pd.DataFrame(rows)
    
    if not df_state.empty:
        # Formato visual
        def mark(val): return "✅ Completo" if val else "—"
        
        view = df_state.copy()
        for c in ["1. Calibrado", "2. Núcleos", "3. Filtrado", "4. Esqueletos/Sholl", "5. Analizado"]:
            view[c] = view[c].apply(mark)
            
        st.dataframe(
            view.set_index("Preparado"),
            use_container_width=True,
            column_config={
                "Grupo": st.column_config.TextColumn("Grupo", width="small"),
            }
        )
    else:
        st.info("No hay datos procesados para mostrar.")
else:
    st.warning("No se encontraron imágenes en `data/raw`. Por favor cargá tus datos para comenzar.")