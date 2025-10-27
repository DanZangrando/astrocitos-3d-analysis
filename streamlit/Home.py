import streamlit as st
from pathlib import Path
from ui.sidebar import render_sidebar

# Welcome page content based on project README
st.set_page_config(page_title="Astrocitos 3D - Análisis", page_icon="🧠", layout="wide")

st.title("Análisis Morfológico 3D de Astrocitos")
render_sidebar(show_calibration=True)

st.markdown(
    """
    **Proyecto:** Reconstrucción, segmentación y análisis cuantitativo de astrocitos (GFAP), otros marcadores y núcleos (DAPI) a partir de imágenes `.lif`.

    **Tecnologías:** Python, Napari, Cellpose, Scikit-image, Scipy, Skan

    ---

    ### Objetivo
    Desarrollar un flujo de trabajo para segmentar y analizar la morfología tridimensional de astrocitos en imágenes de microscopía confocal multi-canal, con énfasis en la identificación de astrocitos completos para cuantificar su complejidad estructural.

    ### Flujo de trabajo (resumen)
    - Segmentación de núcleos (DAPI) con Cellpose.
    - Filtrado por co-localización con GFAP y exclusión por microglía.
    - Limpieza por tamaño físico (µm³) y análisis morfológico avanzado (skeletonización 3D, ramas, Sholl, territorio celular).

    ---
    Esta aplicación Streamlit proporciona:
    1) Una interfaz para definir la calibración física (tamaño de vóxel en µm) por preparado/imágen.
    2) Visualización rápida de los preparados y apertura de Napari con la calibración seleccionada.
    3) Persistencia de parámetros en `data/calibration_overrides.json` para usarlos en los notebooks y en el pipeline.
    """
)

# Links útiles
root = Path(__file__).resolve().parents[1]
notebooks_dir = root / "notebooks"
results_dir = root / "results"

# La navegación lateral de Streamlit Multi-Page se genera automáticamente.
# No añadimos enlaces manuales para evitar duplicados.

st.info("Usa el menú lateral para ir a 'Calibración y Visualización'. Ahí podés fijar la escala (Z/Y/X) y abrir Napari.")
