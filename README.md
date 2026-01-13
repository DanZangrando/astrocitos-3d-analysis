# Análisis Morfológico 3D de Astrocitos

Plataforma interactiva para la reconstrucción, segmentación y análisis morfológico automatizado de astrocitos a partir de imágenes de microscopía confocal.

Este proyecto implementa un pipeline robusto para caracterizar la reactividad astrocitaria bajo condiciones experimentales (e.g., Hipoxia vs Control) utilizando biomarcadores topológicos validados.

## 🚀 Características Principales

*   **Pipeline Unificado**: Desde la imagen cruda (`.tif`/`.lif`) hasta la estadística, en 4 pasos estandarizados.
*   **Segmentación Híbrida**: Combina **Cellpose** (Deep Learning) para núcleos y algoritmos de umbralización adaptativa para procesos gliales.
*   **Análisis Topológico 2D Nativo**: Proyección inteligente y esqueletización basada en grafos para máxima resolución espacial.
*   **Biomarcadores Robustos**: Cuantificación precisa de complejidad, volumen y ramificación.
*   **Interfaz Interactiva**: Construida en **Streamlit** y **Napari** para validación visual inmediata.
*   **Análisis Estadístico**: Comparación automatizada entre grupos evitando pseudoreplicación.

## 📊 Pipeline de Análisis

1.  **Calibración y Visualización**: Definición de escala física (µm/px) y control de calidad.
2.  **Segmentación Nuclear**: Identificación de núcleos (DAPI) con Cellpose.
3.  **Filtrado de Candidatos**: Selección de astrocitos mediante colocalización de GFAP y volumen físico.
4.  **Esqueletización y Sholl**:
    *   Definición de territorios celulares (Voronoi).
    *   Análisis de Sholl (intersecciones radiales).
    *   Topología de grafos (ramas y uniones).

## 🔬 Biomarcadores Clave

El análisis se centra en 4 métricas validadas para describir la morfología astrocitaria:

| Métrica | Significado Biológico |
| :--- | :--- |
| **Radio Crítico (Sholl)** | Distancia de máxima arborización (expansión espacial). |
| **Índice de Ramificación** | Relación Ramas/Uniones (complejidad topológica). |
| **Longitud Total** | Suma de todas las ramas (volumen de exploración). |
| **Número de Terminaciones** | Puntos finales de las ramas (división terminal). |

## 💻 Instalación y Uso

### Prerrequisitos
*   **Python 3.10+**
*   **CUDA Toolkit** (Recomendado para aceleración de GPU con Cellpose)

### Instalación

1.  Clonar el repositorio:
    ```bash
    git clone https://github.com/DanZangrando/astrocitos-3d-analysis.git
    cd astrocitos-3d-analysis
    ```

2.  Instalar dependencias:
    ```bash
    pip install -r requirements.txt
    ```

### Ejecución

Iniciar la interfaz web de Streamlit:

```bash
cd streamlit
streamlit run Home.py
```

El dashboard se abrirá automáticamente en tu navegador (usualmente `http://localhost:8501`). Sigue las instrucciones paso a paso en el menú lateral.

## 📂 Estructura del Proyecto

*   `data/raw`: Imágenes crudas (entrada).
*   `data/processed`: Resultados intermedios y finales organizados por preparado.
*   `streamlit/`: Código de la aplicación web y páginas del pipeline.
*   `paper/methodology/`: Documentación detallada de los métodos científicos.
*   `scripts/`: Scripts auxiliares de pre-procesamiento batch.

## 📄 Metodología

Para una descripción detallada de los algoritmos y librerías utilizadas, consultar [Methodology](paper/methodology/methodology.md).
