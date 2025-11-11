import streamlit as st
from pathlib import Path
from ui.sidebar import render_sidebar
from ui.utils import detect_group
import json
from datetime import datetime
import pandas as pd
import subprocess # <--- CORRECCIÓN
import os # <--- CORRECCIÓN

# Welcome page content based on project README
st.set_page_config(page_title="Astrocitos 3D - Análisis", page_icon="🧠", layout="wide")

st.title("Análisis Morfológico 3D de Astrocitos")
render_sidebar(show_calibration=True)

st.markdown(
    """
    **Proyecto:** Reconstrucción, segmentación y análisis morfológico 3D de astrocitos (GFAP) y núcleos (DAPI) a partir de imágenes `.lif`/`.tif`.

    **Tecnologías:** Python, Streamlit, Napari, Cellpose, scikit-image, SciPy, Skan

    ---

    ### Flujo del pipeline
    El flujo de procesamiento unificado (4 pasos):
    1. **Calibración (01):** Lectura de `calibration.json` para µm/px → Máscara Otsu de fondo DAPI
    2. **Segmentación (02):** Cellpose 3D → `02_cellpose_mask.tif`
    3. **Filtrado (03):** 
       - Filtrado GFAP relativo (StdDev sobre fondo) → `03_gfap_filtered_mask.tif`
       - Filtrado por tamaño físico → `04_final_astrocytes_mask.tif`
    4. **Esqueletización + Sholl (04):** **Pipeline 2D unificado**
       - Proyección 3D→2D (max projection)
       - Territorios Voronoi con zona de exclusión
       - Esqueletización 2D por territorio con conexión de fragmentos
       - Análisis de Sholl 2D nativo integrado (SKAN)
       - Genera `05_skeleton_labels_2d.tif`, `sholl_2d_native.csv`, `sholl_summary.csv`, `sholl_rings_2d_native.json`
    
    **Ventajas del flujo 2D nativo:**
    - ✅ Resolución XY completa (0.38 µm) sin degradación
    - ✅ Sholl 2D más preciso y eficiente
    - ✅ Territorios astrocitarios bien definidos para preparados planos
    - ✅ Más rápido (~10x) que esqueletización 3D con remuestreo isotrópico
    """
)

st.info("Usá el menú lateral para navegar el pipeline. La calibración y parámetros globales viven en streamlit/calibration.json; los resultados por preparado quedan en data/processed/<preparado>/.")

# Inspector/Editor de configuración unificada
root = Path(__file__).resolve().parents[1]
calib_path = root / "streamlit" / "calibration.json"

st.markdown("---")
st.subheader("Inspector de configuración global (calibration.json)")

def _load_calib() -> dict:
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())
        except Exception:
            return {}
    return {}

def _save_calib(data: dict):
    calib_path.parent.mkdir(parents=True, exist_ok=True)
    calib_path.write_text(json.dumps(data, indent=2))

cfg = _load_calib()

colv1, colv2 = st.columns(2)
with colv1:
    st.caption("Vista actual (solo lectura)")
    st.json(cfg or {})
with colv2:
    st.caption("Editar como JSON (avanzado)")
    raw = st.text_area("Contenido de calibration.json", value=json.dumps(cfg, indent=2), height=260)
    bcol1, bcol2 = st.columns(2)
    with bcol1:
        if st.button("💾 Guardar cambios"):
            try:
                new_cfg = json.loads(raw) if raw.strip() else {}
            except Exception as e:
                st.error(f"JSON inválido: {e}")
            else:
                _save_calib(new_cfg)
                st.success("Cambios guardados. Recargá la página (F5) para aplicar.")
                st.rerun()
    with bcol2:
        if st.button("🧰 Exportar backup con timestamp"):
            ts = datetime.now().strftime("%Y%m%d-%H%M%S")
            backup_path = calib_path.with_name(f"calibration.backup-{ts}.json")
            backup_path.write_text(json.dumps(cfg, indent=2))
            st.info(f"Backup exportado: {backup_path.relative_to(root)}")

# --- Dashboard del pipeline por preparado (ACTUALIZADO) ---
st.markdown("---")
st.subheader("Estado del pipeline por preparado")

raw_dir = root / "data" / "raw"
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])

if files:
    stem_to_group_all = {p.stem: detect_group(p, root) for p in files}

    group_filter = st.session_state.get("group_filter", "Todos")
    files_filtered = files if group_filter == "Todos" else [p for p in files if stem_to_group_all.get(p.stem, "CTL") == group_filter]

    rows = []
    for p in files_filtered:
        od = root / "data" / "processed" / p.stem
        rows.append({
            "prepared": p.stem,
            "group": stem_to_group_all.get(p.stem, "CTL"),
            "02_Nucleos": (od/"02_cellpose_mask.tif").exists(),
            "03_Metricas_Nucleo": (od/"03_nucleus_metrics.csv").exists(),
            "04_Astrocitos": (od/"04_final_astrocytes_mask.tif").exists(),
            "04_Skeleton_2D": (od/"05_skeleton_labels_2d.tif").exists(),
            "04_Sholl_2D": (od/"sholl_2d_native.csv").exists(),
        })
    
    df_state = pd.DataFrame(rows)
    
    if not df_state.empty:
        def mark(col):
            return df_state[col].map(lambda v: "✅" if bool(v) else "—")
        
        view = df_state.copy()
        # Actualizar los nombres de las columnas para el dashboard
        cols_to_check = [
            "02_Nucleos", "03_Metricas_Nucleo", "04_Astrocitos", 
            "04_Skeleton_2D", "04_Sholl_2D"
        ]
        for c in cols_to_check:
            if c in view.columns:
                view[c] = mark(c)
                
        st.dataframe(view.set_index(["prepared","group"]), use_container_width=True)
    else:
        st.info("No hay preparados para mostrar con el filtro actual.")
else:
    st.info("No se encontraron archivos en data/raw. Cargá tu dataset para ver el dashboard.")