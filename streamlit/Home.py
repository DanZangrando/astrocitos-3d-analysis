import streamlit as st
from pathlib import Path
from ui.sidebar import render_sidebar
import json
from datetime import datetime
import numpy as np
import pandas as pd
import tifffile
import altair as alt
from scipy import stats

# Welcome page content based on project README
st.set_page_config(page_title="Astrocitos 3D - Análisis", page_icon="🧠", layout="wide")

st.title("Análisis Morfológico 3D de Astrocitos")
render_sidebar(show_calibration=True)

st.markdown(
    """
    **Proyecto:** Reconstrucción, segmentación y análisis morfológico 3D de astrocitos (GFAP) y núcleos (DAPI) a partir de imágenes `.lif`/`.tif`.

    **Tecnologías:** Python, Streamlit, Napari, Cellpose, scikit-image, SciPy, Skan

    ---

    ### Estado actual del pipeline
    Hemos implementado la esqueletización 3D previa al análisis de Sholl. El flujo es:
    1. Calibración física global (µm) — Z, Y, X en `streamlit/calibration.json`.
    2. Segmentación de núcleos (DAPI) — Otsu opcional + Cellpose 3D.
    3. Filtrado por GFAP/Microglía y limpieza por tamaño físico (µm³).
    4. Esqueletización 3D por célula con re-muestreo isotrópico, umbral Otsu/Manual, cierre morfológico, dilatación de semilla, conectividad, radio máximo y:
       - Resolución de solapamientos por cercanía al núcleo.
       - Territorios Voronoi con zona de exclusión en fronteras para evitar entrelazados ambiguos.
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
                st.success("Cambios guardados en streamlit/calibration.json. Recargá la página para verlos reflejados en el sidebar.")
    with bcol2:
        if st.button("🧰 Exportar backup con timestamp"):
            ts = datetime.now().strftime("%Y%m%d-%H%M%S")
            backup_path = calib_path.with_name(f"calibration.backup-{ts}.json")
            backup_path.write_text(json.dumps(cfg, indent=2))
            st.info(f"Backup exportado: {backup_path.relative_to(root)}")

# --- Dashboard del pipeline por preparado (solo tabla de checkeo) ---
st.markdown("---")
st.subheader("Estado del pipeline por preparado")

def _detect_group(p: Path, root: Path) -> str:
    try:
        rel = str(p.relative_to(root)).lower()
    except Exception:
        rel = str(p).lower()
    if "/hip/" in rel:
        return "Hipoxia"
    return "CTL"

raw_dir = root / "data" / "raw"
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if files:
    stem_to_group_all = {p.stem: _detect_group(p, root) for p in files}

    group_filter = st.session_state.get("group_filter", "Todos")
    files_filtered = files if group_filter == "Todos" else [p for p in files if stem_to_group_all.get(p.stem, "CTL") == group_filter]

    import pandas as _pd
    rows = []
    for p in files_filtered:
        od = root / "data" / "processed" / p.stem
        rows.append({
            "prepared": p.stem,
            "group": stem_to_group_all.get(p.stem, "CTL"),
            "01_otsu": (od/"01_otsu_mask.tif").exists(),
            "02_cellpose": (od/"02_cellpose_mask.tif").exists(),
            "03_gfap": (od/"03_gfap_microglia_filtered_mask.tif").exists(),
            "04_final": (od/"04_final_astrocytes_mask.tif").exists(),
            "05_skeleton": (od/"05_skeleton_labels.tif").exists(),
            "06_sholl": (od/"sholl.csv").exists(),
        })
    df_state = _pd.DataFrame(rows)
    if not df_state.empty:
        def mark(col):
            return df_state[col].map(lambda v: "✅" if bool(v) else "—")
        view = df_state.copy()
        for c in ["01_otsu","02_cellpose","03_gfap","04_final","05_skeleton","06_sholl"]:
            view[c] = mark(c)
        st.dataframe(view.set_index(["prepared","group"]), use_container_width=True)
else:
    st.info("No se encontraron archivos en data/raw. Cargá tu dataset para ver el dashboard.")
