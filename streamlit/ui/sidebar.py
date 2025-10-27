from __future__ import annotations
from pathlib import Path
import json
import streamlit as st


def _inject_css():
    st.markdown(
        """
        <style>
        [data-testid="stSidebar"] > div:first-child {
            background: linear-gradient(180deg, #0f172a 0%, #111827 100%);
        }
        [data-testid="stSidebar"] * { color: #e5e7eb; }
        [data-testid="stSidebar"] h1, [data-testid="stSidebar"] h2, [data-testid="stSidebar"] h3 {
            color: #f3f4f6;
        }
        .sidebar-badge { display:inline-block; padding:4px 8px; border-radius:6px; background:#1f2937; color:#d1d5db; font-size:12px; }
        .sidebar-small { font-size:12px; color:#9ca3af; }
        .sidebar-sep { border-bottom:1px solid #374151; margin:8px 0 12px 0; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def _read_global_calibration(calib_path: Path) -> dict:
    try:
        if calib_path.exists():
            return json.loads(calib_path.read_text())
    except Exception:
        pass
    return {}


def render_sidebar(show_calibration: bool = True):
    """
    Renders a styled sidebar with app header and optional global calibration status.
    This function is safe to call from any page, the style is injected once per run.
    """
    _inject_css()

    with st.sidebar:
        st.markdown("### 🧠 Astrocitos 3D")
        st.markdown("<span class='sidebar-badge'>Análisis morfológico</span>", unsafe_allow_html=True)
        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
        st.caption("Usá el selector de páginas arriba para navegar.")

        if show_calibration:
            root = Path(__file__).resolve().parents[1]
            calib_path = root / "calibration.json"
            glob = _read_global_calibration(calib_path)
            if glob:
                z = glob.get("z")
                y = glob.get("y")
                x = glob.get("x")
                st.markdown("**Calibración global activa (µm):**")
                st.write({"z": z, "y": y, "x": x})
            else:
                st.markdown("<span class='sidebar-small'>No hay calibración global guardada aún.</span>", unsafe_allow_html=True)

            # Mostrar parámetros generales del experimento si existen
            try:
                exp_params_path = root / "experiment_params.json"
                if exp_params_path.exists():
                    st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
                    st.markdown("**Parámetros generales del experimento:**")
                    with exp_params_path.open("r") as f:
                        exp_params = json.load(f)
                    # Render simple dict; mantener orden de claves comunes si están
                    preferred_order = [
                        "NUCLEUS_DIAMETER",
                        "MAX_DILATION_ITERATIONS",
                        "GFAP_INTENSITY_THRESHOLD",
                        "MICROGLIA_INTENSITY_THRESHOLD",
                        "MIN_VOLUME_UM3",
                    ]
                    # Reordenar si aplica
                    ordered = {k: exp_params[k] for k in preferred_order if k in exp_params}
                    for k, v in exp_params.items():
                        if k not in ordered:
                            ordered[k] = v
                    st.write(ordered)
                else:
                    st.markdown("<span class='sidebar-small'>No hay parámetros generales guardados aún.</span>", unsafe_allow_html=True)
            except Exception:
                st.markdown("<span class='sidebar-small'>No se pudieron leer los parámetros generales.</span>", unsafe_allow_html=True)
