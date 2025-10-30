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
    Styled sidebar with:
    - App header and unified group filter.
    - Control panel for global parameters (from calibration.json) with tooltips.
    Use this from any page; style is injected once per run.
    """
    _inject_css()

    with st.sidebar:
        st.markdown("### 🧠 Astrocitos 3D")
        st.markdown("<span class='sidebar-badge'>Análisis morfológico</span>", unsafe_allow_html=True)
        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
        st.caption("Usá el selector de páginas arriba para navegar.")

        # Selector unificado de grupo (aplica a todas las páginas)
        st.markdown("**Filtro global de grupo**")
        grp = st.radio(
            "Grupo",
            options=["Todos", "CTL", "Hipoxia"],
            index={"Todos":0, "CTL":1, "Hipoxia":2}.get(st.session_state.get("group_filter", "Todos"), 0),
            horizontal=True,
            key="__sidebar_group_radio__",
        )
        # Guardar selección en session_state para que todas las páginas la lean
        st.session_state["group_filter"] = grp
        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)

        if show_calibration:
            root = Path(__file__).resolve().parents[1]
            calib_path = root / "calibration.json"
            glob = _read_global_calibration(calib_path)
            st.markdown("**Panel de parámetros globales**")
            with st.form("__control_panel__"):
                tz, ty, tx = st.columns(3)
                with tz:
                    glob["z"] = st.number_input("Voxel Z (µm)", value=float(glob.get("z", 1.0)), min_value=0.001, step=0.001, format="%.3f", help="Espaciado físico del vóxel en Z en micrómetros.")
                with ty:
                    glob["y"] = st.number_input("Voxel Y (µm)", value=float(glob.get("y", 1.0)), min_value=0.001, step=0.001, format="%.3f", help="Espaciado físico del vóxel en Y en micrómetros.")
                with tx:
                    glob["x"] = st.number_input("Voxel X (µm)", value=float(glob.get("x", 1.0)), min_value=0.001, step=0.001, format="%.3f", help="Espaciado físico del vóxel en X en micrómetros.")

                st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
                st.markdown("**Canales y Cellpose**")
                c1, c2 = st.columns(2)
                glob["DAPI_CHANNEL_INDEX"] = c1.number_input("Canal DAPI", value=int(glob.get("DAPI_CHANNEL_INDEX", 0)), min_value=0, step=1, help="Índice de canal para DAPI en el stack Z,C,Y,X.")
                glob["GFAP_CHANNEL_INDEX"] = c2.number_input("Canal GFAP", value=int(glob.get("GFAP_CHANNEL_INDEX", 1)), min_value=0, step=1, help="Índice de canal para GFAP en el stack Z,C,Y,X.")
                c3, c4 = st.columns(2)
                glob["MICROGLIA_CHANNEL_INDEX"] = c3.number_input("Canal Microglía", value=int(glob.get("MICROGLIA_CHANNEL_INDEX", 2)), min_value=0, step=1, help="Índice de canal para Microglía en el stack Z,C,Y,X.")
                glob["NUCLEUS_DIAMETER"] = c4.number_input("Diám. núcleo (px)", value=int(glob.get("NUCLEUS_DIAMETER", 30)), min_value=1, step=1, help="Diámetro estimado de núcleos para Cellpose.")
                glob["CELLPOSE_USE_GPU"] = st.checkbox("Usar GPU en Cellpose", value=bool(glob.get("CELLPOSE_USE_GPU", True)), help="Intentar usar GPU si está disponible durante Cellpose.")

                st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
                st.markdown("**Filtrado GFAP/Microglía**")
                f1, f2 = st.columns(2)
                glob["MAX_DILATION_ITERATIONS"] = f1.number_input("Iter. máx. anillo", value=int(glob.get("MAX_DILATION_ITERATIONS", 10)), min_value=1, step=1, help="Cuántas dilataciones evaluar alrededor del núcleo para estimar señal en cascarón.")
                glob["GFAP_INTENSITY_THRESHOLD"] = f2.number_input("Umbral GFAP", value=int(glob.get("GFAP_INTENSITY_THRESHOLD", 40)), min_value=0, step=1, help="Intensidad mínima en cascarón para retener.")
                f3, f4 = st.columns(2)
                glob["MICROGLIA_INTENSITY_THRESHOLD"] = f3.number_input("Umbral Microglía", value=int(glob.get("MICROGLIA_INTENSITY_THRESHOLD", 200)), min_value=0, step=1, help="Si supera el umbral, se descarta por Microglía.")
                glob["MIN_VOLUME_UM3"] = f4.number_input("Vol. mín. (µm³)", value=float(glob.get("MIN_VOLUME_UM3", 75)), min_value=0.0, step=1.0, help="Volumen mínimo del astrocito final.")

                st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
                st.markdown("**Esqueleto (3D)**")
                s1, s2 = st.columns(2)
                glob["SKELETON_TARGET_ISO_UM"] = s1.number_input("Vóxel iso (µm)", value=float(glob.get("SKELETON_TARGET_ISO_UM", min(float(glob.get("z",1.0)), float(glob.get("y",1.0)), float(glob.get("x",1.0))))), min_value=0.05, step=0.01, format="%.2f", help="Tamaño de vóxel al que se remuestrea la ROI para esqueletización.")
                glob["SKELETON_PADDING_UM"] = s2.number_input("Padding ROI (µm)", value=float(glob.get("SKELETON_PADDING_UM", 0.5)), min_value=0.0, step=0.5, help="Margen alrededor de la ROI del núcleo.")
                s3, s4 = st.columns(2)
                glob["SKELETON_SEED_DILATE_UM"] = s3.number_input("Dilatación semilla (µm)", value=float(glob.get("SKELETON_SEED_DILATE_UM", 1.0)), min_value=0.0, step=0.5, help="Dilatación del núcleo en isotrópico para seleccionar GFAP conectado.")
                glob["SKELETON_CONNECTIVITY"] = s4.selectbox("Conectividad 3D", options=[6,26], index=1 if int(glob.get("SKELETON_CONNECTIVITY",26))==26 else 0, help="Conectividad para componente conectado.")
                s5, s6 = st.columns(2)
                glob["SKELETON_CLOSING_UM"] = s5.number_input("Cierre morfológico (µm)", value=float(glob.get("SKELETON_CLOSING_UM", 0.5)), min_value=0.0, step=0.1, help="Cierre para unir hebras antes de esqueletizar.")
                glob["SKELETON_MAX_RADIUS_UM"] = s6.number_input("Radio máx. (µm)", value=float(glob.get("SKELETON_MAX_RADIUS_UM", 100.0)), min_value=0.0, step=1.0, help="Limita el dominio por distancia desde el núcleo.")
                s7, s8 = st.columns(2)
                glob["SKELETON_THRESHOLD_MODE"] = s7.selectbox("Umbral GFAP", options=["otsu","manual"], index=(0 if str(glob.get("SKELETON_THRESHOLD_MODE","otsu")).lower().startswith("otsu") else 1), help="Modo de umbral en la ROI.")
                glob["SKELETON_MANUAL_THRESHOLD"] = s8.number_input("Umbral manual", value=float(glob.get("SKELETON_MANUAL_THRESHOLD", 30.0)), min_value=0.0, step=1.0, help="Valor si el modo es 'manual'.")
                s9, s10 = st.columns(2)
                glob["SKELETON_TERRITORY_VORONOI"] = s9.checkbox("Territorios Voronoi", value=bool(glob.get("SKELETON_TERRITORY_VORONOI", True)), help="Restringir por territorio más cercano al núcleo")
                glob["SKELETON_TERRITORY_EXCLUSION_UM"] = s10.number_input("Exclusión frontera (µm)", value=float(glob.get("SKELETON_TERRITORY_EXCLUSION_UM", 2.0)), min_value=0.0, step=0.1, help="‘Gap’ en la frontera de territorios.")
                s11, s12 = st.columns(2)
                glob["SKELETON_DOMAIN_VOLUME_SOURCE"] = s11.selectbox("Volumen de dominio", options=["gfap","voronoi"], index=(0 if str(glob.get("SKELETON_DOMAIN_VOLUME_SOURCE","gfap")).lower().startswith("gfap") else 1), help="Fuente para volumen del dominio.")
                glob["SKELETON_TUBE_RADIUS_UM"] = s12.number_input("Radio del tubo (µm)", value=float(glob.get("SKELETON_TUBE_RADIUS_UM", 1.5)), min_value=0.0, step=0.1, help="Radio para medir señal alrededor del esqueleto.")
                s13, s14 = st.columns(2)
                glob["SKELETON_PRUNE_ENABLE"] = s13.checkbox("Pruning espículas", value=bool(glob.get("SKELETON_PRUNE_ENABLE", True)), help="Eliminar ramas terminales cortas.")
                glob["SKELETON_PRUNE_MIN_LEN_UM"] = s14.number_input("Long. mín. espícula (µm)", value=float(glob.get("SKELETON_PRUNE_MIN_LEN_UM", 2.0)), min_value=0.0, step=0.5, help="Umbral para pruning.")

                st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
                st.markdown("**Parámetros de Sholl**")
                sh1, sh2 = st.columns(2)
                glob["SHOLL_MIN_RADIUS_UM"] = sh1.number_input("Radio mín. (µm)", value=float(glob.get("SHOLL_MIN_RADIUS_UM", 0.0)), min_value=0.0, step=0.5, help="Radio desde el centroide para el primer anillo.")
                glob["SHOLL_MAX_RADIUS_UM"] = sh2.number_input("Radio máx. (µm)", value=float(glob.get("SHOLL_MAX_RADIUS_UM", 100.0)), min_value=0.0, step=1.0, help="Radio máximo a evaluar.")
                sh3, _ = st.columns([1,1])
                glob["SHOLL_STEP_UM"] = sh3.number_input("Paso entre anillos (µm)", value=float(glob.get("SHOLL_STEP_UM", 20.0)), min_value=0.1, step=0.1, help="Separación entre anillos y semiancho del cascarón.")

                saved = st.form_submit_button("💾 Guardar parámetros globales")
                if saved:
                    try:
                        # limpiar claves no deseadas
                        if "GROUP" in glob:
                            glob.pop("GROUP", None)
                        calib_path.write_text(json.dumps(glob, indent=2))
                        st.success("Parámetros guardados en calibration.json")
                    except Exception as e:
                        st.error(f"No se pudieron guardar parámetros: {e}")

        # Ejecución global del experimento
        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
        st.markdown("**Ejecución global**")
        with st.form("__global_run__"):
            scope_opt = st.selectbox("Ámbito", options=["Grupo CTL", "Grupo Hipoxia", "Todos"], index=2)
            start_from = st.selectbox("Correr desde paso", options=["01", "02", "03", "04", "05"], index=0, help="Recomputa desde el paso elegido y continúa hasta el final")
            overwrite = st.checkbox("Sobrescribir desde el paso", value=True)
            submit = st.form_submit_button("⚙️ Correr experimento (selección)")
            if submit:
                try:
                    from .runner import run_scope, read_calibration
                    root_repo = Path(__file__).resolve().parents[2]  # repo root (para data/raw)
                    cal_path = Path(__file__).resolve().parents[1] / "calibration.json"  # streamlit/calibration.json
                    cal = read_calibration(cal_path)
                    scope = "all"
                    group = None
                    if scope_opt.startswith("Grupo CTL"):
                        scope = "group"; group = "CTL"
                    elif scope_opt.startswith("Grupo Hipoxia"):
                        scope = "group"; group = "Hipoxia"
                    res = run_scope(root_repo, scope=scope, start_step=start_from, cal=cal, selected=None, group=group, overwrite_from_step=overwrite)
                    ok = sum(1 for _, stt in res if not stt.get("error"))
                    st.success(f"Ejecución completada: {ok}/{len(res)} preparados procesados.")
                except Exception as e:
                    st.error(f"No se pudo ejecutar: {e}")
