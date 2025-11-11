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
    Sidebar con:
    - Header de la app y filtro de grupo unificado
    - Panel de parámetros globales del pipeline 2D nativo
    - Herramienta de ejecución batch (individual/grupo/todos)
    """
    _inject_css()

    with st.sidebar:
        st.markdown("### 🧠 Astrocitos 3D")
        st.markdown("<span class='sidebar-badge'>Pipeline 2D Nativo</span>", unsafe_allow_html=True)
        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
        
        # Info rápida del flujo
        with st.expander("ℹ️ Flujo del Pipeline (4 pasos)", expanded=False):
            st.markdown("""
            **01. Calibración + Visualización**
            - Detección de escala física (µm)
            
            **02. Segmentación Nuclear**
            - Umbral Otsu + Cellpose 3D
            
            **03. Filtrado de Astrocitos**
            - Clasificación GFAP+ (relativa)
            - Filtro por volumen mínimo
            
            **04. Skeleton + Sholl 2D**
            - Proyección 3D→2D
            - Territorios Voronoi
            - Esqueletización + Sholl integrado
            """)

        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)

        # Selector unificado de grupo
        st.markdown("**🔍 Filtro de Grupo**")
        grp = st.radio(
            "Grupo",
            options=["Todos", "CTL", "Hipoxia"],
            index={"Todos":0, "CTL":1, "Hipoxia":2}.get(st.session_state.get("group_filter", "Todos"), 0),
            horizontal=True,
            key="__sidebar_group_radio__",
        )
        st.session_state["group_filter"] = grp
        st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)

        if show_calibration:
            root = Path(__file__).resolve().parents[1]
            calib_path = root / "calibration.json"
            glob = _read_global_calibration(calib_path)
            
            st.markdown("**⚙️ Configuración Global**")
            
            with st.expander("📏 Calibración Espacial", expanded=False):
                with st.form("__calib_form__"):
                    st.caption("Espaciado físico de vóxeles")
                    tz, ty, tx = st.columns(3)
                    glob["z"] = tz.number_input("Z (µm)", value=float(glob.get("z", 0.38)), min_value=0.001, step=0.01, format="%.3f")
                    glob["y"] = ty.number_input("Y (µm)", value=float(glob.get("y", 0.38)), min_value=0.001, step=0.01, format="%.3f")
                    glob["x"] = tx.number_input("X (µm)", value=float(glob.get("x", 0.38)), min_value=0.001, step=0.01, format="%.3f")
                    
                    st.caption("Índices de canales (base 0)")
                    c1, c2 = st.columns(2)
                    glob["DAPI_CHANNEL_INDEX"] = c1.number_input("DAPI", value=int(glob.get("DAPI_CHANNEL_INDEX", 0)), min_value=0, step=1)
                    glob["GFAP_CHANNEL_INDEX"] = c2.number_input("GFAP", value=int(glob.get("GFAP_CHANNEL_INDEX", 1)), min_value=0, step=1)
                    
                    if st.form_submit_button("💾 Guardar calibración"):
                        try:
                            calib_path.write_text(json.dumps(glob, indent=2))
                            st.success("✅ Calibración guardada")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error: {e}")

            with st.expander("🔬 Segmentación Nuclear (Paso 02)", expanded=False):
                with st.form("__cellpose_form__"):
                    glob["NUCLEUS_DIAMETER"] = st.number_input(
                        "Diámetro núcleo (px)", 
                        value=int(glob.get("NUCLEUS_DIAMETER", 30)), 
                        min_value=5, step=1,
                        help="Diámetro esperado de núcleos para Cellpose"
                    )
                    glob["CELLPOSE_USE_GPU"] = st.checkbox(
                        "Usar GPU", 
                        value=bool(glob.get("CELLPOSE_USE_GPU", True)),
                        help="Acelera Cellpose si hay GPU disponible"
                    )
                    
                    if st.form_submit_button("💾 Guardar segmentación"):
                        try:
                            calib_path.write_text(json.dumps(glob, indent=2))
                            st.success("✅ Parámetros guardados")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error: {e}")

            with st.expander("🧪 Filtrado de Astrocitos (Paso 03)", expanded=False):
                with st.form("__filter_form__"):
                    st.caption("Clasificación GFAP+ (relativa a fondo)")
                    glob["GFAP_STD_DEV_THRESHOLD"] = st.number_input(
                        "Umbral GFAP (N° StdDev)", 
                        value=float(glob.get("GFAP_STD_DEV_THRESHOLD", 3.0)), 
                        min_value=0.0, step=0.1,
                        help="N° de desviaciones estándar sobre el fondo para GFAP+"
                    )
                    
                    f1, f2 = st.columns(2)
                    glob["SHELL_RADIUS_UM"] = f1.number_input(
                        "Radio shell (µm)", 
                        value=float(glob.get("SHELL_RADIUS_UM", 2.0)), 
                        min_value=0.1, step=0.1,
                        help="Radio del anillo perinuclear para medir GFAP"
                    )
                    glob["MAX_DILATION_ITERATIONS"] = f2.number_input(
                        "Máx. iteraciones", 
                        value=int(glob.get("MAX_DILATION_ITERATIONS", 10)), 
                        min_value=1, step=1,
                        help="Iteraciones de dilatación para formar el shell"
                    )
                    
                    st.caption("Filtro por volumen")
                    glob["MIN_VOLUME_UM3"] = st.number_input(
                        "Volumen mínimo (µm³)", 
                        value=float(glob.get("MIN_VOLUME_UM3", 75)), 
                        min_value=0.0, step=5.0,
                        help="Volumen mínimo del núcleo para ser considerado válido"
                    )
                    
                    st.caption("Avanzado (opcional)")
                    glob["GFAP_INTENSITY_THRESHOLD"] = st.number_input(
                        "Umbral intensidad GFAP (fallback)", 
                        value=float(glob.get("GFAP_INTENSITY_THRESHOLD", 50.0)), 
                        min_value=0.0, step=1.0,
                        help="Umbral absoluto usado si clasificación relativa falla"
                    )
                    
                    if st.form_submit_button("💾 Guardar filtrado"):
                        try:
                            calib_path.write_text(json.dumps(glob, indent=2))
                            st.success("✅ Filtros guardados")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error: {e}")

            with st.expander("🗺️ Pipeline 2D: Skeleton + Sholl (Paso 04)", expanded=False):
                with st.form("__pipeline_2d_form__"):
                    st.caption("Proyección 3D → 2D")
                    glob["PROJECTION_2D_METHOD"] = st.selectbox(
                        "Método de proyección", 
                        options=["max", "mean", "sum"],
                        index=["max", "mean", "sum"].index(str(glob.get("PROJECTION_2D_METHOD", "max")).lower()),
                        help="max=preserva procesos débiles, mean=reduce ruido, sum=maximiza señal"
                    )
                    
                    st.caption("Territorios Voronoi")
                    glob["TERRITORY_EXCLUSION_UM"] = st.number_input(
                        "Gap de exclusión (µm)", 
                        value=float(glob.get("TERRITORY_EXCLUSION_UM", 1.0)), 
                        min_value=0.0, step=0.5,
                        help="Zona de exclusión entre territorios Voronoi"
                    )
                    
                    st.caption("Dominio de análisis")
                    glob["MAX_RADIUS_FROM_NUCLEUS_UM"] = st.number_input(
                        "Radio máximo desde núcleo (µm)", 
                        value=float(glob.get("MAX_RADIUS_FROM_NUCLEUS_UM", 100.0)), 
                        min_value=10.0, max_value=200.0, step=10.0,
                        help="Radio máximo para esqueletizar. Procesos más allá son ignorados (evita señal de fondo)"
                    )
                    
                    st.caption("Conexión de fragmentos GFAP")
                    glob["CONNECT_SKELETON_FRAGMENTS"] = st.checkbox(
                        "Conectar fragmentos", 
                        value=bool(glob.get("CONNECT_SKELETON_FRAGMENTS", True)),
                        help="Conecta señal GFAP discontinua cercana"
                    )
                    if glob["CONNECT_SKELETON_FRAGMENTS"]:
                        glob["CONNECTION_RADIUS_UM"] = st.number_input(
                            "Radio de conexión (µm)", 
                            value=float(glob.get("CONNECTION_RADIUS_UM", 0.5)), 
                            min_value=0.0, step=0.1,
                            help="Distancia máxima para conectar fragmentos"
                        )
                    
                    st.caption("Análisis de Sholl 2D")
                    sh1, sh2, sh3 = st.columns(3)
                    glob["SHOLL_MIN_RADIUS_UM"] = sh1.number_input(
                        "Mín (µm)", 
                        value=float(glob.get("SHOLL_MIN_RADIUS_UM", 0.0)), 
                        min_value=0.0, step=0.5
                    )
                    glob["SHOLL_MAX_RADIUS_UM"] = sh2.number_input(
                        "Máx (µm)", 
                        value=float(glob.get("SHOLL_MAX_RADIUS_UM", 50.0)), 
                        min_value=0.0, step=1.0
                    )
                    glob["SHOLL_STEP_UM"] = sh3.number_input(
                        "Paso (µm)", 
                        value=float(glob.get("SHOLL_STEP_UM", 2.0)), 
                        min_value=0.1, step=0.1
                    )
                    
                    if st.form_submit_button("💾 Guardar pipeline 2D"):
                        try:
                            calib_path.write_text(json.dumps(glob, indent=2))
                            st.success("✅ Pipeline 2D guardado")
                            st.rerun()
                        except Exception as e:
                            st.error(f"Error: {e}")

            # Ejecución global
            st.markdown("<div class='sidebar-sep'></div>", unsafe_allow_html=True)
            st.markdown("**▶️ Ejecución Batch**")
            
            with st.form("__batch_run__"):
                st.caption("Ejecutar pipeline completo sobre múltiples preparados")
                
                scope_opt = st.selectbox(
                    "Ámbito de ejecución",
                    options=["📁 Todos los preparados", "🔵 Solo grupo CTL", "🔴 Solo grupo Hipoxia"],
                    index=0,
                    help="Selecciona qué preparados procesar"
                )
                
                start_from = st.selectbox(
                    "Comenzar desde paso",
                    options=["01", "02", "03", "04"],
                    index=0,
                    help="Paso inicial del pipeline (recomputa desde aquí)"
                )
                
                overwrite = st.checkbox(
                    "Sobrescribir archivos existentes",
                    value=True,
                    help="Si está desmarcado, solo procesa preparados sin salidas"
                )
                
                submit = st.form_submit_button("⚙️ Ejecutar pipeline", use_container_width=True)
                
                if submit:
                    try:
                        from ui import runner
                        
                        root_repo = Path(__file__).resolve().parents[2]
                        cal = _read_global_calibration(calib_path)
                        
                        # Parsear ámbito
                        if "CTL" in scope_opt:
                            scope, group = "group", "CTL"
                        elif "Hipoxia" in scope_opt:
                            scope, group = "group", "Hipoxia"
                        else:
                            scope, group = "all", None
                        
                        # Ejecutar con barra de progreso
                        res = runner.run_scope(
                            root_repo, 
                            scope=scope, 
                            start_step=start_from, 
                            cal=cal, 
                            selected=None, 
                            group=group, 
                            overwrite_from_step=overwrite
                        )
                        
                        # Resumen de resultados
                        ok = sum(1 for _, status in res if not status.get("error"))
                        total = len(res)
                        
                        if ok == total:
                            st.success(f"✅ {ok}/{total} preparados procesados exitosamente")
                        else:
                            st.warning(f"⚠️ {ok}/{total} preparados completados ({total-ok} con errores)")
                            
                    except Exception as e:
                        st.error(f"❌ Error en ejecución: {e}")
                        st.exception(e)
