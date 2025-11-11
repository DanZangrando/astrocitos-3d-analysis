import json
from pathlib import Path
import altair as alt
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
import streamlit as st
from ui.sidebar import render_sidebar
from ui.plots import GROUP_SCALE, PALETTE, boxplot_with_ticks, apply_theme

# Raíz del repo
ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"

st.title("Comparación entre Grupos (CTL vs Hipoxia)")
render_sidebar(show_calibration=True)
apply_theme(st)

# --- Funciones Estadísticas (Sin cambios) ---
def run_comparison_test(a: np.ndarray, b: np.ndarray, alpha: float = 0.05) -> dict:
    """
    Realiza una comparación estadística robusta entre dos grupos (a y b).
    Comprueba Normalidad (Shapiro) -> Elige Welch (normal) o Mann-Whitney (no normal).
    """
    results = {
        'test_name': None, 'comparison_target': None, 'stat_label': None, 
        'stat': np.nan, 'p_value': np.nan, 'is_a_normal': False, 
        'is_b_normal': False, 'message': None, 'n_a': a.size, 'n_b': b.size
    }
    # Manejar casos borde
    if a.size < 3 or b.size < 3:
        results['message'] = f"No hay suficientes datos (CTL n={a.size}, Hipoxia n={b.size}) para un test robusto."
        return results
    a_var, b_var = np.var(a, ddof=1), np.var(b, ddof=1)

    # Test de Normalidad
    if a_var > 1e-10:
        try:
            _, shapiro_a_p = stats.shapiro(a); results['is_a_normal'] = shapiro_a_p > alpha
        except stats.common.ShapiroError: results['is_a_normal'] = False
    if b_var > 1e-10:
        try:
            _, shapiro_b_p = stats.shapiro(b); results['is_b_normal'] = shapiro_b_p > alpha
        except stats.common.ShapiroError: results['is_b_normal'] = False
    
    # Elegir Test
    if results['is_a_normal'] and results['is_b_normal']:
        # Paramétrico: Welch's t-test
        results.update({'test_name': "Test t de Welch", 'comparison_target': "medias", 'stat_label': 't'})
        stat, p_value = stats.ttest_ind(a, b, equal_var=False, nan_policy='omit')
    else:
        # No Paramétrico: Mann-Whitney U
        results.update({'test_name': "U de Mann-Whitney", 'comparison_target': "medianas", 'stat_label': 'U'})
        try:
            if a_var < 1e-10 and b_var < 1e-10 and np.mean(a) == np.mean(b):
                 stat, p_value = np.nan, 1.0; results['message'] = "Datos idénticos."
            else:
                stat, p_value = stats.mannwhitneyu(a, b, alternative='two-sided')
        except ValueError as e:
            results['message'] = f"Error: {e}"; stat, p_value = np.nan, np.nan
    
    results['stat'] = float(stat); results['p_value'] = float(p_value)
    return results

def format_test_results(results: dict) -> str:
    """Formatea los resultados del test en un string markdown legible."""
    if results.get('message'): return f"*{results['message']}*"
    if pd.isna(results['p_value']): return "*Error al calcular estadísticas.*"
    
    p_value = results['p_value']
    p_display = f"p={p_value:.3g}" if p_value >= 0.001 else "p < 0.001"
    
    report_str = (
        f"**Test (sobre N por preparado):** {results['test_name']} (comparando {results['comparison_target']}). "
        f"**{results['stat_label']}={results['stat']:.2f}, {p_display}** "
        f"*(N(CTL)={results['n_a']}, N(Hip)={results['n_b']})*"
    )
    return report_str

def _download_df_button(df: pd.DataFrame, filename_base: str, label: str = "Descargar CSV"):
    """Crea un botón de descarga para un DataFrame."""
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(label=label, data=csv, file_name=f"{filename_base}.csv", mime="text/csv")

def list_raw_images(raw_dir: Path) -> list[Path]:
    """Lista todas las imágenes .tif/.tiff en los subdirectorios CTL y hip."""
    files = []
    for sub in ["CTL", "hip"]:
        d = raw_dir / sub
        if not d.exists(): continue
        for ext in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
            files.extend(sorted(d.rglob(ext)))
    return files

def detect_group_from_path(p: Path) -> str:
    """Detecta el grupo (CTL o Hipoxia) desde la ruta del archivo."""
    return 'Hipoxia' if '/hip/' in str(p).replace('\\', '/').lower() else 'CTL'

# --- Carga de Datos (REFACTORIZADA) ---
@st.cache_data(ttl=60) # Cachear por 1 min (reducido para reflejar cambios más rápido)
def load_all_metrics(raw_dir: Path, proc_dir: Path):
    """
    Carga y fusiona todas las métricas (núcleo, esqueleto, sholl) de todos 
    los preparados procesados.
    """
    files = list_raw_images(raw_dir)
    rows_skel = []
    rows_sholl = []
    rows_nuc = []
    
    for p in files:
        od = proc_dir / p.stem
        group = detect_group_from_path(p)
        
        # 1. Cargar métricas topológicas de esqueleto (skeletons/summary.csv)
        ps = od / "skeletons" / "summary.csv"
        if ps.exists():
            try:
                df = pd.read_csv(ps); df['prepared'] = p.name; df['group'] = group
                rows_skel.append(df)
            except Exception: pass
        
        # 2. Cargar métricas de Sholl 2D (sholl_summary.csv)
        pss = od / "sholl_summary.csv"
        if pss.exists():
            try:
                df = pd.read_csv(pss); df['prepared'] = p.name; df['group'] = group
                rows_sholl.append(df)
            except Exception: pass
            
        # 3. Cargar métricas de núcleo (03_nucleus_metrics.csv)
        pn = od / "03_nucleus_metrics.csv"
        if pn.exists():
            try:
                df = pd.read_csv(pn); df['prepared'] = p.name; df['group'] = group
                rows_nuc.append(df)
            except Exception: pass

    if not rows_skel and not rows_sholl and not rows_nuc:
        return None, None

    # --- Fusión y Agregación ---
    df_skel_all = pd.concat(rows_skel, ignore_index=True) if rows_skel else pd.DataFrame()
    df_nuc_all = pd.concat(rows_nuc, ignore_index=True) if rows_nuc else pd.DataFrame()
    df_sholl_all = pd.concat(rows_sholl, ignore_index=True) if rows_sholl else pd.DataFrame()
    
    # Estrategia de fusión: empezar con el DataFrame que tenga más datos
    # y hacer left joins para preservar todas las filas
    
    if not df_skel_all.empty:
        # Comenzar con skeleton (más completo)
        df_por_celula = df_skel_all.copy()
        
        # Agregar Sholl
        if not df_sholl_all.empty:
            df_por_celula = pd.merge(df_por_celula, df_sholl_all, 
                                     on=['label', 'prepared', 'group'], how='left')
    elif not df_sholl_all.empty:
        # Si no hay skeleton pero sí Sholl, comenzar con Sholl
        df_por_celula = df_sholl_all.copy()
    else:
        # Si no hay ni skeleton ni Sholl, usar DataFrame vacío
        df_por_celula = pd.DataFrame()
    
    # Agregar métricas de núcleo
    if not df_nuc_all.empty:
        # Filtrar solo candidatos a astrocitos
        df_nuc_astros = df_nuc_all[df_nuc_all['is_astrocyte_candidate'] == True].copy()
        
        if not df_nuc_astros.empty and 'label' in df_nuc_astros.columns:
            df_nuc_astros['label'] = df_nuc_astros['label'].astype(int)
            
            if not df_por_celula.empty and 'label' in df_por_celula.columns:
                # Merge con datos existentes
                df_por_celula = pd.merge(df_por_celula, df_nuc_astros, 
                                        on=['label', 'prepared', 'group'], how='left')
            elif df_por_celula.empty:
                # Si no hay skeleton ni sholl, al menos usar núcleos
                df_por_celula = df_nuc_astros
        
    if df_por_celula.empty:
        return None, None

    # --- ¡CRÍTICO! Solución a Pseudoreplicación ---
    # 1. DataFrame para gráficos (todas las células)
    df_plot = df_por_celula.copy()
    
    # 2. DataFrame para estadística (mediana por preparado)
    cols_to_agg = [
        # Métricas de Sholl
        'critical_radius_um', 'peak_intersections', 'auc',
        # Métricas nucleares
        'nucleus_volume_um3', 'nucleus_sphericity',
        # Métricas topológicas básicas
        'n_branches', 'total_branch_length_um', 'mean_branch_length_um',
        'n_endpoints', 'n_junctions',
        # Métricas de tortuosidad
        'tortuosity_mean', 'tortuosity_max', 'tortuosity_std',
        # Índices de complejidad
        'ramification_index', 'termination_index',
        # Distribución de longitudes
        'branch_length_median_um', 'branch_length_p25_um', 'branch_length_p75_um',
        'branch_length_std_um', 'branch_length_cv',
        # Fragmentación
        'n_connected_components'
    ]
    valid_cols_to_agg = [col for col in cols_to_agg if col in df_plot.columns]
    
    if not valid_cols_to_agg or 'prepared' not in df_plot.columns or 'group' not in df_plot.columns:
        return df_plot, pd.DataFrame() # No hay suficientes datos

    df_stats = df_plot.groupby(['prepared', 'group'])[valid_cols_to_agg].median().reset_index()
    
    return df_plot, df_stats

# --- Cargar Datos ---
df_plot, df_stats = load_all_metrics(RAW_DIR, PROC_DIR)

if df_plot is None or df_stats is None:
    st.error("No se encontraron métricas procesadas (summary.csv, sholl_summary.csv, etc.) en `data/processed`.")
    st.stop()

# --- Estado del Pipeline ---
st.markdown("### 📊 Estado del Procesamiento")
total_preps = df_stats['prepared'].nunique()

col_status1, col_status2, col_status3 = st.columns(3)

n_ctl = (df_stats['group'] == 'CTL').sum()
n_hip = (df_stats['group'] == 'Hipoxia').sum()

col_status1.metric("Total Preparados", total_preps)
col_status2.metric("CTL", n_ctl)
col_status3.metric("Hipoxia", n_hip)

# Verificar completitud de datos por tipo
numeric_cols_check = [col for col in df_stats.columns 
                      if col not in ['prepared', 'group'] and df_stats[col].dtype in ['float64', 'int64']]

# Sholl
sholl_cols = [c for c in numeric_cols_check if 'auc' in c or 'critical_radius' in c or 'peak' in c]
if sholl_cols:
    n_with_sholl = df_stats[sholl_cols[0]].notna().sum()
    pct_sholl = (n_with_sholl / total_preps) * 100
    col_status1.metric("Con datos Sholl", f"{n_with_sholl} ({pct_sholl:.0f}%)")

# Topología
topo_cols = [c for c in numeric_cols_check if 'branch' in c or 'tortuosity' in c]
if topo_cols:
    n_with_topo = df_stats[topo_cols[0]].notna().sum()
    pct_topo = (n_with_topo / total_preps) * 100
    col_status2.metric("Con datos topológicos", f"{n_with_topo} ({pct_topo:.0f}%)")

# Núcleo
nuc_cols = [c for c in numeric_cols_check if 'nucleus' in c]
if nuc_cols:
    n_with_nuc = df_stats[nuc_cols[0]].notna().sum()
    pct_nuc = (n_with_nuc / total_preps) * 100
    col_status3.metric("Con datos nucleares", f"{n_with_nuc} ({pct_nuc:.0f}%)")

# Botón para refrescar datos
if st.button("🔄 Refrescar datos (limpiar caché)", help="Útil después de re-ejecutar el pipeline"):
    st.cache_data.clear()
    st.rerun()

st.markdown("---")

# --- Definir métricas disponibles para comparar ---
METRIC_OPTIONS = {
    # Sholl
    'auc': ("Sholl: AUC (Complejidad Total)", "µm·intersec"),
    'peak_intersections': ("Sholl: Pico Intersecciones", "count"),
    'critical_radius_um': ("Sholl: Radio Crítico", "µm"),
    # Núcleo
    'nucleus_volume_um3': ("Volumen del Núcleo", "µm³"),
    'nucleus_sphericity': ("Circularidad Nuclear (2D)", "0-1"),
    # Topología básica
    'total_branch_length_um': ("Longitud Total Esqueleto", "µm"),
    'n_branches': ("Número de Ramas", "count"),
    'mean_branch_length_um': ("Longitud Media de Rama", "µm"),
    'n_endpoints': ("Número de Terminaciones", "count"),
    'n_junctions': ("Número de Bifurcaciones", "count"),
    # Tortuosidad
    'tortuosity_mean': ("Tortuosidad Media (1=recto)", "ratio"),
    'tortuosity_max': ("Tortuosidad Máxima", "ratio"),
    'tortuosity_std': ("Desv. Std. Tortuosidad", "ratio"),
    # Complejidad
    'ramification_index': ("Índice de Ramificación (ramas/junctions)", "ratio"),
    'termination_index': ("Índice de Terminalización (endpoints/junctions)", "ratio"),
    # Distribución longitudes
    'branch_length_median_um': ("Mediana Longitud Ramas", "µm"),
    'branch_length_cv': ("Coef. Variación Longitudes (heterogeneidad)", "ratio"),
    'branch_length_std_um': ("Desv. Std. Longitud Ramas", "µm"),
    # Fragmentación
    'n_connected_components': ("Componentes Conectadas (1=esperado)", "count"),
}

# Filtrar opciones basadas en las columnas que realmente existen
available_metrics = {k: v for k, v in METRIC_OPTIONS.items() if k in df_plot.columns}

st.markdown("### Comparación Global de Métricas")
selected_metric = st.selectbox(
    "Seleccioná la métrica a comparar:",
    options=available_metrics.keys(),
    format_func=lambda k: f"{available_metrics[k][0]} ({available_metrics[k][1]})" if available_metrics[k][1] else available_metrics[k][0]
)

metric_label = f"{available_metrics[selected_metric][0]} ({available_metrics[selected_metric][1]})" if available_metrics[selected_metric][1] else available_metrics[selected_metric][0]

st.markdown(f"#### {metric_label}")

# --- Información de Disponibilidad de Datos ---
if selected_metric in df_stats.columns:
    n_ctl = df_stats.loc[df_stats['group']=='CTL', selected_metric].dropna().shape[0]
    n_hip = df_stats.loc[df_stats['group']=='Hipoxia', selected_metric].dropna().shape[0]
    total_preps = n_ctl + n_hip
    
    col_info1, col_info2, col_info3 = st.columns(3)
    col_info1.metric("📊 Preparados con datos", total_preps)
    col_info2.metric("CTL", n_ctl)
    col_info3.metric("Hipoxia", n_hip)
    
    if total_preps < 6:
        st.warning(f"⚠️ Pocos datos disponibles para esta métrica. Se recomienda re-ejecutar el paso 04 (Sholl) o 05 (topología) en más preparados.")

# --- 1. Gráfico de Distribución (usando df_plot) ---
st.markdown("**Distribución por Célula (para visualización)**")
chart = boxplot_with_ticks(df_plot.dropna(subset=[selected_metric]), selected_metric, 'group', title_x=metric_label)
st.altair_chart(chart, use_container_width=True)

# --- 2. Test Estadístico (usando df_stats) ---
st.markdown("**Test Estadístico (sobre medianas por preparado)**")
if selected_metric not in df_stats.columns:
    st.warning(f"La métrica '{selected_metric}' no se pudo agregar por preparado.")
else:
    try:
        a = df_stats.loc[df_stats['group']=='CTL', selected_metric].dropna().to_numpy()
        b = df_stats.loc[df_stats['group']=='Hipoxia', selected_metric].dropna().to_numpy()
        
        if a.size > 0 and b.size > 0:
            test_results = run_comparison_test(a, b)
            
            # Verificar si hay datos suficientes
            if test_results.get('message'):
                # Datos insuficientes o problema
                st.warning(test_results['message'])
                col_test1, col_test2, col_test3, col_test4 = st.columns(4)
                col_test1.metric("Test Usado", "—")
                col_test2.metric("Estadístico", "nan")
                col_test3.metric("P-valor", "p = nan ns")
                col_test4.metric("Significancia", "⚪")
            elif pd.isna(test_results['p_value']):
                # Error en cálculo
                st.error("No se pudo calcular el test estadístico.")
                col_test1, col_test2, col_test3, col_test4 = st.columns(4)
                col_test1.metric("Test Usado", "—")
                col_test2.metric("Estadístico", "nan")
                col_test3.metric("P-valor", "p = nan ns")
                col_test4.metric("Significancia", "⚪")
            else:
                # Test exitoso
                reporte_stats = format_test_results(test_results)
                
                # Mostrar estadística en formato profesional
                col_test1, col_test2, col_test3, col_test4 = st.columns(4)
                col_test1.metric("Test Usado", test_results['test_name'])
                col_test2.metric("Estadístico", f"{test_results['stat_label']}={test_results['stat']:.3f}")
                
                # Color para p-valor
                p_val = test_results['p_value']
                if p_val < 0.001:
                    p_display = "p < 0.001"
                    sig_level = "***"
                    p_color = "🔴"
                elif p_val < 0.01:
                    p_display = f"p = {p_val:.3f}"
                    sig_level = "**"
                    p_color = "🟠"
                elif p_val < 0.05:
                    p_display = f"p = {p_val:.3f}"
                    sig_level = "*"
                    p_color = "🟡"
                else:
                    p_display = f"p = {p_val:.3f}"
                    sig_level = "ns"
                    p_color = "⚪"
                
                col_test3.metric("P-valor", f"{p_display} {sig_level}")
                col_test4.metric("Significancia", p_color)
                
                st.markdown(reporte_stats)
        else:
            st.caption("No hay suficientes datos por preparado para ambos grupos para el test.")
    except Exception as e:
        st.caption(f"Error al calcular estadísticas: {e}")
        st.exception(e)

st.markdown("---")
st.markdown("### 📊 Datos Agregados y Exportación")

st.markdown("""
Los datos se presentan en dos formatos:
- **Por Célula**: Todas las células individuales (útil para gráficos de distribución)
- **Por Preparado**: Mediana de cada preparado (correcto para análisis estadístico, evita pseudoreplicación)
""")

col1, col2 = st.columns(2)
with col1:
    st.markdown("**Datos por Célula (para gráficos)**")
    st.caption(f"Total: {df_plot.shape[0]} células de {df_plot['prepared'].nunique()} preparados")
    st.dataframe(
        df_plot.round(3),
        use_container_width=True,
        height=400,
        column_config={
            "label": st.column_config.NumberColumn("ID", format="%d"),
            "group": st.column_config.TextColumn("Grupo"),
            "prepared": st.column_config.TextColumn("Preparado")
        }
    )
    _download_df_button(df_plot, "metricas_por_celula", "⬇️ Descargar CSV (por célula)")
    
with col2:
    st.markdown("**Datos por Preparado (para estadística)**")
    st.caption(f"Total: {df_stats.shape[0]} preparados (mediana de cada preparado)")
    st.dataframe(
        df_stats.round(3),
        use_container_width=True,
        height=400,
        column_config={
            "group": st.column_config.TextColumn("Grupo"),
            "prepared": st.column_config.TextColumn("Preparado")
        }
    )
    _download_df_button(df_stats, "metricas_por_preparado", "⬇️ Descargar CSV (por preparado)")

# ============================================================================
# ANÁLISIS DE COMPONENTES PRINCIPALES (PCA)
# ============================================================================
st.markdown("---")
st.markdown("### 🔬 Análisis de Componentes Principales (PCA)")
st.markdown("""
El PCA reduce la dimensionalidad de las métricas morfológicas, permitiendo:
- **Visualizar separación** entre grupos CTL vs Hipoxia en espacio reducido
- **Identificar métricas clave** que más contribuyen a la variabilidad
- **Detectar patrones** y subgrupos morfológicos
""")

# Importar sklearn para PCA
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    # Seleccionar métricas numéricas
    numeric_cols = [col for col in df_stats.columns 
                    if col not in ['prepared', 'group'] and df_stats[col].dtype in ['float64', 'int64']]
    
    if len(numeric_cols) < 2:
        st.warning("No hay suficientes métricas numéricas para PCA.")
    else:
        # Diagnóstico de datos faltantes
        st.markdown("#### 🔍 Diagnóstico de Datos")
        
        # Contar NaN por métrica
        missing_counts = df_stats[numeric_cols].isna().sum()
        missing_pct = (missing_counts / len(df_stats)) * 100
        
        # Filtrar métricas con menos del 30% de datos faltantes
        valid_metrics = missing_pct[missing_pct < 30].index.tolist()
        
        col_diag1, col_diag2 = st.columns(2)
        with col_diag1:
            st.metric("Total Preparados", len(df_stats))
            st.metric("Métricas Totales", len(numeric_cols))
        with col_diag2:
            st.metric("Métricas Válidas (<30% NaN)", len(valid_metrics))
            
        # Mostrar métricas con más datos faltantes
        if (missing_pct > 0).any():
            with st.expander("⚠️ Ver métricas con datos faltantes y recomendaciones"):
                df_missing = pd.DataFrame({
                    'Métrica': missing_counts.index,
                    'N° Faltantes': missing_counts.values,
                    '% Faltante': missing_pct.values
                }).sort_values('% Faltante', ascending=False)
                
                df_missing_display = df_missing[df_missing['% Faltante'] > 0].copy()
                
                # Agregar recomendaciones por tipo de métrica
                def get_recommendation(metric_name):
                    if any(x in metric_name for x in ['auc', 'critical_radius', 'peak_intersections']):
                        return "Paso 05: Re-ejecutar análisis de Sholl"
                    elif 'nucleus' in metric_name:
                        return "Paso 03: Re-ejecutar filtrado de núcleos"
                    else:
                        return "Paso 04: Re-ejecutar esqueletización"
                
                df_missing_display['Recomendación'] = df_missing_display['Métrica'].apply(get_recommendation)
                
                st.dataframe(df_missing_display.round(2), use_container_width=True)
                
                # Resumen de preparados faltantes por tipo
                st.markdown("**📋 Preparados que necesitan re-procesamiento:**")
                
                # Identificar qué preparados faltan para cada tipo
                all_preps = set(df_stats['prepared'].unique())
                
                # Sholl
                sholl_cols = [c for c in numeric_cols if any(x in c for x in ['auc', 'critical_radius', 'peak'])]
                if sholl_cols:
                    preps_with_sholl = set(df_stats[df_stats[sholl_cols[0]].notna()]['prepared'])
                    preps_without_sholl = all_preps - preps_with_sholl
                    if preps_without_sholl:
                        st.markdown(f"- **Sholl (Paso 05)**: {len(preps_without_sholl)} preparados sin datos")
                        # Mostrar lista colapsable usando details/summary HTML
                        preps_list = sorted(list(preps_without_sholl))
                        st.markdown(f"  <details><summary>Ver lista de {len(preps_without_sholl)} preparados</summary><ul>" + 
                                   "".join([f"<li>{p}</li>" for p in preps_list]) + 
                                   "</ul></details>", unsafe_allow_html=True)
                
                # Topología
                topo_cols = [c for c in numeric_cols if any(x in c for x in ['branch', 'tortuosity', 'ramification'])]
                if topo_cols:
                    preps_with_topo = set(df_stats[df_stats[topo_cols[0]].notna()]['prepared'])
                    preps_without_topo = all_preps - preps_with_topo
                    if preps_without_topo:
                        st.markdown(f"- **Topología (Paso 04)**: {len(preps_without_topo)} preparados sin datos")
                        # Mostrar lista colapsable usando details/summary HTML
                        preps_list = sorted(list(preps_without_topo))
                        st.markdown(f"  <details><summary>Ver lista de {len(preps_without_topo)} preparados</summary><ul>" + 
                                   "".join([f"<li>{p}</li>" for p in preps_list]) + 
                                   "</ul></details>", unsafe_allow_html=True)
        
        if len(valid_metrics) < 2:
            st.error("No hay suficientes métricas válidas para PCA (se necesitan al menos 2).")
        else:
            # Preparar datos usando solo métricas válidas
            df_pca_input = df_stats[['prepared', 'group'] + valid_metrics].dropna()
            
            st.info(f"**PCA usando {len(valid_metrics)} métricas en {df_pca_input.shape[0]} preparados**")
            
            if df_pca_input.shape[0] < 3:
                st.warning(f"⚠️ Solo hay {df_pca_input.shape[0]} preparados con datos completos. Se necesitan al menos 3 para PCA.")
                st.markdown("**Preparados con datos completos:**")
                st.write(df_pca_input[['prepared', 'group']].to_dict('records'))
            else:
                X = df_pca_input[valid_metrics].values
                groups = df_pca_input['group'].values
                preparados = df_pca_input['prepared'].values
            
            # Normalizar datos (crítico para PCA)
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            
            # Aplicar PCA
            pca = PCA()
            X_pca = pca.fit_transform(X_scaled)
            
            # --- Varianza Explicada ---
            st.markdown("#### 📊 Varianza Explicada por Componente")
            var_exp = pca.explained_variance_ratio_ * 100
            cum_var = np.cumsum(var_exp)
            
            # Mostrar solo los primeros 10 componentes
            n_components_to_show = min(10, len(var_exp))
            
            df_var = pd.DataFrame({
                'Componente': [f'PC{i+1}' for i in range(n_components_to_show)],
                'Varianza (%)': var_exp[:n_components_to_show],
                'Acumulada (%)': cum_var[:n_components_to_show]
            })
            
            col_var1, col_var2 = st.columns([1, 2])
            with col_var1:
                st.dataframe(df_var.round(2), use_container_width=True, height=300)
            
            with col_var2:
                # Gráfico de varianza explicada (barras + línea acumulada)
                base = alt.Chart(df_var).encode(
                    x=alt.X('Componente:N', sort=None, title='Componente Principal', axis=alt.Axis(labelAngle=0))
                )
                
                chart_var = base.mark_bar().encode(
                    y=alt.Y('Varianza (%):Q', title='Varianza Explicada (%)', axis=alt.Axis(titleColor='#1f77b4')),
                    color=alt.value('#1f77b4')
                ).properties(height=250)
                
                line_var = base.mark_line(point=True, color='red', strokeWidth=2).encode(
                    y=alt.Y('Acumulada (%):Q', title='Varianza Acumulada (%)', axis=alt.Axis(titleColor='red'))
                )
                
                chart_combined = alt.layer(chart_var, line_var).resolve_scale(y='independent')
                st.altair_chart(chart_combined, use_container_width=True)
            
            st.info(f"**PC1 + PC2 explican {cum_var[1]:.1f}% de la varianza total**")
            
            # --- Biplot (PC1 vs PC2) ---
            st.markdown("#### 🎯 Biplot: PC1 vs PC2")
            
            df_biplot = pd.DataFrame({
                'PC1': X_pca[:, 0],
                'PC2': X_pca[:, 1],
                'Grupo': groups,
                'Preparado': preparados
            })
            
            # Gráfico de puntos coloreado por grupo
            scatter = alt.Chart(df_biplot).mark_circle(size=100, opacity=0.8).encode(
                x=alt.X('PC1:Q', title=f'PC1 ({var_exp[0]:.1f}%)'),
                y=alt.Y('PC2:Q', title=f'PC2 ({var_exp[1]:.1f}%)'),
                color=alt.Color('Grupo:N', scale=alt.Scale(domain=['CTL', 'Hipoxia'], 
                                                            range=[PALETTE['CTL'], PALETTE['Hipoxia']])),
                tooltip=['Preparado:N', 'Grupo:N', 'PC1:Q', 'PC2:Q']
            ).properties(width=600, height=400).interactive()
            
            st.altair_chart(scatter, use_container_width=True)
            
            # --- Loadings (Contribución de métricas) ---
            st.markdown("#### 🔑 Loadings: Contribución de Métricas a PC1 y PC2")
            
            loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
            df_loadings = pd.DataFrame(
                loadings[:, :2],  # Solo PC1 y PC2
                columns=['PC1', 'PC2'],
                index=valid_metrics  # Usar solo las métricas válidas usadas en PCA
            )
            df_loadings['Métrica'] = df_loadings.index
            df_loadings['Magnitud_PC1'] = np.abs(df_loadings['PC1'])
            df_loadings['Magnitud_PC2'] = np.abs(df_loadings['PC2'])
            
            # Top 10 métricas más importantes para PC1
            st.markdown("**Top 10 métricas con mayor contribución a PC1:**")
            top_pc1 = df_loadings.nlargest(10, 'Magnitud_PC1')[['Métrica', 'PC1']].reset_index(drop=True)
            
            chart_pc1 = alt.Chart(top_pc1).mark_bar().encode(
                x=alt.X('PC1:Q', title='Loading en PC1'),
                y=alt.Y('Métrica:N', sort='-x', title=''),
                color=alt.condition(
                    alt.datum.PC1 > 0,
                    alt.value('#1f77b4'),
                    alt.value('#d62728')
                )
            ).properties(height=300)
            st.altair_chart(chart_pc1, use_container_width=True)
            
            # Top 10 métricas más importantes para PC2
            st.markdown("**Top 10 métricas con mayor contribución a PC2:**")
            top_pc2 = df_loadings.nlargest(10, 'Magnitud_PC2')[['Métrica', 'PC2']].reset_index(drop=True)
            
            chart_pc2 = alt.Chart(top_pc2).mark_bar().encode(
                x=alt.X('PC2:Q', title='Loading en PC2'),
                y=alt.Y('Métrica:N', sort='-x', title=''),
                color=alt.condition(
                    alt.datum.PC2 > 0,
                    alt.value('#1f77b4'),
                    alt.value('#d62728')
                )
            ).properties(height=300)
            st.altair_chart(chart_pc2, use_container_width=True)
            
            # Tabla completa de loadings
            with st.expander("📋 Ver todos los loadings"):
                df_loadings_display = df_loadings.sort_values('Magnitud_PC1', ascending=False)[['Métrica', 'PC1', 'PC2']]
                st.dataframe(df_loadings_display.round(4), use_container_width=True, height=400)
            
            # --- Interpretación Biológica y Estadísticas de Separación ---
            st.markdown("#### 🧬 Interpretación Biológica y Separación de Grupos")
            
            # Calcular centroides por grupo
            ctl_pc1 = df_biplot[df_biplot['Grupo'] == 'CTL']['PC1'].values
            ctl_pc2 = df_biplot[df_biplot['Grupo'] == 'CTL']['PC2'].values
            hip_pc1 = df_biplot[df_biplot['Grupo'] == 'Hipoxia']['PC1'].values
            hip_pc2 = df_biplot[df_biplot['Grupo'] == 'Hipoxia']['PC2'].values
            
            ctl_mean_pc1 = ctl_pc1.mean()
            hip_mean_pc1 = hip_pc1.mean()
            ctl_mean_pc2 = ctl_pc2.mean()
            hip_mean_pc2 = hip_pc2.mean()
            
            # Distancia euclidiana entre centroides
            dist_centroid = np.sqrt((ctl_mean_pc1 - hip_mean_pc1)**2 + (ctl_mean_pc2 - hip_mean_pc2)**2)
            
            # Test estadístico de separación (PERMANOVA simplificado usando Mann-Whitney en PC1 y PC2)
            _, p_pc1 = stats.mannwhitneyu(ctl_pc1, hip_pc1, alternative='two-sided')
            _, p_pc2 = stats.mannwhitneyu(ctl_pc2, hip_pc2, alternative='two-sided')
            
            # Métricas de clusterización
            st.markdown("##### 📊 Métricas de Separación entre Grupos")
            
            col_sep1, col_sep2, col_sep3, col_sep4 = st.columns(4)
            with col_sep1:
                st.metric("N° CTL", len(ctl_pc1))
            with col_sep2:
                st.metric("N° Hipoxia", len(hip_pc1))
            with col_sep3:
                st.metric("Distancia Centroides", f"{dist_centroid:.2f}", 
                         help="Distancia euclidiana en espacio PC1-PC2")
            with col_sep4:
                # Significancia general (mínimo p-valor)
                p_min = min(p_pc1, p_pc2)
                sig_label = "***" if p_min < 0.001 else "**" if p_min < 0.01 else "*" if p_min < 0.05 else "ns"
                st.metric("Significancia", sig_label, 
                         help=f"Basado en Mann-Whitney (p={p_min:.3g})")
            
            # Detalles de separación por componente
            st.markdown("##### 🎯 Separación por Componente Principal")
            col_pc1_sep, col_pc2_sep = st.columns(2)
            
            with col_pc1_sep:
                st.markdown(f"**PC1 ({var_exp[0]:.1f}% varianza)**")
                st.write(f"- CTL: {ctl_mean_pc1:.2f} ± {ctl_pc1.std():.2f}")
                st.write(f"- Hip: {hip_mean_pc1:.2f} ± {hip_pc1.std():.2f}")
                st.write(f"- **p-valor:** {p_pc1:.3g} {'✓' if p_pc1 < 0.05 else '✗'}")
                
            with col_pc2_sep:
                st.markdown(f"**PC2 ({var_exp[1]:.1f}% varianza)**")
                st.write(f"- CTL: {ctl_mean_pc2:.2f} ± {ctl_pc2.std():.2f}")
                st.write(f"- Hip: {hip_mean_pc2:.2f} ± {hip_pc2.std():.2f}")
                st.write(f"- **p-valor:** {p_pc2:.3g} {'✓' if p_pc2 < 0.05 else '✗'}")
            
            # Interpretación de la separación
            if dist_centroid > 3 and p_min < 0.05:
                interpretation = "🟢 **Separación significativa:** Los grupos muestran diferencias morfológicas claras."
            elif dist_centroid > 3:
                interpretation = "🟡 **Separación notable pero no significativa:** Diferencias visibles pero con variabilidad alta."
            elif p_min < 0.05:
                interpretation = "🟡 **Diferencias estadísticamente significativas:** Grupos separados pero con solapamiento."
            else:
                interpretation = "🔴 **Sin separación clara:** Los grupos se solapan morfológicamente."
            
            st.info(interpretation)
            
            # Interpretación automática basada en loadings
            st.markdown("**Interpretación de componentes:**")
            
            # PC1 interpretación
            top3_pc1_pos = df_loadings[df_loadings['PC1'] > 0].nlargest(3, 'Magnitud_PC1')['Métrica'].tolist()
            top3_pc1_neg = df_loadings[df_loadings['PC1'] < 0].nlargest(3, 'Magnitud_PC1')['Métrica'].tolist()
            
            st.markdown(f"""
            **PC1 ({var_exp[0]:.1f}% varianza):**
            - ➕ Dirección positiva asociada con: *{', '.join(top3_pc1_pos[:3]) if top3_pc1_pos else 'N/A'}*
            - ➖ Dirección negativa asociada con: *{', '.join(top3_pc1_neg[:3]) if top3_pc1_neg else 'N/A'}*
            """)
            
            # PC2 interpretación
            top3_pc2_pos = df_loadings[df_loadings['PC2'] > 0].nlargest(3, 'Magnitud_PC2')['Métrica'].tolist()
            top3_pc2_neg = df_loadings[df_loadings['PC2'] < 0].nlargest(3, 'Magnitud_PC2')['Métrica'].tolist()
            
            st.markdown(f"""
            **PC2 ({var_exp[1]:.1f}% varianza):**
            - ➕ Dirección positiva asociada con: *{', '.join(top3_pc2_pos[:3]) if top3_pc2_pos else 'N/A'}*
            - ➖ Dirección negativa asociada con: *{', '.join(top3_pc2_neg[:3]) if top3_pc2_neg else 'N/A'}*
            """)
            
except ImportError:
    st.error("sklearn no está instalado. Para usar PCA, instalá: `pip install scikit-learn`")