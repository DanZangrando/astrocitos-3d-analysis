import json
from pathlib import Path
import altair as alt
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
import pingouin as pg
import streamlit as st
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
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
    rows_sholl_raw = [] # New list for raw profiles
    
    for p in files:
        od = proc_dir / p.stem
        group = detect_group_from_path(p)
        
        # 1. Cargar métricas topológicas de esqueleto (skeletons/summary.csv)
        ps = od / "skeletons" / "summary.csv"
        if ps.exists():
            try:
                df = pd.read_csv(ps); df['prepared'] = f"{group}-{p.stem}"; df['group'] = group
                rows_skel.append(df)
            except Exception: pass
        
        # 2. Cargar métricas de Sholl 2D (sholl_summary.csv)
        pss = od / "sholl_summary.csv"
        if pss.exists():
            try:
                df = pd.read_csv(pss); df['prepared'] = f"{group}-{p.stem}"; df['group'] = group
                rows_sholl.append(df)
            except Exception: pass

        # 2b. Cargar Perfiles Raw de Sholl (sholl_2d_native.csv)
        pr = od / "sholl_2d_native.csv"
        if pr.exists():
            try:
                # Cargar solo columnas necesarias para ahorrar memoria
                df = pd.read_csv(pr, usecols=['label', 'radius_um', 'intersections']) 
                df['prepared'] = f"{group}-{p.stem}"
                df['group'] = group
                rows_sholl_raw.append(df)
            except Exception: pass
            
        # 3. Cargar métricas de núcleo (03_nucleus_metrics.csv)
        pn = od / "03_nucleus_metrics.csv"
        if pn.exists():
            try:
                df = pd.read_csv(pn); df['prepared'] = f"{group}-{p.stem}"; df['group'] = group
                rows_nuc.append(df)
            except Exception: pass

    if not rows_skel and not rows_sholl and not rows_nuc and not rows_sholl_raw:
        return None, None, None

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
        return df_plot, pd.DataFrame(), pd.DataFrame() # No hay suficientes datos

    df_stats = df_plot.groupby(['prepared', 'group'])[valid_cols_to_agg].median().reset_index()
    
    # --- Fusión de Curvas Sholl (Raw) ---
    # Necesitamos cargar los perfiles raw para la gráfica de curvas
    df_sholl_curves = pd.DataFrame()
    if rows_sholl_raw:
        df_sholl_curves = pd.concat(rows_sholl_raw, ignore_index=True)
    
    return df_plot, df_stats, df_sholl_curves

# --- Cargar Datos ---
df_plot, df_stats, df_sholl_curves = load_all_metrics(RAW_DIR, PROC_DIR)

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
    'critical_radius_um': ("Sholl: Radio Crítico", "µm"),
    'peak_intersections': ("Sholl: Pico Máximo de Intersecciones", "count"),
    
    # Topología Básica
    'total_branch_length_um': ("Longitud Total Esqueleto", "µm"),
    'n_endpoints': ("Número de Terminaciones", "count"),
    'n_branches': ("Número Total de Ramas", "count"),
    'n_junctions': ("Número de Nodos/Bifurcaciones", "count"),
    
    # Análisis de Longitud de Ramas
    'mean_branch_length_um': ("Longitud Media de Rama", "µm"),
    'branch_length_std_um': ("Desviación Estándar de Longitud de Rama", "µm"),
    'branch_length_p75_um': ("Percentil 75 Longitud de Rama", "µm"),
    
    # Índices y Tortuosidad
    'ramification_index': ("Índice de Ramificación (ramas/junctions)", "ratio"),
    'tortuosity_mean': ("Tortuosidad Media", "ratio"),
    'tortuosity_std': ("Desviación Estándar de Tortuosidad", "ratio"),
    
    # Núcleo
    'nucleus_volume_um3': ("Volumen Nuclear", "µm³")
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
st.markdown("### 📈 Análisis de Curvas de Sholl")

if df_sholl_curves.empty:
    st.info("No hay datos de curvas Sholl disponibles para comparar.")
else:
    # 1. Agregación Robusta: Promedio por preparado primero (evita pseudoreplicación)
    # Unidad experimental = Preparado
    df_sholl_prep = df_sholl_curves.groupby(['group', 'prepared', 'radius_um'])['intersections'].mean().reset_index()
    
    # 2. Visualización (Media ± SEM de los preparados)
    st.markdown("**Perfiles de Sholl (Media ± Error Estándar de los Preparados)**")
    
    base = alt.Chart(df_sholl_prep).encode(
        x=alt.X('radius_um:Q', title='Radio (µm)'),
        color=alt.Color('group:N', scale=alt.Scale(domain=['CTL', 'Hipoxia'], range=['#377eb8', '#e41a1c']))
    )
    
    line = base.mark_line(point=False).encode(
        y=alt.Y('mean(intersections):Q', title='Intersecciones (Media)')
    )
    
    band = base.mark_errorband(extent='stderr').encode(
        y=alt.Y('intersections:Q', title='Intersecciones')
    )
    
    chart_sholl = (band + line).properties(height=400)
    st.altair_chart(chart_sholl, use_container_width=True)
    
    # 3. Estadística (Área Bajo la Curva - AUC)
    st.markdown("**Análisis Estadístico (Área Bajo la Curva - AUC)**")
    
    # Calcular la integral (AUC) para cada preparado
    def compute_auc(group_df):
        x = group_df['radius_um'].values
        y = group_df['intersections'].values
        # Ordenar por radio para integrar correctamente
        sort_idx = np.argsort(x)
        return np.trapz(y[sort_idx], x[sort_idx])

    # Aplicar cálculo de AUC
    auc_records = []
    for (group, prep), df_group in df_sholl_prep.groupby(['group', 'prepared']):
        auc = compute_auc(df_group)
        auc_records.append({'group': group, 'prepared': prep, 'sholl_auc': auc})
        
    df_auc = pd.DataFrame(auc_records)

    col_auc_plot, col_auc_stats = st.columns([1, 1])
    
    with col_auc_plot:
        # Gráfico de Boxplot del AUC
        st.markdown("**Distribución del AUC por Grupo**")
        chart_auc = boxplot_with_ticks(df_auc, 'sholl_auc', 'group', title_x="Área Bajo la Curva (AUC)")
        st.altair_chart(chart_auc, use_container_width=True)

    with col_auc_stats:
        # Test estadístico del AUC
        st.markdown("**Test sobre el AUC Total**")
        try:
            a_auc = df_auc.loc[df_auc['group'] == 'CTL', 'sholl_auc'].dropna().to_numpy()
            b_auc = df_auc.loc[df_auc['group'] == 'Hipoxia', 'sholl_auc'].dropna().to_numpy()
            
            if a_auc.size > 0 and b_auc.size > 0:
                test_auc = run_comparison_test(a_auc, b_auc)
                
                if test_auc.get('message') and "No hay suficientes datos" in test_auc['message']:
                     st.warning(test_auc['message'])
                elif pd.isna(test_auc['p_value']):
                     st.error("No se pudo calcular el test sobre AUC.")
                else:
                    st.metric("Test Usado", test_auc['test_name'])
                    st.metric("Estadístico", f"{test_auc['stat_label']}={test_auc['stat']:.3f}")
                    
                    p_val_auc = test_auc['p_value']
                    if p_val_auc < 0.001:
                        p_display = "p < 0.001 *** 🔴"
                    elif p_val_auc < 0.01:
                        p_display = f"p = {p_val_auc:.3f} ** 🟠"
                    elif p_val_auc < 0.05:
                        p_display = f"p = {p_val_auc:.3f} * 🟡"
                    else:
                        p_display = f"p = {p_val_auc:.3f} ns ⚪"
                        
                    st.metric("Significancia (p-valor)", p_display)
                    
                    st.info(f"*(N(CTL)={test_auc['n_a']}, N(Hip)={test_auc['n_b']} preparados)*")
            else:
                st.caption("No hay suficientes datos por preparado para comparar AUC.")
        except Exception as e:
            st.error(f"Error procesando estadística AUC: {e}")
        st.caption("Verificá que los radios sean consistentes entre preparados.")

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

st.markdown("---")
st.markdown("### 🧠 Análisis de Componentes Principales (PCA)")

st.markdown("""
El PCA permite visualizar si el conjunto global de métricas morfológicas y topológicas es suficiente para 
segregar los astrocitos de **Control** versus los de **Hipoxia** en el plano de máxima variación.
""")

pca_data_source = st.radio(
    "Seleccioná la fuente de datos para el PCA:",
    options=["Por Preparado (Mediana)", "Por Célula (Individual)"],
    index=0,
    horizontal=True,
    help="Por preparado evita la pseudoreplicación y muestra la tendencia de cada individuo. Por célula muestra la heterogeneidad de toda la población."
)

df_pca_source = df_stats if pca_data_source == "Por Preparado (Mediana)" else df_plot

# Filtrar solo columnas numéricas que sean métricas (excluyendo IDs como label)
pca_features = [col for col in df_pca_source.select_dtypes(include=['float64', 'int64']).columns 
               if col not in ['label', 'prepared', 'group'] and not col.startswith('Unnamed')]

if len(pca_features) < 2:
    st.warning("No hay suficientes variables numéricas para realizar un PCA.")
else:
    # Eliminar filas con NaNs en las features seleccionadas
    df_pca_clean = df_pca_source.dropna(subset=pca_features).copy()
    
    if df_pca_clean.empty or len(df_pca_clean) < 3:
        st.warning("No hay suficientes datos válidos (sin NaNs) para realizar el PCA.")
    else:
        # Pre-procesamiento: Escalar datos (Z-score)
        X = df_pca_clean[pca_features].values
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # PCA
        pca = PCA(n_components=2)
        components = pca.fit_transform(X_scaled)
        
        # Varianza explicada
        var_ratio = pca.explained_variance_ratio_ * 100
        
        df_pca_clean['PC1'] = components[:, 0]
        df_pca_clean['PC2'] = components[:, 1]
        
        # Visualización Scatter PCA
        st.markdown("**Proyección PCA (PC1 vs PC2)**")
        
        # Tooltip interactivo
        tooltip_cols = ['group', 'prepared']
        if 'label' in df_pca_clean.columns and pca_data_source == "Por Célula (Individual)":
            tooltip_cols.append('label')
            
        scatter = alt.Chart(df_pca_clean).mark_circle(size=60, opacity=0.8).encode(
            x=alt.X('PC1:Q', title=f'PC1 ({var_ratio[0]:.1f}%)'),
            y=alt.Y('PC2:Q', title=f'PC2 ({var_ratio[1]:.1f}%)'),
            color=alt.Color('group:N', scale=alt.Scale(domain=['CTL', 'Hipoxia'], range=['#377eb8', '#e41a1c'])),
            tooltip=tooltip_cols
        ).interactive().properties(height=400)
        
        st.altair_chart(scatter, use_container_width=True)
        
        # Visualización Factor Loadings
        with st.expander("🔍 Ver contribución de las variables originales (Loadings)"):
            st.markdown("Muestra qué variables empujan a las células hacia valores positivos o negativos en cada Componente Principal.")
            
            # Cargar loadings
            loadings = pca.components_.T
            df_loadings = pd.DataFrame(loadings, columns=['PC1', 'PC2'], index=pca_features)
            
            col_l1, col_l2 = st.columns(2)
            
            df_loadings_pc1 = df_loadings.reset_index().rename(columns={'index':'Métrica'}).sort_values('PC1', ascending=False)
            chart_pc1 = alt.Chart(df_loadings_pc1).mark_bar().encode(
                y=alt.Y('Métrica:N', sort='-x', title=None),
                x=alt.X('PC1:Q', title='Peso en PC1'),
                color=alt.condition(alt.datum.PC1 > 0, alt.value("steelblue"), alt.value("orange"))
            ).properties(title="Contribución a PC1")
            
            df_loadings_pc2 = df_loadings.reset_index().rename(columns={'index':'Métrica'}).sort_values('PC2', ascending=False)
            chart_pc2 = alt.Chart(df_loadings_pc2).mark_bar().encode(
                y=alt.Y('Métrica:N', sort='-x', title=None),
                x=alt.X('PC2:Q', title='Peso en PC2'),
                color=alt.condition(alt.datum.PC2 > 0, alt.value("steelblue"), alt.value("orange"))
            ).properties(title="Contribución a PC2")
            
            col_l1.altair_chart(chart_pc1, use_container_width=True)
            col_l2.altair_chart(chart_pc2, use_container_width=True)

