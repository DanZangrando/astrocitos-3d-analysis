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
from ui.plots import GROUP_SCALE, PALETTE, boxplot_with_ticks, boxplot_vertical, apply_theme

# Raíz del repo
ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"

st.title("Group Comparison (CTL vs Hypoxia)")
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
        results['message'] = f"No hay suficientes datos (CTL n={a.size}, Hypoxia n={b.size}) para un test robusto."
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
    """Detecta el grupo (CTL o Hypoxia) desde la ruta del archivo."""
    return 'Hypoxia' if '/hip/' in str(p).replace('\\', '/').lower() else 'CTL'

# --- Carga de Datos (REFACTORIZADA) ---
@st.cache_data(ttl=60) # Cachear por 1 min (reducido para reflejar cambios más rápido)
def load_all_metrics(raw_dir: Path, proc_dir: Path):
    """
    Carga y fusiona todas las métricas (núcleo, esqueleto, sholl) de todos 
    los individuos procesados.
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
                df = pd.read_csv(ps); df['prepared'] = f"{group}-{p.stem}"; df['group'] = group; df['individual'] = p.parent.name
                rows_skel.append(df)
            except Exception: pass
        
        # 2. Cargar métricas de Sholl 2D (sholl_summary.csv)
        pss = od / "sholl_summary.csv"
        if pss.exists():
            try:
                df = pd.read_csv(pss); df['prepared'] = f"{group}-{p.stem}"; df['group'] = group; df['individual'] = p.parent.name
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
                df['individual'] = p.parent.name
                rows_sholl_raw.append(df)
            except Exception: pass
            
        # 3. Cargar métricas de núcleo (03_nucleus_metrics.csv)
        pn = od / "03_nucleus_metrics.csv"
        if pn.exists():
            try:
                df = pd.read_csv(pn); df['prepared'] = f"{group}-{p.stem}"; df['group'] = group; df['individual'] = p.parent.name
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
                                     on=['label', 'prepared', 'group', 'individual'], how='left')
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
                                        on=['label', 'prepared', 'group', 'individual'], how='left')
            elif df_por_celula.empty:
                # Si no hay skeleton ni sholl, al menos usar núcleos
                df_por_celula = df_nuc_astros
        
    if df_por_celula.empty:
        return None, None

    # --- ¡CRÍTICO! Solución a Pseudoreplicación ---
    # 1. DataFrame para gráficos (todas las células)
    df_plot = df_por_celula.copy()
    
    # 2. DataFrame para estadística (media por individuo para evitar pseudoreplicación)
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
    
    if not valid_cols_to_agg or 'individual' not in df_plot.columns or 'group' not in df_plot.columns:
        return df_plot, pd.DataFrame(), pd.DataFrame() # No hay suficientes datos

    df_stats = df_plot.groupby(['individual', 'group'])[valid_cols_to_agg].mean().reset_index()
    
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
total_individuals = df_stats['individual'].nunique()

col_status1, col_status2, col_status3 = st.columns(3)

n_ctl = (df_stats['group'] == 'CTL').sum()
n_hip = (df_stats['group'] == 'Hypoxia').sum()

col_status1.metric("Total Individuos", total_individuals)
col_status2.metric("CTL", n_ctl)
col_status3.metric("Hypoxia", n_hip)

# Verificar completitud de datos por tipo
numeric_cols_check = [col for col in df_stats.columns 
                      if col not in ['individual', 'group'] and df_stats[col].dtype in ['float64', 'int64']]

# Sholl
sholl_cols = [c for c in numeric_cols_check if 'auc' in c or 'critical_radius' in c or 'peak' in c]
if sholl_cols:
    n_with_sholl = df_stats[sholl_cols[0]].notna().sum()
    pct_sholl = (n_with_sholl / total_individuals) * 100
    col_status1.metric("Con datos Sholl", f"{n_with_sholl} ({pct_sholl:.0f}%)")

# Topología
topo_cols = [c for c in numeric_cols_check if 'branch' in c or 'tortuosity' in c]
if topo_cols:
    n_with_topo = df_stats[topo_cols[0]].notna().sum()
    pct_topo = (n_with_topo / total_individuals) * 100
    col_status2.metric("Con datos topológicos", f"{n_with_topo} ({pct_topo:.0f}%)")

# Núcleo
nuc_cols = [c for c in numeric_cols_check if 'nucleus' in c]
if nuc_cols:
    n_with_nuc = df_stats[nuc_cols[0]].notna().sum()
    pct_nuc = (n_with_nuc / total_individuals) * 100
    col_status3.metric("Con datos nucleares", f"{n_with_nuc} ({pct_nuc:.0f}%)")

# Botón para refrescar datos
if st.button("🔄 Refrescar datos (limpiar caché)", help="Útil después de re-ejecutar el pipeline"):
    st.cache_data.clear()
    st.rerun()

st.markdown("---")

# --- Definir métricas disponibles para comparar ---
METRIC_OPTIONS = {
    # Sholl
    'critical_radius_um': ("Sholl: Critical Radius", "µm"),
    'peak_intersections': ("Sholl: Peak Max Intersections", "count"),
    
    # Topología Básica
    'total_branch_length_um': ("Total Skeleton Length", "µm"),
    'n_endpoints': ("Number of Endpoints", "count"),
    'n_branches': ("Total Number of Branches", "count"),
    'n_junctions': ("Number of Junctions", "count"),
    
    # Análisis de Longitud de Ramas
    'mean_branch_length_um': ("Mean Branch Length", "µm"),
    'branch_length_std_um': ("Branch Length Std", "µm"),
    'branch_length_p75_um': ("Branch Length 75th Percentile", "µm"),
    
    # Índices y Tortuosidad
    'ramification_index': ("Ramification Index (branches/junctions)", "ratio"),
    'tortuosity_mean': ("Mean Tortuosity", "ratio"),
    'tortuosity_std': ("Tortuosity Std", "ratio"),
    
    # Núcleo
    'nucleus_volume_um3': ("Nuclear Volume", "µm³")
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
    n_hip = df_stats.loc[df_stats['group']=='Hypoxia', selected_metric].dropna().shape[0]
    total_individuals = n_ctl + n_hip
    
    col_info1, col_info2, col_info3 = st.columns(3)
    col_info1.metric("📊 Individuos con datos", total_individuals)
    col_info2.metric("CTL", n_ctl)
    col_info3.metric("Hypoxia", n_hip)
    
    if total_individuals < 6:
        st.warning(f"⚠️ Pocos datos disponibles para esta métrica. Se recomienda re-ejecutar el paso 04 (Sholl) o 05 (topología) en más individuos.")

# --- 1. Gráfico de Distribución (usando df_plot) ---
col_plot1, col_plot2 = st.columns(2)
with col_plot1:
    distribution_plot_type = st.radio(
        "Tipo de gráfico de distribución:",
        options=["Boxplot Vertical (Recomendado)", "Distribución Horizontal (Barcode)"],
        index=0,
        horizontal=True
    )
with col_plot2:
    data_source_level = st.radio(
        "Nivel de análisis para el gráfico:",
        options=["Por Individuo (Animal)", "Por Célula (Individual)"],
        index=0,
        horizontal=True,
        help="Por individuo muestra un punto por animal (media). Por célula muestra todas las células medidas."
    )

df_to_plot = df_stats if data_source_level == "Por Individuo (Animal)" else df_plot
n_points = df_to_plot[selected_metric].dropna().shape[0]
st.markdown(f"**Distribución {data_source_level} (N={n_points} puntos)**")

if distribution_plot_type == "Boxplot Vertical (Recomendado)":
    chart = boxplot_vertical(df_to_plot.dropna(subset=[selected_metric]), selected_metric, 'group', title_y=metric_label)
else:
    chart = boxplot_with_ticks(df_to_plot.dropna(subset=[selected_metric]), selected_metric, 'group', title_x=metric_label)

st.altair_chart(chart, use_container_width=True)
_download_df_button(df_to_plot.dropna(subset=[selected_metric]), filename_base=f"datos_grafico_{selected_metric}", label="⬇️ Descargar datos de este gráfico (CSV)")

# --- 2. Test Estadístico (usando df_stats) ---
st.markdown("**Test Estadístico (sobre promedios por individuo)**")
if selected_metric not in df_stats.columns:
    st.warning(f"La métrica '{selected_metric}' no se pudo agregar por individuo.")
else:
    try:
        a = df_stats.loc[df_stats['group']=='CTL', selected_metric].dropna().to_numpy()
        b = df_stats.loc[df_stats['group']=='Hypoxia', selected_metric].dropna().to_numpy()
        
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
            st.caption("No hay suficientes datos por individuo para ambos grupos para el test.")
    except Exception as e:
        st.caption(f"Error al calcular estadísticas: {e}")
        st.exception(e)


st.markdown("---")
st.markdown("### 📈 Análisis de Curvas de Sholl")

if df_sholl_curves.empty:
    st.info("No hay datos de curvas Sholl disponibles para comparar.")
else:
    # 1. Agregación Robusta: Promedio por individuo primero (evita pseudoreplicación)
    # Unidad experimental = Individuo
    df_sholl_prep = df_sholl_curves.groupby(['group', 'individual', 'radius_um'])['intersections'].mean().reset_index()
    
    # 2. Visualización (Media ± SEM de los individuos)
    st.markdown("**Perfiles de Sholl (Media ± Error Estándar de los Individuos)**")
    
    base = alt.Chart(df_sholl_prep).encode(
        x=alt.X('radius_um:Q', title='Radius (µm)'),
        color=alt.Color('group:N', scale=alt.Scale(domain=['CTL', 'Hypoxia'], range=['#377eb8', '#e41a1c']))
    )
    
    line = base.mark_line(point=False).encode(
        y=alt.Y('mean(intersections):Q', title='Intersections (Mean)')
    )
    
    band = base.mark_errorband(extent='stderr').encode(
        y=alt.Y('intersections:Q', title='Intersections')
    )
    
    chart_sholl = (band + line).properties(height=400)
    st.altair_chart(chart_sholl, use_container_width=True)
    _download_df_button(df_sholl_prep, filename_base="datos_grafico_curvas_sholl", label="⬇️ Download Data for this Chart (CSV)")
    
    # 3. Estadística (Área Bajo la Curva - AUC)
    st.markdown("**Análisis Estadístico (Área Bajo la Curva - AUC)**")
    
    # Calcular la integral (AUC) para cada individuo
    def compute_auc(group_df):
        x = group_df['radius_um'].values
        y = group_df['intersections'].values
        # Ordenar por radio para integrar correctamente
        sort_idx = np.argsort(x)
        return np.trapz(y[sort_idx], x[sort_idx])

    # Aplicar cálculo de AUC
    auc_records = []
    for (group, individual), df_group in df_sholl_prep.groupby(['group', 'individual']):
        auc = compute_auc(df_group)
        auc_records.append({'group': group, 'individual': individual, 'sholl_auc': auc})
        
    df_auc = pd.DataFrame(auc_records)

    col_auc_plot, col_auc_stats = st.columns([1, 1])
    
    with col_auc_plot:
        # Gráfico de Boxplot del AUC
        st.markdown("**Distribución del AUC por Grupo**")
        chart_auc = boxplot_vertical(df_auc, 'sholl_auc', 'group', title_y="Area Under the Curve (AUC)")
        st.altair_chart(chart_auc, use_container_width=True)
        _download_df_button(df_auc, filename_base="datos_grafico_auc_sholl", label="⬇️ Download Data for this Chart (CSV)")

    with col_auc_stats:
        # Test estadístico del AUC
        st.markdown("**Test sobre el AUC Total**")
        try:
            a_auc = df_auc.loc[df_auc['group'] == 'CTL', 'sholl_auc'].dropna().to_numpy()
            b_auc = df_auc.loc[df_auc['group'] == 'Hypoxia', 'sholl_auc'].dropna().to_numpy()
            
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
                    
                    st.info(f"*(N(CTL)={test_auc['n_a']}, N(Hip)={test_auc['n_b']} individuos)*")
            else:
                st.caption("No hay suficientes datos por individuo para comparar AUC.")
        except Exception as e:
            st.error(f"Error procesando estadística AUC: {e}")
        st.caption("Verificá que los radios sean consistentes entre individuos.")

st.markdown("---")
st.markdown("### 📊 Datos Agregados y Exportación")

st.markdown("""
Los datos se presentan en dos formatos:
- **Por Célula**: Todas las células individuales (útil para gráficos de distribución)
- **Por Individuo**: Mediana/Media de cada animal (correcto para análisis estadístico, evita pseudoreplicación)
""")

col1, col2 = st.columns(2)
with col1:
    st.markdown("**Datos por Célula (para gráficos)**")
    st.caption(f"Total: {df_plot.shape[0]} células de {df_plot['individual'].nunique()} individuos")
    st.dataframe(
        df_plot.round(3),
        use_container_width=True,
        height=400,
        column_config={
            "label": st.column_config.NumberColumn("ID", format="%d"),
            "group": st.column_config.TextColumn("Grupo"),
            "prepared": st.column_config.TextColumn("Preparado"),
            "individual": st.column_config.TextColumn("Individuo")
        }
    )
    _download_df_button(df_plot, "metricas_por_celula", "⬇️ Descargar CSV (por célula)")
    
with col2:
    st.markdown("**Datos por Individuo (para estadística)**")
    st.caption(f"Total: {df_stats.shape[0]} individuos (promedio de cada animal)")
    st.dataframe(
        df_stats.round(3),
        use_container_width=True,
        height=400,
        column_config={
            "group": st.column_config.TextColumn("Grupo"),
            "individual": st.column_config.TextColumn("Individuo")
        }
    )
    _download_df_button(df_stats, "metricas_por_individuo", "⬇️ Descargar CSV (por individuo)")

st.markdown("---")
st.markdown("### 🧠 Análisis de Componentes Principales (PCA)")

st.markdown("""
segregar los astrocitos de **Control** versus los de **Hypoxia** en el plano de máxima variación.
""")

pca_data_source = st.radio(
    "Seleccioná la fuente de datos para el PCA:",
    options=["Por Individuo (Animal)", "Por Célula (Individual)"],
    index=0,
    horizontal=True,
    help="Por individuo evita la pseudoreplicación y muestra la tendencia de cada animal. Por célula muestra la heterogeneidad de toda la población."
)

df_pca_source = df_stats if pca_data_source == "Por Individuo (Animal)" else df_plot

# --- Selección Curada de Variables para PCA ---
with st.expander("🛠️ Configuración de Variables PCA", expanded=True):
    st.markdown("""
    Seleccioná las métricas que querés incluir en el análisis. 
    **Nota:** Incluir demasiadas variables redundantes (como contar ramas y nodos a la vez) puede sesgar el PCA.
    """)
    
    # Diccionario de descripciones para ayudar al usuario
    FEATURE_INFO = {
        'total_branch_length_um': "Longitud total de todos los procesos (Tamaño)",
        'n_endpoints': "Número de terminaciones (Complejidad de ramificación)",
        'mean_branch_length_um': "Promedio de longitud de ramas individuales (Escala)",
        'tortuosity_mean': "Sinuosidad promedio de las ramas (Geometría)",
        'critical_radius_um': "Radio donde ocurre el pico de intersecciones (Extensión Sholl)",
        'peak_intersections': "Máximo número de intersecciones (Densidad Sholl)",
        'auc': "Área bajo la curva Sholl (Complejidad acumulada)",
        'nucleus_volume_um3': "Volumen del núcleo (Tamaño del soma)",
        'n_branches': "Número total de segmentos (Altamente redundante con terminaciones)",
        'n_junctions': "Número de bifurcaciones (Altamente redundante con terminaciones)",
        'ramification_index': "Ratio ramas/nodos (Poco discriminativo en árboles binarios)"
    }
    
    all_numeric = [col for col in df_pca_source.select_dtypes(include=['float64', 'int64']).columns 
                  if col not in ['label', 'prepared', 'group'] and not col.startswith('Unnamed')]
    
    default_curated = [
        'total_branch_length_um', 'n_endpoints', 'mean_branch_length_um', 
        'tortuosity_mean', 'critical_radius_um', 'peak_intersections', 'nucleus_volume_um3'
    ]
    # Filtrar solo las que existen en el dataframe
    default_curated = [c for c in default_curated if c in all_numeric]
    
    selected_pca_features = st.multiselect(
        "Variables a incluir:",
        options=all_numeric,
        default=default_curated,
        format_func=lambda x: f"{x} ({FEATURE_INFO.get(x, 'Métrica')})"
    )
    
if len(selected_pca_features) < 2:
    st.warning("Seleccioná al menos 2 variables para realizar el PCA.")
else:
    # Eliminar filas con NaNs en las features seleccionadas
    df_pca_clean = df_pca_source.dropna(subset=selected_pca_features).copy()
    
    if df_pca_clean.empty or len(df_pca_clean) < 3:
        st.warning("No hay suficientes datos válidos (sin NaNs) para las variables seleccionadas.")
    else:
        # Pre-procesamiento: Escalar datos (Z-score)
        X = df_pca_clean[selected_pca_features].values
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X)
        
        # PCA
        pca = PCA(n_components=2)
        components = pca.fit_transform(X_scaled)
        
        # Varianza explicada
        var_ratio = pca.explained_variance_ratio_ * 100
        
        df_pca_clean['PC1'] = components[:, 0]
        df_pca_clean['PC2'] = components[:, 1]
        
        # --- Visualización Scatter PCA ---
        st.markdown("**Proyección PCA (PC1 vs PC2)**")
        
        # Tooltip interactivo
        tooltip_cols = ['group', 'prepared']
        if 'label' in df_pca_clean.columns and pca_data_source == "Por Célula (Individual)":
            tooltip_cols.append('label')
            
        # --- Análisis Estadístico Adicional ---
        df_centroids = df_pca_clean.groupby('group')[['PC1', 'PC2']].mean().reset_index()
        
        # Test Multivariado (Hotelling T2) sobre PC1 y PC2
        ctl_pc = df_pca_clean[df_pca_clean['group'] == 'CTL'][['PC1', 'PC2']]
        hip_pc = df_pca_clean[df_pca_clean['group'] == 'Hypoxia'][['PC1', 'PC2']]
        
        
        has_stats = False
        stat_error = None
        if len(ctl_pc) >= 3 and len(hip_pc) >= 3:
            try:
                res_mv = pg.multivariate_ttest(ctl_pc, hip_pc)
                # OJO: Pingouin usa 'pval' como nombre de columna
                p_mv = res_mv['pval'].values[0]
                t2_mv = res_mv['T2'].values[0]
                has_stats = True
            except Exception as e:
                stat_error = str(e)
        else:
            stat_error = f"Insuficientes datos (CTL: {len(ctl_pc)}, Hypoxia: {len(hip_pc)}). Se necesitan al menos 3 por grupo."
        
        # Calcular Elipses de Confianza (95%)
        def get_ellipse_df(data, group_name):
            if len(data) < 5: return pd.DataFrame()
            try:
                mu = data.mean(axis=0).values
                cov = data.cov().values
                if np.any(np.isnan(cov)) or np.any(np.isinf(cov)): return pd.DataFrame()
                
                vals, vecs = np.linalg.eigh(cov)
                vals = np.maximum(vals, 1e-10) # Estabilidad numérica
                order = vals.argsort()[::-1]
                vals, vecs = vals[order], vecs[:, order]
                
                # 95% Confidence interval for 2D normal distribution (sqrt(chi2.ppf(0.95, 2)) approx 2.447)
                n_std = 2.447 
                
                t = np.linspace(0, 2*np.pi, 100)
                a, b = n_std * np.sqrt(vals[0]), n_std * np.sqrt(vals[1])
                ellipse_x = mu[0] + a * np.cos(t) * vecs[0, 0] + b * np.sin(t) * vecs[0, 1]
                ellipse_y = mu[1] + a * np.cos(t) * vecs[1, 0] + b * np.sin(t) * vecs[1, 1]
                
                # Crear DataFrame con índice de orden explícito
                df_e = pd.DataFrame({
                    'PC1': ellipse_x, 
                    'PC2': ellipse_y, 
                    'group': group_name,
                    'order_idx': np.arange(len(ellipse_x))
                })
                # Cerrar el loop agregando el primer punto al final
                first_point = df_e.iloc[[0]].copy()
                first_point['order_idx'] = len(ellipse_x)
                return pd.concat([df_e, first_point], ignore_index=True)
            except:
                return pd.DataFrame()

        df_ellipses = pd.concat([
            get_ellipse_df(ctl_pc, 'CTL'),
            get_ellipse_df(hip_pc, 'Hypoxia')
        ], ignore_index=True)

        # --- Visualización Scatter PCA Mejorada ---
        st.markdown("**Proyección PCA con Centroides y Elipses de Confianza (95%)**")
        
        # Base chart
        base = alt.Chart(df_pca_clean).encode(
            x=alt.X('PC1:Q', title=f'PC1 ({var_ratio[0]:.1f}%)'),
            y=alt.Y('PC2:Q', title=f'PC2 ({var_ratio[1]:.1f}%)'),
            color=alt.Color('group:N', scale=alt.Scale(domain=['CTL', 'Hypoxia'], range=['#377eb8', '#e41a1c']))
        )

        # Puntos individuales discretos
        points = base.mark_circle(size=30, opacity=0.3).encode(
            tooltip=tooltip_cols
        )

        # Elipses como contornos limpios (Sin líneas internas fantasma)
        ellipses = alt.Chart(df_ellipses).mark_line(
            strokeWidth=2, opacity=0.9
        ).encode(
            x='PC1:Q', y='PC2:Q', 
            color='group:N',
            detail='group:N', # Separa líneas por grupo
            order='order_idx:Q' # FUERZA el orden de conexión para evitar zigzags
        )

        # Centroides destacados
        centroids = alt.Chart(df_centroids).mark_point(
            size=250, filled=False, shape='cross', strokeWidth=3
        ).encode(
            x='PC1:Q', y='PC2:Q', color='group:N'
        )

        final_pca_chart = (points + ellipses + centroids).properties(height=500).interactive()
        st.altair_chart(final_pca_chart, use_container_width=True)
        _download_df_button(df_pca_clean, filename_base="datos_grafico_pca", label="⬇️ Descargar datos de este gráfico (CSV)")

        # --- Resultados del Test Multivariado (Siempre visible) ---
        st.markdown("#### Estadística de Segregación de Grupos")
        if has_stats:
            c1_pos = df_centroids[df_centroids['group'] == 'CTL'][['PC1', 'PC2']].values[0]
            c2_pos = df_centroids[df_centroids['group'] == 'Hypoxia'][['PC1', 'PC2']].values[0]
            dist_c = np.linalg.norm(c1_pos - c2_pos)
            
            scol1, scol2, scol3 = st.columns(3)
            scol1.metric("Distancia entre Centroides", f"{dist_c:.3f}")
            scol2.metric("Hotelling T²", f"{t2_mv:.2f}")
            
            sig_icon = "🔴" if p_mv < 0.05 else "⚪"
            p_text = f"p < 0.001" if p_mv < 0.001 else f"p = {p_mv:.3f}"
            scol3.metric("P-valor (Multivariado)", f"{p_text} {sig_icon}")
            
            if p_mv < 0.05:
                st.success(f"✅ **Diferencia Significativa:** Existe una segregación global sólida entre grupos en el espacio PCA (p={p_mv:.3g}).")
            else:
                st.info(f"ℹ️ **Sin Diferencia Significativa:** No hay evidencia estadística de segregación clara (p={p_mv:.3g}).")
        else:
            st.warning(f"⚠️ **No se pudo realizar el Test de Hotelling:** {stat_error}")
            st.caption("Asegúrate de tener suficientes datos cargados y procesados para ambos grupos.")

        # --- Guía de Interpretación ---
        with st.expander("❓ ¿Cómo leer este análisis?", expanded=False):
            st.markdown("""
            ### Guía rápida de interpretación del PCA
            
            1. **¿Qué son los ejes PC1 y PC2?**
               Son las variaciones morfológicas más importantes. El porcentaje (%) indica cuánta variabilidad total capturan.
            
            2. **Los Puntos**: Cada punto es un individuo (o célula). Puntos cercanos tienen morfologías similares.
            
            3. **Los Centroides (Cruces ✚)**: Representan el "promedio" de cada grupo. 
               - **Distancia**: Mayor distancia indica mayor diferencia morfológica global.
            
            4. **Las Elipses**: Área de confianza del 95%.
               - **Separadas**: Grupos morfológicamente distintos.
               - **Solapadas**: Gran variabilidad compartida entre grupos.
            
            5. **Test de Hotelling (Multivariado)**:
               - Es el test estadístico para comparar los centroides considerando juntas PC1 y PC2.
               - **p < 0.05 (🔴)**: La separación es estadísticamente real.
            """)
        
        # Visualización Factor Loadings
        with st.expander("🔍 Ver contribución de las variables originales (Loadings)"):
            st.markdown("Muestra qué variables empujan a las células hacia valores positivos o negativos en cada Componente Principal.")
            
            # Cargar loadings
            loadings = pca.components_.T
            df_loadings = pd.DataFrame(loadings, columns=['PC1', 'PC2'], index=selected_pca_features)
            
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

