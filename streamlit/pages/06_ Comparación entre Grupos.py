import json
from pathlib import Path
import altair as alt
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
import streamlit as st
from ui.sidebar import render_sidebar
from ui.plots import GROUP_SCALE, boxplot_with_ticks, apply_theme

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
@st.cache_data(ttl=300) # Cachear por 5 min
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
        
        # 1. Cargar métricas de Sholl 2D (sholl_2d_native.csv + sholl_summary.csv)
        pss = od / "sholl_summary.csv"
        if pss.exists():
            try:
                df = pd.read_csv(pss); df['prepared'] = p.name; df['group'] = group
                rows_sholl.append(df)
            except Exception: pass
            
        # 2. Cargar métricas de núcleo (03_nucleus_metrics.csv)
        pn = od / "03_nucleus_metrics.csv"
        if pn.exists():
            try:
                df = pd.read_csv(pn); df['prepared'] = p.name; df['group'] = group
                rows_nuc.append(df)
            except Exception: pass

    if not rows_sholl and not rows_nuc:
        return None, None

    # --- Fusión y Agregación ---
    df_nuc_all = pd.concat(rows_nuc, ignore_index=True) if rows_nuc else pd.DataFrame()
    df_sholl_all = pd.concat(rows_sholl, ignore_index=True) if rows_sholl else pd.DataFrame()
    
    # Unir métricas por célula (label + prepared)
    df_por_celula = pd.DataFrame()
    if not df_sholl_all.empty:
        df_por_celula = df_sholl_all
    
    # Unir métricas de núcleo
    if not df_nuc_all.empty:
        df_nuc_astros = df_nuc_all[df_nuc_all['is_astrocyte_candidate'] == True].copy()
        if 'label' in df_nuc_astros.columns:
             df_nuc_astros['label'] = df_nuc_astros['label'].astype(int)
        
        if not df_por_celula.empty and 'label' in df_por_celula.columns:
             df_por_celula = pd.merge(df_por_celula, df_nuc_astros, on=['label', 'prepared', 'group'], how='left')
        elif df_por_celula.empty: # Si no hay Sholl, al menos usar los núcleos
             df_por_celula = df_nuc_astros
        
    if df_por_celula.empty:
        return None, None

    # --- ¡CRÍTICO! Solución a Pseudoreplicación ---
    # 1. DataFrame para gráficos (todas las células)
    df_plot = df_por_celula.copy()
    
    # 2. DataFrame para estadística (mediana por preparado)
    cols_to_agg = [
        'critical_radius_um', 'peak_intersections', 'auc',
        'nucleus_volume_um3', 'nucleus_sphericity'
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

# --- Definir métricas disponibles para comparar ---
METRIC_OPTIONS = {
    # Métrica: (Nombre Bonito, Unidad)
    'auc': ("Sholl (AUC)", "µm·intersec"),
    'peak_intersections': ("Sholl (Pico)", "intersec"),
    'critical_radius_um': ("Sholl (Radio Crítico)", "µm"),
    'nucleus_volume_um3': ("Volumen del Núcleo", "µm³"),
    'nucleus_sphericity': ("Esfericidad del Núcleo", "0-1"),
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
            reporte_stats = format_test_results(test_results)
            st.markdown(reporte_stats)
        else:
            st.caption("No hay suficientes datos por preparado para ambos grupos para el test.")
    except Exception as e:
        st.caption(f"Error al calcular estadísticas: {e}")
        st.exception(e)

st.markdown("---")
st.markdown("### Datos Agregados")
col1, col2 = st.columns(2)
with col1:
    st.markdown("**Datos por Célula (para gráficos)**")
    st.dataframe(df_plot.round(3), use_container_width=True, height=400)
    _download_df_button(df_plot, "metricas_por_celula")
with col2:
    st.markdown("**Datos por Preparado (para estadística)**")
    st.dataframe(df_stats.round(3), use_container_width=True, height=400)
    _download_df_button(df_stats, "metricas_por_preparado")