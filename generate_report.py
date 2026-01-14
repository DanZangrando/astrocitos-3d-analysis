
import sys
import pandas as pd
import numpy as np
from scipy import stats
import pingouin as pg
from pathlib import Path

# Add project root to path to import from streamlit.ui
sys.path.append('/home/daniel/Proyectos/astrocitos-3d-analysis')

# Define paths
ROOT = Path('/home/daniel/Proyectos/astrocitos-3d-analysis')
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"

# Import logic from Page 06 (we'll just copy the necessary bits to avoid Streamlit dependency issues if we run this as a script)
# Actually, loading load_all_metrics might fail because of st.cache_data. 
# Better to copy the loading logic or mock st.cache_data.
# Let's copy the essential loading logic for robustness.

def list_raw_images(raw_dir: Path) -> list[Path]:
    files = []
    for sub in ["CTL", "hip"]:
        d = raw_dir / sub
        if not d.exists(): continue
        for ext in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
            files.extend(sorted(d.rglob(ext)))
    return files

def detect_group_from_path(p: Path) -> str:
    return 'Hipoxia' if '/hip/' in str(p).replace('\\', '/').lower() else 'CTL'

def load_metrics():
    files = list_raw_images(RAW_DIR)
    rows_data = [] # Combined flat list for scalar metrics
    rows_sholl_curves = []
    
    for p in files:
        od = PROC_DIR / p.stem
        group = detect_group_from_path(p)
        prep_id = f"{group}-{p.stem}"
        
        # Load Skeleton Summary
        ps = od / "skeletons" / "summary.csv"
        skel_data = pd.DataFrame()
        if ps.exists():
            try: skel_data = pd.read_csv(ps)
            except: pass
            
        # Load Sholl Summary (Scalar metrics like Critical Radius)
        pss = od / "sholl_summary.csv"
        sholl_data = pd.DataFrame()
        if pss.exists():
            try: sholl_data = pd.read_csv(pss)
            except: pass
            
        # Merge scalar cell data
        if not skel_data.empty or not sholl_data.empty:
            # Assume label match or just concat if indices align (usually they do or we rely on label)
            # For simplicity in this reporting script, we'll assume standard merging logic
            # Just taking mean of the prep directly here to simplify matching issues? 
            # No, let's merge properly if possible.
            
            # Simplified: just collect all rows and agg by prep later.
            if not skel_data.empty:
                skel_data['prepared'] = prep_id
                skel_data['group'] = group
                rows_data.append(skel_data)
                
            if not sholl_data.empty:
                sholl_data['prepared'] = prep_id
                sholl_data['group'] = group
                rows_data.append(sholl_data)

        # Load Raw Sholl Curves
        pr = od / "sholl_2d_native.csv"
        if pr.exists():
            try:
                df = pd.read_csv(pr, usecols=['label', 'radius_um', 'intersections']) 
                df['prepared'] = prep_id
                df['group'] = group
                rows_sholl_curves.append(df)
            except: pass

    # Combine
    df_all_scalar = pd.concat(rows_data, ignore_index=True) if rows_data else pd.DataFrame()
    df_curves = pd.concat(rows_sholl_curves, ignore_index=True) if rows_sholl_curves else pd.DataFrame()
    
    return df_all_scalar, df_curves

def run_stat_test(df, metric_col):
    # Aggregation per prep
    df_agg = df.groupby(['prepared', 'group'])[metric_col].median().reset_index()
    
    ctl = df_agg[df_agg['group'] == 'CTL'][metric_col].dropna()
    hip = df_agg[df_agg['group'] == 'Hipoxia'][metric_col].dropna()
    
    if len(ctl) < 3 or len(hip) < 3:
        return f"Insufficient data (N_CTL={len(ctl)}, N_HIP={len(hip)})"

    # Normality
    norm_ctl = stats.shapiro(ctl).pvalue > 0.05 if len(ctl) >= 3 else False
    norm_hip = stats.shapiro(hip).pvalue > 0.05 if len(hip) >= 3 else False
    
    # Selection
    if norm_ctl and norm_hip:
        stat, p = stats.ttest_ind(ctl, hip, equal_var=False)
        test_name = "Welch's t-test"
        stat_name = "t"
    else:
        stat, p = stats.mannwhitneyu(ctl, hip)
        test_name = "Mann-Whitney U"
        stat_name = "U"
        
    return {
        "metric": metric_col,
        "mean_ctl": ctl.mean(),
        "mean_hip": hip.mean(),
        "sem_ctl": ctl.sem(),
        "sem_hip": hip.sem(),
        "test": test_name,
        "stat_val": stat,
        "p_val": p
    }

# Main Execution
print("loading data...")
df_scalar, df_curves = load_metrics()

metrics_of_interest = [
    'ramification_index',
    'total_branch_length_um',
    'n_endpoints',
    'critical_radius_um'
]

print("\n--- SCALAR RESULTS ---")
results_text = []
for m in metrics_of_interest:
    if m in df_scalar.columns:
        res = run_stat_test(df_scalar, m)
        print(res)
        results_text.append(res)
    else:
        print(f"Metric {m} not found in dataframe columns: {df_scalar.columns}")

print("\n--- SHOLL ANOVA ---")
if not df_curves.empty:
    # Aggregation per prep for ANOVA
    df_sholl_prep = df_curves.groupby(['group', 'prepared', 'radius_um'])['intersections'].mean().reset_index()
    
    try:
        aov = pg.mixed_anova(
            data=df_sholl_prep, 
            dv='intersections', 
            within='radius_um', 
            between='group', 
            subject='prepared'
        )
        print(aov)
    except Exception as e:
        print(f"ANOVA Failed: {e}")
