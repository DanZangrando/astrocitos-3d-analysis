import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import sys

# Simulando lo que hace streamlit/pages/06_ Comparación entre Grupos.py
sys.path.append('streamlit')

from ui.utils import list_raw_images, detect_group

raw_dir = Path("data/raw")
proc_dir = Path("data/processed")

def detect_group_from_path(p: Path) -> str:
    return 'Hypoxia' if '/hip/' in str(p).replace('\\', '/').lower() else 'CTL'

def load_data():
    files = list_raw_images(raw_dir)
    rows_skel = []
    
    for p in files:
        od = proc_dir / p.stem
        group = detect_group_from_path(p)
        ps = od / "skeletons" / "summary.csv"
        if ps.exists():
            df = pd.read_csv(ps)
            df['prepared'] = f"{group}-{p.stem}"; df['group'] = group; df['individual'] = p.parent.name
            rows_skel.append(df)
            
    df_plot = pd.concat(rows_skel, ignore_index=True)
    return df_plot

df_plot = load_data()
df_pca_source = df_plot

selected_pca_features = ['total_branch_length_um', 'n_endpoints', 'mean_branch_length_um', 'tortuosity_mean']

df_pca_clean = df_pca_source.dropna(subset=selected_pca_features).copy()
print("df_pca_clean columns:", df_pca_clean.columns)

tooltip_cols = ['group', 'individual']
tooltip_cols.append('prepared')
if 'label' in df_pca_clean.columns:
    tooltip_cols.append('label')
    
print("Tooltip cols:", tooltip_cols)

