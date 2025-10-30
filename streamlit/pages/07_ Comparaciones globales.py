import json
from pathlib import Path

import altair as alt
import numpy as np
import pandas as pd
import tifffile
from scipy import stats
import streamlit as st
from ui.sidebar import render_sidebar

# Raíz del repo
ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PROC_DIR = ROOT / "data" / "processed"

st.title("Comparaciones globales entre grupos")
render_sidebar(show_calibration=True)

# Estética: fondo y tarjetas
st.markdown(
    """
    <style>
    .stApp { background: linear-gradient(180deg, #0b1220 0%, #101826 100%) !important; }
    .block-container { padding-top: 1rem; }
    .section-card { background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.08); border-radius: 12px; padding: 1rem 1.25rem; margin: 1rem 0 1.25rem 0; }
    .section-card h4 { margin-top: 0; }
    </style>
    """,
    unsafe_allow_html=True,
)

# Paleta consistente
GROUP_SCALE = alt.Scale(domain=['CTL', 'Hipoxia'], range=['#1f77b4', '#d62728'])


def boxplot_with_ticks(df: pd.DataFrame, x_col: str, color_col: str, title_x: str = ""):
    """Horizontal boxplot with semi-transparent ticks for distribution."""
    base = alt.Chart(df)
    box = base.mark_boxplot(outliers=True).encode(
        x=alt.X(f"{x_col}:Q", title=title_x),
        y=alt.Y(f"{color_col}:N", title=None),
        color=alt.Color(f"{color_col}:N", scale=GROUP_SCALE, legend=None),
    )
    ticks = base.mark_tick(opacity=0.35, thickness=1).encode(
        x=alt.X(f"{x_col}:Q"),
        y=alt.Y(f"{color_col}:N", title=None),
        color=alt.Color(f"{color_col}:N", scale=GROUP_SCALE, legend=None),
    )
    return box + ticks


def violin_density(df: pd.DataFrame, x_col: str, color_col: str, title_x: str = "", height: int = 160):
    return (
        alt.Chart(df)
        .transform_density(x_col, as_=[x_col, 'density'], groupby=[color_col])
        .mark_area(opacity=0.25)
        .encode(
            x=alt.X(f"{x_col}:Q", title=title_x),
            y=alt.Y('density:Q', title='densidad'),
            color=alt.Color(f"{color_col}:N", scale=GROUP_SCALE),
        )
        .properties(height=height)
    )


def histogram(df: pd.DataFrame, x_col: str, color_col: str, title_x: str, maxbins: int = 30, height: int = 160):
    return (
        alt.Chart(df)
        .mark_bar(opacity=0.45)
        .encode(
            x=alt.X(f"{x_col}:Q", bin=alt.Bin(maxbins=maxbins), title=title_x),
            y=alt.Y('count():Q', title='conteo'),
            color=alt.Color(f"{color_col}:N", scale=GROUP_SCALE),
        )
        .properties(height=height)
    )


def filter_outliers_iqr(df: pd.DataFrame, col: str, by_group: bool = True, k: float = 1.5) -> pd.DataFrame:
    if df.empty or col not in df.columns:
        return df
    if by_group and 'group' in df.columns:
        def _f(g):
            q1 = g[col].quantile(0.25)
            q3 = g[col].quantile(0.75)
            iqr = q3 - q1
            lo, hi = q1 - k*iqr, q3 + k*iqr
            return g[(g[col] >= lo) & (g[col] <= hi)]
        return df.groupby('group', group_keys=False).apply(_f)
    q1 = df[col].quantile(0.25)
    q3 = df[col].quantile(0.75)
    iqr = q3 - q1
    lo, hi = q1 - k*iqr, q3 + k*iqr
    return df[(df[col] >= lo) & (df[col] <= hi)]


def list_raw_images(raw_dir: Path) -> list[Path]:
    files = []
    for sub in ["CTL", "hip"]:
        d = raw_dir / sub
        if not d.exists():
            continue
        for ext in ("*.tif", "*.tiff", "*.TIF", "*.TIFF"):
            files.extend(sorted(d.rglob(ext)))
    return files


def detect_group_from_path(p: Path) -> str:
    return 'Hipoxia' if '/hip/' in str(p).replace('\\', '/').lower() else 'CTL'


# Recolectar datos agregados
files = list_raw_images(RAW_DIR)
rows_sum = []
rows_sh = []
rows_sh_full = []

for p in files:
    od = PROC_DIR / p.stem
    pg = od / "skeletons" / "summary.csv"
    if pg.exists():
        try:
            dfi = pd.read_csv(pg)
            dfi['prepared'] = p.name
            dfi['group'] = detect_group_from_path(p)
            rows_sum.append(dfi)
        except Exception:
            pass
    ps = od / "sholl.csv"
    if ps.exists():
        try:
            dff = pd.read_csv(ps)
            dff['prepared'] = p.name
            dff['group'] = detect_group_from_path(p)
            rows_sh_full.append(dff)
            # Pico por célula
            dfi = dff.sort_values(['label','intersections'], ascending=[True, False]).groupby('label', as_index=False).first()
            rows_sh.append(dfi)
        except Exception:
            pass

st.markdown("---")

# Opciones de visualización y filtros
with st.expander("Opciones de visualización y filtros", expanded=True):
    show_violin = st.checkbox("Mostrar violines", value=True)
    show_hist = st.checkbox("Mostrar histogramas", value=True)
    apply_outliers = st.checkbox("Filtrar outliers (IQR)", value=True)
    by_group = st.checkbox("Filtrar por grupo (IQR por grupo)", value=True)
    iqr_k = st.slider("Factor IQR", 1.0, 3.0, 1.5, 0.5)
    maxbins = st.slider("Bins del histograma", 10, 60, 30, step=5)

# 1) N° astrocitos por preparado
count_rows = []
for p in files:
    od = PROC_DIR / p.stem
    fn = od / "04_final_astrocytes_mask.tif"
    if fn.exists():
        try:
            lab = tifffile.imread(fn)
            count_rows.append({
                "prepared": p.name,
                "group": detect_group_from_path(p),
                "n_astro": int((np.unique(lab) > 0).sum()),
            })
        except Exception:
            pass
if count_rows:
    df_counts = pd.DataFrame(count_rows)
    st.markdown("#### N° de astrocitos por preparado")
    st.altair_chart(boxplot_with_ticks(df_counts, 'n_astro', 'group', title_x='N° astrocitos por preparado'), use_container_width=True)
    try:
        if set(df_counts['group'].unique()) >= {'CTL','Hipoxia'}:
            a = df_counts.loc[df_counts['group']=='CTL','n_astro'].to_numpy()
            b = df_counts.loc[df_counts['group']=='Hipoxia','n_astro'].to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (n_astro): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
    except Exception:
        pass

# 2) Longitud y volumen de dominio por célula
if rows_sum:
    df_all = pd.concat(rows_sum, ignore_index=True)
    
    def _render_metric(df: pd.DataFrame, col: str, title: str):
        if col not in df.columns or not df[col].notna().any():
            return df
        base_df = df.dropna(subset=[col])
        before = base_df.shape[0]
        fdf = filter_outliers_iqr(base_df, col, by_group=by_group, k=iqr_k) if apply_outliers else base_df
        removed = before - fdf.shape[0]
        st.altair_chart(boxplot_with_ticks(fdf, col, 'group', title_x=title), use_container_width=True)
        if removed > 0:
            st.caption(f"Filtrados {removed} outliers ({removed/before:.1%}) por IQR{' por grupo' if by_group else ''}.")
        # resumen por grupo
        try:
            grp = fdf.groupby('group')[col]
            stats_df = grp.agg(['count','mean','median']).reset_index().rename(columns={'mean':'media','median':'mediana'})
            st.dataframe(stats_df, use_container_width=True)
        except Exception:
            pass
        if show_violin:
            st.altair_chart(violin_density(fdf, col, 'group', title_x=title), use_container_width=True)
        if show_hist:
            st.altair_chart(histogram(fdf, col, 'group', title_x=title, maxbins=maxbins), use_container_width=True)
        return fdf
    # --- Por célula ---
    st.markdown("#### Por célula: Longitud total")
    if 'total_length_um' in df_all.columns:
        f_len = _render_metric(df_all, 'total_length_um', 'Longitud total (µm)')
        try:
            a = f_len.loc[f_len['group']=='CTL','total_length_um'].dropna().to_numpy()
            b = f_len.loc[f_len['group']=='Hipoxia','total_length_um'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (longitud total): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass
    st.markdown("#### Por célula: Volumen de dominio")
    if 'domain_volume_um3' in df_all.columns and df_all['domain_volume_um3'].notna().any():
        f_dom = _render_metric(df_all, 'domain_volume_um3', 'Volumen dominio (µm³)')
        try:
            a = f_dom.loc[f_dom['group']=='CTL','domain_volume_um3'].dropna().to_numpy()
            b = f_dom.loc[f_dom['group']=='Hipoxia','domain_volume_um3'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (volumen dominio): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass

    # --- Promedio por preparado ---
    if 'total_length_um' in df_all.columns:
        mean_len = df_all.groupby(['prepared','group'], as_index=False)['total_length_um'].mean()
        st.markdown("#### Longitud promedio por preparado")
        st.altair_chart(
            boxplot_with_ticks(mean_len, 'total_length_um', 'group', title_x='Longitud promedio (µm)'),
            use_container_width=True,
        )
        try:
            a = mean_len.loc[mean_len['group']=='CTL','total_length_um'].dropna().to_numpy()
            b = mean_len.loc[mean_len['group']=='Hipoxia','total_length_um'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (longitud promedio por preparado): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass
    if 'domain_volume_um3' in df_all.columns and df_all['domain_volume_um3'].notna().any():
        mean_dom = df_all.groupby(['prepared','group'], as_index=False)['domain_volume_um3'].mean()
        st.markdown("#### Volumen de dominio promedio por preparado")
        st.altair_chart(
            boxplot_with_ticks(mean_dom, 'domain_volume_um3', 'group', title_x='Volumen de dominio promedio (µm³)'),
            use_container_width=True,
        )
        try:
            a = mean_dom.loc[mean_dom['group']=='CTL','domain_volume_um3'].dropna().to_numpy()
            b = mean_dom.loc[mean_dom['group']=='Hipoxia','domain_volume_um3'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (volumen dominio promedio por preparado): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass

    # 3) Volumen en tubo alrededor del esqueleto
    if 'tube_volume_um3' in df_all.columns and df_all['tube_volume_um3'].notna().any():
        st.markdown("#### Volumen en tubo alrededor del esqueleto")
        f_tube = _render_metric(df_all, 'tube_volume_um3', 'Volumen en tubo (µm³)')
        try:
            a = f_tube.loc[f_tube['group']=='CTL','tube_volume_um3'].dropna().to_numpy()
            b = f_tube.loc[f_tube['group']=='Hipoxia','tube_volume_um3'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (volumen tubo por célula): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass
        # Suma por preparado
        try:
            sums = df_all.groupby(['prepared','group'])['tube_volume_um3'].sum().reset_index()
            st.altair_chart(boxplot_with_ticks(sums, 'tube_volume_um3', 'group', title_x='Volumen en tubo — suma por preparado (µm³)'), use_container_width=True)
            if set(sums['group'].unique()) >= {'CTL','Hipoxia'}:
                a = sums.loc[sums['group']=='CTL','tube_volume_um3'].to_numpy()
                b = sums.loc[sums['group']=='Hipoxia','tube_volume_um3'].to_numpy()
                if a.size and b.size:
                    u = stats.mannwhitneyu(a,b, alternative='two-sided')
                    st.caption(f"Mann–Whitney U (suma por preparado): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass

    # 4) GFAP y grosor
    if 'tube_sum_intensity_per_um' in df_all.columns and df_all['tube_sum_intensity_per_um'].notna().any():
        st.markdown("#### GFAP — Intensidad en el tubo por µm (proxy de expresión)")
        f_gfap = _render_metric(df_all, 'tube_sum_intensity_per_um', 'Intensidad/µm (tubo)')
        try:
            a = f_gfap.loc[f_gfap['group']=='CTL','tube_sum_intensity_per_um'].dropna().to_numpy()
            b = f_gfap.loc[f_gfap['group']=='Hipoxia','tube_sum_intensity_per_um'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (GFAP intensidad/µm): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass
    if 'local_radius_um_median' in df_all.columns and df_all['local_radius_um_median'].notna().any():
        st.markdown("#### Grosor de prolongaciones (radio local mediano)")
        f_rad = _render_metric(df_all, 'local_radius_um_median', 'Radio local mediano (µm)')
        try:
            a = f_rad.loc[f_rad['group']=='CTL','local_radius_um_median'].dropna().to_numpy()
            b = f_rad.loc[f_rad['group']=='Hipoxia','local_radius_um_median'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (radio local mediano): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass

    # 5) Ramificación
    if 'n_branches' in df_all.columns and df_all['n_branches'].notna().any():
        st.markdown("#### Ramificación (número de ramas)")
        f_nb = _render_metric(df_all, 'n_branches', 'N° ramas por célula')
        try:
            a = f_nb.loc[f_nb['group']=='CTL','n_branches'].dropna().to_numpy()
            b = f_nb.loc[f_nb['group']=='Hipoxia','n_branches'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (n° ramas): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass
    if 'branch_density_per_100um' in df_all.columns and df_all['branch_density_per_100um'].notna().any():
        st.markdown("#### Densidad de ramificación (ramas por 100 µm)")
        f_bd = _render_metric(df_all, 'branch_density_per_100um', 'Ramas por 100 µm')
        try:
            a = f_bd.loc[f_bd['group']=='CTL','branch_density_per_100um'].dropna().to_numpy()
            b = f_bd.loc[f_bd['group']=='Hipoxia','branch_density_per_100um'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (densidad ramas): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass

    # Interpretación
    try:
        st.markdown("### Interpretación (Hipoxia vs CTL)")
        def med(df, col):
            s = df[col].dropna()
            return float(s.median()) if not s.empty else np.nan
        def arrow(delta):
            return '↑' if delta > 0 else ('↓' if delta < 0 else '→')
        g = {}
        for col in [
            'tube_sum_intensity_per_um', 'local_radius_um_median', 'total_length_um',
            'n_branches', 'branch_density_per_100um', 'domain_volume_um3']:
            if col in df_all.columns:
                g[col] = {
                    'CTL': med(df_all.loc[df_all['group']=='CTL'], col),
                    'Hipoxia': med(df_all.loc[df_all['group']=='Hipoxia'], col),
                }
        lines = []
        def add_line(name, key, expected):
            if key not in g: return
            ctl = g[key]['CTL']; hip = g[key]['Hipoxia']
            if np.isnan(ctl) or np.isnan(hip): return
            delta = hip - ctl
            ok = (delta>0 and expected=="up") or (delta<0 and expected=="down")
            lines.append(f"- {name}: Hipoxia {arrow(delta)} CTL (medianas: {hip:.2f} vs {ctl:.2f}) — {'concordante' if ok else 'discordante'}")
        add_line("GFAP intensidad/µm (tubo)", 'tube_sum_intensity_per_um', "up")
        add_line("Grosor (radio local mediano)", 'local_radius_um_median', "up")
        add_line("Longitud total", 'total_length_um', "down")
        add_line("N° ramas", 'n_branches', "down")
        add_line("Densidad de ramas/100µm", 'branch_density_per_100um', "down")
        add_line("Volumen de dominio", 'domain_volume_um3', "down")
        if rows_sh_full:
            df_sh_full = pd.concat(rows_sh_full, ignore_index=True)
            band = df_sh_full[(df_sh_full['radius_um'] >= 5.0) & (df_sh_full['radius_um'] <= 20.0)]
            if not band.empty:
                ctl = float(band.loc[band['group']=='CTL','intersections'].mean()) if (band['group']=='CTL').any() else np.nan
                hip = float(band.loc[band['group']=='Hipoxia','intersections'].mean()) if (band['group']=='Hipoxia').any() else np.nan
                if not (np.isnan(ctl) or np.isnan(hip)):
                    delta = hip - ctl
                    ok = delta < 0
                    lines.append(f"- Complejidad Sholl 5–20µm: Hipoxia {arrow(delta)} CTL (promedios: {hip:.2f} vs {ctl:.2f}) — {'concordante' if ok else 'discordante'}")
        if lines:
            st.markdown("\n".join(lines))
            st.caption("‘Concordante’ indica dirección esperada por literatura para astrogliosis reactiva.")
    except Exception:
        pass

# Sholl: curvas por grupo y tests
if rows_sh or rows_sh_full:
    if rows_sh_full:
        df_curves = pd.concat(rows_sh_full, ignore_index=True)
        st.markdown("#### Sholl — Curvas por grupo")
        try:
            med = df_curves.groupby(['group','radius_um'])['intersections'].median().reset_index()
            q1 = df_curves.groupby(['group','radius_um'])['intersections'].quantile(0.25).reset_index().rename(columns={'intersections':'q1'})
            q3 = df_curves.groupby(['group','radius_um'])['intersections'].quantile(0.75).reset_index().rename(columns={'intersections':'q3'})
            band = med.merge(q1,on=['group','radius_um']).merge(q3,on=['group','radius_um'])
            max_r = float(df_curves['radius_um'].max()) if not df_curves.empty else 50.0
            line = alt.Chart(med).mark_line().encode(
                x=alt.X('radius_um:Q', scale=alt.Scale(domain=[0, max_r])),
                y='intersections:Q',
                color=alt.Color('group:N', scale=GROUP_SCALE)
            )
            ribbon = alt.Chart(band).mark_area(opacity=0.15).encode(
                x=alt.X('radius_um:Q', scale=alt.Scale(domain=[0, max_r])),
                y='q1:Q', y2='q3:Q',
                color=alt.Color('group:N', scale=GROUP_SCALE)
            )
            st.altair_chart(ribbon + line, use_container_width=True)
        except Exception:
            pass

        # Tests por radio con BH
        try:
            radii = sorted(df_curves['radius_um'].unique())
            pvals = []
            for r in radii:
                a = df_curves.loc[(df_curves['group']=='CTL') & (df_curves['radius_um']==r), 'intersections'].dropna().to_numpy()
                b = df_curves.loc[(df_curves['group']=='Hipoxia') & (df_curves['radius_um']==r), 'intersections'].dropna().to_numpy()
                if a.size and b.size:
                    u = stats.mannwhitneyu(a,b, alternative='two-sided')
                    pvals.append((r, float(u.pvalue)))
            if pvals:
                m = len(pvals)
                ranked = sorted(pvals, key=lambda x: x[1])
                adj = []
                for i,(r,p) in enumerate(ranked, start=1):
                    adj_p = min(p * m / i, 1.0)
                    adj.append((r, p, adj_p))
                sig = [(r,p,ap) for (r,p,ap) in adj if ap <= 0.05]
                if sig:
                    st.caption(f"Rangos con diferencias (FDR≤0.05): {len(sig)}/{m}")
                    st.dataframe(pd.DataFrame(sig, columns=['radius_um','p','p_adj']).round(4), use_container_width=True)
        except Exception:
            pass

        # Cobertura por radio
        try:
            totals = df_curves.groupby('group')['label'].nunique().to_dict()
            df_pos = df_curves.assign(pos=(df_curves['intersections'] > 0).astype(int))
            cov = df_pos.groupby(['group','radius_um'])['pos'].sum().reset_index()
            cov['fraction'] = cov.apply(lambda r: (r['pos'] / max(1, totals.get(r['group'], 1))), axis=1)
            max_r = float(df_curves['radius_um'].max()) if not df_curves.empty else 50.0
            cov_chart = alt.Chart(cov).mark_line().encode(
                x=alt.X('radius_um:Q', title='radio (µm)', scale=alt.Scale(domain=[0, max_r])),
                y=alt.Y('fraction:Q', title='fracción de células con intersección>0'),
                color=alt.Color('group:N', scale=GROUP_SCALE)
            ).properties(height=220)
            st.altair_chart(cov_chart, use_container_width=True)
        except Exception:
            pass

        # Curvas normalizadas por célula
        with st.expander("Curvas normalizadas por célula (radio 0–1)", expanded=False):
            try:
                d = df_curves.copy()
                dpos = d[d['intersections'] > 0]
                rmax = dpos.groupby(['prepared','label'], as_index=False)['radius_um'].max().rename(columns={'radius_um':'r_max'})
                dn = d.merge(rmax, on=['prepared','label'], how='left')
                dn = dn[dn['r_max'].notna() & (dn['r_max'] > 0)]
                dn['norm_r'] = (dn['radius_um'] / dn['r_max']).clip(0, 1.0)
                step = float(st.selectbox("Resolución de bins", options=[0.02, 0.05, 0.10], index=1))
                dn['bin'] = (np.round(dn['norm_r'] / step) * step).clip(0, 1.0)
                medn = dn.groupby(['group','bin'])['intersections'].median().reset_index()
                q1n = dn.groupby(['group','bin'])['intersections'].quantile(0.25).reset_index().rename(columns={'intersections':'q1'})
                q3n = dn.groupby(['group','bin'])['intersections'].quantile(0.75).reset_index().rename(columns={'intersections':'q3'})
                bandn = medn.merge(q1n, on=['group','bin']).merge(q3n, on=['group','bin'])
                line_n = alt.Chart(medn).mark_line().encode(
                    x=alt.X('bin:Q', title='radio normalizado (0–1)', scale=alt.Scale(domain=[0,1])),
                    y='intersections:Q',
                    color=alt.Color('group:N', scale=GROUP_SCALE)
                )
                rib_n = alt.Chart(bandn).mark_area(opacity=0.15).encode(
                    x=alt.X('bin:Q', title='radio normalizado (0–1)', scale=alt.Scale(domain=[0,1])),
                    y='q1:Q', y2='q3:Q', color=alt.Color('group:N', scale=GROUP_SCALE)
                )
                st.altair_chart(rib_n + line_n, use_container_width=True)
                # Cobertura normalizada
                totals = dn.groupby('group')['label'].nunique().to_dict()
                dn_pos = dn[dn['intersections'] > 0].drop_duplicates(['group','label','bin'])
                covn = dn_pos.groupby(['group','bin'])['label'].nunique().reset_index().rename(columns={'label':'n_labels'})
                covn['fraction'] = covn.apply(lambda r: (r['n_labels'] / max(1, totals.get(r['group'], 1))), axis=1)
                covn_chart = alt.Chart(covn).mark_line().encode(
                    x=alt.X('bin:Q', title='radio normalizado (0–1)', scale=alt.Scale(domain=[0,1])),
                    y=alt.Y('fraction:Q', title='fracción de células con intersección>0'),
                    color=alt.Color('group:N', scale=GROUP_SCALE)
                ).properties(height=220)
                st.altair_chart(covn_chart, use_container_width=True)

                # Tests por ventanas con BH
                win_edges = [(0.10,0.30), (0.30,0.50), (0.50,0.70), (0.70,0.90)]
                rows = []
                for a,b in win_edges:
                    band = dn[(dn['norm_r'] >= a) & (dn['norm_r'] < b)]
                    if band.empty:
                        continue
                    per_cell = band.groupby(['group','prepared','label'], as_index=False)['intersections'].mean().rename(columns={'intersections':'mean_intersections_norm'})
                    a_vals = per_cell.loc[per_cell['group']=='CTL','mean_intersections_norm'].dropna().to_numpy()
                    b_vals = per_cell.loc[per_cell['group']=='Hipoxia','mean_intersections_norm'].dropna().to_numpy()
                    if a_vals.size and b_vals.size:
                        u = stats.mannwhitneyu(a_vals, b_vals, alternative='two-sided')
                        rows.append({'win': f"{a:.2f}-{b:.2f}", 'U': float(u.statistic), 'p': float(u.pvalue), 'ctl_mean': float(np.mean(a_vals)), 'hip_mean': float(np.mean(b_vals))})
                if rows:
                    dfw = pd.DataFrame(rows)
                    m = dfw.shape[0]
                    order = np.argsort(dfw['p'].to_numpy())
                    p_sorted = dfw['p'].to_numpy()[order]
                    p_adj_sorted = np.minimum.accumulate(p_sorted * m / (np.arange(m)+1)[::-1])[::-1]
                    dfw.loc[order, 'p_adj'] = p_adj_sorted
                    dfw['dir'] = dfw.apply(lambda r: '↑' if r['hip_mean']>r['ctl_mean'] else ('↓' if r['hip_mean']<r['ctl_mean'] else '→'), axis=1)
                    st.markdown("##### Tests por ventanas de radio normalizado (BH)")
                    st.dataframe(dfw[['win','U','p','p_adj','ctl_mean','hip_mean','dir']].round(4), use_container_width=True)
            except Exception:
                pass

    # Pico por grupo
    if rows_sh:
        df_peak = pd.concat(rows_sh, ignore_index=True)
        st.markdown("#### Sholl — Pico de intersecciones por grupo")
        st.altair_chart(boxplot_with_ticks(df_peak, 'intersections', 'group', title_x='Pico de intersecciones'), use_container_width=True)
        try:
            a = df_peak.loc[df_peak['group']=='CTL','intersections'].dropna().to_numpy()
            b = df_peak.loc[df_peak['group']=='Hipoxia','intersections'].dropna().to_numpy()
            if a.size and b.size:
                u = stats.mannwhitneyu(a,b, alternative='two-sided')
                st.caption(f"Mann–Whitney U (pico intersecciones): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
        except Exception:
            pass

# Complejidad 5–20 µm (curvas completas)
if rows_sh_full:
    try:
        df_sh_full = pd.concat(rows_sh_full, ignore_index=True)
        band = df_sh_full[(df_sh_full['radius_um'] >= 5.0) & (df_sh_full['radius_um'] <= 20.0)]
        if not band.empty:
            band_cell = band.groupby(['prepared','group','label'], as_index=False)['intersections'].mean().rename(columns={'intersections':'mean_intersections_5_20um'})
            st.markdown("#### Complejidad de ramificación (Sholl 5–20 µm)")
            st.altair_chart(boxplot_with_ticks(band_cell, 'mean_intersections_5_20um', 'group', title_x='Intersecciones promedio (5–20 µm)'), use_container_width=True)
            if set(band_cell['group'].unique()) >= {'CTL','Hipoxia'}:
                a = band_cell.loc[band_cell['group']=='CTL','mean_intersections_5_20um'].dropna().to_numpy()
                b = band_cell.loc[band_cell['group']=='Hipoxia','mean_intersections_5_20um'].dropna().to_numpy()
                if a.size and b.size:
                    u = stats.mannwhitneyu(a,b, alternative='two-sided')
                    st.caption(f"Mann–Whitney U (Sholl 5–20 µm): U={float(u.statistic):.1f}, p={float(u.pvalue):.3g}")
    except Exception:
        pass
