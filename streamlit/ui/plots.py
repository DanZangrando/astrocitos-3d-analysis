from __future__ import annotations
# Se eliminaron: BytesIO, Optional, Union, vl_convert
import altair as alt

# Shared color scale for groups
PALETTE = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}
GROUP_SCALE = alt.Scale(domain=["CTL", "Hipoxia"], range=[PALETTE["CTL"], PALETTE["Hipoxia"]])


def configure_altair() -> None:
    """Ensure Altair uses inline data and no row limits for export.

    This improves compatibility with vl-convert and export backends.
    """
    try:
        alt.data_transformers.enable('default', max_rows=None)
    except Exception:
        # Safe to ignore if environment doesn't support this call
        pass


def apply_theme(st_module) -> None:
    """Inject a subtle gradient background and card-like sections."""
    st_module.markdown(
        """
        <style>
        .stApp { background: linear-gradient(180deg, #0b1220 0%, #101826 100%) !important; }
        .block-container { padding-top: 2.5rem; }
        .block-container h1 { margin-top: 0.25rem !important; }
        .section-card { background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.08);
                         border-radius: 12px; padding: 1rem 1.25rem; margin: 1rem 0 1.25rem 0; }
        .section-card h4 { margin-top: 0; }
        </style>
        """,
        unsafe_allow_html=True,
    )


def boxplot_with_ticks(df, x_col: str, color_col: str, title_x: str = ""):
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


def violin_density(df, x_col: str, color_col: str, title_x: str = "", height: int = 160):
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


def histogram(df, x_col: str, color_col: str, title_x: str, maxbins: int = 30, height: int = 160):
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
