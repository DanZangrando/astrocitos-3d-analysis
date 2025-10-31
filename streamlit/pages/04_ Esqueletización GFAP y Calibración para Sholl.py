import sys
import os
import json
import subprocess # <--- CORRECCIÓN
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
import pandas as pd
import altair as alt
from ui.sidebar import render_sidebar
# --- Importar la lógica unificada ---
from ui import pipeline, runner

st.title("Esqueletización GFAP y Calibración para Sholl")
render_sidebar(show_calibration=True)

root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
calib_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer.py"

# --- Utilidades ---
def _read_global_calibration():
    if calib_path.exists():
        try:
            return json.loads(calib_path.read_text())
        except Exception:
            pass
    return {}

def _save_global_calibration(cal_data: dict):
    calib_path.parent.mkdir(parents=True, exist_ok=True)
    calib_path.write_text(json.dumps(cal_data, indent=2))

def get_output_dir_for_image(img_path: Path) -> Path:
    base_name = img_path.stem
    out_dir = root / "data" / "processed" / base_name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir

def _detect_group(p: Path, root: Path) -> str:
    try:
        rel = str(p.relative_to(root)).lower()
    except Exception:
        rel = str(p).lower()
    if "/hip/" in rel:
        return "Hipoxia"
    if "/ctl/" in rel:
        return "CTL"
    return "CTL"

# 1) Selección de imagen
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
if not files:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

group_filter = st.session_state.get("group_filter", "Todos")
if group_filter == "Todos":
    files_avail = files
else:
    files_avail = [p for p in files if _detect_group(p, root) == group_filter]
if not files_avail:
    st.info("No hay preparados para el grupo seleccionado.")
    st.stop()

labels = [str(p.relative_to(root)) for p in files_avail]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files_avail))), format_func=lambda i: labels[i])
img_path = files_avail[idx]
out_dir = get_output_dir_for_image(img_path)
group = _detect_group(img_path, root)
glob_calib = _read_global_calibration()

mask_final = out_dir / "04_final_astrocytes_mask.tif"
if not mask_final.exists():
    st.error("No se encontró 04_final_astrocytes_mask.tif. Ejecutá el 'Paso 03: Filtrado' primero.")
    st.stop()
def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"
st.markdown(_group_badge_html(group), unsafe_allow_html=True)
st.caption(f"Usando máscara de entrada: {mask_final.relative_to(root)}")


# --------- Recalcular por ámbito (Batch) ---------
with st.expander("Recalcular por ámbito (Batch)", expanded=False):
    scope = st.radio("Ámbito (Batch)", options=["Grupo", "Todos"], horizontal=True, key="p04_scope")
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p04_scope_group")
    if st.button("▶️ Recalcular (desde 04)", key="p04_recalc"):
        cal = _read_global_calibration()
        sc = "group" if scope=="Grupo" else "all"
        res = runner.run_scope(root, scope=sc, start_step="04", cal=cal, selected=None, group=scope_group, overwrite_from_step=True)
        ok = sum(1 for _, stt in res if not stt.get("error"))
        st.success(f"Listo: {ok}/{len(res)} preparados procesados desde 04.")

# 2) Parámetros de esqueletización (UI)
n_channels = 1
try:
    arr_preview, axes_prev = pipeline.load_image_any(img_path)
    vol_prev = pipeline.reorder_to_zcyx(arr_preview, axes_prev)
    n_channels = vol_prev.shape[1]
except Exception as e:
    st.warning(f"No se pudo precargar la imagen: {e}")

gcol1, gcol2, gcol3 = st.columns(3)
with gcol1:
    gfap_default = int(glob_calib.get("GFAP_CHANNEL_INDEX", 1 if n_channels > 1 else 0))
    gfap_default = min(max(0, gfap_default), max(0, n_channels-1))
    gfap_idx = st.number_input("Índice de canal GFAP", value=int(gfap_default), min_value=0, max_value=max(0, n_channels-1), step=1)
    
z_um = float(glob_calib.get('z', 1.0)); y_um = float(glob_calib.get('y', 0.3)); x_um = float(glob_calib.get('x', 0.3))
def_iso = float(min(z_um, y_um, x_um))

def _get_param(name: str, fallback):
    return glob_calib.get(name, fallback)

# Defaults
default_target_iso_um = float(_get_param("SKELETON_TARGET_ISO_UM", def_iso))
default_padding_um = float(_get_param("SKELETON_PADDING_UM", 2.0))
default_seed_dilate_um = float(_get_param("SKELETON_SEED_DILATE_UM", 2.0))
default_connectivity = int(_get_param("SKELETON_CONNECTIVITY", 26))
default_closing_um = float(_get_param("SKELETON_CLOSING_UM", 0.8))
default_territory_voronoi = bool(_get_param("SKELETON_TERRITORY_VORONOI", False))
default_territory_excl_um = float(_get_param("SKELETON_TERRITORY_EXCLUSION_UM", 1.0))
default_domain_source = str(_get_param("SKELETON_DOMAIN_VOLUME_SOURCE", "gfap")).lower()
default_prune_enable = bool(_get_param("SKELETON_PRUNE_ENABLE", False))
default_prune_min_len_um = float(_get_param("SKELETON_PRUNE_MIN_LEN_UM", 2.0))
default_tube_radius_um = float(_get_param("SKELETON_TUBE_RADIUS_UM", 1.5))
default_max_radius_um = float(_get_param("SKELETON_MAX_RADIUS_UM", 50.0))

with gcol2:
    target_iso_um = st.number_input("Vóxel isotrópico objetivo (µm)", value=float(default_target_iso_um), min_value=0.05, step=0.01, format="%.2f")
with gcol3:
    padding_um = st.number_input("Padding ROI (µm)", value=float(default_padding_um), min_value=0.0, step=0.5)

mcol1, mcol2, mcol3 = st.columns(3)
thresh_options = ["Otsu (ROI)", "Manual"]
_def_idx = 1 if str(_get_param("SKELETON_THRESHOLD_MODE", "otsu")).lower().startswith("manual") else 0
with mcol1:
    thresh_mode = st.selectbox("Umbral GFAP", options=thresh_options, index=_def_idx)
with mcol2:
    manual_thr = st.number_input("Umbral manual", value=float(_get_param("SKELETON_MANUAL_THRESHOLD", 50)), min_value=0.0, step=1.0, format="%.0f")
with mcol3:
    closing_um = st.number_input("Cierre morfológico (µm)", value=float(default_closing_um), min_value=0.0, step=0.2, format="%.2f")

seedcol1, seedcol2 = st.columns(2)
with seedcol1:
    seed_dilate_um = st.number_input("Dilatación de semilla (µm)", value=float(default_seed_dilate_um), min_value=0.0, step=0.5)
with seedcol2:
    connectivity = st.selectbox("Conectividad 3D", options=[6, 26], index=(1 if default_connectivity == 26 else 0))

radiuscol = st.columns(1)[0]
with radiuscol:
    max_radius_um = st.number_input("Radio máximo desde núcleo (µm)", value=float(round(default_max_radius_um, 2)), min_value=0.0, step=0.5)

tcol1, tcol2 = st.columns(2)
with tcol1:
    territory_voronoi = st.checkbox("Territorios por proximidad (Voronoi)", value=bool(default_territory_voronoi))
with tcol2:
    territory_excl_um = st.number_input("Zona de exclusión en frontera (µm)", value=float(default_territory_excl_um), min_value=0.0, step=0.2, format="%.2f")

vcol1, vcol2 = st.columns(2)
with vcol1:
    domain_volume_source = st.selectbox(
        "Volumen de dominio por célula",
        options=["GFAP conectado", "Territorio Voronoi"],
        index=(0 if default_domain_source.startswith("gfap") else 1),
    )
with vcol2:
    with st.expander("Opciones avanzadas: pruning topológico"):
        prune_enable = st.checkbox("Pruning: eliminar espículas cortas", value=bool(default_prune_enable))
        prune_min_len_um = st.number_input(
            "Longitud mínima de espícula (µm)",
            value=float(default_prune_min_len_um), min_value=0.0, step=0.5
        )

tube_radius_um = st.number_input(
    "Radio del tubo alrededor del esqueleto (µm)",
    value=float(default_tube_radius_um), min_value=0.0, step=0.2, format="%.2f"
)

conflict_resolve = st.checkbox(
    "Resolver solapamientos por cercanía al núcleo (µm)",
    value=bool(glob_calib.get("SKELETON_CONFLICT_RESOLVE", True)),
    help="Si dos esqueletos se superponen, se asigna el vóxel al astro más cercano al núcleo."
)

save_skel_params = st.button("💾 Guardar parámetros del skeleton (global)")
if save_skel_params:
    exp_params_path = root / "streamlit" / "calibration.json"
    exp = glob_calib.copy()
    _mode_str = "manual" if str(thresh_mode).lower().startswith("manual") else "otsu"
    exp.update({
        "GFAP_CHANNEL_INDEX": int(gfap_idx),
        "SKELETON_TARGET_ISO_UM": float(target_iso_um),
        "SKELETON_PADDING_UM": float(padding_um),
        "SKELETON_SEED_DILATE_UM": float(seed_dilate_um),
        "SKELETON_CONNECTIVITY": int(connectivity),
        "SKELETON_CLOSING_UM": float(closing_um),
        "SKELETON_MAX_RADIUS_UM": float(max_radius_um),
        "SKELETON_THRESHOLD_MODE": _mode_str,
        "SKELETON_MANUAL_THRESHOLD": float(manual_thr),
        "SKELETON_TERRITORY_VORONOI": bool(territory_voronoi),
        "SKELETON_TERRITORY_EXCLUSION_UM": float(territory_excl_um),
        "SKELETON_DOMAIN_VOLUME_SOURCE": ("gfap" if domain_volume_source.lower().startswith("gfap") else "voronoi"),
        "SKELETON_PRUNE_ENABLE": bool(prune_enable),
        "SKELETON_PRUNE_MIN_LEN_UM": float(prune_min_len_um),
        "SKELETON_TUBE_RADIUS_UM": float(tube_radius_um),
        "SKELETON_CONFLICT_RESOLVE": bool(conflict_resolve),
    })
    _save_global_calibration(exp)
    st.success(f"Parámetros del skeleton guardados en {exp_params_path.relative_to(root)}")

# --- CORRECCIÓN: Lógica de 'run_skel' ---
run_skel = st.button("🕸️ Esqueletizar y Guardar (Preparado actual)")
open_napari = st.button("👁️ Abrir en Napari (con skeleton)")

if run_skel:
    try:
        # 1. Guardar parámetros de la UI
        exp = glob_calib.copy()
        _mode_str = "manual" if str(thresh_mode).lower().startswith("manual") else "otsu"
        exp.update({
            "GFAP_CHANNEL_INDEX": int(gfap_idx), "SKELETON_TARGET_ISO_UM": float(target_iso_um),
            "SKELETON_PADDING_UM": float(padding_um), "SKELETON_SEED_DILATE_UM": float(seed_dilate_um),
            "SKELETON_CONNECTIVITY": int(connectivity), "SKELETON_CLOSING_UM": float(closing_um),
            "SKELETON_MAX_RADIUS_UM": float(max_radius_um), "SKELETON_THRESHOLD_MODE": _mode_str,
            "SKELETON_MANUAL_THRESHOLD": float(manual_thr), "SKELETON_TERRITORY_VORONOI": bool(territory_voronoi),
            "SKELETON_TERRITORY_EXCLUSION_UM": float(territory_excl_um),
            "SKELETON_DOMAIN_VOLUME_SOURCE": ("gfap" if domain_volume_source.lower().startswith("gfap") else "voronoi"),
            "SKELETON_PRUNE_ENABLE": bool(prune_enable), "SKELETON_PRUNE_MIN_LEN_UM": float(prune_min_len_um),
            "SKELETON_TUBE_RADIUS_UM": float(tube_radius_um), "SKELETON_CONFLICT_RESOLVE": bool(conflict_resolve),
        })
        _save_global_calibration(exp)
        st.info("Parámetros de UI guardados en calibration.json")
        
        # 2. Ejecutar SOLO el Paso 04
        with st.spinner(f"Ejecutando Paso 04 (Esqueletización) para {img_path.stem}..."):
            pipeline.run_advanced_skeleton_and_save(
                img_path=img_path,
                mask_path=mask_final, # input
                out_dir=out_dir,
                cal=exp,
                conflict_resolve=bool(conflict_resolve)
            )
        st.success("Esqueletización completada.")
        st.session_state["__refresh_metrics_skel"] = True
    except Exception as e:
        st.error(f"Error en esqueletización: {e}")
        st.exception(e)

if open_napari:
    try:
        z = float(glob_calib.get('z', 1.0)); y = float(glob_calib.get('y', 0.3)); x = float(glob_calib.get('x', 0.3))
        out_skel_labels = out_dir / "05_skeleton_labels.tif"
        cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
        if mask_final.exists():
            cmd += ["--final", str(mask_final)]
        if out_skel_labels.exists():
            cmd += ["--skeleton", str(out_skel_labels)]
        env = os.environ.copy()
        subprocess.Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")

# --------- 4) Métricas ---------
st.markdown("---")
st.markdown("### 📈 Métricas del esqueleto (desde `skeletons/summary.csv`)")

def _render_metrics_table_and_charts():
    summary_path = out_dir / "skeletons" / "summary.csv"
    if not summary_path.exists():
        st.info("No hay resultados de esqueletización (summary.csv) aún. Ejecutá 'Esqueletizar y Guardar'.")
        return

    try:
        df = pd.read_csv(summary_path)
    except Exception as e:
        st.error(f"No se pudo leer summary.csv: {e}")
        return

    st.dataframe(df.round(3), use_container_width=True)

    # Resumen moderno: tarjetas con métricas clave
    n_cells = int(df.shape[0])
    mlen = float(df["total_length_um"].median()) if "total_length_um" in df.columns and df["total_length_um"].notna().any() else 0.0
    m_n_branches = float(df["n_branches"].median()) if "n_branches" in df.columns and df["n_branches"].notna().any() else 0.0
    m_tortuosity = float(df["mean_tortuosity"].median()) if "mean_tortuosity" in df.columns and df["mean_tortuosity"].notna().any() else 0.0
    mvol = float(df["domain_volume_um3"].median()) if "domain_volume_um3" in df.columns and df["domain_volume_um3"].notna().any() else 0.0
    
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Astrocitos (n)", n_cells)
    c2.metric("Longitud Mediana (µm)", f"{mlen:.1f}")
    c3.metric("N° Ramas Mediano", f"{m_n_branches:.1f}")
    c4.metric("Tortuosidad Mediana", f"{m_tortuosity:.2f}")
    c5.metric("Volumen Dominio Mediano (µm³)", f"{mvol:.0f}")

    # Gráficas
    base = alt.Chart(df)
    charts = []
    if "total_length_um" in df.columns:
        charts.append(base.mark_bar().encode(x=alt.X("total_length_um:Q", bin=alt.Bin(maxbins=30), title="Longitud total (µm)"), y="count()").properties(height=200))
    if "n_branches" in df.columns:
        charts.append(base.mark_bar().encode(x=alt.X("n_branches:Q", bin=alt.Bin(maxbins=20), title="# ramas"), y="count()").properties(height=200))
    if "mean_tortuosity" in df.columns:
        charts.append(base.mark_bar().encode(x=alt.X("mean_tortuosity:Q", bin=alt.Bin(maxbins=30), title="Tortuosidad media"), y="count()").properties(height=200))
    if "domain_volume_um3" in df.columns:
        charts.append(base.mark_bar().encode(x=alt.X("domain_volume_um3:Q", bin=alt.Bin(maxbins=30), title="Volumen de dominio (µm³)"), y="count()").properties(height=200))
    
    if charts:
        st.altair_chart(alt.vconcat(*charts).resolve_scale(x='independent'), use_container_width=True)

if st.button("🔄 Actualizar métricas") or st.session_state.get("__refresh_metrics_skel", False):
    _render_metrics_table_and_charts()
    st.session_state["__refresh_metrics_skel"] = False