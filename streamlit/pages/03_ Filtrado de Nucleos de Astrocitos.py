import sys
import os
import json
import subprocess # <--- CORRECCIÓN
from pathlib import Path
import numpy as np
import streamlit as st
import tifffile
from skimage.measure import regionprops
import pandas as pd
import altair as alt
from ui.sidebar import render_sidebar
# --- Importar la lógica unificada ---
from ui import pipeline, runner

st.title("Filtrado de Núcleos de Astrocitos")
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

def existing_results(out_dir: Path) -> dict:
    return {
        "cellpose": (out_dir / "02_cellpose_mask.tif").exists(),
        "gfap_filtered": (out_dir / "03_gfap_microglia_filtered_mask.tif").exists(),
        "final_mask": (out_dir / "04_final_astrocytes_mask.tif").exists(),
        "nucleus_metrics": (out_dir / "03_nucleus_metrics.csv").exists(),
    }

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

# --------- 1) Selección de imagen ---------
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
status = existing_results(out_dir)
group = _detect_group(img_path, root)
glob_calib = _read_global_calibration()

st.markdown("### Estado de resultados guardados")
st.write({
    "cellpose_mask": status["cellpose"],
    "nucleus_metrics": status["nucleus_metrics"],
    "gfap_filtered_mask": status["gfap_filtered"],
    "final_mask": status["final_mask"],
    "output_dir": str(out_dir.relative_to(root)),
})
def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"
st.markdown(_group_badge_html(group), unsafe_allow_html=True)

# --- CORRECCIÓN: "Recalcular por ámbito" es solo para BATCH ---
with st.expander("Recalcular por ámbito (Batch)", expanded=False):
    scope = st.radio("Ámbito (Batch)", options=["Grupo", "Todos"], horizontal=True, key="p03_scope")
    
    scope_group = None
    if scope == "Grupo":
        scope_group = st.selectbox("Grupo", options=["CTL","Hipoxia"], index=0, key="p03_scope_group")
        
    if st.button("▶️ Recalcular (desde 03)", key="p03_recalc"):
        cal = _read_global_calibration()
        sc = "group" if scope == "Grupo" else "all"
        
        res = runner.run_scope(root, scope=sc, start_step="03", cal=cal, selected=None, group=scope_group, overwrite_from_step=True)
        ok = sum(1 for _, stt in res if not stt.get("error"))
        st.success(f"Listo: {ok}/{len(res)} preparados procesados desde 03.")

# --------- 2) Canales y parámetros ---------
n_channels = 1
try:
    arr_preview, axes_prev = pipeline.load_image_any(img_path)
    vol_prev = pipeline.reorder_to_zcyx(arr_preview, axes_prev)
    n_channels = vol_prev.shape[1]
except Exception as e:
    st.warning(f"No se pudo precargar la imagen: {e}")

st.markdown("### Selección de canales")
cc1, cc2 = st.columns(2)
with cc1:
    gfap_default = int(glob_calib.get("GFAP_CHANNEL_INDEX", 1 if n_channels > 1 else 0))
    gfap_default = min(max(0, gfap_default), max(0, n_channels-1))
    gfap_idx = st.number_input("Índice de canal GFAP", value=int(gfap_default), min_value=0, max_value=max(0, n_channels-1), step=1)
with cc2:
    micro_default = int(glob_calib.get("MICROGLIA_CHANNEL_INDEX", 2 if n_channels > 2 else min(1, n_channels-1)))
    micro_default = min(max(0, micro_default), max(0, n_channels-1))
    micro_idx = st.number_input("Índice de canal Microglía", value=int(micro_default), min_value=0, max_value=max(0, n_channels-1), step=1)

st.markdown("### Parámetros de filtrado (GFAP / Microglía)")
st.info("Usamos umbrales **relativos** (desviaciones estándar sobre el fondo) para mayor robustez científica.")
fcol1, fcol2, fcol3 = st.columns(3)
with fcol1:
    shell_radius_um = st.number_input("Radio del 'Shell' (µm)", value=float(glob_calib.get("SHELL_RADIUS_UM", 2.0)), min_value=0.1, step=0.1)
with fcol2:
    gfap_std_thr = st.number_input("Umbral GFAP (N° StdDev)", value=float(glob_calib.get("GFAP_STD_DEV_THRESHOLD", 3.0)), min_value=0.0, step=0.1)
with fcol3:
    micro_std_thr = st.number_input("Umbral Microglía (N° StdDev)", value=float(glob_calib.get("MICROGLIA_STD_DEV_THRESHOLD", 5.0)), min_value=0.0, step=0.1)

st.markdown("### Filtro por tamaño (volumen físico)")
min_volume_um3 = st.number_input("Volumen mínimo (µm³)", value=int(glob_calib.get("MIN_VOLUME_UM3", 75)), min_value=0, step=1)

save_experiment_params = st.button("💾 Guardar parámetros del experimento (global)")
if save_experiment_params:
    glob_calib.update({
        "GFAP_CHANNEL_INDEX": int(gfap_idx),
        "MICROGLIA_CHANNEL_INDEX": int(micro_idx),
        "SHELL_RADIUS_UM": float(shell_radius_um),
        "GFAP_STD_DEV_THRESHOLD": float(gfap_std_thr),
        "MICROGLIA_STD_DEV_THRESHOLD": float(micro_std_thr),
        "MIN_VOLUME_UM3": int(min_volume_um3),
    })
    _save_global_calibration(glob_calib)
    st.success(f"Parámetros del experimento guardados en {calib_path.relative_to(root)}")

# --------- 3) Lógica de filtrado (REFACTORIZADA) ---------
st.markdown("### Ejecutar Filtrado (Para el preparado actual)")
run_filter = st.button("🔬 Ejecutar filtrado GFAP/Microglía y Tamaño")

@st.cache_data(show_spinner="Cargando canales y máscaras...")
def load_data_for_filtering(img_path, out_dir, gfap_idx, micro_idx):
    arr, axes = pipeline.load_image_any(img_path)
    vol = pipeline.reorder_to_zcyx(arr, axes)
    gfap_channel = vol[:, int(gfap_idx), :, :]
    microglia_channel = vol[:, int(micro_idx), :, :]
    
    otsu_path = out_dir / "01_otsu_mask.tif"
    cellpose_path = out_dir / "02_cellpose_mask.tif"
    
    if not otsu_path.exists() or not cellpose_path.exists():
        missing = "Otsu (01)" if not otsu_path.exists() else "Cellpose (02)"
        st.error(f"No se encontró {missing}. Ejecutá primero la página de Otsu/Cellpose.")
        return None
        
    otsu_mask = tifffile.imread(otsu_path)
    cellpose_masks = tifffile.imread(cellpose_path)
    return cellpose_masks, gfap_channel, microglia_channel, otsu_mask

if run_filter:
    try:
        # 1. Guardar los parámetros de la UI en calibration.json
        glob_calib.update({
            "GFAP_CHANNEL_INDEX": int(gfap_idx),
            "MICROGLIA_CHANNEL_INDEX": int(micro_idx),
            "SHELL_RADIUS_UM": float(shell_radius_um),
            "GFAP_STD_DEV_THRESHOLD": float(gfap_std_thr),
            "MICROGLIA_STD_DEV_THRESHOLD": float(micro_std_thr),
            "MIN_VOLUME_UM3": int(min_volume_um3),
        })
        _save_global_calibration(glob_calib)
        st.info("Parámetros de UI guardados en calibration.json")
        
        # 2. Cargar los datos necesarios
        data = load_data_for_filtering(img_path, out_dir, gfap_idx, micro_idx)
        
        if data:
            cellpose_masks, gfap_channel, microglia_channel, otsu_mask = data
            
            # 3. Ejecutar SOLO el Paso 03 (Filtrado)
            with st.spinner("Ejecutando Paso 03 (Filtrado relativo)..."):
                gfap_filtered_mask, df_metrics = pipeline.run_filter_and_save(
                    cellpose_masks, gfap_channel, microglia_channel, 
                    otsu_mask, glob_calib, out_dir
                )
            n_kept = int((df_metrics["is_astrocyte_candidate"]).sum())
            st.success(f"Filtro GFAP/Microglía guardado. {n_kept} candidatos retenidos.")
            
            # 4. Ejecutar SOLO el Paso 04 (Tamaño)
            with st.spinner("Ejecutando Paso 04 (Filtro por tamaño)…"):
                final_mask = pipeline.run_size_filter_and_save(gfap_filtered_mask, glob_calib, out_dir)
            n_final = len(np.unique(final_mask)) - 1
            st.success(f"Filtro de tamaño guardado. {n_final} astrocitos finales.")
            
            st.session_state["__refresh_metrics"] = True # Refrescar métricas
            
    except Exception as e:
        st.error(f"Error en filtrado: {e}")
        st.exception(e)

# --------- 4) Ver en Napari ---------
open_napari_image = st.button("👁️ Abrir en Napari (solo imagen)")
open_napari_with_masks = st.button("🧪 Abrir en Napari con máscaras disponibles")

def _launch_napari(include_masks: bool):
    env = os.environ.copy()
    z = float(glob_calib.get('z', 1.0)); y = float(glob_calib.get('y', 0.3)); x = float(glob_calib.get('x', 0.3))
    cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z), "--y", str(y), "--x", str(x)]
    if include_masks:
        paths = {
            "--otsu": out_dir / "01_otsu_mask.tif",
            "--cellpose": out_dir / "02_cellpose_mask.tif",
            "--gfap": out_dir / "03_gfap_microglia_filtered_mask.tif",
            "--final": out_dir / "04_final_astrocytes_mask.tif"
        }
        for flag, p in paths.items():
            if p.exists():
                cmd += [flag, str(p)]
    try:
        subprocess.Popen(cmd, env=env)
        st.info("Napari lanzado en una ventana separada.")
    except Exception as e:
        st.error(f"No se pudo lanzar Napari: {e}")

if open_napari_image: _launch_napari(include_masks=False)
if open_napari_with_masks: _launch_napari(include_masks=True)


# --------- 5) Métricas y resumen (Actualizado) ---------
st.markdown("### 📈 Métricas del pipeline y resumen")

def render_metrics_section():
    df_metrics_path = out_dir / "03_nucleus_metrics.csv"
    if not df_metrics_path.exists():
        st.info("No hay métricas de núcleo (03_nucleus_metrics.csv). Ejecutá el filtrado primero.")
        return
        
    df_metrics = pd.read_csv(df_metrics_path)
    
    n_cellpose = df_metrics.shape[0]
    n_gfap = int(df_metrics["is_astrocyte_candidate"].sum())
    
    final_mask_path = out_dir / "04_final_astrocytes_mask.tif"
    n_final = 0
    if final_mask_path.exists():
        final_mask = tifffile.imread(final_mask_path)
        n_final = len(np.unique(final_mask)) - 1

    ret_gfap = (n_gfap / n_cellpose * 100.0) if n_cellpose else 0.0
    drop_final = n_gfap - n_final # Caída de GFAP a Final

    mcol1, mcol2, mcol3 = st.columns(3)
    mcol1.metric("Núcleos (Cellpose)", n_cellpose)
    mcol2.metric("Candidatos (GFAP/Micro)", n_gfap, delta=f"{ret_gfap:.1f}% retenidos")
    mcol3.metric("Astrocitos (Tamaño)", n_final, delta=f"−{drop_final}" if drop_final > 0 else None)

    chart_df = pd.DataFrame({
        "Etapa": ["Cellpose", "GFAP/Microglía", "Final (Volumen)"],
        "Cantidad": [n_cellpose, n_gfap, n_final]
    })
    bar = alt.Chart(chart_df).mark_bar(cornerRadiusTopLeft=3, cornerRadiusTopRight=3).encode(
        x=alt.X("Etapa:N", sort=None),
        y=alt.Y("Cantidad:Q"),
        color=alt.Color("Etapa:N", legend=None)
    ).properties(height=220)
    st.altair_chart(bar, use_container_width=True)

    st.markdown("#### Distribución de Métricas de Núcleos (Todos los núcleos)")
    
    chart_vol = alt.Chart(df_metrics).mark_bar().encode(
        x=alt.X("nucleus_volume_um3:Q", bin=alt.Bin(maxbins=40), title="Volumen del Núcleo (µm³)"),
        y=alt.Y("count()", title="Conteo"),
        color=alt.Color("is_astrocyte_candidate:N", title="Candidato a Astrocito")
    ).properties(height=200)
    st.altair_chart(chart_vol, use_container_width=True)
    
    chart_sph = alt.Chart(df_metrics.dropna(subset=['nucleus_sphericity'])).mark_bar().encode(
        x=alt.X("nucleus_sphericity:Q", bin=alt.Bin(maxbins=40), title="Esfericidad del Núcleo (0-1)"),
        y=alt.Y("count()", title="Conteo"),
        color=alt.Color("is_astrocyte_candidate:N", title="Candidato a Astrocito")
    ).properties(height=200)
    st.altair_chart(chart_sph, use_container_width=True)
    
    st.markdown("#### Espacio de Características del Filtrado")
    scatter_filt = alt.Chart(df_metrics).mark_circle(opacity=0.7).encode(
        x=alt.X("shell_gfap_mean:Q", title="GFAP Medio (Shell)"),
        y=alt.Y("shell_microglia_mean:Q", title="Microglía Media (Shell)"),
        color=alt.Color("is_astrocyte_candidate:N", title="Candidato"),
        tooltip=["label", "shell_gfap_mean", "shell_microglia_mean", "is_astrocyte_candidate"]
    ).interactive().properties(height=300)
    st.altair_chart(scatter_filt, use_container_width=True)

    st.markdown("#### Detalle por núcleo")
    st.dataframe(df_metrics.round(3), use_container_width=True)

refresh = st.button("🔄 Actualizar métricas")
if refresh or st.session_state.get("__refresh_metrics", False):
    render_metrics_section()
    st.session_state["__refresh_metrics"] = False