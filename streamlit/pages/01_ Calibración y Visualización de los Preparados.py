import streamlit as st
import sys
import json
from pathlib import Path
import numpy as np
import tifffile
import subprocess
import os
from ui.sidebar import render_sidebar
from ui.utils import detect_group

# set_page_config debe definirse solo una vez en la app (en Home.py)
st.title("Calibración y Visualización de los Preparados")
render_sidebar(show_calibration=True)

# Desde streamlit/pages necesitamos subir dos niveles hasta la raíz del repo
root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
# Archivo unificado para calibración y parámetros globales
overrides_path = root / "streamlit" / "calibration.json"
napari_script = root / "streamlit" / "napari_viewer.py"

# Utils

def _to_um(value: float, unit: str | None):
    if value is None:
        return None
    if not unit:
        return float(value)
    unit = unit.strip().lower()
    if unit in {"µm","um","micrometer","micrometre"}:
        return float(value)
    if unit in {"nm","nanometer","nanometre"}:
        return float(value) / 1000.0
    if unit in {"mm","millimeter","millimetre"}:
        return float(value) * 1000.0
    if unit in {"m","meter","metre"}:
        return float(value) * 1_000_000.0
    return float(value)


def read_calibration_from_tiff(tif_path: Path):
    z_um = y_um = x_um = None
    with tifffile.TiffFile(str(tif_path)) as tf:
        ome_xml = tf.ome_metadata
        ij_meta = getattr(tf, 'imagej_metadata', None)
        try:
            page0 = tf.pages[0]
            xres_tag = page0.tags.get('XResolution')
            yres_tag = page0.tags.get('YResolution')
            res_unit_tag = page0.tags.get('ResolutionUnit')
            XRES = float(xres_tag.value[0]) / float(xres_tag.value[1]) if xres_tag else None
            YRES = float(yres_tag.value[0]) / float(yres_tag.value[1]) if yres_tag else None
            RES_UNIT = res_unit_tag.value if res_unit_tag else None
        except Exception:
            XRES = YRES = RES_UNIT = None
    # Parse OME-XML
    if ome_xml:
        import xml.etree.ElementTree as ET
        try:
            root = ET.fromstring(ome_xml)
            ns = [
                "http://www.openmicroscopy.org/Schemas/OME/2016-06",
                "http://www.openmicroscopy.org/Schemas/OME/2015-01",
                "http://www.openmicroscopy.org/Schemas/OME/2013-06",
            ]
            pixels = None
            for n in ns:
                pixels = root.find(f".//{{{n}}}Pixels")
                if pixels is not None:
                    break
            if pixels is None:
                pixels = root.find('.//Pixels')
            if pixels is not None:
                z = pixels.attrib.get('PhysicalSizeZ')
                y = pixels.attrib.get('PhysicalSizeY')
                x = pixels.attrib.get('PhysicalSizeX')
                zu = pixels.attrib.get('PhysicalSizeZUnit') or pixels.attrib.get('PhysicalSizeZUnitSymbol')
                yu = pixels.attrib.get('PhysicalSizeYUnit') or pixels.attrib.get('PhysicalSizeYUnitSymbol')
                xu = pixels.attrib.get('PhysicalSizeXUnit') or pixels.attrib.get('PhysicalSizeXUnitSymbol')
                z_um = _to_um(float(z), zu) if z is not None else None
                y_um = _to_um(float(y), yu) if y is not None else None
                x_um = _to_um(float(x), xu) if x is not None else None
        except Exception:
            pass
    # Fallback ImageJ
    if ij_meta and z_um is None:
        try:
            unit = ij_meta.get('unit')
            spacing = ij_meta.get('spacing')
            z_um = _to_um(float(spacing), unit) if spacing is not None else None
        except Exception:
            pass
    # Fallback Resolution tags
    if (x_um is None or y_um is None) and RES_UNIT and (XRES or YRES):
        if RES_UNIT == 2:  # inch
            if XRES:
                x_um = 25400.0 / XRES
            if YRES:
                y_um = 25400.0 / YRES
        elif RES_UNIT == 3:  # centimeter
            if XRES:
                x_um = 10000.0 / XRES
            if YRES:
                y_um = 10000.0 / YRES
    return z_um, y_um, x_um


def load_overrides():
    if overrides_path.exists():
        try:
            return json.loads(overrides_path.read_text())  # Unificado: {"z","y","x", ...otros parámetros}
        except Exception:
            return {}
    return {}


def save_override_global(z: float, y: float, x: float):
    # Mantener cualquier otro parámetro ya guardado
    cur = load_overrides() or {}
    cur.update({"z": float(z), "y": float(y), "x": float(x)})
    overrides_path.parent.mkdir(parents=True, exist_ok=True)
    overrides_path.write_text(json.dumps(cur, indent=2))


# UI
st.subheader("1) Selección de imagen")

all_tifs = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.tiff")])
all_files = all_tifs
if not all_files:
    st.warning("No se encontraron archivos .tif/.tiff en data/raw.")
    st.stop()

groups = [detect_group(p, root) for p in all_files]
g_ctl = sum(1 for g in groups if g == "CTL")
g_hip = sum(1 for g in groups if g == "Hipoxia")

cc1, cc2 = st.columns(2)
cc1.metric("Total de preparados (.tif/.tiff)", len(all_files))
cc2.metric("CTL / Hipoxia", f"{g_ctl} / {g_hip}", help="Conteo por grupo principal")

# Filtro por grupo
# Filtro por grupo (unificado desde la barra lateral)
group_filter = st.session_state.get("group_filter", "Todos")
if group_filter == "Todos":
    files_avail = all_files
else:
    files_avail = [p for p in all_files if detect_group(p, root) == group_filter]
if not files_avail:
    st.info("No hay preparados para el grupo seleccionado.")
    st.stop()

# Buscar imágenes en data/raw (.tif/.tiff)
labels = [str(p.relative_to(root)) for p in files_avail]
idx = st.selectbox("Elegí un preparado", options=list(range(len(files_avail))), format_func=lambda i: labels[i])
img_path = files_avail[idx]
sel_group = detect_group(img_path, root)
def _group_badge_html(group: str) -> str:
    color = {"CTL": "#1f77b4", "Hipoxia": "#d62728"}.get(group, "#7f7f7f")
    return f"<span style='background:{color};color:white;padding:3px 8px;border-radius:999px;font-weight:600;font-size:0.85rem;'>{group}</span>"
st.markdown(_group_badge_html(sel_group), unsafe_allow_html=True)

st.subheader("2) Calibración física global del experimento (µm)")

# Valores detectados
z_det, y_det, x_det = read_calibration_from_tiff(img_path)
ov = load_overrides() or {}

col1, col2 = st.columns(2)
with col1:
    st.write("Valores detectados (solo lectura) en la imagen seleccionada:")
    st.json({"detected": {"z": z_det, "y": y_det, "x": x_det}})
    # Indicadores rápidos
    mc1, mc2, mc3 = st.columns(3)
    mc1.metric("Z detectado", f"{z_det:.4f} µm" if z_det else "—")
    mc2.metric("Y detectado", f"{y_det:.4f} µm" if y_det else "—")
    mc3.metric("X detectado", f"{x_det:.4f} µm" if x_det else "—")
    # Sugerencia si Y/X parecen 1.0: recomendar 0.3 µm
    if (y_det is None or (isinstance(y_det, (int,float)) and abs(y_det - 1.0) < 1e-6)) or (x_det is None or (isinstance(x_det, (int,float)) and abs(x_det - 1.0) < 1e-6)):
        st.warning("No se detectó calibración fiable en Y/X (≈1.0 µm). Para este experimento suele ser ~0.3 µm. Ajustá abajo y guardá.")
with col2:
    st.write("Ingresá/ajustá la calibración global a usar en todo el experimento:")
    z_in = st.number_input("Z (µm)", value=float(ov.get('z', z_det or 1.0)), min_value=0.0001, step=0.01, format="%0.4f")
    y_in = st.number_input("Y (µm)", value=float(ov.get('y', y_det or 0.3)), min_value=0.0001, step=0.01, format="%0.4f")
    x_in = st.number_input("X (µm)", value=float(ov.get('x', x_det or 0.3)), min_value=0.0001, step=0.01, format="%0.4f")

save_col, open_col = st.columns(2)
with save_col:
    if st.button("💾 Guardar calibración global"):
        save_override_global(z_in, y_in, x_in)
        st.success(f"Calibración global guardada en {overrides_path.relative_to(root)}.")

with open_col:
    if st.button("🧪 Abrir en Napari con esta calibración"):
        # Abrimos Napari en un proceso separado para no bloquear Streamlit
        env = os.environ.copy()
        # Obtener índices de canales desde calibración
        cal = json.loads(overrides_path.read_text()) if overrides_path.exists() else {}
        dapi_idx = int(cal.get("DAPI_CHANNEL_INDEX", 0))
        gfap_idx = int(cal.get("GFAP_CHANNEL_INDEX", 1))
        cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z_in), "--y", str(y_in), "--x", str(x_in),
               "--dapi_idx", str(dapi_idx), "--gfap_idx", str(gfap_idx)]
        try:
            env["NAPARI_DISABLE_PLUGIN_AUTOLOAD"] = "1"
            subprocess.Popen(cmd, env=env)
            st.info("Napari lanzado en una ventana separada.")
        except Exception as e:
            st.error(f"No se pudo lanzar Napari: {e}")

st.markdown("---")
st.caption("Tip: Si tus archivos no traen metadatos físicos confiables, usá esta página para fijar Y y X (por ejemplo 0.3 µm) además de Z. Los indicadores arriba te muestran si los metadatos fueron detectados correctamente.")