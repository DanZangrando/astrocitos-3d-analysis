import streamlit as st
import sys
import json
from pathlib import Path
import numpy as np
import tifffile
import subprocess
import os
from ui.sidebar import render_sidebar

# set_page_config debe definirse solo una vez en la app (en Home.py)
st.title("Calibración y Visualización de los Preparados")
render_sidebar(show_calibration=True)

# Desde streamlit/pages necesitamos subir dos niveles hasta la raíz del repo
root = Path(__file__).resolve().parents[2]
raw_dir = root / "data" / "raw"
# Guardamos calibración global del experimento dentro de la app de Streamlit
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


def read_calibration_from_lif(lif_path: Path):
    """Lee metadatos de un archivo Leica .lif usando readlif.
    Devuelve (z_um, y_um, x_um) o None si no se encuentran.
    """
    z_um = y_um = x_um = None
    try:
        import readlif
        rdr = readlif.Reader(str(lif_path))
        # Tomamos la primera serie / imagen
        try:
            series_list = rdr.getSeries() if hasattr(rdr, 'getSeries') else None
        except Exception:
            series_list = None
        series0 = series_list[0] if series_list else None
        scale = None
        if series0 is not None:
            scale = getattr(series0, 'scale', None)
            if scale is None:
                info = getattr(series0, 'info', None)
                if isinstance(info, dict):
                    scale = info.get('scale')
        if scale is None:
            # Fallback: algunas versiones exponen scale en la imagen directamente
            try:
                img0 = rdr.get_image(0)
                scale = getattr(img0, 'scale', None)
            except Exception:
                pass
        if scale is not None and len(scale) >= 3:
            # Se asume (Z, Y, X) en micrómetros
            z_um = float(scale[0]) if scale[0] is not None else None
            y_um = float(scale[1]) if scale[1] is not None else None
            x_um = float(scale[2]) if scale[2] is not None else None
    except Exception:
        pass
    return z_um, y_um, x_um


def load_overrides():
    if overrides_path.exists():
        try:
            return json.loads(overrides_path.read_text())  # formato: {"z":..,"y":..,"x":..}
        except Exception:
            return {}
    return {}


def save_override_global(z: float, y: float, x: float):
    data = {"z": float(z), "y": float(y), "x": float(x)}
    overrides_path.parent.mkdir(parents=True, exist_ok=True)
    overrides_path.write_text(json.dumps(data, indent=2))


# UI
st.subheader("1) Selección de imagen")

# Buscar imágenes en data/raw (.tif y .lif)
files = sorted([p for p in raw_dir.rglob("*.tif")] + [p for p in raw_dir.rglob("*.lif")])
if not files:
    st.warning("No se encontraron archivos .tif en data/raw.")
    st.stop()

labels = [str(p.relative_to(root)) for p in files]
idx = st.selectbox("Elegí un preparado (.tif)", options=list(range(len(files))), format_func=lambda i: labels[i])
img_path = files[idx]

st.subheader("2) Calibración física global del experimento (µm)")

# Valores detectados
suffix = img_path.suffix.lower()
if suffix == '.lif':
    z_det, y_det, x_det = read_calibration_from_lif(img_path)
else:
    z_det, y_det, x_det = read_calibration_from_tiff(img_path)
ov = load_overrides() or {}

col1, col2 = st.columns(2)
with col1:
    st.write("Valores detectados (solo lectura) en la imagen seleccionada:")
    st.json({"detected": {"z": z_det, "y": y_det, "x": x_det}})
    # Sugerencia si Y/X parecen 1.0: recomendar 0.3 µm
    if (y_det is None or abs(y_det - 1.0) < 1e-6) or (x_det is None or abs(x_det - 1.0) < 1e-6):
        st.warning("No se detectó calibración fiable en Y/X (1.0 µm). Para este experimento suele ser ~0.3 µm. Ajustá abajo y guardá.")
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
        cmd = [sys.executable, str(napari_script), "--path", str(img_path), "--z", str(z_in), "--y", str(y_in), "--x", str(x_in)]
        try:
            subprocess.Popen(cmd, env=env)
            st.info("Napari lanzado en una ventana separada.")
        except Exception as e:
            st.error(f"No se pudo lanzar Napari: {e}")

st.markdown("---")
st.caption("Tip: Si tus TIFFs no traen metadatos físicos, usá esta página para fijar Y y X (por ejemplo 0.3 µm) además de Z.")
