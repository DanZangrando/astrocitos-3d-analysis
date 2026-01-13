"""
Pipeline 2D Unificado para Análisis de Astrocitos

Este módulo implementa un flujo de trabajo completamente 2D usando funcionalidades nativas de SKAN:
1. Proyección de máscaras nucleares y señal GFAP (max projection)
2. Partición territorial usando diagramas de Voronoi con zona de exclusión
3. Esqueletización 2D por territorio con conexión opcional de fragmentos
4. Análisis de Sholl 2D usando skan.sholl_analysis() (cuenta cruces de branches)
5. Visualización con skan.draw.sholl_shells()

Ventajas sobre el enfoque 3D y implementaciones custom:
- Mayor resolución espacial (mantiene 0.38 µm en XY vs 1 µm isotrópico)
- Sholl nativo de SKAN: cuenta cruces de BRANCHES (método estándar en neurociencia)
- Voronoi 2D simple y robusto
- Territorios astrocitarios bien definidos para preparados planos
- Máxima robustez: usa implementaciones validadas de SKAN (no reinventamos la rueda)
- Coherencia con literatura: sholl_analysis() implementa el método clásico de Sholl (1953)

Nota técnica sobre Sholl:
- SKAN cuenta cruces de BRANCHES/PATHS en cada capa concéntrica
- No cuenta píxeles individuales (eso inflaría artificialmente el conteo)
- Un branch puede cruzar múltiples capas → múltiples intersecciones
- Método estándar: mide complejidad de ramificación real
"""

from __future__ import annotations
from pathlib import Path
from typing import Tuple
import json

import numpy as np
import pandas as pd
import tifffile
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.ndimage import distance_transform_edt, binary_dilation, generate_binary_structure
from skimage.morphology import skeletonize, disk, binary_closing
from skimage.filters import threshold_otsu
from skimage.measure import regionprops, label as cc_label
import skan
from skan import Skeleton, summarize, sholl_analysis, draw


def project_mask_labels_2d(mask_3d: np.ndarray, method: str = 'max') -> np.ndarray:
    """
    Proyecta máscara 3D de labels a 2D preservando etiquetas.
    
    Para cada posición (y,x), toma el label del slice donde ese label
    tiene mayor presencia según el método (max o sum de voxels).
    """
    if mask_3d.ndim != 3:
        raise ValueError(f"Esperaba volumen 3D, recibido shape={mask_3d.shape}")
    
    labels = np.unique(mask_3d)
    labels = labels[labels > 0]
    
    result = np.zeros(mask_3d.shape[1:], dtype=mask_3d.dtype)
    
    for lab in labels:
        binary_3d = (mask_3d == lab)
        if method == 'max':
            # Presencia en cualquier slice
            projection = np.any(binary_3d, axis=0)
        elif method == 'sum':
            # Suma de voxels por columna
            projection = np.sum(binary_3d, axis=0) > 0
        else:
            projection = np.max(binary_3d, axis=0)
        
        result[projection] = lab
    
    return result


def project_3d_to_2d(vol_3d: np.ndarray, method: str = 'max') -> np.ndarray:
    """Proyecta volumen 3D de intensidad a 2D."""
    if vol_3d.ndim != 3:
        raise ValueError(f"Esperaba volumen 3D, recibido shape={vol_3d.shape}")
    
    if method == 'max':
        return np.max(vol_3d, axis=0)
    elif method == 'mean':
        return np.mean(vol_3d, axis=0)
    elif method == 'sum':
        return np.sum(vol_3d, axis=0)
    else:
        return np.max(vol_3d, axis=0)


def create_circular_mask_from_nucleus(
    nuclear_mask_2d: np.ndarray,
    max_radius_um: float,
    spacing: Tuple[float, float]
) -> np.ndarray:
    """
    Crea una máscara circular de radio máximo desde el centroide nuclear.
    
    Args:
        nuclear_mask_2d: Máscara binaria del núcleo en 2D
        max_radius_um: Radio máximo en µm
        spacing: (y_um, x_um) espaciado de píxeles
        
    Returns:
        Máscara binaria 2D con círculo de radio máximo desde núcleo
    """
    if not np.any(nuclear_mask_2d):
        return np.zeros_like(nuclear_mask_2d, dtype=bool)
    
    y_um, x_um = spacing
    
    # Calcular centroide del núcleo
    props = regionprops(nuclear_mask_2d.astype(np.uint8))
    if len(props) == 0:
        return np.zeros_like(nuclear_mask_2d, dtype=bool)
    
    centroid_y, centroid_x = props[0].centroid
    
    # Crear grid de coordenadas
    Y, X = np.ogrid[:nuclear_mask_2d.shape[0], :nuclear_mask_2d.shape[1]]
    
    # Calcular distancia euclidiana en µm desde centroide
    dist_y_um = (Y - centroid_y) * y_um
    dist_x_um = (X - centroid_x) * x_um
    dist_um = np.sqrt(dist_y_um**2 + dist_x_um**2)
    
    # Máscara circular
    circular_mask = (dist_um <= max_radius_um)
    
    return circular_mask


def filter_skeleton_by_nucleus_connectivity(
    skeleton: np.ndarray,
    nuclear_mask: np.ndarray
) -> np.ndarray:
    """
    Filtra el esqueleto para quedarse SOLO con la componente conectada que toca el núcleo.
    
    Esto elimina fragmentos aislados dentro del territorio que no están conectados
    al soma del astrocito, incluso después del intento de conexión de fragmentos.
    
    Args:
        skeleton: Esqueleto binario (puede tener múltiples componentes)
        nuclear_mask: Máscara del núcleo
        
    Returns:
        Esqueleto filtrado con solo la componente conectada al núcleo
    """
    if not np.any(skeleton):
        return skeleton
    
    # Etiquetar componentes conectadas del esqueleto
    labeled_skeleton, n_components = cc_label(skeleton, return_num=True, connectivity=2)
    
    if n_components == 1:
        # Solo hay una componente, no hay nada que filtrar
        return skeleton
    
    # Encontrar qué componente(s) tocan el núcleo
    nuclear_overlap = labeled_skeleton * nuclear_mask
    components_touching_nucleus = np.unique(nuclear_overlap)
    components_touching_nucleus = components_touching_nucleus[components_touching_nucleus > 0]
    
    if len(components_touching_nucleus) == 0:
        # Ninguna componente toca el núcleo directamente
        # Buscar la componente más cercana al núcleo
        nucleus_coords = np.argwhere(nuclear_mask)
        if len(nucleus_coords) == 0:
            return skeleton
        
        nucleus_centroid = nucleus_coords.mean(axis=0)
        
        min_dist = np.inf
        closest_component = 1
        
        for comp_id in range(1, n_components + 1):
            comp_coords = np.argwhere(labeled_skeleton == comp_id)
            if len(comp_coords) == 0:
                continue
            
            # Distancia mínima entre esta componente y el centroide nuclear
            dists = np.linalg.norm(comp_coords - nucleus_centroid, axis=1)
            min_comp_dist = np.min(dists)
            
            if min_comp_dist < min_dist:
                min_dist = min_comp_dist
                closest_component = comp_id
        
        # Quedarse solo con la componente más cercana
        filtered_skeleton = (labeled_skeleton == closest_component)
        
        n_removed = n_components - 1
        print(f"      → Removidas {n_removed} componentes desconectadas (no tocaban núcleo)")
        
        return filtered_skeleton
    
    # Si hay componentes tocando el núcleo, quedarse con la(s) que lo tocan
    filtered_skeleton = np.isin(labeled_skeleton, components_touching_nucleus)
    
    n_removed = n_components - len(components_touching_nucleus)
    if n_removed > 0:
        print(f"      → Removidas {n_removed} componentes desconectadas del núcleo")
    
    return filtered_skeleton


def compute_voronoi_territories_2d(
    mask_2d: np.ndarray,
    exclusion_gap_um: float,
    spacing: Tuple[float, float]
) -> dict:
    """
    Calcula territorios Voronoi para cada célula con zona de exclusión.
    
    Args:
        mask_2d: Máscara 2D con labels nucleares
        exclusion_gap_um: Distancia de exclusión entre territorios (µm)
        spacing: (y_um, x_um) espaciado de píxeles
        
    Returns:
        Dict {label: territory_mask_2d}
    """
    y_um, x_um = spacing
    exclusion_gap_px = int(exclusion_gap_um / max(y_um, x_um))
    
    labels = np.unique(mask_2d)
    labels = labels[labels > 0]
    
    if len(labels) == 0:
        return {}
    
    # Calcular centroides
    props = regionprops(mask_2d.astype(np.int32))
    centroids = {p.label: p.centroid for p in props}
    
    # Si solo hay una célula, todo es su territorio
    if len(labels) == 1:
        return {labels[0]: np.ones_like(mask_2d, dtype=bool)}
    
    # Crear mapa de distancias al centroide más cercano
    dist_map = np.full(mask_2d.shape, np.inf, dtype=np.float32)
    label_map = np.zeros(mask_2d.shape, dtype=np.uint16)
    
    for lab, (cy, cx) in centroids.items():
        # Crear grilla de coordenadas
        yy, xx = np.mgrid[0:mask_2d.shape[0], 0:mask_2d.shape[1]]
        
        # Distancia euclidiana ponderada por espaciado
        distances = np.sqrt(((yy - cy) * y_um)**2 + ((xx - cx) * x_um)**2)
        
        # Actualizar donde esta célula es más cercana
        closer = distances < dist_map
        dist_map[closer] = distances[closer]
        label_map[closer] = lab
    
    # Crear zona de exclusión en bordes entre territorios
    territories = {}
    for lab in labels:
        territory = (label_map == lab)
        
        if exclusion_gap_px > 0:
            # Encontrar bordes con otros territorios
            from scipy.ndimage import binary_erosion
            eroded = binary_erosion(territory, iterations=exclusion_gap_px)
            territories[lab] = eroded
        else:
            territories[lab] = territory
    
    return territories


def connect_nearby_fragments(
    binary_mask: np.ndarray,
    radius_um: float,
    spacing: Tuple[float, float]
) -> np.ndarray:
    """
    Conecta fragmentos binarios cercanos mediante dilatación controlada.
    
    Útil para unir procesos astrocitarios fragmentados en la señal de GFAP.
    """
    y_um, x_um = spacing
    radius_px = int(radius_um / max(y_um, x_um))
    
    if radius_px < 1:
        return binary_mask
    
    # Dilatación + cierre
    selem = disk(radius_px)
    dilated = binary_dilation(binary_mask, selem)
    closed = binary_closing(dilated, selem)
    
    return closed


def run_unified_2d_skeleton_and_sholl(
    mask_3d: np.ndarray,
    gfap_3d: np.ndarray,
    cal: dict,
    out_dir: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Pipeline 2D unificado: proyección → Voronoi → esqueleto → Sholl.
    
    Returns:
        (skeleton_metrics_df, sholl_df)
    """
    print("\n" + "="*70)
    print("PIPELINE 2D UNIFICADO: ESQUELETIZACIÓN + SHOLL")
    print("="*70)
    
    # Parámetros
    y_um, x_um = float(cal.get('y', 0.3)), float(cal.get('x', 0.3))
    spacing = (y_um, x_um)
    
    proj_method = str(cal.get("PROJECTION_2D_METHOD", "max")).lower()
    exclusion_gap_um = float(cal.get("TERRITORY_EXCLUSION_UM", 1.0))
    connect_fragments = bool(cal.get("CONNECT_SKELETON_FRAGMENTS", True))
    connection_radius_um = float(cal.get("CONNECTION_RADIUS_UM", 0.5))
    
    # Parámetro de radio máximo (NUEVO)
    max_radius_um = float(cal.get("MAX_RADIUS_FROM_NUCLEUS_UM", 100.0))
    
    # Sholl parameters
    sholl_min_r = float(cal.get("SHOLL_MIN_RADIUS_UM", 5.0))
    sholl_max_r = float(cal.get("SHOLL_MAX_RADIUS_UM", 100.0))
    sholl_step = float(cal.get("SHOLL_STEP_UM", 2.0))
    
    print(f"\n📊 Parámetros:")
    print(f"  - Proyección: {proj_method}")
    print(f"  - Resolución XY: {y_um:.3f} × {x_um:.3f} µm")
    print(f"  - Gap Voronoi: {exclusion_gap_um:.2f} µm")
    print(f"  - Radio máximo desde núcleo: {max_radius_um:.1f} µm")
    print(f"  - Conectar fragmentos: {connect_fragments}")
    if connect_fragments:
        print(f"  - Radio de conexión: {connection_radius_um:.2f} µm")
    print(f"  - Sholl: {sholl_min_r:.1f} - {sholl_max_r:.1f} µm (step={sholl_step:.1f})")
    
    # PASO 1: Proyecciones
    print("\n[1/5] Proyectando volúmenes 3D → 2D...")
    mask_2d = project_mask_labels_2d(mask_3d, method=proj_method)
    gfap_2d = project_3d_to_2d(gfap_3d, method=proj_method)
    
    # Guardar proyecciones
    tifffile.imwrite(out_dir / "04_final_astrocytes_mask_2d.tif", mask_2d.astype(np.uint16))
    tifffile.imwrite(out_dir / "gfap_projection_2d.tif", gfap_2d.astype(np.float32))
    
    labels = np.unique(mask_2d)
    labels = labels[labels > 0]
    print(f"  ✓ {len(labels)} astrocitos detectados")
    
    # PASO 2: Territorios Voronoi
    print("\n[2/5] Calculando territorios Voronoi con exclusión...")
    territories = compute_voronoi_territories_2d(mask_2d, exclusion_gap_um, spacing)
    print(f"  ✓ {len(territories)} territorios definidos")
    
    # PASO 3: Esqueletización por territorio
    print("\n[3/5] Esqueletizando GFAP por territorio...")
    skeleton_combined = np.zeros_like(mask_2d, dtype=np.uint16)
    skeleton_metrics = []
    
    skel_dir = out_dir / "skeletons"
    skel_dir.mkdir(parents=True, exist_ok=True)
    
    for label in labels:
        try:
            # 3a. Extraer territorio y núcleo
            territory_mask = territories.get(label, mask_2d == label)
            nuclear_mask = (mask_2d == label)
            
            if not np.any(territory_mask):
                print(f"  ⚠ Label {label}: territorio vacío")
                continue
            
            # 3b. Crear máscara de radio máximo desde núcleo
            circular_mask = create_circular_mask_from_nucleus(nuclear_mask, max_radius_um, spacing)
            
            # Combinar: territorio Voronoi ∩ círculo de radio máximo
            domain_mask = territory_mask & circular_mask
            
            # Debug: mostrar reducción de dominio
            territory_area_px = np.sum(territory_mask)
            domain_area_px = np.sum(domain_mask)
            reduction_pct = 100.0 * (1 - domain_area_px / territory_area_px) if territory_area_px > 0 else 0
            
            if reduction_pct > 5:  # Solo mostrar si hay reducción significativa
                print(f"    ↳ Filtro de radio máximo: {territory_area_px} → {domain_area_px} px ({reduction_pct:.1f}% reducción)")
            
            if not np.any(domain_mask):
                print(f"  ⚠ Label {label}: dominio vacío tras aplicar radio máximo")
                continue
            
            # 3c. GFAP dentro del dominio
            gfap_domain = gfap_2d * domain_mask
            gfap_region = gfap_domain[domain_mask]
            
            if len(gfap_region) == 0 or np.std(gfap_region) == 0:
                print(f"  ⚠ Label {label}: sin señal GFAP en dominio")
                continue
            
            # 3d. Umbralización (Otsu local)
            try:
                threshold_val = threshold_otsu(gfap_region)
            except:
                threshold_val = np.percentile(gfap_region, 75)
            
            gfap_binary = (gfap_domain > threshold_val) & domain_mask
            
            # 3e. Conexión de fragmentos (opcional)
            if connect_fragments:
                gfap_binary = connect_nearby_fragments(gfap_binary, connection_radius_um, spacing)
                gfap_binary = gfap_binary & domain_mask  # Re-aplicar máscara de dominio
            
            # 3f. Esqueletización
            skeleton = skeletonize(gfap_binary)
            
            # 3g. FILTRAR COMPONENTES DESCONECTADAS DEL NÚCLEO
            # Esto elimina fragmentos aislados que no están conectados al soma
            skeleton = filter_skeleton_by_nucleus_connectivity(skeleton, nuclear_mask)
            
            if not np.any(skeleton):
                print(f"  ⚠ Label {label}: esqueleto vacío tras filtrar componentes desconectadas")
                continue
            
            skeleton_combined[skeleton] = label
            
            # Métricas básicas
            n_skeleton_px = np.sum(skeleton)
            skeleton_length_um = n_skeleton_px * np.sqrt(y_um * x_um)  # Aproximación
            
            # Propiedades nucleares
            props = regionprops((nuclear_mask).astype(np.int32))
            if len(props) > 0:
                p = props[0]
                nuclear_area_um2 = p.area * y_um * x_um
                centroid_y, centroid_x = p.centroid
            else:
                nuclear_area_um2 = 0
                centroid_y, centroid_x = 0, 0
            
            skeleton_metrics.append({
                'label': int(label),
                'skeleton_length_um': float(skeleton_length_um),
                'skeleton_pixels': int(n_skeleton_px),
                'nuclear_area_um2': float(nuclear_area_um2),
                'centroid_y_um': float(centroid_y * y_um),
                'centroid_x_um': float(centroid_x * x_um),
                'territory_area_um2': float(np.sum(territory_mask) * y_um * x_um),
                'gfap_threshold': float(threshold_val)
            })
            
        except Exception as e:
            print(f"  ✗ Label {label}: error - {e}")
            continue
    
    # Guardar esqueleto combinado
    out_skeleton_2d = out_dir / "05_skeleton_labels_2d.tif"
    tifffile.imwrite(out_skeleton_2d, skeleton_combined.astype(np.uint16))
    print(f"  ✓ Esqueletos guardados: {out_skeleton_2d.name}")
    
    # PASO 4: Análisis SKAN nativo
    print("\n[4/5] Analizando topología con SKAN...")
    
    for i, row in enumerate(skeleton_metrics):
        label = row['label']
        skel_binary = (skeleton_combined == label)
        
        if not np.any(skel_binary):
            continue
        
        try:
            # Crear objeto Skeleton de SKAN
            skel_obj = Skeleton(skel_binary, spacing=spacing)
            
            # Métricas topológicas
            branch_data = summarize(skel_obj)
            
            if not branch_data.empty:
                n_branches = len(branch_data)
                total_branch_length = branch_data['branch-distance'].sum()
                
                # Robust counting using degrees from the graph
                # degrees: degree of each node index
                degrees = skel_obj.degrees
                n_endpoints = int(np.sum(degrees == 1))
                n_junctions = int(np.sum(degrees > 2))
                
                # NUEVAS MÉTRICAS DE TORTUOSIDAD
                tortuosity = branch_data['branch-distance'] / branch_data['euclidean-distance'].replace(0, np.nan)
                tortuosity_mean = float(tortuosity.mean()) if not tortuosity.isna().all() else 1.0
                tortuosity_max = float(tortuosity.max()) if not tortuosity.isna().all() else 1.0
                tortuosity_std = float(tortuosity.std()) if not tortuosity.isna().all() else 0.0
                
                # MÉTRICAS DE COMPLEJIDAD RAMIFICACIONAL
                ramification_index = float(n_branches / max(n_junctions, 1))
                termination_index = float(n_endpoints / max(n_junctions, 1)) if n_junctions > 0 else np.nan
                
                # DISTRIBUCIÓN DE LONGITUDES
                branch_length_median = float(branch_data['branch-distance'].median())
                branch_length_p25 = float(branch_data['branch-distance'].quantile(0.25))
                branch_length_p75 = float(branch_data['branch-distance'].quantile(0.75))
                branch_length_std = float(branch_data['branch-distance'].std())
                branch_length_cv = float(branch_length_std / branch_data['branch-distance'].mean()) if branch_data['branch-distance'].mean() > 0 else 0.0
                
                # FRAGMENTACIÓN
                n_components = int(branch_data['skeleton-id'].nunique())
                
                skeleton_metrics[i].update({
                    'n_branches': int(n_branches),
                    'total_branch_length_um': float(total_branch_length),
                    'n_endpoints': int(n_endpoints),
                    'n_junctions': int(n_junctions),
                    'mean_branch_length_um': float(total_branch_length / n_branches) if n_branches > 0 else 0.0,
                    # Nuevas métricas
                    'tortuosity_mean': tortuosity_mean,
                    'tortuosity_max': tortuosity_max,
                    'tortuosity_std': tortuosity_std,
                    'ramification_index': ramification_index,
                    'termination_index': termination_index,
                    'branch_length_median_um': branch_length_median,
                    'branch_length_p25_um': branch_length_p25,
                    'branch_length_p75_um': branch_length_p75,
                    'branch_length_std_um': branch_length_std,
                    'branch_length_cv': branch_length_cv,
                    'n_connected_components': n_components
                })
            else:
                skeleton_metrics[i].update({
                    'n_branches': 0,
                    'total_branch_length_um': 0.0,
                    'n_endpoints': 0,
                    'n_junctions': 0,
                    'mean_branch_length_um': 0.0,
                    'tortuosity_mean': 1.0,
                    'tortuosity_max': 1.0,
                    'tortuosity_std': 0.0,
                    'ramification_index': 0.0,
                    'termination_index': np.nan,
                    'branch_length_median_um': 0.0,
                    'branch_length_p25_um': 0.0,
                    'branch_length_p75_um': 0.0,
                    'branch_length_std_um': 0.0,
                    'branch_length_cv': 0.0,
                    'n_connected_components': 0
                })
        except Exception as e:
            print(f"  ⚠ SKAN falló para label {label}: {e}")
            skeleton_metrics[i].update({
                'n_branches': 0,
                'total_branch_length_um': 0.0,
                'n_endpoints': 0,
                'n_junctions': 0,
                'mean_branch_length_um': 0.0,
                'tortuosity_mean': 1.0,
                'tortuosity_max': 1.0,
                'tortuosity_std': 0.0,
                'ramification_index': 0.0,
                'termination_index': np.nan,
                'branch_length_median_um': 0.0,
                'branch_length_p25_um': 0.0,
                'branch_length_p75_um': 0.0,
                'branch_length_std_um': 0.0,
                'branch_length_cv': 0.0,
                'n_connected_components': 0
            })
    
    df_skeleton = pd.DataFrame(skeleton_metrics)
    df_skeleton.to_csv(skel_dir / "summary.csv", index=False)
    print(f"  ✓ Métricas guardadas: skeletons/summary.csv")
    
    # PASO 5: Análisis de Sholl 2D usando SKAN nativo
    print("\n[5/5] Calculando análisis de Sholl 2D con SKAN nativo...")
    print("  ℹ️  SKAN sholl_analysis() cuenta cruces de branches (método estándar)")
    sholl_results = []
    radii = np.arange(sholl_min_r, sholl_max_r + sholl_step, sholl_step)
    
    for label in labels:
        skel_binary = (skeleton_combined == label)
        nuclear_mask = (mask_2d == label)
        
        if not np.any(skel_binary) or not np.any(nuclear_mask):
            continue
        
        # Centroide nuclear
        props = regionprops(nuclear_mask.astype(np.int32))
        if len(props) == 0:
            continue
        
        centroid_y_px, centroid_x_px = props[0].centroid
        
        try:
            # Crear objeto Skeleton para Sholl
            skel_obj = Skeleton(skel_binary, spacing=spacing)
            
            # SKAN sholl_analysis():
            # - Centro en UNIDADES FÍSICAS (µm), no píxeles
            # - Radios en UNIDADES FÍSICAS (µm)
            # - Cuenta cruces de BRANCHES (paths), no píxeles individuales
            
            # Convertir centro de píxeles a µm
            center_um = np.array([centroid_y_px * y_um, centroid_x_px * x_um])
            
            # Radios ya están en µm
            _, _, counts = sholl_analysis(skel_obj, center=center_um, shells=radii)
            
            # Guardar resultados
            for radius_um, count in zip(radii, counts):
                sholl_results.append({
                    'label': int(label),
                    'radius_um': float(radius_um),
                    'intersections': int(count)
                })
                
        except Exception as e:
            print(f"  ⚠ Sholl SKAN falló para label {label}: {e}")
            # Sin fallback - si SKAN falla, hay un problema con el esqueleto
            continue
    
    df_sholl = pd.DataFrame(sholl_results)
    
    # Guardar curvas completas
    df_sholl.to_csv(out_dir / "sholl_2d.csv", index=False)
    df_sholl.to_csv(out_dir / "sholl_2d_native.csv", index=False)
    print(f"  ✓ Resultados Sholl guardados: sholl_2d.csv (+ sholl_2d_native.csv para compatibilidad)")
    
    # GENERAR SHOLL_SUMMARY.CSV con métricas escalares (AUC, peak, critical_radius)
    # Esto es CRÍTICO para la página 06 de comparación entre grupos
    if not df_sholl.empty:
        from scipy.integrate import trapezoid
        
        summary_records = []
        for label in df_sholl['label'].unique():
            label_data = df_sholl[df_sholl['label'] == label].sort_values('radius_um')
            radii = label_data['radius_um'].to_numpy()
            intersections = label_data['intersections'].to_numpy()
            
            if len(radii) < 2:
                # Sin suficientes puntos para calcular métricas
                summary_records.append({
                    'label': int(label),
                    'critical_radius_um': np.nan,
                    'peak_intersections': np.nan,
                    'auc': np.nan
                })
                continue
            
            # Calcular métricas escalares
            peak_idx = np.argmax(intersections)
            critical_radius_um = radii[peak_idx]
            peak_intersections = intersections[peak_idx]
            auc = trapezoid(intersections, radii)
            
            summary_records.append({
                'label': int(label),
                'critical_radius_um': float(critical_radius_um),
                'peak_intersections': float(peak_intersections),
                'auc': float(auc)
            })
        
        df_sholl_summary = pd.DataFrame(summary_records)
        df_sholl_summary.to_csv(out_dir / "sholl_summary.csv", index=False)
        print(f"  ✓ Resumen Sholl guardado: sholl_summary.csv (AUC, peak, critical_radius)")
    else:
        print(f"  ⚠ No se pudo generar sholl_summary.csv (sin datos Sholl)")
    
    # Guardar configuración de anillos para visualización con draw.sholl_shells()
    # Formato: diccionario {label: {centroid_um, radii_um}} para compatibilidad con visualizador
    rings_dict = {}
    for label in labels:
        nuclear_mask = (mask_2d == label)
        if not np.any(nuclear_mask):
            continue
        
        props = regionprops(nuclear_mask.astype(np.int32))
        if len(props) == 0:
            continue
        
        cy_px, cx_px = props[0].centroid
        
        # Formato compatible con napari_viewer_2d.py (espera diccionario)
        rings_dict[str(int(label))] = {
            'center_y_px': float(cy_px),
            'center_x_px': float(cx_px),
            'centroid_um': [0.0, float(cy_px * y_um), float(cx_px * x_um)],  # [Z, Y, X] para compatibilidad
            'radii_um': radii.tolist(),
            'radii_px': (radii / np.mean([y_um, x_um])).tolist()  # Para draw.sholl_shells()
        }
    
    # Guardar con ambos nombres para compatibilidad
    with open(out_dir / "sholl_rings.json", 'w') as f:
        json.dump(rings_dict, f, indent=2)
    with open(out_dir / "sholl_rings_2d_native.json", 'w') as f:
        json.dump(rings_dict, f, indent=2)
    print(f"  ✓ Anillos Sholl guardados: sholl_rings.json (+ sholl_rings_2d_native.json para compatibilidad)")
    
    print("\n" + "="*70)
    print(f"✓ PIPELINE COMPLETADO: {len(skeleton_metrics)} astrocitos procesados")
    print("="*70 + "\n")
    
    return df_skeleton, df_sholl
