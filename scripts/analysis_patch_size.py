"""
Analysis: Clearing Patch Size vs Neighbor Dry-Down Rate
========================================================
Paper 2 complementary analysis.

Tests hypothesis: larger contiguous clearings produce faster dry-down
in adjacent intact forest than smaller clearings.

Method:
1. From Hansen GFC, identify contiguous clearing patches and their areas
2. For each SMAP forest pixel within 36km of a clearing, assign the size
   of the nearest clearing patch
3. Bin by clearing size and compute dry-down rates per bin
4. Test whether clearing size predicts neighbor dry-down rate

Uses the edge_distance_9km.nc grid and adds a "nearest clearing size" field.
Because the connected-component labeling was done at 300m in the grid builder,
we re-derive patch sizes here from the saved grid + re-running labels on the
lossyear tiles (at downsampled resolution for memory).
"""

import numpy as np
import netCDF4 as nc
from scipy import stats
from scipy.ndimage import label as ndlabel
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')
import os, sys
import rasterio

# ── Paths ──────────────────────────────────────────────────────────
HANSEN_DIR = r'D:\amazon paper 2\data\hansen'
SMAP_PATH = r'D:\amazon paper\data\smap\smap_l4_rootzone_sm_monthly_amazon_2015_2024.nc'
EDGE_PATH = r'D:\amazon paper 2\data\edge_distance_9km.nc'
FIG_DIR = r'D:\amazon paper 2\figures'
OUT_NC = r'D:\amazon paper 2\data\patch_size_9km.nc'

# ── Parameters ─────────────────────────────────────────────────────
TC_THRESHOLD = 50
LOSS_CUTOFF = 15
DS = 10  # downsample factor (30m -> 300m)
PIXEL_AREA_300M_KM2 = (300 * 300) / 1e6  # 0.09 km2

# Patch size bins (km2) — fine-grained
SIZE_BINS = [1, 3, 10, 30, 100, 300, 1000, 100000]
SIZE_LABELS = ['1-3', '3-10', '10-30', '30-100', '100-300', '300-1000', '1000+']
SIZE_COLORS = ['#ffffcc', '#fee08b', '#fdae61', '#fc8d59', '#f46d43', '#d73027', '#7f0000']

# Only analyze forest pixels within this distance of an edge
MAX_DIST_KM = 36  # focus on near-edge forest

# SMAP grid params
SMAP_RES_DEG = 9 / 111
LAT_MIN, LAT_MAX = -14, 2
LON_MIN, LON_MAX = -70, -44

TREECOVER_TILES = [
    ("00N", "060W"), ("00N", "070W"),
    ("10S", "060W"), ("10S", "070W"),
]


def build_patch_size_map():
    """
    For each tile, compute connected-component clearing patches at 300m
    and create a map of patch sizes. Then for each SMAP grid cell,
    record the size of the nearest substantial clearing.
    """
    print("Building patch-size map from Hansen tiles...")

    # We'll accumulate results at the SMAP grid level
    smap_lats = np.arange(LAT_MIN, LAT_MAX, SMAP_RES_DEG) + SMAP_RES_DEG / 2
    smap_lons = np.arange(LON_MIN, LON_MAX, SMAP_RES_DEG) + SMAP_RES_DEG / 2

    nearest_patch_size = np.full((len(smap_lats), len(smap_lons)), np.nan)

    for lat_str, lon_str in TREECOVER_TILES:
        tc_file = os.path.join(HANSEN_DIR,
            f"Hansen_GFC-2023-v1.11_treecover2000_{lat_str}_{lon_str}.tif")
        ly_file = os.path.join(HANSEN_DIR,
            f"Hansen_GFC-2023-v1.11_lossyear_{lat_str}_{lon_str}.tif")

        if not os.path.exists(tc_file) or not os.path.exists(ly_file):
            print(f"  Skipping tile {lat_str}_{lon_str}")
            continue

        print(f"\n  Tile {lat_str}_{lon_str}:")

        with rasterio.open(tc_file) as ds:
            tc = ds.read(1)
            transform = ds.transform
        with rasterio.open(ly_file) as ds:
            ly = ds.read(1)

        # Build cleared mask
        was_forest = tc >= TC_THRESHOLD
        cleared = was_forest & (ly > 0) & (ly <= LOSS_CUTOFF)
        del tc, ly, was_forest

        # Downsample to 300m
        h, w = cleared.shape
        h_ds, w_ds = h // DS, w // DS
        cleared_ds = cleared[:h_ds*DS, :w_ds*DS].reshape(
            h_ds, DS, w_ds, DS).mean(axis=(1, 3)) > 0.5
        del cleared

        # Connected components
        labels, n_comp = ndlabel(cleared_ds)
        print(f"    Components: {n_comp:,}")

        # Compute area of each component (in km2)
        comp_sizes_px = np.bincount(labels.ravel())
        comp_sizes_km2 = comp_sizes_px * PIXEL_AREA_300M_KM2
        comp_sizes_km2[0] = 0  # background

        # Filter: only patches >= 1 km2
        min_px = int(1.0 / PIXEL_AREA_300M_KM2)
        substantial = comp_sizes_px >= min_px
        substantial[0] = False
        n_subst = np.sum(substantial)
        print(f"    Substantial patches (>=1 km2): {n_subst:,}")

        # Create patch-size array at 300m: each pixel gets the area of its patch
        patch_size_map_ds = comp_sizes_km2[labels]
        patch_size_map_ds[~substantial[labels]] = 0  # zero out small patches

        # Tile geo info
        west = transform.c
        north = transform.f
        px_w = transform.a * DS  # 300m pixel width in degrees
        px_h = -transform.e * DS

        # For each SMAP cell, find the size of the nearest substantial clearing
        for i, slat in enumerate(smap_lats):
            for j, slon in enumerate(smap_lons):
                # SMAP cell bounds
                cell_s = slat - SMAP_RES_DEG / 2
                cell_n = slat + SMAP_RES_DEG / 2
                cell_w = slon - SMAP_RES_DEG / 2
                cell_e = slon + SMAP_RES_DEG / 2

                # Check if cell is within this tile
                tile_south = north - h_ds * px_h
                tile_east = west + w_ds * px_w
                if cell_n < tile_south or cell_s > north:
                    continue
                if cell_e < west or cell_w > tile_east:
                    continue

                # Convert to 300m pixel indices
                col_s = int((cell_w - west) / px_w)
                col_e = int((cell_e - west) / px_w)
                row_s = int((north - cell_n) / px_h)
                row_e = int((north - cell_s) / px_h)

                col_s = max(0, col_s)
                col_e = min(w_ds, col_e)
                row_s = max(0, row_s)
                row_e = min(h_ds, row_e)

                if row_e <= row_s or col_e <= col_s:
                    continue

                # Find the largest clearing patch within/near this cell
                sub = patch_size_map_ds[row_s:row_e, col_s:col_e]
                if np.any(sub > 0):
                    # Take the max patch size touching this cell
                    max_size = np.max(sub)
                    if np.isnan(nearest_patch_size[i, j]) or max_size > nearest_patch_size[i, j]:
                        nearest_patch_size[i, j] = max_size

        del labels, comp_sizes_px, comp_sizes_km2, patch_size_map_ds, cleared_ds

    return smap_lats, smap_lons, nearest_patch_size


def load_smap():
    """Load SMAP L4 root-zone SM."""
    print('\nLoading SMAP data...')
    ds = nc.Dataset(SMAP_PATH, 'r')
    time_raw = ds.variables['time'][:]
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    sm = ds.variables['sm_rootzone'][:]
    base = datetime(2000, 1, 1)
    dates = np.array([base + timedelta(days=float(t)) for t in time_raw])
    years = np.array([d.year for d in dates])
    months = np.array([d.month for d in dates])
    sm = np.ma.masked_less_equal(sm, 0)
    sm = np.ma.masked_greater(sm, 1)
    sm = np.ma.filled(sm, np.nan)
    ds.close()
    return lat, lon, sm, years, months


def load_edge_grid():
    """Load edge-distance and forest fraction."""
    ds = nc.Dataset(EDGE_PATH, 'r')
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    dist = ds.variables['mean_dist_km'][:]
    frac = ds.variables['frac_forest'][:]
    ds.close()
    return lat, lon, dist, frac


def match_and_classify(smap_lat, smap_lon, ps_lat, ps_lon, patch_sizes,
                       edge_lat, edge_lon, edge_dist, edge_frac):
    """
    Match patch-size grid to SMAP grid.
    Only include forest pixels within MAX_DIST_KM of an edge.
    Classify by patch size bin.
    """
    print('\nMatching and classifying...')
    bin_map = np.full((len(smap_lat), len(smap_lon)), -1, dtype=np.int32)
    tol = 0.05

    for i, slat in enumerate(smap_lat):
        for j, slon in enumerate(smap_lon):
            # Match to edge grid
            ei = np.argmin(np.abs(edge_lat - slat))
            ej = np.argmin(np.abs(edge_lon - slon))
            if abs(edge_lat[ei] - slat) > tol or abs(edge_lon[ej] - slon) > tol:
                continue
            dist = edge_dist[ei, ej]
            frac = edge_frac[ei, ej]
            if not np.isfinite(dist) or dist > MAX_DIST_KM or frac < 0.5:
                continue

            # Match to patch-size grid
            pi = np.argmin(np.abs(ps_lat - slat))
            pj = np.argmin(np.abs(ps_lon - slon))
            if abs(ps_lat[pi] - slat) > tol or abs(ps_lon[pj] - slon) > tol:
                continue
            ps = patch_sizes[pi, pj]
            if not np.isfinite(ps) or ps < 1:
                continue

            # Classify by patch size
            for b in range(len(SIZE_LABELS)):
                if SIZE_BINS[b] <= ps < SIZE_BINS[b + 1]:
                    bin_map[i, j] = b
                    break

    for b, label in enumerate(SIZE_LABELS):
        print(f'  Bin {label} km2: {np.sum(bin_map == b)} cells')

    return bin_map


def compute_drydown_by_bin(sm, years, months, bin_map):
    """Compute dry-down rates per patch-size bin."""
    unique_years = sorted(set(years))
    results = {label: [] for label in SIZE_LABELS}

    for b, label in enumerate(SIZE_LABELS):
        cells = np.argwhere(bin_map == b)
        if len(cells) == 0:
            continue

        ts = np.full(len(years), np.nan)
        for t in range(len(years)):
            vals = [sm[t, ci, cj] for ci, cj in cells if np.isfinite(sm[t, ci, cj])]
            if vals:
                ts[t] = np.mean(vals)

        for yr in unique_years:
            idx = [(months == m) & (years == yr) for m in range(3, 9)]
            vals = []
            for m_mask in idx:
                if np.any(m_mask):
                    v = ts[m_mask]
                    v = v[np.isfinite(v)]
                    if len(v) > 0:
                        vals.append(v[0])
            if len(vals) >= 4:
                slope, _, _, _, _ = stats.linregress(np.arange(len(vals)), vals)
                results[label].append(slope)

    return results


def compute_efold_by_bin(sm, years, months, bin_map):
    """Compute e-folding times per patch-size bin."""
    unique_years = sorted(set(years))
    results = {label: [] for label in SIZE_LABELS}

    for b, label in enumerate(SIZE_LABELS):
        cells = np.argwhere(bin_map == b)
        if len(cells) == 0:
            continue

        ts = np.full(len(years), np.nan)
        for t in range(len(years)):
            vals = [sm[t, ci, cj] for ci, cj in cells if np.isfinite(sm[t, ci, cj])]
            if vals:
                ts[t] = np.mean(vals)

        for yr in unique_years:
            idx = [(months == m) & (years == yr) for m in range(3, 9)]
            vals = []
            for m_mask in idx:
                if np.any(m_mask):
                    v = ts[m_mask]
                    v = v[np.isfinite(v)]
                    if len(v) > 0:
                        vals.append(v[0])
            if len(vals) >= 4:
                arr = np.array(vals)
                if np.all(arr > 0):
                    slope, _, _, _, _ = stats.linregress(np.arange(len(arr)), np.log(arr))
                    if slope < 0:
                        tau = -1.0 / slope
                        if 1 < tau < 100:
                            results[label].append(tau)

    return results


def plot_results(drydown, efold, bin_map):
    """Create figure."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # (a) Dry-down rate by clearing size
    ax = axes[0]
    box_data = [drydown.get(label, []) for label in SIZE_LABELS]
    valid_labels = [l for l, d in zip(SIZE_LABELS, box_data) if d]
    valid_data = [d for d in box_data if d]
    valid_colors = [c for c, d in zip(SIZE_COLORS, box_data) if d]

    if valid_data:
        bp = ax.boxplot(valid_data, widths=0.5, patch_artist=True,
                       medianprops=dict(color='black', lw=1.5))
        for patch, c in zip(bp['boxes'], valid_colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        ax.set_xticklabels(valid_labels, fontsize=9)

    ax.set_xlabel('Nearest clearing size (km$^2$)')
    ax.set_ylabel('Dry-down rate (m$^3$ m$^{-3}$ month$^{-1}$)')
    ax.set_title('(a) Neighbor dry-down by clearing size', fontweight='bold', loc='left')
    ax.axhline(0, color='grey', ls='--', lw=0.7)
    ax.grid(True, axis='y', alpha=0.3)

    # (b) E-folding time by clearing size
    ax = axes[1]
    efold_data = [efold.get(label, []) for label in SIZE_LABELS]
    valid_labels_e = [l for l, d in zip(SIZE_LABELS, efold_data) if d]
    valid_data_e = [d for d in efold_data if d]
    valid_colors_e = [c for c, d in zip(SIZE_COLORS, efold_data) if d]

    if valid_data_e:
        bp2 = ax.boxplot(valid_data_e, widths=0.5, patch_artist=True,
                        medianprops=dict(color='black', lw=1.5))
        for patch, c in zip(bp2['boxes'], valid_colors_e):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        ax.set_xticklabels(valid_labels_e, fontsize=9)

    ax.set_xlabel('Nearest clearing size (km$^2$)')
    ax.set_ylabel('E-folding time (months)')
    ax.set_title('(b) Moisture residence time by clearing size', fontweight='bold', loc='left')
    ax.grid(True, axis='y', alpha=0.3)

    # (c) Sample size
    ax = axes[2]
    counts = [np.sum(bin_map == b) for b in range(len(SIZE_LABELS))]
    ax.bar(range(len(SIZE_LABELS)), counts, color=SIZE_COLORS, alpha=0.7)
    ax.set_xticks(range(len(SIZE_LABELS)))
    ax.set_xticklabels(SIZE_LABELS, fontsize=9)
    ax.set_xlabel('Nearest clearing size (km$^2$)')
    ax.set_ylabel('Number of 9km forest cells')
    ax.set_title('(c) Sample size', fontweight='bold', loc='left')
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, 'fig_patch_size_drydown.png')
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'Figure saved: {out}')


def run_statistics(drydown, efold):
    """Statistical tests."""
    print('\n' + '='*60)
    print('STATISTICAL TESTS')
    print('='*60)

    groups = [np.array(drydown[l]) for l in SIZE_LABELS if drydown[l]]
    if len(groups) >= 2:
        F, p = stats.f_oneway(*groups)
        print(f'\nANOVA on dry-down rates: F={F:.3f}, p={p:.6f}')

    # Smallest vs largest clearing
    smallest = SIZE_LABELS[0]
    largest = SIZE_LABELS[-1]
    if drydown[smallest] and drydown[largest]:
        t, p = stats.ttest_ind(drydown[smallest], drydown[largest], equal_var=False)
        print(f'\n{smallest} km2 vs {largest} km2:')
        print(f'  Welch t = {t:.3f}, p = {p:.6f}')
        print(f'  Small clearing mean: {np.mean(drydown[smallest]):.5f}')
        print(f'  Large clearing mean: {np.mean(drydown[largest]):.5f}')

    print('\nE-folding times:')
    for label in SIZE_LABELS:
        vals = efold[label]
        if vals:
            print(f'  {label} km2: {np.mean(vals):.1f} +/- {np.std(vals):.1f} months (n={len(vals)})')

    # Spearman
    bin_centers = [np.sqrt(SIZE_BINS[i] * SIZE_BINS[i+1]) for i in range(len(SIZE_LABELS))]
    means = [np.mean(drydown[l]) if drydown[l] else np.nan for l in SIZE_LABELS]
    valid = [(x, m) for x, m in zip(bin_centers, means) if np.isfinite(m)]
    if len(valid) >= 3:
        xs, ms = zip(*valid)
        rho, p = stats.spearmanr(xs, ms)
        print(f'\nSpearman: clearing size vs dry-down rate: rho={rho:.3f}, p={p:.4f}')


# ── Main ───────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*60)
    print("Patch-Size Analysis: Clearing Size vs Neighbor Dry-Down")
    print(f"Max distance from edge: {MAX_DIST_KM} km")
    print(f"Size bins: {SIZE_BINS} km2")
    print("="*60)

    # Build patch-size map from Hansen tiles
    ps_lat, ps_lon, patch_sizes = build_patch_size_map()

    n_valid = np.sum(np.isfinite(patch_sizes) & (patch_sizes > 0))
    print(f"\nPatch-size grid: {len(ps_lat)} x {len(ps_lon)}, {n_valid} cells with data")

    # Load SMAP and edge grid
    smap_lat, smap_lon, sm, years, months = load_smap()
    edge_lat, edge_lon, edge_dist, edge_frac = load_edge_grid()

    # Match and classify
    bin_map = match_and_classify(
        smap_lat, smap_lon, ps_lat, ps_lon, patch_sizes,
        edge_lat, edge_lon, edge_dist, edge_frac
    )

    # Compute dry-down rates
    print('\nComputing dry-down rates by patch size...')
    drydown = compute_drydown_by_bin(sm, years, months, bin_map)

    # Compute e-folding times
    print('Computing e-folding times by patch size...')
    efold = compute_efold_by_bin(sm, years, months, bin_map)

    # Statistics
    run_statistics(drydown, efold)

    # Plot
    print('\nPlotting...')
    plot_results(drydown, efold, bin_map)

    print('\nDone.')
