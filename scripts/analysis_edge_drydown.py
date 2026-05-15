"""
Analysis: Edge-Distance Stratified Dry-Down Rates
===================================================
Paper 2 core analysis.

Combines:
1. SMAP L4 root-zone soil moisture (9km, 2015-2024) from Paper 1
2. Edge-distance grid (9km, from build_edge_distance_grid.py)

Tests hypothesis: intact forest pixels closer to deforestation edges
dry faster than deep-interior forest pixels.

Method:
- Classify SMAP pixels into edge-distance bins using the distance grid
- For each bin, compute post-wet-season dry-down rates (Mar-Aug slope)
  exactly as in Paper 1's analysis7_smap_residence.py
- Compare rates across bins using ANOVA + pairwise tests
- Plot distance-decay curve: dry-down rate vs distance to edge

Also computes:
- Lag-1 autocorrelation (persistence) by distance bin
- Seasonal cycle by distance bin
- E-folding dry-down time by distance bin
"""

import numpy as np
import netCDF4 as nc
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')
import os, sys

# ── Paths ──────────────────────────────────────────────────────────
SMAP_PATH = r'D:\amazon paper\data\smap\smap_l4_rootzone_sm_monthly_amazon_2015_2024.nc'
EDGE_PATH = r'D:\amazon paper 2\data\edge_distance_9km.nc'
FIG_DIR = r'D:\amazon paper 2\figures'
os.makedirs(FIG_DIR, exist_ok=True)

# ── Edge-distance bins ─────────────────────────────────────────────
# Bin edges in km — chosen to give roughly logarithmic spacing
BIN_EDGES = [0, 9, 18, 36, 72, 144, 500]
BIN_LABELS = ['0-9', '9-18', '18-36', '36-72', '72-144', '144+']
BIN_COLORS = ['#d73027', '#fc8d59', '#fee08b', '#d9ef8b', '#91cf60', '#1a9850']

# Minimum forest fraction to include a cell (avoid mixed cells)
MIN_FOREST_FRAC = 0.5


def load_smap():
    """Load SMAP L4 root-zone SM — same loader as Paper 1."""
    print('Loading SMAP data ...')
    ds = nc.Dataset(SMAP_PATH, 'r')
    time_raw = ds.variables['time'][:]
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    sm = ds.variables['sm_rootzone'][:]

    base = datetime(2000, 1, 1)
    dates = np.array([base + timedelta(days=float(t)) for t in time_raw])
    years = np.array([d.year for d in dates])
    months = np.array([d.month for d in dates])
    print(f'  Time: {dates[0]:%Y-%m} to {dates[-1]:%Y-%m}, {len(dates)} months')
    print(f'  Grid: {len(lat)} x {len(lon)}')

    sm = np.ma.masked_less_equal(sm, 0)
    sm = np.ma.masked_greater(sm, 1)
    sm = np.ma.filled(sm, np.nan)

    ds.close()
    return lat, lon, sm, dates, years, months


def load_edge_grid():
    """Load edge-distance grid from build_edge_distance_grid.py."""
    print('Loading edge-distance grid ...')
    ds = nc.Dataset(EDGE_PATH, 'r')
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    mean_dist = ds.variables['mean_dist_km'][:]
    frac_forest = ds.variables['frac_forest'][:]
    frac_cleared = ds.variables['frac_cleared'][:]
    ds.close()
    print(f'  Grid: {len(lat)} x {len(lon)}')
    print(f'  Valid cells: {np.sum(np.isfinite(mean_dist))}')
    return lat, lon, mean_dist, frac_forest, frac_cleared


def match_grids(smap_lat, smap_lon, edge_lat, edge_lon, edge_dist, edge_frac_f):
    """
    Map edge-distance values onto SMAP grid cells.
    Uses nearest-neighbor matching within tolerance.
    """
    print('Matching grids ...')
    dist_on_smap = np.full((len(smap_lat), len(smap_lon)), np.nan)
    frac_on_smap = np.full((len(smap_lat), len(smap_lon)), np.nan)

    tol = 0.05  # degrees (~5km tolerance for matching)

    matched = 0
    for i, slat in enumerate(smap_lat):
        for j, slon in enumerate(smap_lon):
            # Find nearest edge-grid cell
            lat_idx = np.argmin(np.abs(edge_lat - slat))
            lon_idx = np.argmin(np.abs(edge_lon - slon))

            if (abs(edge_lat[lat_idx] - slat) < tol and
                abs(edge_lon[lon_idx] - slon) < tol):
                dist_on_smap[i, j] = edge_dist[lat_idx, lon_idx]
                frac_on_smap[i, j] = edge_frac_f[lat_idx, lon_idx]
                if np.isfinite(edge_dist[lat_idx, lon_idx]):
                    matched += 1

    print(f'  Matched {matched} cells')
    return dist_on_smap, frac_on_smap


def classify_pixels(dist_grid, frac_grid):
    """
    Assign each SMAP pixel to an edge-distance bin.
    Returns bin_map (same shape as dist_grid) with bin indices 0..N-1, or -1 for excluded.
    """
    bin_map = np.full(dist_grid.shape, -1, dtype=np.int32)

    for b in range(len(BIN_LABELS)):
        lo = BIN_EDGES[b]
        hi = BIN_EDGES[b + 1]

        mask = (dist_grid >= lo) & (dist_grid < hi) & (frac_grid >= MIN_FOREST_FRAC)
        bin_map[mask] = b

    for b in range(len(BIN_LABELS)):
        n = np.sum(bin_map == b)
        print(f'  Bin {BIN_LABELS[b]} km: {n} cells')

    return bin_map


def compute_drydown_by_bin(sm, years, months, bin_map, smap_lat, smap_lon):
    """
    Compute post-wet-season dry-down rate (Mar-Aug linear slope)
    for each distance bin, exactly as Paper 1 does for arc vs interior.

    Returns dict: bin_label -> list of slopes (m3/m3/month).
    """
    unique_years = sorted(set(years))
    results = {label: [] for label in BIN_LABELS}

    for b, label in enumerate(BIN_LABELS):
        # Get all SMAP cells in this bin
        cells = np.argwhere(bin_map == b)
        if len(cells) == 0:
            continue

        # Spatial mean SM for this bin at each time step
        ts = np.full(len(years), np.nan)
        for t in range(len(years)):
            vals = []
            for ci, cj in cells:
                v = sm[t, ci, cj]
                if np.isfinite(v):
                    vals.append(v)
            if vals:
                ts[t] = np.mean(vals)

        # Compute dry-down slope for each year
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
                x = np.arange(len(vals))
                slope, _, _, _, _ = stats.linregress(x, vals)
                results[label].append(slope)

    return results


def compute_efold_by_bin(sm, years, months, bin_map):
    """
    Compute e-folding dry-down time for each bin.
    Fit exponential decay to Mar-Aug SM for each year.
    """
    unique_years = sorted(set(years))
    results = {label: [] for label in BIN_LABELS}

    for b, label in enumerate(BIN_LABELS):
        cells = np.argwhere(bin_map == b)
        if len(cells) == 0:
            continue

        ts = np.full(len(years), np.nan)
        for t in range(len(years)):
            vals = [sm[t, ci, cj] for ci, cj in cells
                    if np.isfinite(sm[t, ci, cj])]
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
                # Fit: SM(t) = SM0 * exp(-t/tau)
                # log(SM) = log(SM0) - t/tau
                vals_arr = np.array(vals)
                if np.all(vals_arr > 0):
                    x = np.arange(len(vals_arr))
                    log_vals = np.log(vals_arr)
                    slope, _, _, _, _ = stats.linregress(x, log_vals)
                    if slope < 0:
                        tau = -1.0 / slope  # e-folding time in months
                        if 1 < tau < 100:
                            results[label].append(tau)

    return results


def plot_results(drydown, efold, bin_map):
    """Create multi-panel figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # ── (a) Dry-down rate by distance bin ─────────────────────────
    ax = axes[0, 0]
    box_data = [drydown.get(label, []) for label in BIN_LABELS]
    valid_labels = [label for label, d in zip(BIN_LABELS, box_data) if len(d) > 0]
    valid_data = [d for d in box_data if len(d) > 0]
    valid_colors = [c for c, d in zip(BIN_COLORS, box_data) if len(d) > 0]

    if valid_data:
        bp = ax.boxplot(valid_data, widths=0.5, patch_artist=True,
                       medianprops=dict(color='black', lw=1.5))
        for patch, c in zip(bp['boxes'], valid_colors):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        ax.set_xticklabels(valid_labels, fontsize=8)

    ax.set_xlabel('Distance to deforestation edge (km)')
    ax.set_ylabel('Dry-down rate (m$^3$ m$^{-3}$ month$^{-1}$)')
    ax.set_title('(a) Post-wet-season dry-down by edge distance', fontweight='bold', loc='left')
    ax.axhline(0, color='grey', ls='--', lw=0.7)
    ax.grid(True, axis='y', alpha=0.3)

    # ── (b) Mean dry-down rate vs distance (the decay curve) ──────
    ax = axes[0, 1]
    bin_centers = [(BIN_EDGES[i] + BIN_EDGES[i+1])/2 for i in range(len(BIN_LABELS))]
    means = [np.mean(drydown[label]) if drydown[label] else np.nan for label in BIN_LABELS]
    sems = [stats.sem(drydown[label]) if len(drydown[label]) > 1 else np.nan for label in BIN_LABELS]

    valid_x = [x for x, m in zip(bin_centers, means) if np.isfinite(m)]
    valid_m = [m for m in means if np.isfinite(m)]
    valid_s = [s for s, m in zip(sems, means) if np.isfinite(m)]

    if valid_x:
        ax.errorbar(valid_x, valid_m, yerr=valid_s, fmt='o-', color='#d73027',
                    lw=2, ms=8, capsize=4, label='Dry-down rate')
        ax.set_xlabel('Distance to deforestation edge (km)')
        ax.set_ylabel('Mean dry-down rate (m$^3$ m$^{-3}$ month$^{-1}$)')
        ax.set_title('(b) Edge-distance decay curve', fontweight='bold', loc='left')
        ax.axhline(0, color='grey', ls='--', lw=0.7)
        ax.grid(True, alpha=0.3)

    # ── (c) E-folding time by distance bin ────────────────────────
    ax = axes[1, 0]
    efold_data = [efold.get(label, []) for label in BIN_LABELS]
    valid_labels_e = [label for label, d in zip(BIN_LABELS, efold_data) if len(d) > 0]
    valid_data_e = [d for d in efold_data if len(d) > 0]
    valid_colors_e = [c for c, d in zip(BIN_COLORS, efold_data) if len(d) > 0]

    if valid_data_e:
        bp2 = ax.boxplot(valid_data_e, widths=0.5, patch_artist=True,
                        medianprops=dict(color='black', lw=1.5))
        for patch, c in zip(bp2['boxes'], valid_colors_e):
            patch.set_facecolor(c)
            patch.set_alpha(0.6)
        ax.set_xticklabels(valid_labels_e, fontsize=8)

    ax.set_xlabel('Distance to deforestation edge (km)')
    ax.set_ylabel('E-folding time (months)')
    ax.set_title('(c) Moisture residence time by edge distance', fontweight='bold', loc='left')
    ax.grid(True, axis='y', alpha=0.3)

    # Paper 1 reference lines
    ax.axhline(13.7, color='#1B7837', ls='--', lw=1, alpha=0.7, label='Paper 1: intact (13.7 mo)')
    ax.axhline(8.9, color='#D95F02', ls='--', lw=1, alpha=0.7, label='Paper 1: arc (8.9 mo)')
    ax.legend(frameon=False, fontsize=8)

    # ── (d) Number of cells per bin ───────────────────────────────
    ax = axes[1, 1]
    counts = [np.sum(bin_map == b) for b in range(len(BIN_LABELS))]
    bars = ax.bar(range(len(BIN_LABELS)), counts, color=BIN_COLORS, alpha=0.7)
    ax.set_xticks(range(len(BIN_LABELS)))
    ax.set_xticklabels(BIN_LABELS, fontsize=8)
    ax.set_xlabel('Distance to deforestation edge (km)')
    ax.set_ylabel('Number of 9km SMAP cells')
    ax.set_title('(d) Sample size per distance bin', fontweight='bold', loc='left')
    ax.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, 'fig_edge_drydown.png')
    plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'Figure saved: {out}')


def run_statistics(drydown, efold):
    """Run ANOVA and pairwise tests."""
    print('\n' + '='*60)
    print('STATISTICAL TESTS')
    print('='*60)

    # One-way ANOVA on dry-down rates
    groups = [np.array(drydown[label]) for label in BIN_LABELS if drydown[label]]
    if len(groups) >= 2:
        F, p = stats.f_oneway(*groups)
        print(f'\nOne-way ANOVA on dry-down rates: F={F:.3f}, p={p:.6f}')

    # Pairwise: nearest bin vs farthest bin
    nearest_label = BIN_LABELS[0]
    farthest_label = BIN_LABELS[-1]
    if drydown[nearest_label] and drydown[farthest_label]:
        t, p = stats.ttest_ind(drydown[nearest_label], drydown[farthest_label], equal_var=False)
        print(f'\nNearest ({nearest_label} km) vs farthest ({farthest_label} km):')
        print(f'  Welch t = {t:.3f}, p = {p:.6f}')
        print(f'  Nearest mean: {np.mean(drydown[nearest_label]):.5f} m3/m3/month')
        print(f'  Farthest mean: {np.mean(drydown[farthest_label]):.5f} m3/m3/month')

    # E-folding time comparison
    print('\nE-folding times by bin:')
    for label in BIN_LABELS:
        vals = efold[label]
        if vals:
            print(f'  {label} km: {np.mean(vals):.1f} +/- {np.std(vals):.1f} months (n={len(vals)})')

    # Spearman correlation: bin center vs mean dry-down rate
    bin_centers = [(BIN_EDGES[i] + BIN_EDGES[i+1])/2 for i in range(len(BIN_LABELS))]
    means = [np.mean(drydown[label]) if drydown[label] else np.nan for label in BIN_LABELS]
    valid = [(x, m) for x, m in zip(bin_centers, means) if np.isfinite(m)]
    if len(valid) >= 3:
        xs, ms = zip(*valid)
        rho, p = stats.spearmanr(xs, ms)
        print(f'\nSpearman: distance vs dry-down rate: rho={rho:.3f}, p={p:.4f}')


# ── Main ───────────────────────────────────────────────────────────
if __name__ == '__main__':
    # Load data
    smap_lat, smap_lon, sm, dates, years, months = load_smap()
    edge_lat, edge_lon, edge_dist, edge_frac_f, edge_frac_c = load_edge_grid()

    # Match grids
    dist_on_smap, frac_on_smap = match_grids(
        smap_lat, smap_lon, edge_lat, edge_lon, edge_dist, edge_frac_f
    )

    # Classify into distance bins
    bin_map = classify_pixels(dist_on_smap, frac_on_smap)

    # Compute dry-down rates by bin
    print('\nComputing dry-down rates by distance bin...')
    drydown = compute_drydown_by_bin(sm, years, months, bin_map, smap_lat, smap_lon)

    # Compute e-folding times by bin
    print('\nComputing e-folding times by distance bin...')
    efold = compute_efold_by_bin(sm, years, months, bin_map)

    # Statistics
    run_statistics(drydown, efold)

    # Plot
    print('\nPlotting...')
    plot_results(drydown, efold, bin_map)

    print('\nDone.')
