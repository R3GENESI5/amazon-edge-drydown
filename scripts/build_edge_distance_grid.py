"""
Build deforestation edge-distance grid for Paper 2.
====================================================
1. Load Hansen treecover2000 + lossyear tiles
2. Create binary forest/cleared mask (forest = treecover2000 >= 50% AND no loss by 2015)
3. Compute Euclidean distance from each forest pixel to nearest cleared pixel (at 30m)
4. Aggregate to SMAP L4 9km grid:
   - Mean distance-to-edge per 9km cell
   - Fraction of cell that is forest
   - Fraction of cell that is deforested
5. Save as NetCDF for use by dry-down analysis

Uses Hansen GFC v1.11:
- treecover2000: % canopy cover in 2000 (0-100)
- lossyear: year of loss (1-23 = 2001-2023, 0 = no loss)

Processing strategy: tile-by-tile to manage memory (~600MB per treecover tile).
"""

import numpy as np
import os
import sys
from datetime import datetime

# Lazy imports — check availability up front
missing = []
try:
    import rasterio
    from rasterio.merge import merge as rasterio_merge
    from rasterio.transform import from_bounds
    from rasterio.warp import reproject, Resampling
except ImportError:
    missing.append('rasterio')
try:
    from scipy.ndimage import distance_transform_edt
except ImportError:
    missing.append('scipy')
try:
    import netCDF4 as nc
except ImportError:
    missing.append('netCDF4')

if missing:
    print(f"Missing packages: {', '.join(missing)}")
    print("Install with: pip install rasterio scipy netCDF4")
    sys.exit(1)

# ── Paths ──────────────────────────────────────────────────────────
HANSEN_DIR = r"D:\amazon paper 2\data\hansen"
OUT_DIR = r"D:\amazon paper 2\data"
OUT_NC = os.path.join(OUT_DIR, "edge_distance_9km.nc")

# ── Parameters ─────────────────────────────────────────────────────
TC_THRESHOLD = 50       # Minimum treecover2000 (%) to count as "was forest"
LOSS_CUTOFF = 15        # lossyear <= 15 means cleared by 2015 (start of SMAP)
MIN_CLEARING_KM2 = 1.0  # Minimum contiguous clearing area to count as "edge"
SMAP_RES_DEG = 9 / 111  # ~9km in degrees (~0.081 deg)

# Study domain (matches Paper 1 extended)
LAT_MIN, LAT_MAX = -14, 2
LON_MIN, LON_MAX = -70, -44

# Hansen tile definitions we have
TREECOVER_TILES = [
    ("00N", "060W"), ("00N", "070W"),
    ("10S", "060W"), ("10S", "070W"),
]
# These cover 10N to 20S, 80W to 60W — the core Amazon

def tile_bounds(lat_str, lon_str):
    """Convert Hansen tile name to (west, south, east, north)."""
    lat_val = int(lat_str[:-1])
    lon_val = int(lon_str[:-1])
    if lat_str.endswith('S'):
        lat_val = -lat_val
    if lon_str.endswith('W'):
        lon_val = -lon_val
    # Each tile is 10x10 degrees
    return (lon_val, lat_val - 10, lon_val + 10, lat_val)


def process_tile(lat_str, lon_str):
    """
    Process one 10x10 degree Hansen tile.
    Returns (forest_mask, cleared_mask, transform, shape) at 30m.

    To keep memory manageable, we work in blocks.
    """
    tc_file = os.path.join(
        HANSEN_DIR,
        f"Hansen_GFC-2023-v1.11_treecover2000_{lat_str}_{lon_str}.tif"
    )
    ly_file = os.path.join(
        HANSEN_DIR,
        f"Hansen_GFC-2023-v1.11_lossyear_{lat_str}_{lon_str}.tif"
    )

    if not os.path.exists(tc_file):
        print(f"  WARNING: {tc_file} not found, skipping tile")
        return None
    if not os.path.exists(ly_file):
        print(f"  WARNING: {ly_file} not found, skipping tile")
        return None

    print(f"  Loading treecover: {os.path.basename(tc_file)}")
    with rasterio.open(tc_file) as tc_ds:
        tc = tc_ds.read(1)          # uint8, 0-100
        transform = tc_ds.transform
        crs = tc_ds.crs

    print(f"  Loading lossyear: {os.path.basename(ly_file)}")
    with rasterio.open(ly_file) as ly_ds:
        ly = ly_ds.read(1)          # uint8, 0-23

    print(f"  Tile shape: {tc.shape}, dtype: {tc.dtype}")

    # Was-forest: treecover >= threshold in 2000
    was_forest = tc >= TC_THRESHOLD

    # Cleared by 2015: had forest in 2000 AND lost it by lossyear <= LOSS_CUTOFF
    cleared_raw = was_forest & (ly > 0) & (ly <= LOSS_CUTOFF)

    n_raw = np.sum(cleared_raw)
    print(f"  Raw cleared pixels: {n_raw:,}")

    # Filter to substantial clearings only (>= MIN_CLEARING_KM2)
    # At 30m, 1 km2 ~ 1111 pixels. Working at 30m on 40k x 40k is too large
    # for scipy.ndimage.label, so we downsample 10x to 300m, label, filter, upsample.
    from scipy.ndimage import label as ndlabel, zoom

    DS = 10  # downsample factor (30m -> 300m)
    min_px_300m = int(MIN_CLEARING_KM2 * 1e6 / (300 * 300))  # ~11 pixels at 300m

    # Downsample: a 300m pixel is "cleared" if majority of its 10x10 sub-pixels are
    h, w = cleared_raw.shape
    h_ds, w_ds = h // DS, w // DS
    cleared_ds = cleared_raw[:h_ds*DS, :w_ds*DS].reshape(h_ds, DS, w_ds, DS).mean(axis=(1,3)) > 0.5

    print(f"  Downsampled to {cleared_ds.shape} for connected-component labeling...")
    labels, n_components = ndlabel(cleared_ds)
    print(f"  Found {n_components:,} connected components")

    # Count pixels per component
    component_sizes = np.bincount(labels.ravel())
    # Keep only components >= min_px_300m
    keep = component_sizes >= min_px_300m
    keep[0] = False  # background
    n_kept = np.sum(keep) - (1 if keep[0] else 0)
    print(f"  Components >= {MIN_CLEARING_KM2} km2: {n_kept:,}")

    # Build filtered mask at 300m
    cleared_ds_filtered = keep[labels]

    # Upsample back to 30m using nearest-neighbor
    cleared_filtered = np.zeros_like(cleared_raw)
    cleared_filtered[:h_ds*DS, :w_ds*DS] = np.repeat(
        np.repeat(cleared_ds_filtered, DS, axis=0), DS, axis=1
    )
    # Intersect with original cleared mask (only keep pixels that were actually cleared)
    cleared = cleared_raw & cleared_filtered

    n_filtered = np.sum(cleared)
    print(f"  Filtered cleared pixels: {n_filtered:,} ({100*n_filtered/max(n_raw,1):.1f}% of raw)")

    # Still-forest in 2015: had forest and NOT cleared
    forest = was_forest & ~cleared_raw  # forest = not cleared by ANY loss, not just large

    n_forest = np.sum(forest)
    n_total = tc.size
    print(f"  Forest (2015): {n_forest:,} px ({100*n_forest/n_total:.1f}%)")
    print(f"  Substantial clearings: {n_filtered:,} px ({100*n_filtered/n_total:.1f}%)")

    # Free memory
    del tc, ly, was_forest, cleared_raw, cleared_ds, labels, component_sizes
    del keep, cleared_ds_filtered, cleared_filtered

    return forest, cleared, transform, crs


def compute_distance_to_edge(forest, cleared, pixel_size_m=30):
    """
    Compute distance from each forest pixel to nearest cleared pixel.
    Uses scipy EDT on the inverse of the cleared mask.

    Returns distance array (same shape) in km. Non-forest pixels = NaN.
    """
    print("  Computing distance to deforestation edge (this may take a while)...")

    # EDT needs a binary array where 0 = the features we measure distance TO
    # We want distance from forest pixels to nearest cleared pixel
    # So: 0 where cleared, 1 elsewhere
    mask = ~cleared  # True where NOT cleared

    # EDT gives distance in pixel units from each True pixel to nearest False pixel
    dist_px = distance_transform_edt(mask)

    # Convert to km
    dist_km = dist_px * pixel_size_m / 1000.0

    # Only keep values for forest pixels
    dist_km[~forest] = np.nan

    print(f"  Distance range: {np.nanmin(dist_km):.1f} to {np.nanmax(dist_km):.1f} km")
    print(f"  Median distance (forest pixels): {np.nanmedian(dist_km):.1f} km")

    del mask, dist_px
    return dist_km


def aggregate_to_smap_grid(forest, cleared, dist_km, transform):
    """
    Aggregate 30m products to SMAP 9km grid.

    For each 9km cell, compute:
    - frac_forest: fraction of cell that is still forest
    - frac_cleared: fraction cleared by 2015
    - mean_dist_km: mean distance-to-edge for forest pixels
    - median_dist_km: median distance-to-edge for forest pixels
    - n_forest_px: number of forest pixels (for weighting)
    """
    print("  Aggregating to 9km SMAP grid...")

    nrows, ncols = forest.shape

    # Get lat/lon bounds from transform
    # rasterio transform: (west, pixel_width, 0, north, 0, -pixel_height)
    west = transform.c
    north = transform.f
    px_w = transform.a
    px_h = -transform.e  # positive
    east = west + ncols * px_w
    south = north - nrows * px_h

    # Clip to study domain
    lat_start = max(south, LAT_MIN)
    lat_end = min(north, LAT_MAX)
    lon_start = max(west, LON_MIN)
    lon_end = min(east, LON_MAX)

    # SMAP grid edges
    smap_lats = np.arange(lat_start, lat_end, SMAP_RES_DEG)
    smap_lons = np.arange(lon_start, lon_end, SMAP_RES_DEG)
    n_slat = len(smap_lats)
    n_slon = len(smap_lons)

    print(f"  SMAP grid: {n_slat} x {n_slon} cells")

    frac_forest = np.full((n_slat, n_slon), np.nan)
    frac_cleared = np.full((n_slat, n_slon), np.nan)
    mean_dist = np.full((n_slat, n_slon), np.nan)
    median_dist = np.full((n_slat, n_slon), np.nan)
    n_forest_px = np.zeros((n_slat, n_slon), dtype=np.int32)

    for i, slat in enumerate(smap_lats):
        for j, slon in enumerate(smap_lons):
            # 9km cell bounds
            cell_s = slat
            cell_n = slat + SMAP_RES_DEG
            cell_w = slon
            cell_e = slon + SMAP_RES_DEG

            # Convert to pixel indices
            col_start = int((cell_w - west) / px_w)
            col_end = int((cell_e - west) / px_w)
            row_start = int((north - cell_n) / px_h)
            row_end = int((north - cell_s) / px_h)

            # Clip to array bounds
            col_start = max(0, col_start)
            col_end = min(ncols, col_end)
            row_start = max(0, row_start)
            row_end = min(nrows, row_end)

            if row_end <= row_start or col_end <= col_start:
                continue

            # Extract sub-arrays
            f_sub = forest[row_start:row_end, col_start:col_end]
            c_sub = cleared[row_start:row_end, col_start:col_end]
            d_sub = dist_km[row_start:row_end, col_start:col_end]

            n_total = f_sub.size
            if n_total == 0:
                continue

            n_f = np.sum(f_sub)
            n_c = np.sum(c_sub)

            frac_forest[i, j] = n_f / n_total
            frac_cleared[i, j] = n_c / n_total
            n_forest_px[i, j] = n_f

            if n_f > 0:
                forest_dists = d_sub[f_sub]
                valid = forest_dists[np.isfinite(forest_dists)]
                if len(valid) > 0:
                    mean_dist[i, j] = np.mean(valid)
                    median_dist[i, j] = np.median(valid)

    # Cell centers (for NetCDF)
    lat_centers = smap_lats + SMAP_RES_DEG / 2
    lon_centers = smap_lons + SMAP_RES_DEG / 2

    return {
        'lat': lat_centers,
        'lon': lon_centers,
        'frac_forest': frac_forest,
        'frac_cleared': frac_cleared,
        'mean_dist_km': mean_dist,
        'median_dist_km': median_dist,
        'n_forest_px': n_forest_px,
    }


def save_netcdf(results, out_path):
    """Save aggregated results as NetCDF."""
    print(f"\nSaving to {out_path}...")

    ds = nc.Dataset(out_path, 'w', format='NETCDF4')
    ds.title = "Deforestation edge-distance grid for Amazon basin (9km)"
    ds.source = "Hansen GFC v1.11 (treecover2000, lossyear) aggregated to ~9km"
    ds.history = f"Created {datetime.now().isoformat()} by build_edge_distance_grid.py"
    ds.treecover_threshold = f"{TC_THRESHOLD}%"
    ds.loss_cutoff_year = f"2001-20{LOSS_CUTOFF:02d}"

    lat = results['lat']
    lon = results['lon']

    ds.createDimension('lat', len(lat))
    ds.createDimension('lon', len(lon))

    lat_var = ds.createVariable('lat', 'f4', ('lat',))
    lat_var.units = 'degrees_north'
    lat_var[:] = lat

    lon_var = ds.createVariable('lon', 'f4', ('lon',))
    lon_var.units = 'degrees_east'
    lon_var[:] = lon

    for name, data, units, desc in [
        ('frac_forest', results['frac_forest'], '1',
         'Fraction of 9km cell that was forest in 2015 (treecover2000>=50% and no loss by 2015)'),
        ('frac_cleared', results['frac_cleared'], '1',
         'Fraction of 9km cell cleared by 2015 (had forest in 2000, lost by 2015)'),
        ('mean_dist_km', results['mean_dist_km'], 'km',
         'Mean distance from forest pixels to nearest cleared pixel within 9km cell'),
        ('median_dist_km', results['median_dist_km'], 'km',
         'Median distance from forest pixels to nearest cleared pixel within 9km cell'),
        ('n_forest_px', results['n_forest_px'], '1',
         'Number of 30m forest pixels in 9km cell'),
    ]:
        var = ds.createVariable(name, 'f4', ('lat', 'lon'),
                               fill_value=np.nan if name != 'n_forest_px' else -1)
        var.units = units
        var.long_name = desc
        var[:] = data

    ds.close()
    print(f"  Saved: {os.path.getsize(out_path)/1024:.1f} KB")


# ── Main ───────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("="*60)
    print("Build Edge-Distance Grid for Paper 2")
    print(f"Treecover threshold: {TC_THRESHOLD}%")
    print(f"Loss cutoff: 2001-20{LOSS_CUTOFF:02d}")
    print(f"SMAP resolution: ~{SMAP_RES_DEG:.4f} deg (~9km)")
    print(f"Study domain: {LAT_MIN} to {LAT_MAX}N, {LON_MIN} to {LON_MAX}E")
    print("="*60)

    all_results = []

    for lat_str, lon_str in TREECOVER_TILES:
        print(f"\n--- Tile {lat_str}_{lon_str} ---")

        result = process_tile(lat_str, lon_str)
        if result is None:
            continue

        forest, cleared, transform, crs = result

        # Compute distance to deforestation edge
        dist_km = compute_distance_to_edge(forest, cleared)

        # Aggregate to SMAP grid
        agg = aggregate_to_smap_grid(forest, cleared, dist_km, transform)
        all_results.append(agg)

        # Free memory before next tile
        del forest, cleared, dist_km
        print(f"  Tile done. SMAP cells: {len(agg['lat'])} x {len(agg['lon'])}")

    if not all_results:
        print("ERROR: No tiles processed!")
        sys.exit(1)

    # Merge tile results into single grid
    print("\nMerging tiles...")
    # Combine all lat/lon values
    all_lats = np.unique(np.concatenate([r['lat'] for r in all_results]))
    all_lons = np.unique(np.concatenate([r['lon'] for r in all_results]))
    all_lats = np.sort(all_lats)
    all_lons = np.sort(all_lons)

    merged = {
        'lat': all_lats,
        'lon': all_lons,
        'frac_forest': np.full((len(all_lats), len(all_lons)), np.nan),
        'frac_cleared': np.full((len(all_lats), len(all_lons)), np.nan),
        'mean_dist_km': np.full((len(all_lats), len(all_lons)), np.nan),
        'median_dist_km': np.full((len(all_lats), len(all_lons)), np.nan),
        'n_forest_px': np.zeros((len(all_lats), len(all_lons)), dtype=np.int32),
    }

    for r in all_results:
        for i, lat_val in enumerate(r['lat']):
            for j, lon_val in enumerate(r['lon']):
                gi = np.argmin(np.abs(all_lats - lat_val))
                gj = np.argmin(np.abs(all_lons - lon_val))

                # If this cell has more forest pixels, prefer its values
                if r['n_forest_px'][i, j] > merged['n_forest_px'][gi, gj]:
                    for key in ['frac_forest', 'frac_cleared', 'mean_dist_km', 'median_dist_km']:
                        merged[key][gi, gj] = r[key][i, j]
                    merged['n_forest_px'][gi, gj] = r['n_forest_px'][i, j]

    n_valid = np.sum(np.isfinite(merged['frac_forest']))
    print(f"Merged grid: {len(all_lats)} x {len(all_lons)}, {n_valid} valid cells")

    # Save
    save_netcdf(merged, OUT_NC)

    # Summary stats
    ff = merged['frac_forest']
    fc = merged['frac_cleared']
    md = merged['mean_dist_km']
    valid = np.isfinite(md)

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Grid cells with data: {np.sum(valid)}")
    print(f"Forest fraction: {np.nanmean(ff):.3f} (mean), {np.nanmedian(ff):.3f} (median)")
    print(f"Cleared fraction: {np.nanmean(fc):.3f} (mean)")
    print(f"Distance to edge: {np.nanmean(md):.1f} km (mean), {np.nanmedian(md):.1f} km (median)")
    print(f"  Min: {np.nanmin(md):.1f} km")
    print(f"  Max: {np.nanmax(md):.1f} km")
    print(f"  P25: {np.nanpercentile(md[valid], 25):.1f} km")
    print(f"  P75: {np.nanpercentile(md[valid], 75):.1f} km")
    print("\nDone.")
