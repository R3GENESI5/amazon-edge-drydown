"""
Extract Amazon basin boundary from HydroSHEDS HydroBASINS data.
Also creates a simple bounding-box shapefile as fallback.
Logs to create_boundary_log.txt.
"""
import os, sys, zipfile
LOG = open(r"D:\amazon paper 2\scripts\create_boundary_log.txt", "w", buffering=1)
def log(msg):
    LOG.write(msg + "\n"); LOG.flush()

BOUNDARY = r"D:\amazon paper 2\data\boundaries"
hybas_zip = os.path.join(BOUNDARY, "hybas_sa_lev01-06.zip")

# Step 1: Extract HydroBASINS zip
log("Step 1: Extracting HydroBASINS...")
hybas_dir = os.path.join(BOUNDARY, "hybas_sa")
if os.path.exists(hybas_zip):
    os.makedirs(hybas_dir, exist_ok=True)
    with zipfile.ZipFile(hybas_zip, 'r') as z:
        names = z.namelist()
        log(f"  Contents: {len(names)} files")
        for n in names:
            log(f"    {n}")
        z.extractall(hybas_dir)
    log("  Extracted OK")
else:
    log(f"  {hybas_zip} not found!")

# Step 2: Try to extract Amazon basin using geopandas
log("\nStep 2: Extract Amazon basin polygon...")
try:
    import geopandas as gpd
    import glob
    
    # Find level-1 shapefile (largest basins)
    shps = glob.glob(os.path.join(hybas_dir, "**", "*lev01*.shp"), recursive=True)
    if not shps:
        shps = glob.glob(os.path.join(hybas_dir, "**", "*.shp"), recursive=True)
    log(f"  Found shapefiles: {shps}")
    
    if shps:
        # Level 1 has the major basins; Amazon is HYBAS_ID 6010000010
        # or we can find it by area (largest in South America)
        for shp in shps:
            gdf = gpd.read_file(shp)
            log(f"  {os.path.basename(shp)}: {len(gdf)} features, cols={list(gdf.columns)}")
            if 'HYBAS_ID' in gdf.columns:
                # Amazon basin at level 1
                amazon = gdf[gdf['HYBAS_ID'] == 6010000010]
                if len(amazon) == 0:
                    # Try finding by largest area
                    gdf['area_calc'] = gdf.geometry.area
                    amazon = gdf.nlargest(1, 'area_calc')
                    log(f"  Selected largest basin: HYBAS_ID={amazon.iloc[0].get('HYBAS_ID','?')}")
                
                if len(amazon) > 0:
                    out = os.path.join(BOUNDARY, "amazon_basin.shp")
                    amazon.to_file(out)
                    log(f"  Saved Amazon basin to {out}")
                    break
        
except ImportError as e:
    log(f"  geopandas not available: {e}")
    log("  Will create a simple bounding box instead.")

# Step 3: Create a simple Amazon bounding box as fallback
log("\nStep 3: Creating Amazon bounding box shapefile...")
try:
    import geopandas as gpd
    from shapely.geometry import box
    
    # Approximate Amazon basin bounding box
    # Covers most of the basin: ~5N to 20S, 80W to 44W
    amazon_bbox = box(-80, -20, -44, 5)
    gdf_bbox = gpd.GeoDataFrame(
        {'name': ['Amazon_bbox'], 'geometry': [amazon_bbox]},
        crs='EPSG:4326'
    )
    out_bbox = os.path.join(BOUNDARY, "amazon_bbox.shp")
    gdf_bbox.to_file(out_bbox)
    log(f"  Saved bounding box to {out_bbox}")
except Exception as e:
    log(f"  Failed to create bbox shapefile: {e}")

# Step 4: Delete the bogus PRODES file
log("\nStep 4: Cleaning up bogus PRODES file...")
prodes_bogus = r"D:\amazon paper 2\data\prodes\prodes_yearly_biome.zip"
if os.path.exists(prodes_bogus) and os.path.getsize(prodes_bogus) < 100:
    os.remove(prodes_bogus)
    log(f"  Removed {prodes_bogus} (was only {2} bytes)")

log("\nDONE")
LOG.close()
