"""
Download Amazon deforestation boundary data for Paper 2.
Tries PRODES first, then MapBiomas, then Hansen Global Forest Change.
Also downloads Amazon basin boundary from Natural Earth.
"""
import os
import sys
import requests
import time
import json

# Directories
DATA_DIR = r"D:\amazon paper 2\data"
PRODES_DIR = os.path.join(DATA_DIR, "prodes")
MAPBIOMAS_DIR = os.path.join(DATA_DIR, "mapbiomas")
HANSEN_DIR = os.path.join(DATA_DIR, "hansen")
BOUNDARY_DIR = os.path.join(DATA_DIR, "boundaries")

for d in [PRODES_DIR, MAPBIOMAS_DIR, HANSEN_DIR, BOUNDARY_DIR]:
    os.makedirs(d, exist_ok=True)

def download_file(url, dest_path, desc="", timeout=300):
    """Download a file with progress reporting."""
    print(f"  Downloading: {desc or url}")
    print(f"  URL: {url}")
    print(f"  Dest: {dest_path}")
    try:
        resp = requests.get(url, stream=True, timeout=timeout, 
                           headers={"User-Agent": "Mozilla/5.0"}, allow_redirects=True)
        resp.raise_for_status()
        total = int(resp.headers.get('content-length', 0))
        downloaded = 0
        with open(dest_path, 'wb') as f:
            for chunk in resp.iter_content(chunk_size=8192*16):
                f.write(chunk)
                downloaded += len(chunk)
                if total > 0 and downloaded % (1024*1024*10) < 8192*16:
                    pct = downloaded / total * 100
                    print(f"    {downloaded/(1024*1024):.1f} MB / {total/(1024*1024):.1f} MB ({pct:.0f}%)")
        size_mb = os.path.getsize(dest_path) / (1024*1024)
        print(f"  SUCCESS: {size_mb:.1f} MB downloaded")
        return True
    except Exception as e:
        print(f"  FAILED: {e}")
        if os.path.exists(dest_path):
            os.remove(dest_path)
        return False

def try_prodes():
    """OPTION 1: Try PRODES from TerraBrasilis/INPE."""
    print("\n" + "="*60)
    print("OPTION 1: PRODES (INPE/TerraBrasilis)")
    print("="*60)
    
    # TerraBrasilis provides PRODES data via their API and direct downloads
    # The accumulated deforestation shapefile is the key product
    urls = [
        # PRODES accumulated deforestation - shapefile via TerraBrasilis download service
        ("http://terrabrasilis.dpi.inpe.br/download/dataset/legal-amz-prodes/vector/yearly_deforestation.zip",
         os.path.join(PRODES_DIR, "prodes_yearly_deforestation.zip"),
         "PRODES yearly deforestation shapefile"),
        # Alternative: accumulated deforestation 
        ("http://terrabrasilis.dpi.inpe.br/download/dataset/legal-amz-prodes/vector/accumulated_deforestation.zip",
         os.path.join(PRODES_DIR, "prodes_accumulated_deforestation.zip"),
         "PRODES accumulated deforestation shapefile"),
        # Alternative direct link patterns
        ("http://terrabrasilis.dpi.inpe.br/file-delivery/download/prodes-amz/vector/yearly_deforestation_biome.zip",
         os.path.join(PRODES_DIR, "prodes_yearly_biome.zip"),
         "PRODES yearly deforestation biome"),
    ]
    
    success = False
    for url, dest, desc in urls:
        if download_file(url, dest, desc, timeout=120):
            success = True
            break
        time.sleep(2)
    
    if not success:
        # Try the TerraBrasilis GeoServer WFS for a smaller extract
        print("\n  Trying TerraBrasilis GeoServer WFS...")
        wfs_url = (
            "http://terrabrasilis.dpi.inpe.br/geoserver/prodes-amz/ows?"
            "service=WFS&version=1.0.0&request=GetFeature"
            "&typeName=prodes-amz:yearly_deforestation"
            "&outputFormat=SHAPE-ZIP&maxFeatures=50000"
        )
        dest = os.path.join(PRODES_DIR, "prodes_wfs_deforestation.zip")
        if download_file(wfs_url, dest, "PRODES WFS yearly deforestation", timeout=180):
            success = True
    
    return success

def try_mapbiomas():
    """OPTION 2: Try MapBiomas land cover."""
    print("\n" + "="*60)
    print("OPTION 2: MapBiomas")
    print("="*60)
    print("  MapBiomas requires authentication/API access.")
    print("  Checking direct download links...")
    
    # MapBiomas data is typically accessed via Google Earth Engine or their platform
    # Direct downloads are available for some products
    urls = [
        # MapBiomas Collection 9 - Amazon biome
        ("https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_9/lclu/coverage/brasil_coverage_2022.tif",
         os.path.join(MAPBIOMAS_DIR, "mapbiomas_brasil_2022.tif"),
         "MapBiomas Collection 9 Brazil 2022"),
        # Alternative URL patterns
        ("https://storage.googleapis.com/mapbiomas-public/brasil/collection-8/lclu/coverage/brasil_coverage_2022.tif",
         os.path.join(MAPBIOMAS_DIR, "mapbiomas_c8_brasil_2022.tif"),
         "MapBiomas Collection 8 Brazil 2022"),
    ]
    
    for url, dest, desc in urls:
        if download_file(url, dest, desc, timeout=600):
            return True
        time.sleep(2)
    
    return False

def try_hansen():
    """OPTION 3: Hansen Global Forest Change (UMD)."""
    print("\n" + "="*60)
    print("OPTION 3: Hansen Global Forest Change (UMD)")
    print("="*60)
    
    # Hansen data is organized as 10x10 degree tiles
    # Amazon basin spans roughly 10N to 20S, 80W to 40W
    # File naming: Hansen_GFC-2023-v1.11_lossyear_00N_080W.tif
    # We need treecover2000 and lossyear layers
    
    base_url = "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11"
    
    # Define tiles covering the Amazon
    # Latitude bands: 10N, 00N, 10S, 20S
    # Longitude bands: 80W, 70W, 60W, 50W, 40W
    lat_tiles = ["10N", "00N", "10S", "20S"]
    lon_tiles = ["080W", "070W", "060W", "050W"]
    
    # First just download the lossyear layer (smaller, more useful for our purpose)
    # and treecover2000 for baseline forest extent
    layers = ["lossyear"]  # Start with lossyear only
    
    success_count = 0
    total_tiles = len(lat_tiles) * len(lon_tiles) * len(layers)
    
    for layer in layers:
        for lat in lat_tiles:
            for lon in lon_tiles:
                fname = f"Hansen_GFC-2023-v1.11_{layer}_{lat}_{lon}.tif"
                url = f"{base_url}/{fname}"
                dest = os.path.join(HANSEN_DIR, fname)
                
                if os.path.exists(dest) and os.path.getsize(dest) > 1000:
                    print(f"  Already exists: {fname}")
                    success_count += 1
                    continue
                
                if download_file(url, dest, fname, timeout=300):
                    success_count += 1
                else:
                    print(f"  Skipping tile {fname}")
                time.sleep(1)
    
    print(f"\n  Downloaded {success_count}/{total_tiles} tiles")
    
    # Also try to get treecover2000 for at least the core Amazon tiles
    if success_count > 0:
        print("\n  Now downloading treecover2000 for core Amazon tiles...")
        core_tiles = [("00N", "060W"), ("00N", "070W"), ("10S", "060W"), ("10S", "070W")]
        for lat, lon in core_tiles:
            fname = f"Hansen_GFC-2023-v1.11_treecover2000_{lat}_{lon}.tif"
            url = f"{base_url}/{fname}"
            dest = os.path.join(HANSEN_DIR, fname)
            if os.path.exists(dest) and os.path.getsize(dest) > 1000:
                print(f"  Already exists: {fname}")
                continue
            download_file(url, dest, fname, timeout=300)
            time.sleep(1)
    
    return success_count > 0

def download_amazon_boundary():
    """Download Amazon basin boundary shapefile."""
    print("\n" + "="*60)
    print("AMAZON BASIN BOUNDARY")
    print("="*60)
    
    # Try multiple sources for Amazon basin boundary
    urls = [
        # HydroBASINS from HydroSHEDS - Amazon basin level 1
        # Natural Earth large-scale rivers/lakes
        ("https://naciscdn.org/naturalearth/10m/physical/ne_10m_rivers_lake_centerlines.zip",
         os.path.join(BOUNDARY_DIR, "ne_10m_rivers_lake_centerlines.zip"),
         "Natural Earth rivers (contains Amazon)"),
        # WWF HydroSHEDS basins
        ("https://data.hydrosheds.org/file/HydroBASINS/standard/hybas_sa_lev01-06_v1c.zip",
         os.path.join(BOUNDARY_DIR, "hybas_sa_lev01-06.zip"),
         "HydroSHEDS South America basins L1-6"),
    ]
    
    success = False
    for url, dest, desc in urls:
        if download_file(url, dest, desc, timeout=300):
            success = True
    
    # Also try a simpler approach: download the Amazon basin from a known source
    # The Amazon basin outline from GRDC or similar
    simple_urls = [
        ("https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip",
         os.path.join(BOUNDARY_DIR, "ne_110m_countries.zip"),
         "Natural Earth countries (for clipping)"),
        ("https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip",
         os.path.join(BOUNDARY_DIR, "ne_10m_countries.zip"),
         "Natural Earth 10m countries"),
    ]
    
    for url, dest, desc in simple_urls:
        if download_file(url, dest, desc, timeout=120):
            success = True
    
    return success

def report_downloads():
    """Report what was downloaded."""
    print("\n" + "="*60)
    print("DOWNLOAD REPORT")
    print("="*60)
    
    for dirname in [PRODES_DIR, MAPBIOMAS_DIR, HANSEN_DIR, BOUNDARY_DIR]:
        files = []
        if os.path.exists(dirname):
            for f in os.listdir(dirname):
                fp = os.path.join(dirname, f)
                if os.path.isfile(fp):
                    size_mb = os.path.getsize(fp) / (1024*1024)
                    files.append((f, size_mb))
        
        rel = os.path.relpath(dirname, DATA_DIR)
        print(f"\n  {rel}/")
        if files:
            for fname, size in sorted(files):
                print(f"    {fname}: {size:.1f} MB")
            total = sum(s for _, s in files)
            print(f"    --- Total: {total:.1f} MB")
        else:
            print(f"    (empty)")

if __name__ == "__main__":
    print("Amazon Deforestation Data Downloader")
    print("For Paper 2: Lateral Moisture Extraction")
    print(f"Python: {sys.executable}")
    print(f"Data dir: {DATA_DIR}")
    
    # Try PRODES first
    prodes_ok = try_prodes()
    
    # Try MapBiomas if PRODES failed
    mapbiomas_ok = False
    if not prodes_ok:
        mapbiomas_ok = try_mapbiomas()
    
    # Try Hansen if both failed
    hansen_ok = False
    if not prodes_ok and not mapbiomas_ok:
        hansen_ok = try_hansen()
    elif not prodes_ok:
        # Even if MapBiomas worked, Hansen is very reliable - get it too
        print("\n  Skipping Hansen (MapBiomas succeeded)")
    else:
        print("\n  Skipping Hansen (PRODES succeeded)")
    
    # If nothing worked from options 1-2, definitely try Hansen
    if not prodes_ok and not mapbiomas_ok and not hansen_ok:
        print("\n  WARNING: All primary options failed!")
        print("  Retrying Hansen with individual tile checks...")
        hansen_ok = try_hansen()
    
    # Always try to get Amazon boundary
    boundary_ok = download_amazon_boundary()
    
    # Report
    report_downloads()
    
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"  PRODES:     {'OK' if prodes_ok else 'FAILED'}")
    print(f"  MapBiomas:  {'OK' if mapbiomas_ok else 'SKIPPED' if prodes_ok else 'FAILED'}")
    print(f"  Hansen:     {'OK' if hansen_ok else 'SKIPPED' if (prodes_ok or mapbiomas_ok) else 'FAILED'}")
    print(f"  Boundary:   {'OK' if boundary_ok else 'FAILED'}")
    print("\nDone.")
