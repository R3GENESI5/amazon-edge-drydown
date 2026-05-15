"""
Download Amazon deforestation data for Paper 2.
Logs output to download_log.txt.
"""
import os, sys, time, requests

LOG = open(r"D:\amazon paper 2\scripts\download_log.txt", "w", buffering=1)
def log(msg):
    LOG.write(msg + "\n")
    LOG.flush()

DATA = r"D:\amazon paper 2\data"
PRODES = os.path.join(DATA, "prodes")
HANSEN = os.path.join(DATA, "hansen")
BOUNDARY = os.path.join(DATA, "boundaries")
for d in [PRODES, HANSEN, BOUNDARY]:
    os.makedirs(d, exist_ok=True)

def dl(url, dest, desc="", timeout=300):
    log(f"  DL: {desc}")
    log(f"    URL: {url}")
    try:
        r = requests.get(url, stream=True, timeout=timeout,
                        headers={"User-Agent": "Mozilla/5.0"}, allow_redirects=True)
        r.raise_for_status()
        total = int(r.headers.get('content-length', 0))
        got = 0
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=131072):
                f.write(chunk)
                got += len(chunk)
        sz = os.path.getsize(dest) / (1024*1024)
        log(f"    OK: {sz:.1f} MB")
        return True
    except Exception as e:
        log(f"    FAIL: {e}")
        if os.path.exists(dest):
            os.remove(dest)
        return False

log("=== OPTION 1: PRODES ===")
prodes_ok = False
prodes_urls = [
    ("http://terrabrasilis.dpi.inpe.br/download/dataset/legal-amz-prodes/vector/yearly_deforestation.zip",
     os.path.join(PRODES, "prodes_yearly_deforestation.zip"), "PRODES yearly deforestation shp"),
    ("http://terrabrasilis.dpi.inpe.br/download/dataset/legal-amz-prodes/vector/accumulated_deforestation.zip",
     os.path.join(PRODES, "prodes_accumulated_deforestation.zip"), "PRODES accumulated deforestation"),
    ("http://terrabrasilis.dpi.inpe.br/file-delivery/download/prodes-amz/vector/yearly_deforestation_biome.zip",
     os.path.join(PRODES, "prodes_yearly_biome.zip"), "PRODES yearly biome"),
]
for url, dest, desc in prodes_urls:
    if dl(url, dest, desc, timeout=120):
        prodes_ok = True
        break
    time.sleep(2)

if not prodes_ok:
    log("  Trying WFS...")
    wfs = ("http://terrabrasilis.dpi.inpe.br/geoserver/prodes-amz/ows?"
           "service=WFS&version=1.0.0&request=GetFeature"
           "&typeName=prodes-amz:yearly_deforestation"
           "&outputFormat=SHAPE-ZIP&maxFeatures=50000")
    prodes_ok = dl(wfs, os.path.join(PRODES, "prodes_wfs.zip"), "PRODES WFS", 180)

log(f"PRODES result: {'OK' if prodes_ok else 'FAILED'}")

log("\n=== OPTION 3: HANSEN GFC ===")
log("(Skipping Option 2 MapBiomas - requires GEE auth)")
base = "https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11"
lats = ["10N", "00N", "10S", "20S"]
lons = ["080W", "070W", "060W", "050W"]
hansen_ok = 0
hansen_total = len(lats) * len(lons)

for lat in lats:
    for lon in lons:
        fn = f"Hansen_GFC-2023-v1.11_lossyear_{lat}_{lon}.tif"
        dest = os.path.join(HANSEN, fn)
        if os.path.exists(dest) and os.path.getsize(dest) > 1000:
            log(f"  Already have: {fn}")
            hansen_ok += 1
            continue
        if dl(f"{base}/{fn}", dest, fn, 600):
            hansen_ok += 1
        time.sleep(1)

log(f"Hansen lossyear tiles: {hansen_ok}/{hansen_total}")

# Also get treecover2000 for core tiles
log("\nDownloading treecover2000 for core Amazon tiles...")
core = [("00N","060W"),("00N","070W"),("10S","060W"),("10S","070W")]
tc_ok = 0
for lat, lon in core:
    fn = f"Hansen_GFC-2023-v1.11_treecover2000_{lat}_{lon}.tif"
    dest = os.path.join(HANSEN, fn)
    if os.path.exists(dest) and os.path.getsize(dest) > 1000:
        log(f"  Already have: {fn}")
        tc_ok += 1
        continue
    if dl(f"{base}/{fn}", dest, fn, 600):
        tc_ok += 1
    time.sleep(1)
log(f"Hansen treecover2000 tiles: {tc_ok}/{len(core)}")

log("\n=== AMAZON BASIN BOUNDARY ===")
bnd_ok = False
bnd_urls = [
    ("https://naciscdn.org/naturalearth/10m/physical/ne_10m_rivers_lake_centerlines.zip",
     os.path.join(BOUNDARY, "ne_10m_rivers.zip"), "NE 10m rivers"),
    ("https://naciscdn.org/naturalearth/110m/cultural/ne_110m_admin_0_countries.zip",
     os.path.join(BOUNDARY, "ne_110m_countries.zip"), "NE 110m countries"),
    ("https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip",
     os.path.join(BOUNDARY, "ne_10m_countries.zip"), "NE 10m countries"),
]
for url, dest, desc in bnd_urls:
    if dl(url, dest, desc, 120):
        bnd_ok = True

# Also try HydroSHEDS basins
dl("https://data.hydrosheds.org/file/HydroBASINS/standard/hybas_sa_lev01-06_v1c.zip",
   os.path.join(BOUNDARY, "hybas_sa_lev01-06.zip"), "HydroSHEDS SA basins", 300)

log("\n=== DOWNLOAD REPORT ===")
for dname in [PRODES, HANSEN, BOUNDARY]:
    label = os.path.relpath(dname, DATA)
    files = []
    if os.path.exists(dname):
        for f in os.listdir(dname):
            fp = os.path.join(dname, f)
            if os.path.isfile(fp):
                files.append((f, os.path.getsize(fp)/(1024*1024)))
    log(f"\n  {label}/")
    if files:
        for fn, sz in sorted(files):
            log(f"    {fn}: {sz:.1f} MB")
        log(f"    Total: {sum(s for _,s in files):.1f} MB")
    else:
        log(f"    (empty)")

log(f"\nSUMMARY: PRODES={'OK' if prodes_ok else 'FAIL'}, Hansen={hansen_ok}/{hansen_total} tiles, Boundary={'OK' if bnd_ok else 'FAIL'}")
log("DONE")
LOG.close()
