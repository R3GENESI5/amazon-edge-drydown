# Data Sources and Download Instructions

Raw datasets are not included in this repository (combined size ~3 GB). They are publicly available from the sources below.

---

## SMAP L4 v7 root-zone soil moisture
- **Variable:** Monthly root-zone soil moisture (m³/m³), 9 km EASE-Grid 2.0
- **Period:** April 2015 – December 2024 (117 months)
- **Source:** NASA Earthdata / NSIDC
- **Product:** SPL4SMGP (Geophysical Variables) — monthly average
- **URL:** https://nsidc.org/data/spl4smgp
- **Reference:** Reichle et al. (2017) doi:10.1175/JHM-D-16-0291.1
- **Place in:** `data/smap/`

## Hansen Global Forest Change v1.11
- **Variables:** `treecover2000` (% canopy cover, 2000), `lossyear` (year of forest loss, 1=2001 to 23=2023)
- **Resolution:** 30 m
- **Source:** University of Maryland GLAD, hosted on Google Earth Engine and storage.googleapis.com
- **URL:** https://storage.googleapis.com/earthenginepartners-hansen/GFC-2023-v1.11
- **Tiles needed (core Amazon):** `00N_060W`, `00N_070W`, `10S_060W`, `10S_070W`
- **Reference:** Hansen et al. (2013) doi:10.1126/science.1244693
- **Place in:** `data/hansen/`

## HydroSHEDS HydroBASINS Level 1
- **Variable:** Amazon basin polygon (vector)
- **Source:** WWF / HydroSHEDS
- **URL:** https://www.hydrosheds.org/products/hydrobasins
- **Reference:** Lehner et al. (2008) doi:10.1029/2008EO100001
- **Place in:** `data/boundaries/`

---

## File organisation expected by scripts

```
data/
├── smap/
│   └── SMAP_L4_SM_gph_M.YYYYMM.h5     # 117 monthly files
├── hansen/
│   ├── Hansen_GFC-2023-v1.11_treecover2000_00N_060W.tif
│   ├── Hansen_GFC-2023-v1.11_lossyear_00N_060W.tif
│   └── ... (4 tiles × 2 layers = 8 files)
└── boundaries/
    └── hybas_sa_lev01_v1c.shp
```
