# Standard domain definitions for all analyses
# Amazon ALLJ Paper - Canonical regions
#
# Created to resolve audit finding: inconsistent spatial domains across scripts.
# All scripts should import from this file for consistency.

STUDY_DOMAIN = dict(lat_min=-12, lat_max=6, lon_min=-75, lon_max=-40)
ARC_DEFORESTATION = dict(lat_min=-12, lat_max=-5, lon_min=-55, lon_max=-45)
INTACT_INTERIOR = dict(lat_min=-5, lat_max=0, lon_min=-65, lon_max=-55)
TRANSECT_BAND = dict(lat_min=-5, lat_max=0)  # 0-5S for flux transects
CAPE_REGION = dict(lat_min=-8, lat_max=-2, lon_min=-60, lon_max=-50)
FLUX_TRANSECT_LON = -52  # longitude for flux transect
