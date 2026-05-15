# Paper 2: Edge Effects and Minimum Viable Restoration Scales

## Research gaps addressed

### Gap 2: Minimum viable restoration scale per hydrological lever
- Literature quantifies deforestation thresholds (~50km crossover, ~40% basin tipping point)
- But NO paper quantifies the minimum RESTORATION patch size needed to recover hydrological function
- Staal et al. (2020), Zemp et al. (2017), Spracklen et al. (2012) all work from deforestation side
- The restoration side (how much forest do you need to REBUILD) is unquantified

### Gap 4: Deforestation edge effects on adjacent intact forest dry-down
- Laurance et al. (2002) documented edge effects at 100-300m (microclimate, tree mortality)
- Garcia-Carreras & Parker (2011) modeled mesoscale circulation changes
- But NO observational study has measured whether deforestation edges accelerate soil moisture dry-down in adjacent intact forest at mesoscale (9km SMAP resolution)
- Our SMAP residence time methodology from Paper 1 can test this directly

## Key hypotheses
1. Intact forest pixels adjacent to deforestation edges dry faster than deep interior pixels
2. The edge effect penetration depth is measurable at 9km SMAP resolution
3. Larger contiguous clearings produce faster dry-down in adjacent forest than smaller ones
4. There exists a minimum restoration patch size below which hydrological function does not recover

## Data
- SMAP L4 root-zone SM (9km, 2015-2024) — already have from Paper 1
- MODIS MCD12Q1 land cover (500m, annual) — need to download
- ERA5 moisture flux fields — already have from Paper 1

## Relation to Paper 1
- Paper 1 established the dry-down methodology and the 14 vs 9 month residence time finding
- Paper 2 extends this spatially: instead of two boxes (arc vs intact), we resolve the gradient
- Paper 2 answers the operational questions ARARA raised
