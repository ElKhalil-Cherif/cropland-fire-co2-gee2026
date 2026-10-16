# Open-source crop-specific model of CO₂ emissions from cropland fires (Google Earth Engine)

Code and results for the paper:

> Cherif, E.K.; Oulad Dahman, C. *An Open-Source, Crop-Specific Model for Estimating CO₂ Emissions from Cropland Fires in Google Earth Engine: Application to the Lower Mississippi River Basin.* (submitted, 2026)

The model links every burned cropland pixel to its crop type (USDA Cropland Data Layer) and converts it into CO₂ with the IPCC bottom-up equation

    E = A × Y × RPR × DM × CF × EF × 10⁻⁶

where A is the burned area (ha), Y is the FAOSTAT crop yield (kg ha⁻¹), RPR is the residue-to-product ratio, DM is the dry-matter fraction, CF is the combustion factor and EF is the emission factor (1441 g CO₂ kg⁻¹ dry matter). Uncertainty is propagated with a 10,000-draw Monte Carlo simulation.

## Repository content

| Path | What it does |
|---|---|
| `gee/01_CF2020_all_in_one.js` | Main GEE script. Part 1: whole basin, MODIS MCD64A1, annual and monthly CO₂ by crop, VIIRS consistency check (prints `RESULTS_JSON` and `MAP_THUMBNAIL_URL`). Part 2: Testing Zone (a), Sep–Nov, Sentinel-2 + VIIRS vs MCD64A1 vs FireCCI51 (export task `CF2020_compare_zoneA`). |
| `gee/02_zoneA_comparison_standalone.js` | Part 2 alone, with yields and emission coefficients written in the code (no yield-table asset needed). |
| `gee/03_crop_type_map_2020.js` | Crop type of the MCD64A1 burned cropland (prints `CROP_MAP_URL`, used for Fig. 3b). |
| `gee/tiled_full_basin/` | Two-stage version of the Sentinel-2 method (stage 1 exports monthly burn masks per Sentinel-2 tile, stage 2 computes emissions). Needed for large areas; set the asset paths to your own project. |
| `python/monte_carlo.py` | Monte Carlo uncertainty (run outside GEE because of the free-tier memory limit). Usage: `python monte_carlo.py results.json` |
| `python/figures_maps.py` | Fig. 1 (study area) and Fig. 3 (emission intensity and crop type). |
| `python/figure_flowchart.py` | Fig. 2 (model structure). |
| `results/basin_2020_MCD64A1_results.json` | Basin results for 2020 (annual and monthly area and CO₂ by crop, yields, coefficients, VIIRS check, Monte Carlo). |
| `results/zoneA_Sep-Nov_2020_comparison.csv` | Testing Zone (a) comparison output. Each record is `month|CDL crop id|ha|t CO2`. |

## How to run

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/) and paste a script from `gee/`.
2. Press **Run**. Results are printed in the Console; heavy computations are exported from the **Tasks** tab to Google Drive.
3. Python scripts: `pip install -r requirements.txt`. The map script needs the Natural Earth GeoJSON files `ne_50m_land`, `ne_50m_admin_0_countries`, `ne_50m_admin_1_states_provinces_lakes`, `ne_10m_lakes` and `ne_10m_rivers_lake_centerlines` (from https://github.com/nvkelso/natural-earth-vector, folder `geojson/`) in a folder `python/naturalearth/`.

### Input data

All satellite and crop data are read from the public GEE catalog: `MODIS/061/MCD64A1`, `ESA/CCI/FireCCI/5_1`, `COPERNICUS/S2_SR_HARMONIZED`, `USDA/NASS/CDL/2020`.
Three inputs are user assets of the authors and must be readable (or replaced by your own copies):

- VIIRS 375 m active fires 2020: `projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/NOAA20_2020` and `.../SUIOMI2020` (from NASA FIRMS)
- FAOSTAT yields: `projects/ee-cherifelkhalil-etu/assets/crop_yield_table` (from https://www.fao.org/faostat/)
- Sentinel-2 tile grid: `projects/ee-cherifelkhalil-etu/assets/S2/sentinel_2_ID_Tiles` (script 01 only; script 02 reads the tiles from the Sentinel-2 collection)

### Applying the model elsewhere

Change the area of interest (`AOI_COORDS`), the year, the crop map and the parameter table (`PARAMS`, yields). The burned-area module can use MCD64A1, FireCCI51 or the Sentinel-2 method.

## Known issue

In the runs reported in the paper, the GEE adaptation of the Sentinel-2 + VIIRS method returned no burned pixels in Testing Zone (a), although MCD64A1 and FireCCI51 both mapped burned cropland there. The cause is being investigated. The original, validated Python version of the field-scale method is at https://github.com/ElKhalil-Cherif/burned-area-cropland.

## Related work

- Oulad Dahman, C.; Cherif, E.K.; Boulaassal, H. Estimating cropland fire emissions with remote sensing: A literature review and data strategy. *E3S Web Conf.* 2025, 676, 03002. https://doi.org/10.1051/e3sconf/202567603002
- Cherif, E.K.; Oulad Dahman, C.; Bernardino, A. Automated high-resolution mapping of cropland burned areas in the Lower Mississippi River Basin using VIIRS and Sentinel-2 imagery. IGARSS 2026.

## License

MIT (see `LICENSE`).
