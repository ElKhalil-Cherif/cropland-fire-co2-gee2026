# Open-source crop-specific model of CO₂ emissions from cropland fires (Google Earth Engine)

Code and results for the paper:

> Cherif, E.K.; Oulad Dahman, C. *An Open-Source, Crop-Specific Model for Estimating CO₂ Emissions from Cropland Fires in Google Earth Engine: Application to the Lower Mississippi River Basin.* (submitted, 2026)

The model links every burned cropland pixel to its crop type (USDA Cropland Data Layer) and converts it into CO₂ with the IPCC bottom-up equation

    E = A × Y × RPR × DM × CF × EF × 10⁻⁶

where A is the burned area (ha), Y the FAOSTAT crop yield (kg ha⁻¹), RPR the residue-to-product ratio, DM the dry-matter fraction, CF the combustion factor and EF the emission factor (1441 g CO₂ kg⁻¹ dry matter). Uncertainty is propagated with a 10,000-draw Monte Carlo simulation.

The burned-area module has two tiers:

- **Tier 1** – MODIS MCD64A1 (500 m), whole basin, January–December 2020.
- **Tier 2** – field-scale Sentinel-2 method guided by VIIRS 375 m active fires (30 m), compared with MCD64A1, FireCCI51 and GFED4.1s over the same sub-area (basin north of 33.75° N) and period (September–November 2020).

## Files

| File | Purpose |
|---|---|
| `01_basin_MCD64A1_emissions.js` | GEE script, Tier 1: annual and monthly CO₂ by crop for the whole basin and VIIRS consistency check (Console: `RESULTS_JSON`, `MAP_THUMBNAIL_URL`). |
| `02_subarea_S2_vs_MODIS_FireCCI.js` | GEE script, Tier 2: Sentinel-2 + VIIRS burned area and MCD64A1 and FireCCI51 on the same 30 m crop grid, by crop and month, with spatial overlaps. Exports one CSV per month (Tasks tab). `REGION` selects the sub-area (`NORTH_HALF` in the paper). |
| `03_crop_type_map_2020.js` | GEE script: crop type of the MCD64A1 burned cropland (Fig. 3b, `CROP_MAP_URL`). |
| `monte_carlo.py` | Monte Carlo uncertainty, run outside GEE: `python monte_carlo.py results_basin_2020_MCD64A1.json` |
| `analyse_comparison_gfed.py` | Sums the three Tier 2 CSVs by crop and month and extracts GFED4.1s agricultural CO₂ for the same sub-area and months: `python analyse_comparison_gfed.py results_subarea_S2_MODIS_FireCCI_m09.csv ..._m10.csv ..._m11.csv GFED4.1s_2020_beta.hdf5` |
| `aoi_coords.json` | Study-area polygon (used by the GFED extraction). |
| `figures_maps.py`, `figure_flowchart.py`, `figure_comparison.py` | Figs. 1 and 3, Fig. 2 and Fig. 6. |
| `results_basin_2020_MCD64A1.json` | Tier 1 results (annual and monthly area and CO₂ by crop, yields, coefficients, VIIRS check). |
| `results_subarea_S2_MODIS_FireCCI_m09/m10/m11.csv` | Tier 2 results used in the paper (original workflow order: NDVI vegetation filter, NBR/MIRBI/BAIS2 Otsu consensus, then 500 m VIIRS mask). One row per CDL crop: burned area (ha) of each product, pairwise overlaps and t CO₂. |
| `results_subarea_summary.json` | Totals by crop and month, ratios and GFED4.1s monthly agricultural dry matter and CO₂. |
| `sensitivity_3km_4of5_m09/m10/m11.csv` | Sensitivity run (3 km VIIRS buffer applied before thresholding, no NDVI filter, 4 of 5 indices incl. NDTI and SAVI), reported in Section 3.4; it flags harvested fields as burned. |

## How to run

1. Open the [Google Earth Engine Code Editor](https://code.earthengine.google.com/), paste a `.js` script and press **Run**. Results appear in the Console; script 02 exports its CSV files from the **Tasks** tab to Google Drive.
2. Python: `pip install -r requirements.txt`. The map script needs the Natural Earth GeoJSON files `ne_50m_land`, `ne_50m_admin_0_countries`, `ne_50m_admin_1_states_provinces_lakes`, `ne_10m_lakes` and `ne_10m_rivers_lake_centerlines` (https://github.com/nvkelso/natural-earth-vector, folder `geojson/`) in a folder `naturalearth/`.
3. GFED4.1s monthly files are available at https://www.globalfiredata.org.

### Input data

Public GEE catalog: `MODIS/061/MCD64A1`, `ESA/CCI/FireCCI/5_1`, `COPERNICUS/S2_SR_HARMONIZED`, `USDA/NASS/CDL/2020`.
User assets of the authors (must be readable, or replaced by your own copies):

- VIIRS 375 m active fires 2020 (NASA FIRMS): `projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/NOAA20_2020` and `.../SUIOMI2020`. These assets have no GEE time stamp; script 02 selects them by month with the `ACQ_DATE` property (milliseconds).
- FAOSTAT yields: `projects/ee-cherifelkhalil-etu/assets/crop_yield_table` (script 01; script 02 has the yields written in the code).

### Applying the model elsewhere

Change the area of interest (`AOI_COORDS`), the year, the crop map and the parameter table (yields, RPR, DM, CF, EF). The burned-area module can use MCD64A1, FireCCI51 or the Sentinel-2 method. The Sentinel-2 settings (processing order, `VIIRS_BUFFER_M`, compositing) strongly affect the result and should be calibrated against reference burned fields.

## Related work

- Oulad Dahman, C.; Cherif, E.K.; Boulaassal, H. Estimating cropland fire emissions with remote sensing: A literature review and data strategy. *E3S Web Conf.* 2025, 676, 03002. https://doi.org/10.1051/e3sconf/202567603002
- Cherif, E.K.; Oulad Dahman, C.; Bernardino, A. Automated high-resolution mapping of cropland burned areas in the Lower Mississippi River Basin using VIIRS and Sentinel-2 imagery. IGARSS 2026. Code: https://github.com/ElKhalil-Cherif/burned-area-cropland

## License

MIT (see `LICENSE`).
