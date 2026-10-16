// =====================================================================
// 01  CROPLAND FIRE CO2 2020 - WHOLE BASIN, MODIS MCD64A1 (Tier 1)
//
// Lower Mississippi River Basin, January-December 2020.
// MODIS MCD64A1 burned area x USDA CDL 2020 crop type x crop-specific
// emission coefficients (FAOSTAT yield x RPR x DM x CF x EF).
// Prints annual and monthly CO2 by crop and the VIIRS consistency check
// (RESULTS_JSON) and a link to the emission-intensity map (MAP_THUMBNAIL_URL).
// The Monte Carlo analysis is run outside GEE with monte_carlo.py.
//
// HOW TO USE: press Run; results appear in the Console in 1-3 minutes.
// =====================================================================

// ======================= PART 1: MODIS, WHOLE BASIN =======================
var YEAR = 2020;
var SCALE = 120;
var MC_DRAWS = 10000;
var A_REL_SD = 0.25;   // burned-area uncertainty (1 SD, relative) - assumption, not validated
var EF_REL_SD = 0.15;  // emission-factor uncertainty (1 SD, relative)

var AOI_COORDS = [[[-89.833008,37.308108],[-91.606201,35.715921],[-91.861084,35.019116],[-92.243408,34.664959],[-92.243408,34.266222],[-91.518311,33.655409],[-91.584229,33.215125],[-92.221436,32.430002],[-91.836914,31.834296],[-92.485108,31.084225],[-92.287354,30.792170],[-93.667236,30.109015],[-91.965910,29.738164],[-90.477905,29.661294],[-90.322998,29.920179],[-90.427368,30.108103],[-90.744324,30.152648],[-90.793213,30.705001],[-91.363349,30.721486],[-91.520508,31.040631],[-90.853638,32.144007],[-89.925293,33.279201],[-89.914307,34.229610],[-89.474854,35.104721],[-88.486084,35.599556],[-88.272949,36.469881],[-88.290528,37.049043],[-88.244385,37.358198],[-88.002686,37.749851],[-89.004639,37.913901],[-89.917603,37.904376],[-90.038452,37.541083],[-89.833008,37.308108]]];
var aoi = ee.Geometry.Polygon(AOI_COORDS);

var CROP_IDS = [1, 2, 3, 4, 5, 10, 13, 24, 28, 29, 45, 46];
var CROP_NAMES = ee.Dictionary({'1': 'Corn', '2': 'Cotton', '3': 'Rice', '4': 'Sorghum', '5': 'Soybeans',
  '10': 'Peanuts', '13': 'Popcorn', '24': 'Winter wheat', '28': 'Oats', '29': 'Millet',
  '45': 'Sugarcane', '46': 'Sweet potato'});

// FAOSTAT item names (a list = first one found in the table is used)
var CROP_TO_FAO = {
  '1': ['Maize (corn)', 'Maize'], '13': ['Maize (corn)', 'Maize'],
  '2': ['Seed cotton, unginned', 'Seed cotton'], '3': ['Rice', 'Rice, paddy'],
  '4': ['Sorghum'], '5': ['Soya beans', 'Soybeans'], '29': ['Millet'], '28': ['Oats'],
  '10': ['Groundnuts, excluding shelled', 'Groundnuts, with shell', 'Groundnuts, in shell'],  // FIX 3
  '46': ['Sweet potatoes'], '45': ['Sugar cane'], '24': ['Wheat']
};

// Table 1 of the paper. rpr = central value used for the main estimate;
// rprMin/rprMax = range sampled in the Monte Carlo.
var PARAMS = {
  '1':  {RPR: 1.00, rprMin: 1.00, rprMax: 1.04, DMF: 0.87, CF: 0.80, EF: 1441},
  '2':  {RPR: 2.05, rprMin: 2.00, rprMax: 2.10, DMF: 0.85, CF: 0.90, EF: 1441},  // FIX 2
  '3':  {RPR: 1.40, rprMin: 1.08, rprMax: 1.40, DMF: 0.89, CF: 0.80, EF: 1441},
  '4':  {RPR: 1.40, rprMin: 1.40, rprMax: 1.46, DMF: 0.89, CF: 0.85, EF: 1441},
  '5':  {RPR: 2.10, rprMin: 2.10, rprMax: 2.22, DMF: 0.91, CF: 0.85, EF: 1441},
  '10': {RPR: 1.00, rprMin: 1.00, rprMax: 1.22, DMF: 0.94, CF: 0.85, EF: 1441},
  '13': {RPR: 1.00, rprMin: 1.00, rprMax: 1.00, DMF: 0.87, CF: 0.80, EF: 1441},
  '24': {RPR: 1.30, rprMin: 1.30, rprMax: 1.32, DMF: 0.89, CF: 0.90, EF: 1441},
  '28': {RPR: 1.30, rprMin: 1.30, rprMax: 1.46, DMF: 0.89, CF: 0.85, EF: 1441},
  '29': {RPR: 1.40, rprMin: 1.40, rprMax: 1.50, DMF: 0.90, CF: 0.85, EF: 1441},
  '45': {RPR: 0.47, rprMin: 0.27, rprMax: 0.47, DMF: 0.85, CF: 0.80, EF: 1441},  // FIX 1
  '46': {RPR: 0.40, rprMin: 0.40, rprMax: 0.42, DMF: 0.22, CF: 0.85, EF: 1441}
};

var yields = ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/crop_yield_table')
  .filter(ee.Filter.eq('Area', 'United States of America'))
  .filter(ee.Filter.eq('Year', YEAR));

// Yield in kg/ha from the FAOSTAT table (-1 if the item is missing)
function yieldForCrop(id) {
  var rows = yields.filter(ee.Filter.inList('Item', CROP_TO_FAO[String(id)]));
  var row = rows.first();
  return ee.Number(ee.Algorithms.If(rows.size().gt(0),
    ee.Number.parse(ee.String(ee.Feature(row).get('Value')).replace(',', '', 'g')), -1));
}
var YIELDS = {};
CROP_IDS.forEach(function(id) { YIELDS[id] = yieldForCrop(id); });

// t CO2 per burned hectare = Y(kg/ha) x RPR x DMF x CF x EF(g/kg) / 1e6
var factorValues = CROP_IDS.map(function(id) {
  var p = PARAMS[String(id)];
  return YIELDS[id].max(0).multiply(p.RPR * p.DMF * p.CF * p.EF).divide(1e6);
});
var FACTORS = ee.Dictionary.fromLists(CROP_IDS.map(String), factorValues);

var crop = ee.Image('USDA/NASS/CDL/' + YEAR).select('cropland').clip(aoi);
var cropMask = crop.remap(CROP_IDS, CROP_IDS.map(function() { return 1; }), 0);
var modis = ee.ImageCollection('MODIS/061/MCD64A1')
  .filterDate(YEAR + '-01-01', (YEAR + 1) + '-01-01').filterBounds(aoi);

function burnForMonth(month) {
  var start = ee.Date.fromYMD(YEAR, month, 1);
  return modis.filterDate(start, start.advance(1, 'month')).select('BurnDate')
    .max().gt(0).rename('burn').clip(aoi);
}

// Burned hectares per crop for one burn mask -> ee.Dictionary {'5': ha, ...}
function areaByCrop(burn) {
  var area = ee.Image.pixelArea().divide(1e4).rename('area_ha').updateMask(burn).updateMask(cropMask);
  var groups = ee.List(ee.Dictionary(area.addBands(crop.rename('crop_id')).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'crop_id'}),
    geometry: aoi, scale: SCALE, crs: 'EPSG:5070', tileScale: 16, maxPixels: 1e13   // FIX 4
  })).get('groups', ee.List([])));
  return ee.Dictionary.fromLists(
    groups.map(function(g) { return ee.Number(ee.Dictionary(g).get('crop_id')).format('%.0f'); }),
    groups.map(function(g) { return ee.Dictionary(g).get('sum'); }));
}

function cropTable(areaDict, label) {
  return ee.FeatureCollection(CROP_IDS.map(function(id) {
    var ha = ee.Number(ee.Dictionary(areaDict).get(String(id), 0));
    var f = ee.Number(FACTORS.get(String(id)));
    return ee.Feature(null, {month: label, crop_id: id, crop_name: CROP_NAMES.get(String(id)),
      burned_area_ha: ha, factor_tCO2_ha: f, emission_tCO2: ha.multiply(f)});
  }));
}

// ---------------- Annual (union of the 12 monthly masks) ----------------
var annualBurn = modis.select('BurnDate').max().gt(0).rename('burn').clip(aoi);
var annualArea = areaByCrop(annualBurn);
var annualCrop = cropTable(annualArea, 'annual');
var annualAreaHa = annualCrop.aggregate_sum('burned_area_ha');
var annualCO2 = annualCrop.aggregate_sum('emission_tCO2');

// ---------------- Monthly x crop ----------------
var monthlyCrop = ee.FeatureCollection([1,2,3,4,5,6,7,8,9,10,11,12].map(function(m) {
  return cropTable(areaByCrop(burnForMonth(m)), m);
})).flatten();
var monthlyTotals = ee.FeatureCollection([1,2,3,4,5,6,7,8,9,10,11,12].map(function(m) {
  var s = monthlyCrop.filter(ee.Filter.eq('month', m));
  return ee.Feature(null, {month: m, burned_area_ha: s.aggregate_sum('burned_area_ha'),
    emission_tCO2: s.aggregate_sum('emission_tCO2')});
}));

// ---------------- Monte Carlo (FIX 5) ----------------
// Not printed: too heavy for Earth Engine memory. monte_carlo.py runs the same
// Monte Carlo outside Earth Engine from the per-crop areas and yields.
var draws = ee.FeatureCollection(ee.List.sequence(1, MC_DRAWS).map(function(i) {
  return ee.Feature(null, {draw: i});
}));
CROP_IDS.forEach(function(id, j) {
  draws = draws.randomColumn('zA_' + id, 100 + j, 'normal')
               .randomColumn('zE_' + id, 200 + j, 'normal')
               .randomColumn('uR_' + id, 300 + j, 'uniform');
});
draws = draws.map(function(f) {
  var props = {}; var total = ee.Number(0);
  CROP_IDS.forEach(function(id) {
    var p = PARAMS[String(id)];
    var ha = ee.Number(annualArea.get(String(id), 0))
      .multiply(ee.Number(1).add(ee.Number(f.get('zA_' + id)).multiply(A_REL_SD)).max(0));
    var rpr = ee.Number(f.get('uR_' + id)).multiply(p.rprMax - p.rprMin).add(p.rprMin);
    var ef = ee.Number(1).add(ee.Number(f.get('zE_' + id)).multiply(EF_REL_SD)).max(0).multiply(p.EF);
    var e = ha.multiply(YIELDS[id].max(0)).multiply(rpr).multiply(p.DMF * p.CF).multiply(ef).divide(1e6);
    props['E_' + id] = e; total = total.add(e);
  });
  props.E_total = total;
  return ee.Feature(null, props);
});
function ci(prop) {
  var v = draws.aggregate_array(prop);
  return ee.List([v.reduce(ee.Reducer.percentile([2.5])), v.reduce(ee.Reducer.median()),
                  v.reduce(ee.Reducer.percentile([97.5]))]);
}

// ---------------- VIIRS consistency ----------------
var viirs = ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/NOAA20_2020')
  .merge(ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/SUIOMI2020'))
  .filterBounds(aoi);
var direct = annualBurn.unmask(0).reduceRegions({collection: viirs, reducer: ee.Reducer.first(), scale: 500, tileScale: 8});
var buffered = annualBurn.unmask(0).reduceRegions({
  collection: viirs.map(function(f) { return f.buffer(375); }),
  reducer: ee.Reducer.max(), scale: 500, tileScale: 8});
var viirsN = viirs.size();

// ---------------- Console output ----------------
function rows(fc, keys) {
  return fc.map(function(f) {
    var s = ee.String('');
    keys.forEach(function(k, i) {
      var v = f.get(k);
      s = s.cat(i ? '|' : '').cat(ee.String(ee.Algorithms.If(ee.Algorithms.ObjectType(v).equals('String'),
        v, ee.Number(v).format('%.3f'))));
    });
    return f.set('row', s);
  }).aggregate_array('row');
}
print('YIELDS_USED_kg_ha (-1 = item missing)', ee.Dictionary.fromLists(CROP_IDS.map(String),
  CROP_IDS.map(function(id) { return YIELDS[id]; })));
print('FACTORS_tCO2_per_ha', FACTORS);
print('ANNUAL_BURNED_CROPLAND_HA', annualAreaHa);
print('ANNUAL_CO2_TONNES', annualCO2);
print('MONTHLY_SUM_CO2_TONNES', monthlyTotals.aggregate_sum('emission_tCO2'));
print('CROP_ROWS crop|ha|tCO2', rows(annualCrop.sort('emission_tCO2', false), ['crop_name', 'burned_area_ha', 'emission_tCO2']));
print('MONTHLY_ROWS month|ha|tCO2', rows(monthlyTotals, ['month', 'burned_area_ha', 'emission_tCO2']));
print('MONTHLY_CROP_ROWS month|crop|ha|tCO2',
  rows(monthlyCrop.filter(ee.Filter.gt('burned_area_ha', 0)), ['month', 'crop_name', 'burned_area_ha', 'emission_tCO2']));
print('VIIRS_POINTS_TOTAL', viirsN);
print('VIIRS_DIRECT_MATCH', direct.filter(ee.Filter.eq('first', 1)).size());
print('VIIRS_MATCH_WITHIN_375M', buffered.filter(ee.Filter.eq('max', 1)).size());

// ONE line with every number (text, so nothing is folded in the Console).
// Copy this whole line into results.json for monte_carlo.py.
print('RESULTS_JSON', ee.String.encodeJSON(ee.Dictionary({
  yields: ee.Dictionary.fromLists(CROP_IDS.map(String), CROP_IDS.map(function(id) { return YIELDS[id]; })),
  factors: FACTORS,
  annual_ha: annualAreaHa,
  annual_t: annualCO2,
  monthly_sum_t: monthlyTotals.aggregate_sum('emission_tCO2'),
  crops: rows(annualCrop, ['crop_name', 'burned_area_ha', 'emission_tCO2']),
  months: rows(monthlyTotals, ['month', 'burned_area_ha', 'emission_tCO2']),
  month_crop: rows(monthlyCrop.filter(ee.Filter.gt('burned_area_ha', 0)), ['month', 'crop_name', 'burned_area_ha', 'emission_tCO2']),
  viirs_total: viirsN,
  viirs_direct: direct.filter(ee.Filter.eq('first', 1)).size(),
  viirs_375m: buffered.filter(ee.Filter.eq('max', 1)).size()
})));

// ---------------- Map ----------------
var factorImage = crop.remap(CROP_IDS, factorValues, 0);
var emissionDensity = factorImage.multiply(100).updateMask(annualBurn).updateMask(factorImage.gt(0))
  .rename('tCO2_km2');   // t/ha x 100 = t/km2 of burned cropland
Map.setCenter(-90.8, 33.5, 6);
Map.addLayer(emissionDensity, {min: 0, max: 2000, palette: ['fff7bc', 'fec44f', 'fe9929', 'd95f0e', '993404']},
  '2020 CO2 (t per km2 burned cropland)');
Map.addLayer(annualBurn.selfMask(), {palette: ['ff0000']}, 'MCD64A1 burned area', false);

// Link to the emission-intensity map (Fig. 3a)
print('MAP_THUMBNAIL_URL', emissionDensity.getThumbURL({region: aoi, crs: 'EPSG:4326', min: 0, max: 2000, dimensions: 1600, format: 'png',
                palette: ['fff7bc', 'fec44f', 'fe9929', 'd95f0e', '993404']}));

