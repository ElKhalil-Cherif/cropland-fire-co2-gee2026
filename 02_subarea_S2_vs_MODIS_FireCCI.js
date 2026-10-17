// =====================================================================
// 02  CROPLAND FIRE CO2 2020 - SENTINEL-2 + VIIRS METHOD (Tier 2) vs MODIS MCD64A1
// Same sub-area, same months, same crops, same CO2 coefficients.
// (FireCCI51 is computed as well, as a second coarse product.)
//
// v3: follows the order of the original six-phase workflow (Cherif et al.,
// IGARSS 2026), with three indices only (NBR, MIRBI, BAIS2; CSI, NDTI and
// SAVI removed):
//   1. SCL cloud/shadow/water mask + CDL cropland mask
//   2. pre-fire and post-fire images
//   3. NDVI on the post-fire image, Otsu threshold, remove green vegetation
//   4. NBR, MIRBI, BAIS2 on pre and vegetation-filtered post; dIndex = pre - post
//   5. Otsu threshold for each dIndex over all cropland of the Sentinel-2 tile
//   6. consensus: burned only where all three indices agree
//   7. keep only pixels within 500 m of a VIIRS fire on cropland (same month)
//   8. remove objects smaller than the minimum mapping unit
// VIIRS points are assigned to months with their ACQ_DATE property (ms).
//
// HOW TO USE
//  1. Choose REGION below, press Run, then in the Tasks tab press Run on the
//     three tasks CF2020_v3_<REGION>_m09 / _m10 / _m11.
//  2. Download the three CSV files from Google Drive and analyse them with
//     analyse_comparison_gfed.py.
// =====================================================================

// 'ZONE_A' (Testing Zone a, 6 tiles, quick) | 'SOUTH_HALF' | 'NORTH_HALF' | 'FULL'
var REGION = 'NORTH_HALF';
var MONTHS = [9, 10, 11];
var SPLIT_LAT = 33.75;   // middle of the basin (29.6-37.9 N)

var YEAR = 2020;
var AOI_COORDS = [[[-89.833008,37.308108],[-91.606201,35.715921],[-91.861084,35.019116],[-92.243408,34.664959],[-92.243408,34.266222],[-91.518311,33.655409],[-91.584229,33.215125],[-92.221436,32.430002],[-91.836914,31.834296],[-92.485108,31.084225],[-92.287354,30.792170],[-93.667236,30.109015],[-91.965910,29.738164],[-90.477905,29.661294],[-90.322998,29.920179],[-90.427368,30.108103],[-90.744324,30.152648],[-90.793213,30.705001],[-91.363349,30.721486],[-91.520508,31.040631],[-90.853638,32.144007],[-89.925293,33.279201],[-89.914307,34.229610],[-89.474854,35.104721],[-88.486084,35.599556],[-88.272949,36.469881],[-88.290528,37.049043],[-88.244385,37.358198],[-88.002686,37.749851],[-89.004639,37.913901],[-89.917603,37.904376],[-90.038452,37.541083],[-89.833008,37.308108]]];
var aoi = ee.Geometry.Polygon(AOI_COORDS);

// Crops and CO2 coefficients, t CO2 per burned ha (paper Table 2:
// FAOSTAT 2020 US yield x RPR x DM x CF x 1441 g/kg / 1e6). Peanuts = 0 (no yield).
var CROP_IDS = [1, 2, 3, 4, 5, 10, 13, 24, 28, 29, 45, 46];
var YIELD = {'1': 10760.5, '2': 2750.7, '3': 8540.1, '4': 4594.6, '5': 3432.7, '10': 0,
  '13': 10760.5, '24': 3341.7, '28': 2335.2, '29': 1054.0, '45': 85400.5, '46': 21744.3};
var PARAMS = {
  '1':  {RPR: 1.00, DMF: 0.87, CF: 0.80}, '2':  {RPR: 2.05, DMF: 0.85, CF: 0.90},
  '3':  {RPR: 1.40, DMF: 0.89, CF: 0.80}, '4':  {RPR: 1.40, DMF: 0.89, CF: 0.85},
  '5':  {RPR: 2.10, DMF: 0.91, CF: 0.85}, '10': {RPR: 1.00, DMF: 0.94, CF: 0.85},
  '13': {RPR: 1.00, DMF: 0.87, CF: 0.80}, '24': {RPR: 1.30, DMF: 0.89, CF: 0.90},
  '28': {RPR: 1.30, DMF: 0.89, CF: 0.85}, '29': {RPR: 1.40, DMF: 0.90, CF: 0.85},
  '45': {RPR: 0.47, DMF: 0.85, CF: 0.80}, '46': {RPR: 0.40, DMF: 0.22, CF: 0.85}
};
var FACTORS = ee.Dictionary.fromLists(CROP_IDS.map(String), CROP_IDS.map(function(id) {
  var p = PARAMS[String(id)];
  return YIELD[String(id)] * p.RPR * p.DMF * p.CF * 1441 / 1e6;
}));
print('t CO2 per burned ha (check against Table 2)', FACTORS);

var crop = ee.Image('USDA/NASS/CDL/' + YEAR).select('cropland');
var cropMask = crop.remap(CROP_IDS, CROP_IDS.map(function() { return 1; }), 0);

// VIIRS 375 m active fires 2020 (S-NPP + NOAA-20), shared from the old account
var viirs = ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/NOAA20_2020')
  .merge(ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/SUIOMI2020'))
  .filterBounds(aoi);

var CMP_MONTHS = MONTHS;
var S2_SCALE = 30;
var S2_CRS = 'EPSG:5070';
var MIN_PATCH_PIXELS = 1;         // minimum mapping unit in 30 m pixels (1 pixel = 0.09 ha)
var VIIRS_BUFFER_M = 500;         // buffer around VIIRS fires (original workflow: 500 m)
var CLOUD_PERCENT = 20;
var HIST_SCALE = 60;              // pixel size used only to build the Otsu histograms
var PRE_DAYS = 30;
var POST_EXTRA_DAYS = 15;
var cropMaskS2 = cropMask.eq(1);

var REGIONS = {
  ZONE_A: ee.Geometry.Rectangle([-92.1, 29.6, -91.0, 30.75]),
  SOUTH_HALF: ee.Geometry.Rectangle([-94, 29, -87, SPLIT_LAT]),
  NORTH_HALF: ee.Geometry.Rectangle([-94, SPLIT_LAT, -87, 38.5]),
  FULL: aoi
};
var zone = REGIONS[REGION].intersection(aoi, 1);

// VIIRS points of one month. The VIIRS assets store the date in the
// property ACQ_DATE as milliseconds since 1970 (no Earth Engine time stamp),
// so the month is selected on ACQ_DATE.
function viirsInMonth(month) {
  var start = ee.Date.fromYMD(YEAR, month, 1);
  var end = start.advance(1, 'month');
  return viirs.filter(ee.Filter.or(
    ee.Filter.date(start, end),
    ee.Filter.and(ee.Filter.gte('ACQ_DATE', start.millis()), ee.Filter.lt('ACQ_DATE', end.millis()))));
}

// --------------------------- SENTINEL-2 ------------------------------
function maskS2(image) {
  var scl = image.select('SCL');
  // no data, saturated, cloud shadow, clouds, cirrus, snow (+60 m buffer)
  var bad = scl.eq(0).or(scl.eq(1)).or(scl.eq(3)).or(scl.eq(8)).or(scl.eq(9))
    .or(scl.eq(10)).or(scl.eq(11)).focalMax(60, 'circle', 'meters');
  return image.select(['B4', 'B5', 'B6', 'B8', 'B8A', 'B11', 'B12'])
    .divide(10000).clamp(0, 1)
    .updateMask(bad.not().and(scl.neq(6)))   // also drop water
    .copyProperties(image, ['system:time_start']);
}

function addIndices(image) {
  image = ee.Image(image);
  var nbr = image.normalizedDifference(['B8', 'B12']).rename('NBR');
  var mirbi = image.expression('10*swir2 - 9.8*swir1 + 2', {
    swir2: image.select('B12'), swir1: image.select('B11')}).rename('MIRBI');
  var bais2 = image.expression(
    '(1 - sqrt((re2*re3*nirn)/red)) * (((swir2 - nirn)/sqrt(swir2 + nirn)) + 1)', {
      re2: image.select('B5'), re3: image.select('B6'), nirn: image.select('B8A'),
      red: image.select('B4'), swir2: image.select('B12')}).rename('BAIS2');
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  return image.addBands([nbr, mirbi, bais2, ndvi]);
}

// ------------------------------ OTSU ---------------------------------
// Returns the threshold that best separates the histogram into 2 classes.
function otsu(histogram) {
  histogram = ee.Dictionary(histogram);
  var counts = ee.Array(histogram.get('histogram'));
  var means = ee.Array(histogram.get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(), [0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  var bss = ee.List.sequence(1, ee.Number(size).subtract(1)).map(function(i) {
    i = ee.Number(i);
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMean = means.slice(0, 0, i).multiply(aCounts)
      .reduce(ee.Reducer.sum(), [0]).get([0]).divide(aCount.max(1e-12));
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount.max(1e-12));
    return aCount.multiply(aMean.subtract(mean).pow(2))
      .add(bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  var best = ee.Number(ee.List(ee.Array(bss).argmax()).get(0));
  return ee.Number(means.toList().get(best));
}

// ------------------------- ONE TILE, ONE MONTH -----------------------
// Returns a byte band 'mMM': 1 = burned that month, 0 = not burned.
// VIIRS fires of one month that fall on cropland (Phase 2, step 1).
function viirsCropFires(month) {
  var pts = viirsInMonth(month).filterBounds(zone);
  return cropMask.rename('crop').reduceRegions({collection: pts, reducer: ee.Reducer.first(),
    scale: 30, crs: S2_CRS}).filter(ee.Filter.eq('first', 1));
}
var viirsByMonth = {};
MONTHS.forEach(function(m) {
  var buffers = viirsCropFires(m).map(function(f) { return ee.Feature(f).buffer(VIIRS_BUFFER_M, 50); });
  viirsByMonth[m] = ee.Image(0).byte().paint(buffers, 1).selfMask();
});

function histogram(image, band, geom) {
  return image.select(band).reduceRegion({
    reducer: ee.Reducer.histogram({maxBuckets: 64}), geometry: geom,
    scale: HIST_SCALE, crs: S2_CRS, maxPixels: 1e9, tileScale: 4}).get(band);
}
// Otsu threshold; 'fallback' is used when the tile has no valid pixels.
function threshold(image, band, geom, fallback) {
  var h = histogram(image, band, geom);
  var ok = ee.Algorithms.If(h, ee.List(ee.Dictionary(h).get('bucketMeans')).size().gte(2), false);
  return ee.Number(ee.Algorithms.If(ok, otsu(h), fallback));
}

function monthBurn(month, tileId, tileGeom) {
  var mm = ('0' + month).slice(-2);
  var monthStart = ee.Date.fromYMD(YEAR, month, 1);
  var monthEnd = monthStart.advance(1, 'month');

  // Step 1: Sentinel-2 L2A of this tile, SCL mask (in maskS2) and cropland mask.
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(tileGeom)
    .filter(ee.Filter.eq('MGRS_TILE', tileId))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PERCENT));
  var preCol = s2.filterDate(monthStart.advance(-PRE_DAYS, 'day'), monthStart).map(maskS2);
  var postCol = s2.filterDate(monthStart, monthEnd.advance(POST_EXTRA_DAYS, 'day')).map(maskS2);

  // Step 2: pre-fire image = median of the 30 days before the month;
  // post-fire image = per pixel, the date with the lowest NBR in the month + 15 days
  // (GEE works with composites; this keeps short-lived burn scars).
  var pre = addIndices(preCol.median()).updateMask(cropMaskS2);
  var post = postCol.map(function(img) {
    img = addIndices(img);
    return img.addBands(img.select('NBR').multiply(-1).rename('negNBR'));
  }).qualityMosaic('negNBR').updateMask(cropMaskS2);

  // Step 3: vegetation filter on the post-fire image: NDVI + Otsu.
  var ndviT = threshold(post, 'NDVI', tileGeom, 1);
  var postFiltered = post.updateMask(post.select('NDVI').lte(ndviT));

  // Step 4: dIndex = pre - post for NBR, MIRBI and BAIS2.
  var d = ee.Image.cat([
    pre.select('NBR').subtract(postFiltered.select('NBR')).rename('dNBR'),
    pre.select('MIRBI').subtract(postFiltered.select('MIRBI')).rename('dMIRBI'),
    pre.select('BAIS2').subtract(postFiltered.select('BAIS2')).rename('dBAIS2')
  ]);

  // Step 5: Otsu threshold of each dIndex over all cropland of the tile.
  // NBR drops after burning (dNBR high); MIRBI and BAIS2 rise (dMIRBI, dBAIS2 low).
  var tNBR = threshold(d, 'dNBR', tileGeom, 1e9);
  var tMIRBI = threshold(d, 'dMIRBI', tileGeom, -1e9);
  var tBAIS2 = threshold(d, 'dBAIS2', tileGeom, -1e9);

  // Step 6: consensus = logical intersection of the three masks.
  var burn = d.select('dNBR').gt(tNBR)
    .and(d.select('dMIRBI').lt(tMIRBI))
    .and(d.select('dBAIS2').lt(tBAIS2));

  // Step 7: keep only pixels within 500 m of a VIIRS fire on cropland in the same month.
  burn = burn.updateMask(viirsByMonth[month]).selfMask();

  // Step 8: minimum mapping unit.
  if (MIN_PATCH_PIXELS > 1) {
    var patchSize = burn.connectedPixelCount({maxSize: 256, eightConnected: true});
    burn = burn.updateMask(patchSize.gte(MIN_PATCH_PIXELS));
  }

  var result = burn.unmask(0).byte().rename('m' + mm);
  var empty = ee.Image(0).byte().rename('m' + mm);
  var ok = preCol.size().gt(0).and(postCol.size().gt(0));
  return ee.Image(ee.Algorithms.If(ok, result, empty));
}



// Sentinel-2 tiles over the sub-area; each tile gets its own Otsu thresholds.
var s2Zone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(zone)
  .filterDate(YEAR + '-08-01', YEAR + '-12-16');
var tileIds = s2Zone.aggregate_array('MGRS_TILE').distinct();

var CROP_FACT = FACTORS;
function oneMonth(month) {
  var mm = ('0' + month).slice(-2);
  var s2 = ee.ImageCollection.fromImages(tileIds.map(function(id) {
    var tileGeom = s2Zone.filter(ee.Filter.eq('MGRS_TILE', id)).geometry().dissolve(100).intersection(zone, 100);
    return monthBurn(month, ee.String(id), tileGeom).clip(tileGeom);
  })).max().unmask(0).gt(0);
  var start = ee.Date.fromYMD(YEAR, month, 1);
  function coarse(id) {
    return ee.ImageCollection(id).select('BurnDate').filterDate(start, start.advance(1, 'month'))
      .max().unmask(0).gt(0);
  }
  var modis = coarse('MODIS/061/MCD64A1');
  var cci = coarse('ESA/CCI/FireCCI/5_1');
  var ha = ee.Image.pixelArea().divide(1e4);
  var bands = ee.Image.cat([
    ha.rename('cropland'),
    ha.updateMask(s2).rename('s2'),
    ha.updateMask(modis).rename('modis'),
    ha.updateMask(cci).rename('firecci'),
    ha.updateMask(s2.and(modis)).rename('s2_and_modis'),
    ha.updateMask(s2.and(cci)).rename('s2_and_firecci'),
    ha.updateMask(modis.and(cci)).rename('modis_and_firecci')
  ]).updateMask(cropMaskS2);
  var names = ['cropland', 's2', 'modis', 'firecci', 's2_and_modis', 's2_and_firecci', 'modis_and_firecci'];
  var groups = ee.List(bands.addBands(crop.rename('crop')).reduceRegion({
    reducer: ee.Reducer.sum().repeat(names.length).group({groupField: names.length, groupName: 'crop'}),
    geometry: zone, scale: S2_SCALE, crs: S2_CRS, tileScale: 16, maxPixels: 1e13
  }).get('groups'));
  // one row per crop: hectares of each layer and t CO2 of s2 / modis / firecci
  var rows = groups.map(function(g) {
    g = ee.Dictionary(g);
    var id = ee.Number(g.get('crop')).format('%.0f');
    var sums = ee.List(g.get('sum'));
    var f = ee.Number(CROP_FACT.get(id, 0));
    var props = ee.Dictionary.fromLists(names.map(function(n) { return n + '_ha'; }), sums)
      .set('crop', id).set('month', month).set('region', REGION)
      .set('s2_tCO2', ee.Number(sums.get(1)).multiply(f))
      .set('modis_tCO2', ee.Number(sums.get(2)).multiply(f))
      .set('firecci_tCO2', ee.Number(sums.get(3)).multiply(f));
    return ee.Feature(null, props);
  });
  Export.table.toDrive({
    collection: ee.FeatureCollection(rows),
    description: 'CF2020_v3_' + REGION + '_m' + mm,
    fileFormat: 'CSV',
    selectors: ['region', 'month', 'crop', 'cropland_ha', 's2_ha', 'modis_ha', 'firecci_ha',
      's2_and_modis_ha', 's2_and_firecci_ha', 'modis_and_firecci_ha', 's2_tCO2', 'modis_tCO2', 'firecci_tCO2']
  });
  return s2;
}

var s2Maps = MONTHS.map(oneMonth);
print('Region', REGION, 'Sentinel-2 tiles', tileIds);
MONTHS.forEach(function(m) {
  print('VIIRS points in region, month ' + m, viirsInMonth(m).filterBounds(zone).size(),
        'of which on cropland', viirsCropFires(m).size());
});
print('Now open the Tasks tab and press Run on the ' + MONTHS.length + ' CF2020_v3 tasks.');

Map.centerObject(zone, REGION === 'ZONE_A' ? 9 : 7);
Map.addLayer(zone, {color: 'yellow'}, 'Sub-area ' + REGION);
