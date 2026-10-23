// =====================================================================
// CROPLAND FIRE CO2 2020 - SENTINEL-2 + VIIRS METHOD (Tier 2) vs MODIS MCD64A1 (Tier 1)
// Same area (whole basin), same 12 months, same 30 m grid, same crops,
// same CO2 coefficients for both tiers. (FireCCI51 is computed as well.)
//
// v5: pre/post pairing of the original workflow (Cherif et al., IGARSS 2026):
//  - for each Sentinel-2 tile, the clear acquisition dates form a time series;
//    each date is paired with the NEXT date of the same tile (pre = date i,
//    post = date i+1) when they are less than 30 days apart, so every post image
//    becomes the pre image of the following pair. No compositing, no median.
//  - each pair is processed on its own, tile by tile: NDVI Otsu on the post
//    image, dIndex = pre - post, one Otsu threshold per dIndex computed over all
//    cropland of that tile, consensus of the 3 indices, then the VIIRS mask.
//  - VIIRS mask of a pair = 500 m buffers of VIIRS fires on cropland acquired
//    between the pre date and the post date (inclusive).
//  - month of a burned pixel = month of the VIIRS fire that confirms it (the
//    fire within 500 m in the pair window), not the month of the post image.
//    Example: pre 15 Sep, post 2 Oct, VIIRS fire 25 Sep -> September.
//    If fires of two months confirm the same pixel, the earlier month is used,
//    so no pixel is counted twice. The monthly map is the union over pairs and tiles.
//  - pairs without any VIIRS cropland fire of the month in their window are
//    skipped (the VIIRS mask would remove every pixel anyway): quota only.
//
// Steps of each pair (original order; CSI, NDTI and SAVI removed):
//   1. SCL cloud/shadow/water mask + CDL cropland mask
//   2. pre-fire image (date i) and post-fire image (date i+1), same tile
//   3. NDVI on the post-fire image, Otsu threshold, remove green vegetation
//   4. NBR, MIRBI, BAIS2 on pre and vegetation-filtered post; dIndex = pre - post
//   5. Otsu threshold for each dIndex over all cropland of the tile
//   6. consensus: burned only where all three indices agree
//   7. keep only pixels within 500 m of a VIIRS fire on cropland in the pair window
//      (month of the pixel = month of that fire)
//   8. minimum mapping unit (1 pixel at 30 m)
//
// HOW TO USE
//  1. Press Run, then in the Tasks tab press Run on the 12 tasks
//     CF2020_v5_FULL_m01 ... _m12 (one month per task).
//  2. If a month fails (memory or time limit), run that month again with
//     REGION = 'NORTH_HALF' and then 'SOUTH_HALF' (MONTHS = [that month]);
//     the two CSV files of that month are simply added together.
//  3. Download the CSV files from Google Drive.
//  4. EXPORT_MAPS = true also exports the monthly 30 m burned-area maps
//     (needed for the validation).
// =====================================================================

// 'ZONE_A' (Testing Zone a, 6 tiles, quick) | 'SOUTH_HALF' | 'NORTH_HALF' | 'FULL'
var REGION = 'FULL';
var MONTHS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
var EXPORT_MAPS = true;
var PAIR_MAX_DAYS = 30;   // pre and post images of a pair must be less than 30 days apart
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
var cropMaskS2 = cropMask.eq(1);

var REGIONS = {
  ZONE_A: ee.Geometry.Rectangle([-92.1, 29.6, -91.0, 30.75]),
  SOUTH_HALF: ee.Geometry.Rectangle([-94, 29, -87, SPLIT_LAT]),
  NORTH_HALF: ee.Geometry.Rectangle([-94, SPLIT_LAT, -87, 38.5]),
  FULL: aoi
};
var zone = REGIONS[REGION].intersection(aoi, 1);

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

// ------------------------------ VIIRS --------------------------------
// VIIRS fires of the whole year that fall on cropland (Phase 2, step 1).
// ACQ_DATE is the acquisition day in milliseconds since 1970.
var cropFires = cropMask.rename('crop').reduceRegions({
  collection: viirs.filterBounds(zone.buffer(VIIRS_BUFFER_M, 100)), reducer: ee.Reducer.first(),
  scale: 30, crs: S2_CRS}).filter(ee.Filter.eq('first', 1));
function viirsCropFires(month) {
  var start = ee.Date.fromYMD(YEAR, month, 1);
  return cropFires.filter(ee.Filter.and(ee.Filter.gte('ACQ_DATE', start.millis()),
    ee.Filter.lt('ACQ_DATE', start.advance(1, 'month').millis())));
}

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

// ------------------------- ONE TILE, ONE PAIR ------------------------
// All scenes of one tile acquired on one day, cloud-masked and mosaicked.
function dayImage(col, day) {
  var d = ee.Date(day);
  return addIndices(col.filterDate(d, d.advance(1, 'day')).map(maskS2).mosaic()).updateMask(cropMaskS2);
}

// Returns a byte image: 1 = burned between preDay and postDay and confirmed by a
// VIIRS fire of the month [mStart, mEnd), 0 = not burned.
function pairBurn(col, preDay, postDay, tileGeom, mStart, mEnd) {
  var preD = ee.Date(preDay), postD = ee.Date(postDay);
  // Step 7 input: VIIRS cropland fires acquired between the two images.
  var fires = cropFires.filter(ee.Filter.and(
      ee.Filter.gte('ACQ_DATE', preD.millis()),
      ee.Filter.lt('ACQ_DATE', postD.advance(1, 'day').millis())))
    .filterBounds(tileGeom.buffer(VIIRS_BUFFER_M, 100));

  // Steps 1-2: pre-fire and post-fire images (SCL mask in maskS2, cropland mask).
  var pre = dayImage(col, preDay);
  var post = dayImage(col, postDay);

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

  // Step 7: keep only pixels within 500 m of a VIIRS cropland fire of this pair
  // that was detected in this month, and not near an earlier-month fire of the pair.
  function bufferImage(fc) {
    return ee.Image(0).byte().rename('burned').paint(fc.map(function(f) { return ee.Feature(f).buffer(VIIRS_BUFFER_M, 50); }), 1);
  }
  var firesMonth = fires.filter(ee.Filter.and(ee.Filter.gte('ACQ_DATE', mStart.millis()), ee.Filter.lt('ACQ_DATE', mEnd.millis())));
  var firesBefore = fires.filter(ee.Filter.lt('ACQ_DATE', mStart.millis()));
  var viirsMask = bufferImage(firesMonth).and(bufferImage(firesBefore).not());
  burn = burn.updateMask(viirsMask.selfMask()).selfMask();

  // Step 8: minimum mapping unit.
  if (MIN_PATCH_PIXELS > 1) {
    burn = burn.updateMask(burn.connectedPixelCount({maxSize: 256, eightConnected: true}).gte(MIN_PATCH_PIXELS));
  }
  var empty = ee.Image(0).byte().rename('burned');
  return ee.Image(ee.Algorithms.If(firesMonth.size().gt(0), burn.unmask(0).byte().rename('burned'), empty));
}

// ------------------------- ONE TILE, ONE MONTH -----------------------
// Pairs consecutive acquisition dates of the tile; keeps the pairs that are less
// than PAIR_MAX_DAYS apart and whose window [pre, post] overlaps the month.
function tileMonthBurn(tileId, month) {
  var mStart = ee.Date.fromYMD(YEAR, month, 1);
  var mEnd = mStart.advance(1, 'month');
  var tileAll = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filter(ee.Filter.eq('MGRS_TILE', tileId))
    .filterDate(mStart.advance(-PAIR_MAX_DAYS, 'day'), mEnd.advance(PAIR_MAX_DAYS, 'day'));
  var tileGeom = tileAll.geometry().dissolve(100);   // whole tile: Otsu per tile
  var col = tileAll.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PERCENT));
  var days = col.aggregate_array('system:time_start')
    .map(function(t) { return ee.Date(t).format('YYYY-MM-dd'); }).distinct().sort();
  var pairs = days.slice(0, -1).zip(days.slice(1));   // [[d0, d1], [d1, d2], ...]
  var maps = pairs.map(function(p) {
    p = ee.List(p);
    var preDay = ee.String(p.get(0)), postDay = ee.String(p.get(1));
    var gap = ee.Date(postDay).difference(ee.Date(preDay), 'day');
    var overlaps = ee.Date(postDay).millis().gte(mStart.millis())
      .and(ee.Date(preDay).millis().lt(mEnd.millis()));
    return ee.Algorithms.If(gap.lt(PAIR_MAX_DAYS).and(overlaps),
      pairBurn(col, preDay, postDay, tileGeom, mStart, mEnd), ee.Image(0).byte().rename('burned'));
  });
  return ee.ImageCollection.fromImages(maps.add(ee.Image(0).byte().rename('burned'))).max().clip(tileGeom);
}

// Sentinel-2 tiles over the region; each tile gets its own Otsu thresholds.
var s2Zone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(zone)
  .filterDate(YEAR + '-01-01', (YEAR + 1) + '-01-01');
var tileIds = s2Zone.aggregate_array('MGRS_TILE').distinct();

var CROP_FACT = FACTORS;
function oneMonth(month) {
  var mm = ('0' + month).slice(-2);
  var s2 = ee.ImageCollection.fromImages(tileIds.map(function(id) {
    return tileMonthBurn(ee.String(id), month);
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
    description: 'CF2020_v5_' + REGION + '_m' + mm,
    fileFormat: 'CSV',
    selectors: ['region', 'month', 'crop', 'cropland_ha', 's2_ha', 'modis_ha', 'firecci_ha',
      's2_and_modis_ha', 's2_and_firecci_ha', 'modis_and_firecci_ha', 's2_tCO2', 'modis_tCO2', 'firecci_tCO2']
  });
  if (EXPORT_MAPS) {
    Export.image.toDrive({
      image: s2.updateMask(cropMaskS2).unmask(0).byte().clip(zone),
      description: 'CF2020_v5_map_' + REGION + '_m' + mm,
      region: zone, scale: S2_SCALE, crs: S2_CRS, maxPixels: 1e13
    });
  }
  return s2;
}

var s2Maps = MONTHS.map(oneMonth);
print('Region', REGION, 'Sentinel-2 tiles', tileIds);
MONTHS.forEach(function(m) {
  print('VIIRS fires on cropland, month ' + m, viirsCropFires(m).filterBounds(zone).size());
});
print('Now open the Tasks tab and press Run on the ' + MONTHS.length + ' CF2020_v5 tasks.');

Map.centerObject(zone, REGION === 'ZONE_A' ? 9 : 7);
Map.addLayer(zone, {color: 'yellow'}, 'Region ' + REGION);
