// =====================================================================
// CROPLAND FIRE CO2 2020 - TESTING ZONE (a) COMPARISON, SEP-NOV
// Self-contained version for a NEW Earth Engine account.
// Sentinel-2 + VIIRS method (30 m) vs MODIS MCD64A1 vs FireCCI51,
// same crops and same CO2 coefficients (paper Table 2).
//
// BEFORE RUNNING: the only file still read from the old account is the
// VIIRS fire points. In the OLD account: Assets tab > VIIRS_DATA_2020 >
// click NOAA20_2020 > Share > tick "Anyone can read" > Done.
// Do the same for SUIOMI2020.
//
// HOW TO USE: press Run, then Tasks tab > Run on CF2020_compare_zoneA.
// When COMPLETED, open CF2020_compare_zoneA.csv in this account's
// Google Drive and paste the text of the "json" column to Claude.
// =====================================================================

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

// ============ PART 2: SENTINEL-2 vs MODIS vs FireCCI51, TESTING ZONE (a) ============
var CMP_MONTHS = [9, 10, 11];     // peak crop-residue burning months
var S2_SCALE = 30;
var S2_CRS = 'EPSG:5070';
var MIN_PATCH_M2 = 3000;          // smallest burned patch kept
var VIIRS_BUFFER_M = 3000;        // search distance around each VIIRS fire
var CLOUD_PERCENT = 20;
var MIN_VOTES = 4;                // burned when >= 4 of 5 indices agree
var PRE_DAYS = 30;
var POST_EXTRA_DAYS = 15;

// Testing Zone (a), paper Figure 1: Lafayette to Baton Rouge
var zone = ee.Geometry.Rectangle([-92.1, 29.6, -91.0, 30.75]).intersection(aoi, 1);
var cropMaskS2 = cropMask.eq(1);

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
  var ndti = image.normalizedDifference(['B11', 'B12']).rename('NDTI');
  var savi = image.expression('1.5*(nir - red)/(nir + red + 0.5)', {
    nir: image.select('B8'), red: image.select('B4')}).rename('SAVI');
  return image.addBands([nbr, mirbi, bais2, ndti, savi]);
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

// Safe version: if there are no valid pixels, return a huge threshold
// so nothing is marked as burned (never a false "everything burned").
function safeThreshold(histograms, band) {
  var h = ee.Dictionary(histograms).get(band, null);
  var ok = ee.Algorithms.If(h,
    ee.List(ee.Dictionary(h).get('bucketMeans')).size().gte(2), false);
  return ee.Number(ee.Algorithms.If(ok, otsu(h), 1e9));
}

// ------------------------- ONE TILE, ONE MONTH -----------------------
// Returns a byte band 'mMM': 1 = burned that month, 0 = not burned.
function monthBurn(month, tileId, tileGeom) {
  var mm = ('0' + month).slice(-2);
  var monthStart = ee.Date.fromYMD(YEAR, month, 1);
  var monthEnd = monthStart.advance(1, 'month');

  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(tileGeom)
    .filter(ee.Filter.eq('MGRS_TILE', tileId))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', CLOUD_PERCENT));
  var preCol = s2.filterDate(monthStart.advance(-PRE_DAYS, 'day'), monthStart);
  var postCol = s2.filterDate(monthStart, monthEnd.advance(POST_EXTRA_DAYS, 'day'));

  // Indices computed ONCE per composite.
  // Pre-fire: median. Post-fire: the "most burned looking" date per pixel
  // (lowest NBR), so short-lived crop burn scars are not averaged away.
  var pre = addIndices(preCol.map(maskS2).median());
  var post = postCol.map(maskS2).map(function(img) {
    img = addIndices(img);
    return img.addBands(img.select('NBR').multiply(-1).rename('negNBR'));
  }).qualityMosaic('negNBR');

  // All 5 changes point the same way: bigger value = more likely burned.
  var changes = ee.Image.cat([
    pre.select('NBR').subtract(post.select('NBR')).rename('dNBR'),
    post.select('MIRBI').subtract(pre.select('MIRBI')).rename('dMIRBI'),
    post.select('BAIS2').subtract(pre.select('BAIS2')).rename('dBAIS2'),
    pre.select('NDTI').subtract(post.select('NDTI')).rename('dNDTI'),   // NDTI drops after burning
    pre.select('SAVI').subtract(post.select('SAVI')).rename('dSAVI')
  ]);

  // VIIRS mask FIRST: Otsu only sees cropland pixels near active fires.
  var viirsBuffer = viirs.filterDate(monthStart, monthEnd).filterBounds(tileGeom)
    .map(function(f) { return ee.Feature(f).buffer(VIIRS_BUFFER_M, 100); });
  var viirsMask = ee.Image(0).byte().paint(viirsBuffer, 1).selfMask();
  var valid = changes.mask().reduce(ee.Reducer.min()).and(viirsMask).and(cropMaskS2);
  var analysis = changes.updateMask(valid);

  // One 64-bucket histogram call for all 5 indices.
  var hist = analysis.reduceRegion({
    reducer: ee.Reducer.histogram({maxBuckets: 64}),
    geometry: tileGeom, scale: S2_SCALE, crs: S2_CRS, maxPixels: 1e9, tileScale: 4
  });

  var votes = analysis.select('dNBR').gt(safeThreshold(hist, 'dNBR'))
    .add(analysis.select('dMIRBI').gt(safeThreshold(hist, 'dMIRBI')))
    .add(analysis.select('dBAIS2').gt(safeThreshold(hist, 'dBAIS2')))
    .add(analysis.select('dNDTI').gt(safeThreshold(hist, 'dNDTI')))
    .add(analysis.select('dSAVI').gt(safeThreshold(hist, 'dSAVI')));

  // Remove tiny patches. Only burned pixels are processed (selfMask).
  // No reproject(): exports already run at 30 m in EPSG:5070, and forcing
  // 30 m on the map preview is what caused "Object too large".
  var burn = votes.gte(MIN_VOTES).selfMask();
  var minPixels = Math.ceil(MIN_PATCH_M2 / (S2_SCALE * S2_SCALE));
  var patchSize = burn.connectedPixelCount({maxSize: 256, eightConnected: true});
  burn = burn.updateMask(patchSize.gte(minPixels));

  var result = burn.unmask(0).byte().rename('m' + mm);
  var empty = ee.Image(0).byte().rename('m' + mm);
  var hasImages = preCol.size().gt(0).and(postCol.size().gt(0));
  return ee.Image(ee.Algorithms.If(hasImages, result, empty));
}


// Sentinel-2 burn map for the zone: each Sentinel-2 tile gets its own Otsu thresholds.
// (Tiles are read from the Sentinel-2 collection itself, so no tile asset is needed.)
var s2Zone = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterBounds(zone)
  .filterDate(YEAR + '-08-01', YEAR + '-12-16');
var tileIds = s2Zone.aggregate_array('MGRS_TILE').distinct();
var s2Burn = ee.ImageCollection.fromImages(tileIds.map(function(id) {
  var tileGeom = s2Zone.filter(ee.Filter.eq('MGRS_TILE', id)).geometry().dissolve(100).intersection(zone, 100);
  return ee.Image.cat(CMP_MONTHS.map(function(m) { return monthBurn(m, ee.String(id), tileGeom); }))
    .clip(tileGeom);
})).max().unmask(0);

// Burn month 9/10/11 (0 = not burned). A pixel burned twice counts in its first month.
function firstMonth(bandFor) {
  var out = ee.Image(0).byte();
  for (var k = CMP_MONTHS.length - 1; k >= 0; k--) out = out.where(bandFor(CMP_MONTHS[k]).gt(0), CMP_MONTHS[k]);
  return out.rename('month');
}
var monthS2 = firstMonth(function(m) { return s2Burn.select('m' + ('0' + m).slice(-2)); });
function fromBurnDate(collectionId) {
  var col = ee.ImageCollection(collectionId).select('BurnDate')
    .filterDate(YEAR + '-01-01', (YEAR + 1) + '-01-01');
  return firstMonth(function(m) {
    var start = ee.Date.fromYMD(YEAR, m, 1);
    return col.filterDate(start, start.advance(1, 'month')).max().unmask(0);
  });
}
var monthModis = fromBurnDate('MODIS/061/MCD64A1');
var monthCci = fromBurnDate('ESA/CCI/FireCCI/5_1');

// Burned ha and t CO2 per crop and month, all on the same 30 m CDL grid -> 'month|crop|ha|tCO2'
var cropFull = ee.Image('USDA/NASS/CDL/' + YEAR).select('cropland');
function areaByCropMonth(monthImg) {
  var key = cropFull.toInt32().multiply(100).add(monthImg.toInt32()).rename('key');
  var area = ee.Image.pixelArea().divide(1e4).rename('ha')
    .updateMask(cropMaskS2).updateMask(monthImg.gt(0));
  var groups = ee.List(ee.Dictionary(area.addBands(key).reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'key'}),
    geometry: zone, scale: S2_SCALE, crs: S2_CRS, tileScale: 16, maxPixels: 1e13
  })).get('groups', ee.List([])));
  return groups.map(function(g) {
    g = ee.Dictionary(g);
    var k = ee.Number(g.get('key')).round();
    var cropId = k.divide(100).floor().format('%.0f');
    var ha = ee.Number(g.get('sum'));
    var t = ha.multiply(ee.Number(FACTORS.get(cropId, 0)).max(0));
    return k.mod(100).format('%.0f').cat('|').cat(cropId).cat('|')
      .cat(ha.format('%.1f')).cat('|').cat(t.format('%.1f'));
  });
}

// Spatial agreement between the three maps (burned any time Sep-Nov)
var b2 = monthS2.gt(0), bM = monthModis.gt(0), bC = monthCci.gt(0);
function haWhere(mask, name) {
  return ee.Image.pixelArea().divide(1e4).updateMask(cropMaskS2).updateMask(mask).rename(name);
}
var overlap = ee.Image.cat([haWhere(ee.Image(1), 'cropland'), haWhere(b2, 's2'), haWhere(bM, 'modis'),
    haWhere(bC, 'firecci'), haWhere(b2.and(bM), 's2_and_modis'), haWhere(b2.and(bC), 's2_and_firecci'),
    haWhere(bM.and(bC), 'modis_and_firecci')])
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: zone, scale: S2_SCALE, crs: S2_CRS,
                 tileScale: 16, maxPixels: 1e13});

var compare = ee.Dictionary({
  zone_bbox: [-92.1, 29.6, -91.0, 30.75],
  months: CMP_MONTHS,
  s2_tiles: tileIds,
  s2: areaByCropMonth(monthS2),
  modis: areaByCropMonth(monthModis),
  firecci: areaByCropMonth(monthCci),
  overlap_ha: overlap
});
Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(null, {json: ee.String.encodeJSON(compare)})]),
  description: 'CF2020_compare_zoneA', fileFormat: 'CSV', selectors: ['json']
});
print('PART 2: open the Tasks tab and press Run on CF2020_compare_zoneA.');

// Map layers for Testing Zone (a) (zoom in to see Sentinel-2 detail)
var monthPal = {min: 9, max: 11, palette: ['ffff00', 'ff8800', 'ff0000']};
Map.centerObject(zone, 9);
Map.addLayer(zone, {color: 'yellow'}, 'Testing zone (a)', false);
Map.addLayer(monthModis.selfMask().updateMask(cropMaskS2).clip(zone), monthPal, 'Zone (a) MODIS Sep-Nov', false);
Map.addLayer(monthCci.selfMask().updateMask(cropMaskS2).clip(zone), monthPal, 'Zone (a) FireCCI51 Sep-Nov', false);
