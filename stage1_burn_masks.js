// =====================================================================
// CROPLAND FIRE CO2 - STAGE 1: BURNED-AREA MASKS
// Lower Mississippi River Basin, Sentinel-2 + VIIRS
//
// What this script does:
//   For every Sentinel-2 tile and every month it finds burned cropland
//   pixels (Sentinel-2 change + Otsu thresholds, only near VIIRS fires)
//   and saves them ONCE as an Earth Engine asset.
//   Stage 2 then reads these saved masks to make all tables and maps,
//   so the heavy work is never repeated.
//
// How to use:
//   1. Run with TEST_MODE = true. A burned-area layer for one tile and
//      one month appears on the map in about a minute. Check that it
//      looks right (burned fields in red).
//   2. Set TEST_MODE = false and REGION = 'TEST_ZONE', press Run, then
//      press "Run all" in the Tasks tab. One task per tile (12 months
//      in each), instead of one task per tile per month.
//   3. When those tasks finish, run Stage 2.
//   4. Repeat step 2 with REGION = 'FULL_AOI' for the whole basin.
// =====================================================================

// ----------------------------- SETTINGS ------------------------------
var YEAR = 2020;
var TEST_MODE = true;           // true = quick map preview of one tile, no exports
var TEST_MONTH = 10;            // month shown in TEST_MODE (10 = October)
var REGION = 'TEST_ZONE';       // 'TEST_ZONE' (paper, area a) or 'FULL_AOI'
var MONTHS = [9, 10, 11];       // months exported (Sep-Nov = peak burning). Use all 12 for a full year

var SCALE = 30;                 // output pixel size (m)
var CRS = 'EPSG:5070';          // CONUS Albers equal-area
var MIN_PATCH_M2 = 3000;        // smallest burned patch kept (3000 m2 = 4 pixels of 30 m)
var VIIRS_BUFFER_M = 3000;      // search distance around each VIIRS fire point
var CLOUD_PERCENT = 20;         // scene cloud filter (paper, section 2.2.4: < 20 %)
var MIN_VOTES = 4;              // pixel is burned when >= 4 of 5 indices agree
var PRE_DAYS = 30;              // pre-fire window before the month starts
var POST_EXTRA_DAYS = 15;       // post-fire window runs to 15 days after month end

var OUTPUT_ROOT = 'projects/burnished-edge-468212-i6/assets/';
var COLLECTION = OUTPUT_ROOT + 'CF' + YEAR + '_burn';   // Stage 2 reads this

// ------------------------------ AREAS --------------------------------
var AOI_COORDS = [[[-89.833008,37.308108],[-91.606201,35.715921],[-91.861084,35.019116],[-92.243408,34.664959],[-92.243408,34.266222],[-91.518311,33.655409],[-91.584229,33.215125],[-92.221436,32.430002],[-91.836914,31.834296],[-92.485108,31.084225],[-92.287354,30.79217],[-93.667236,30.109015],[-91.965910,29.738164],[-90.477905,29.661294],[-90.322998,29.920179],[-90.427368,30.108103],[-90.744324,30.152648],[-90.793213,30.705001],[-91.363349,30.721486],[-91.520508,31.040631],[-90.853638,32.144007],[-89.925293,33.279201],[-89.914307,34.229610],[-89.474854,35.104721],[-88.486084,35.599556],[-88.272949,36.469881],[-88.290528,37.049043],[-88.244385,37.358198],[-88.002686,37.749851],[-89.004639,37.913901],[-89.917603,37.904376],[-90.038452,37.541083],[-89.833008,37.308108]]];
var aoi = ee.Geometry.Polygon(AOI_COORDS, null, false);

// Testing Zone (a) from the paper's Figure 1: Lafayette to Baton Rouge, south Louisiana
var testZone = ee.Geometry.Rectangle([-92.1, 29.6, -91.0, 30.75]).intersection(aoi, 1);
var region = REGION === 'FULL_AOI' ? aoi : testZone;
var TEST_POINT = ee.Geometry.Point([-91.75, 30.05]);  // tile used in TEST_MODE

// ------------------------------ INPUTS -------------------------------
// 12 crops of the paper (CDL codes)
var cropIds = [1, 2, 3, 4, 5, 10, 13, 24, 28, 29, 45, 46];
var crop = ee.Image('USDA/NASS/CDL/' + YEAR).select('cropland');
var cropMask = crop.remap(cropIds, cropIds.map(function() { return 1; }), 0).eq(1);

var viirs = ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/NOAA20_2020')
  .merge(ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/VIIRS_DATA_2020/SUIOMI2020'))
  .filterBounds(region);

var tiles = ee.FeatureCollection('projects/ee-cherifelkhalil-etu/assets/S2/sentinel_2_ID_Tiles')
  .filterBounds(region);

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
  var valid = changes.mask().reduce(ee.Reducer.min()).and(viirsMask).and(cropMask);
  var analysis = changes.updateMask(valid);

  // One 64-bucket histogram call for all 5 indices.
  var hist = analysis.reduceRegion({
    reducer: ee.Reducer.histogram({maxBuckets: 64}),
    geometry: tileGeom, scale: SCALE, crs: CRS, maxPixels: 1e9, tileScale: 4
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
  var minPixels = Math.ceil(MIN_PATCH_M2 / (SCALE * SCALE));
  var patchSize = burn.connectedPixelCount({maxSize: 256, eightConnected: true});
  burn = burn.updateMask(patchSize.gte(minPixels));

  var result = burn.unmask(0).byte().rename('m' + mm);
  var empty = ee.Image(0).byte().rename('m' + mm);
  var hasImages = preCol.size().gt(0).and(postCol.size().gt(0));
  return ee.Image(ee.Algorithms.If(hasImages, result, empty));
}

// All months of one tile, as one multi-band image.
function tileImage(tileId, months) {
  var tileGeom = tiles.filter(ee.Filter.eq('Name', tileId)).geometry().intersection(region, 1);
  var bands = months.map(function(m) { return monthBurn(m, tileId, tileGeom); });
  return {
    image: ee.Image.cat(bands).clip(tileGeom)
      .set({year: YEAR, tile_id: tileId, scale_m: SCALE}),
    geometry: tileGeom
  };
}

// ------------------------------- RUN ---------------------------------
// setCenter needs no server call (centerObject on a computed area could fail with a network error)
Map.setCenter(-91.55, 30.2, 9);
Map.addLayer(aoi, {color: 'white'}, 'Full study area', false);
Map.addLayer(testZone, {color: 'yellow'}, 'Testing zone (a)', false);

if (TEST_MODE) {
  // Quick visual check of one tile and one month. No export, no waiting.
  ee.Feature(tiles.filterBounds(TEST_POINT).first()).get('Name').evaluate(function(tileId) {
    print('TEST_MODE: tile', tileId, 'month', TEST_MONTH);
    var t = tileImage(tileId, [TEST_MONTH]);
    // Preview only a 30 km x 30 km box (thresholds still come from the whole tile).
    var box = TEST_POINT.buffer(15000).bounds();
    Map.setCenter(-91.75, 30.05, 12);
    Map.addLayer(box, {color: 'yellow'}, 'Preview box', false);
    Map.addLayer(crop.updateMask(cropMask).clip(box), {}, 'CDL crops', false);
    Map.addLayer(t.image.selfMask().clip(box), {palette: ['ff0000']}, 'Burned ' + YEAR + '-' + TEST_MONTH);
    print('Keep the map zoomed in (zoom 11 or closer) to see the red burned fields.');
    print('Burned area inside the preview box (ha):',
      ee.Image.pixelArea().divide(1e4).updateMask(t.image.select(0))
        .reduceRegion({reducer: ee.Reducer.sum(), geometry: box,
                       scale: SCALE, crs: CRS, maxPixels: 1e9, tileScale: 4}));
  });
} else {
  // Create the output collection once (ignored if it already exists).
  try { ee.data.createAsset({type: 'ImageCollection'}, COLLECTION); } catch (e) {}

  var months = MONTHS;
  // Only ONE request to the server here; everything else runs in the tasks.
  tiles.aggregate_array('Name').distinct().sort().evaluate(function(ids) {
    print('Tiles to export (' + ids.length + ' tasks):', ids);
    ids.forEach(function(tileId) {
      var t = tileImage(tileId, months);
      var name = 'CF' + YEAR + '_' + tileId + (REGION === 'FULL_AOI' ? '' : '_zoneA');
      Export.image.toAsset({
        image: t.image,
        description: name,
        assetId: COLLECTION + '/' + name,
        region: t.geometry,
        crs: CRS,
        scale: SCALE,
        maxPixels: 1e13,
        pyramidingPolicy: {'.default': 'max'}
      });
    });
  });
}
