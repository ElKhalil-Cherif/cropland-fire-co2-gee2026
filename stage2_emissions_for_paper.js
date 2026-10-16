// =====================================================================
// CROPLAND FIRE CO2 - STAGE 2: EMISSIONS, TABLES AND FIGURES FOR THE PAPER
//
// Reads the burned-area masks saved by Stage 1 (fast: no Sentinel-2
// processing here) and produces everything the paper needs:
//
//   Fig. 3  yearly CO2 map (t CO2/ha) + uncertainty (SD) map   -> GeoTIFF
//   Table 2 / Fig. 4  yearly CO2 per crop with 95% CI           -> CSV
//   Fig. 5  crop area (km2 and %) per crop                       -> CSV
//   Fig. 6  monthly CO2 per crop                                 -> CSV
//           monthly total CO2 with 95% CI (Monte Carlo, N=1000)  -> CSV
//
// Equation (paper, Eq. 2):
//   E (t CO2) = A (ha) x Yield (t/ha) x RPR x DM x CF x EF (kg/t) / 1000
//
// How to use: press Run, then "Run all" in the Tasks tab.
// Files go to Google Drive, folder EXPORT_FOLDER.
// =====================================================================

// ----------------------------- SETTINGS ------------------------------
var YEAR = 2020;
var SCALE = 30;
var CRS = 'EPSG:5070';
var COLLECTION = 'projects/burnished-edge-468212-i6/assets/CF' + YEAR + '_burn';
var EXPORT_FOLDER = 'CroplandFire_' + YEAR;

var DO_TEST_ZONE = true;        // paper results, Testing Zone (a)
var DO_FULL_AOI = false;        // whole basin (turn on after Stage 1 ran for FULL_AOI)
var EXPORT_MAPS = true;         // yearly CO2 + SD GeoTIFF
var EXPORT_MONTHLY_MAPS = false;// 12 monthly CO2 GeoTIFFs (usually not needed)
var SHOW_CHARTS = true;         // charts in the Console for Testing Zone (a)

var COUNT_ONCE_PER_YEAR = true; // a field burned in 2 months is counted once (first month)

// Monte Carlo (paper, section 2.2.4)
var MC_DRAWS = 1000;
var A_REL_SD = 0.20;            // burned-area error (20 %); replace with your validation result
var EF_REL_SD = 0.15;           // emission-factor spread (15 %)
                                // RPR: uniform between the min and max values of Table 1

// ---------------------------- PARAMETERS -----------------------------
// RPR, DM, CF, EF come from Table 1 of the paper (min / max where two sources).
// YIELD (t fresh product / ha) = USDA NASS 2020 U.S. averages, converted from
// bu/ac, lb/ac, cwt/ac or tons/ac.  CHECK these against the yields used in
// your original script / NASS QuickStats before publishing.
var CROPS = [
  {id: 1,  name: 'Corn',           yield: 10.80, rpr: [1.00, 1.04], dm: 0.87,  cf: 0.80, ef: [1412, 1441]},
  {id: 2,  name: 'Cotton',         yield: 0.95,  rpr: [2.00, 2.10], dm: 0.85,  cf: 0.90, ef: [1441, 1441]},
  {id: 3,  name: 'Rice',           yield: 8.54,  rpr: [1.08, 1.40], dm: 0.89,  cf: 0.80, ef: [1378, 1441]},
  {id: 4,  name: 'Sorghum',        yield: 4.59,  rpr: [1.40, 1.46], dm: 0.89,  cf: 0.85, ef: [1441, 1441]},
  {id: 5,  name: 'Soybeans',       yield: 3.43,  rpr: [2.10, 2.22], dm: 0.91,  cf: 0.85, ef: [1376, 1441]},
  {id: 10, name: 'Peanuts',        yield: 4.25,  rpr: [1.00, 1.22], dm: 0.94,  cf: 0.85, ef: [1441, 1441]},
  {id: 13, name: 'Popcorn',        yield: 4.80,  rpr: [1.00, 1.00], dm: 0.87,  cf: 0.80, ef: [1441, 1441]},
  {id: 24, name: 'Winter Wheat',   yield: 3.42,  rpr: [1.30, 1.32], dm: 0.89,  cf: 0.90, ef: [1441, 1441]},
  {id: 28, name: 'Oats',           yield: 2.34,  rpr: [1.30, 1.46], dm: 0.89,  cf: 0.85, ef: [1441, 1441]},
  {id: 29, name: 'Millet',         yield: 1.90,  rpr: [1.40, 1.50], dm: 0.90,  cf: 0.85, ef: [1441, 1441]},
  {id: 45, name: 'Sugarcane',      yield: 77.8,  rpr: [0.27, 0.47], dm: 0.85,  cf: 0.80, ef: [1441, 1441]},
  {id: 46, name: 'Sweet Potatoes', yield: 23.0,  rpr: [0.40, 0.42], dm: 0.255, cf: 0.85, ef: [1441, 1441]}
];

function avg(a) { return (a[0] + a[1]) / 2; }
// Central emission factor per hectare burned (t CO2 / ha)
CROPS.forEach(function(c) {
  c.tPerHa = c.yield * avg(c.rpr) * c.dm * c.cf * avg(c.ef) / 1000;
  // relative SD of one pixel's estimate (independent errors)
  var rprRel = (c.rpr[1] - c.rpr[0]) / Math.sqrt(12) / avg(c.rpr);
  c.relSd = Math.sqrt(A_REL_SD * A_REL_SD + EF_REL_SD * EF_REL_SD + rprRel * rprRel);
});
var cropIds = CROPS.map(function(c) { return c.id; });
var MONTHS = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

// ------------------------------ AREAS --------------------------------
var AOI_COORDS = [[[-89.833008,37.308108],[-91.606201,35.715921],[-91.861084,35.019116],[-92.243408,34.664959],[-92.243408,34.266222],[-91.518311,33.655409],[-91.584229,33.215125],[-92.221436,32.430002],[-91.836914,31.834296],[-92.485108,31.084225],[-92.287354,30.79217],[-93.667236,30.109015],[-91.965910,29.738164],[-90.477905,29.661294],[-90.322998,29.920179],[-90.427368,30.108103],[-90.744324,30.152648],[-90.793213,30.705001],[-91.363349,30.721486],[-91.520508,31.040631],[-90.853638,32.144007],[-89.925293,33.279201],[-89.914307,34.229610],[-89.474854,35.104721],[-88.486084,35.599556],[-88.272949,36.469881],[-88.290528,37.049043],[-88.244385,37.358198],[-88.002686,37.749851],[-89.004639,37.913901],[-89.917603,37.904376],[-90.038452,37.541083],[-89.833008,37.308108]]];
var aoi = ee.Geometry.Polygon(AOI_COORDS, null, false);
var testZone = ee.Geometry.Rectangle([-92.1, 29.6, -91.0, 30.75]).intersection(aoi, 1);

// ------------------------------ INPUTS -------------------------------
var crop = ee.Image('USDA/NASS/CDL/' + YEAR).select('cropland');
var cropMask = crop.remap(cropIds, cropIds.map(function() { return 1; }), 0).eq(1);

var stack = ee.ImageCollection(COLLECTION).max().unmask(0);  // max = tile overlaps counted once
var stackBands = stack.bandNames();
function monthBand(m) {
  var name = 'm' + ('0' + m).slice(-2);
  return ee.Image(ee.Algorithms.If(stackBands.contains(name),
    stack.select(name), ee.Image(0).byte())).rename('b');
}

// burnMonth: 1-12 = month the pixel burned (first month if COUNT_ONCE_PER_YEAR), 0 = not burned
var burnMonth = ee.Image(0).byte();
if (COUNT_ONCE_PER_YEAR) {
  for (var m = 12; m >= 1; m--) burnMonth = burnMonth.where(monthBand(m).eq(1), m);
} else {
  for (var k = 1; k <= 12; k++) burnMonth = burnMonth.where(monthBand(k).eq(1), k);
}
burnMonth = burnMonth.updateMask(cropMask).rename('burn_month');
var burned = burnMonth.gt(0);

// Per-pixel CO2 (t/ha) and its SD
var tPerHa = crop.remap(cropIds, CROPS.map(function(c) { return c.tPerHa; }), 0).toFloat();
var relSd = crop.remap(cropIds, CROPS.map(function(c) { return c.relSd; }), 0).toFloat();
var co2Ha = tPerHa.updateMask(burned).rename('co2_t_ha');
var co2SdHa = tPerHa.multiply(relSd).updateMask(burned).rename('co2_sd_t_ha');

// ---------------------------- STATISTICS -----------------------------
// ONE reduction gives burned hectares for every (crop, month) pair.
function burnedHaByCropMonth(geom) {
  var img = ee.Image.pixelArea().divide(1e4).rename('ha')
    .addBands(crop.toInt32().multiply(100).add(burnMonth.toInt32()).rename('key'))
    .updateMask(burned);
  var groups = ee.List(img.reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'key'}),
    geometry: geom, scale: SCALE, crs: CRS, maxPixels: 1e13, tileScale: 4
  }).get('groups'));
  var keys = groups.map(function(g) { return ee.Number(ee.Dictionary(g).get('key')).format('%.0f'); });
  var sums = groups.map(function(g) { return ee.Dictionary(g).get('sum'); });
  return ee.Dictionary.fromLists(keys, sums);   // '510' -> soybeans, October
}

function cropAreaTable(geom) {
  var groups = ee.List(ee.Image.pixelArea().divide(1e6).rename('km2')
    .addBands(crop.rename('crop'))
    .updateMask(cropMask)
    .reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'crop'}),
      geometry: geom, scale: SCALE, crs: CRS, maxPixels: 1e13, tileScale: 4
    }).get('groups'));
  var d = ee.Dictionary.fromLists(
    groups.map(function(g) { return ee.Number(ee.Dictionary(g).get('crop')).format('%.0f'); }),
    groups.map(function(g) { return ee.Dictionary(g).get('sum'); }));
  var total = ee.Number(d.values().reduce(ee.Reducer.sum()));
  return ee.FeatureCollection(CROPS.map(function(c) {
    var km2 = ee.Number(d.get(String(c.id), 0));
    return ee.Feature(null, {crop_id: c.id, crop: c.name,
      area_km2: km2, percent: km2.divide(total.max(1e-9)).multiply(100)});
  }));
}

function runRegion(label, geom) {
  var ha = burnedHaByCropMonth(geom);
  function A(cid, m) { return ee.Number(ha.get(String(cid * 100 + m), 0)); }

  // Monthly x crop (long table, zeros included -> easy plotting, Fig. 6)
  var monthly = [];
  CROPS.forEach(function(c) {
    MONTHS.forEach(function(m) {
      monthly.push(ee.Feature(null, {crop_id: c.id, crop: c.name, month: m,
        burned_ha: A(c.id, m), co2_t: A(c.id, m).multiply(c.tPerHa)}));
    });
  });
  monthly = ee.FeatureCollection(monthly);

  // Yearly burned area per crop
  var yearHa = {};
  CROPS.forEach(function(c) {
    yearHa[c.id] = ee.Number(ee.List(MONTHS.map(function(m) { return A(c.id, m); }))
      .reduce(ee.Reducer.sum()));
  });

  // ---- Monte Carlo: independent draws per crop ----
  var draws = ee.FeatureCollection(ee.List.sequence(1, MC_DRAWS).map(function(i) {
    return ee.Feature(null, {draw: i});
  }));
  CROPS.forEach(function(c, j) {
    draws = draws.randomColumn('zA_' + c.id, 100 + j, 'normal')
                 .randomColumn('zE_' + c.id, 200 + j, 'normal')
                 .randomColumn('uR_' + c.id, 300 + j, 'uniform');
  });
  draws = draws.map(function(f) {
    var out = {};
    var total = ee.Number(0);
    var monthTot = MONTHS.map(function() { return ee.Number(0); });
    CROPS.forEach(function(c) {
      var aF = ee.Number(1).add(ee.Number(f.get('zA_' + c.id)).multiply(A_REL_SD)).max(0);
      var eF = ee.Number(1).add(ee.Number(f.get('zE_' + c.id)).multiply(EF_REL_SD)).max(0);
      var rpr = ee.Number(f.get('uR_' + c.id)).multiply(c.rpr[1] - c.rpr[0]).add(c.rpr[0]);
      var perHa = rpr.multiply(c.yield * c.dm * c.cf * avg(c.ef) / 1000).multiply(eF).multiply(aF);
      var e = yearHa[c.id].multiply(perHa);
      out['E_' + c.id] = e;
      total = total.add(e);
      MONTHS.forEach(function(m, i) { monthTot[i] = monthTot[i].add(A(c.id, m).multiply(perHa)); });
    });
    out.E_total = total;
    MONTHS.forEach(function(m, i) { out['M_' + m] = monthTot[i]; });
    return ee.Feature(null, out);
  });

  function mcStats(prop) {
    var vals = draws.aggregate_array(prop);
    return {
      mean: vals.reduce(ee.Reducer.mean()),
      sd: vals.reduce(ee.Reducer.stdDev()),
      p025: vals.reduce(ee.Reducer.percentile([2.5])),
      p975: vals.reduce(ee.Reducer.percentile([97.5]))
    };
  }

  // Table 2 / Fig. 4: yearly per crop + total, with 95 % CI
  var yearly = CROPS.map(function(c) {
    var s = mcStats('E_' + c.id);
    return ee.Feature(null, {crop_id: c.id, crop: c.name,
      burned_ha: yearHa[c.id], co2_t: yearHa[c.id].multiply(c.tPerHa),
      co2_mc_mean: s.mean, co2_mc_sd: s.sd, co2_ci95_low: s.p025, co2_ci95_high: s.p975});
  });
  var totalHa = ee.Number(ee.List(CROPS.map(function(c) { return yearHa[c.id]; })).reduce(ee.Reducer.sum()));
  var totalCo2 = ee.Number(ee.List(CROPS.map(function(c) {
    return yearHa[c.id].multiply(c.tPerHa); })).reduce(ee.Reducer.sum()));
  var st = mcStats('E_total');
  yearly.push(ee.Feature(null, {crop_id: 0, crop: 'TOTAL', burned_ha: totalHa, co2_t: totalCo2,
    co2_mc_mean: st.mean, co2_mc_sd: st.sd, co2_ci95_low: st.p025, co2_ci95_high: st.p975}));
  yearly = ee.FeatureCollection(yearly);

  // Fig. 6 confidence band: monthly total with 95 % CI
  var monthlyCi = ee.FeatureCollection(MONTHS.map(function(m) {
    var s = mcStats('M_' + m);
    var central = ee.Number(ee.List(CROPS.map(function(c) {
      return A(c.id, m).multiply(c.tPerHa); })).reduce(ee.Reducer.sum()));
    return ee.Feature(null, {month: m, co2_t: central, co2_mc_mean: s.mean,
      co2_ci95_low: s.p025, co2_ci95_high: s.p975});
  }));

  var areas = cropAreaTable(geom);

  // ---- Exports (CSV to Drive) ----
  var tag = label + '_' + YEAR;
  Export.table.toDrive({collection: yearly, description: 'Table2_yearly_CO2_by_crop_' + tag,
    folder: EXPORT_FOLDER, fileFormat: 'CSV',
    selectors: ['crop_id', 'crop', 'burned_ha', 'co2_t', 'co2_mc_mean', 'co2_mc_sd', 'co2_ci95_low', 'co2_ci95_high']});
  Export.table.toDrive({collection: monthly, description: 'Fig6_monthly_CO2_by_crop_' + tag,
    folder: EXPORT_FOLDER, fileFormat: 'CSV',
    selectors: ['crop_id', 'crop', 'month', 'burned_ha', 'co2_t']});
  Export.table.toDrive({collection: monthlyCi, description: 'Fig6_monthly_total_CI95_' + tag,
    folder: EXPORT_FOLDER, fileFormat: 'CSV',
    selectors: ['month', 'co2_t', 'co2_mc_mean', 'co2_ci95_low', 'co2_ci95_high']});
  Export.table.toDrive({collection: areas, description: 'Fig5_crop_area_' + tag,
    folder: EXPORT_FOLDER, fileFormat: 'CSV',
    selectors: ['crop_id', 'crop', 'area_km2', 'percent']});

  // ---- Exports (GeoTIFF to Drive) ----
  if (EXPORT_MAPS) {
    Export.image.toDrive({
      image: co2Ha.addBands(co2SdHa).addBands(burnMonth.toFloat()).clip(geom),
      description: 'Fig3_CO2_map_' + tag, folder: EXPORT_FOLDER,
      region: geom, crs: CRS, scale: SCALE, maxPixels: 1e13,
      fileFormat: 'GeoTIFF', formatOptions: {cloudOptimized: true}
    });
  }
  if (EXPORT_MONTHLY_MAPS) {
    MONTHS.forEach(function(m) {
      Export.image.toDrive({
        image: tPerHa.updateMask(burnMonth.eq(m)).rename('co2_t_ha').clip(geom),
        description: 'CO2_map_' + tag + '_M' + ('0' + m).slice(-2), folder: EXPORT_FOLDER,
        region: geom, crs: CRS, scale: SCALE, maxPixels: 1e13, fileFormat: 'GeoTIFF'
      });
    });
  }

  return {yearly: yearly, monthly: monthly, areas: areas, monthlyCi: monthlyCi};
}

// ------------------------------- RUN ---------------------------------
print('Emission factor per burned hectare (t CO2/ha):',
  CROPS.map(function(c) { return c.name + ': ' + c.tPerHa.toFixed(2); }));

Map.setCenter(-91.55, 30.2, 9);
Map.addLayer(crop.updateMask(cropMask).clip(testZone), {}, 'CDL crops', false);
Map.addLayer(co2Ha.clip(DO_FULL_AOI ? aoi : testZone),
  {min: 0, max: 40, palette: ['ffffb2', 'fecc5c', 'fd8d3c', 'f03b20', 'bd0026']},
  'CO2 ' + YEAR + ' (t/ha)');
Map.addLayer(testZone, {color: 'yellow'}, 'Testing zone (a)', false);

if (DO_TEST_ZONE) {
  var zone = runRegion('TestZoneA', testZone);
  if (SHOW_CHARTS) {
    print(ui.Chart.feature.byFeature(zone.yearly.filter(ee.Filter.neq('crop', 'TOTAL')), 'crop', ['co2_t'])
      .setChartType('ColumnChart')
      .setOptions({title: 'Fig. 4  Yearly CO2 by crop, ' + YEAR + ' (t)', legend: {position: 'none'}}));
    print(ui.Chart.feature.byFeature(zone.areas, 'crop', ['area_km2'])
      .setChartType('ColumnChart')
      .setOptions({title: 'Fig. 5  Crop area in Testing Zone (a) (km2)', legend: {position: 'none'}}));
    print(ui.Chart.feature.groups(zone.monthly, 'month', 'co2_t', 'crop')
      .setChartType('LineChart')
      .setOptions({title: 'Fig. 6  Monthly CO2 by crop, ' + YEAR + ' (t)'}));
    print('Table 2 (Testing Zone a):', zone.yearly);
  }
}
if (DO_FULL_AOI) {
  runRegion('FullAOI', aoi);
}
print('Done building tasks. Open the Tasks tab and press "Run all".');
