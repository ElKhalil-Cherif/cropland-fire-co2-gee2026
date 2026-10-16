// CROP TYPE OF BURNED CROPLAND, 2020 (MCD64A1 x CDL) - for the paper map
// Press Run, then copy the CROP_MAP_URL line from the Console and paste it to Claude.
var AOI_COORDS = [[[-89.833008,37.308108],[-91.606201,35.715921],[-91.861084,35.019116],[-92.243408,34.664959],[-92.243408,34.266222],[-91.518311,33.655409],[-91.584229,33.215125],[-92.221436,32.430002],[-91.836914,31.834296],[-92.485108,31.084225],[-92.287354,30.792170],[-93.667236,30.109015],[-91.965910,29.738164],[-90.477905,29.661294],[-90.322998,29.920179],[-90.427368,30.108103],[-90.744324,30.152648],[-90.793213,30.705001],[-91.363349,30.721486],[-91.520508,31.040631],[-90.853638,32.144007],[-89.925293,33.279201],[-89.914307,34.229610],[-89.474854,35.104721],[-88.486084,35.599556],[-88.272949,36.469881],[-88.290528,37.049043],[-88.244385,37.358198],[-88.002686,37.749851],[-89.004639,37.913901],[-89.917603,37.904376],[-90.038452,37.541083],[-89.833008,37.308108]]];
var aoi = ee.Geometry.Polygon(AOI_COORDS);
// classes: 1 soybeans, 2 rice, 3 sugarcane, 4 corn (incl. popcorn), 5 cotton, 6 other crops
var crop = ee.Image('USDA/NASS/CDL/2020').select('cropland');
var cls = crop.remap([5, 3, 45, 1, 13, 2, 4, 10, 24, 28, 29, 46], [1, 2, 3, 4, 4, 5, 6, 6, 6, 6, 6, 6], 0).selfMask();
var burned = ee.ImageCollection('MODIS/061/MCD64A1').filterDate('2020-01-01', '2021-01-01')
  .select('BurnDate').max().gt(0);
// Crop type of each burned pixel (the map is drawn at about 500 m, so no reprojection is needed)
var burnedCrop = cls.updateMask(burned).clip(aoi).rename('crop');
var vis = {min: 1, max: 6, palette: ['8c2d04', 'ff9f1c', '6a3d9a', 'ffd23f', '1f78b4', '9e9e9e']};
Map.setCenter(-90.8, 33.5, 6);
Map.addLayer(burnedCrop, vis, 'Burned cropland by crop type, 2020');
print('CROP_MAP_URL', burnedCrop.visualize(vis).getThumbURL({region: aoi, crs: 'EPSG:4326', dimensions: 1600, format: 'png'}));
