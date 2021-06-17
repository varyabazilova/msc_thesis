// https://code.earthengine.google.com/66b11e6161403f6a030f1b1d190ccc34


// prototype time-series inspector App: sourse code
// draw a polygon around the lakes, see the time-series chart 

// App link: 
//https://varyabazilova.users.earthengine.app/view/glacierlakestestapp

/// find water on Landsat images, recode to 1/0, summ the water occurance, make time-series plots (based on Landsat 8 data)
/// water detection is based on tresholds on NDWI, mNDWI, NDVI, band ratios (to exclude snow), vis Blue band (to exclude shadows) and slope
/// GLIMS dataset is used to construct a buffer around glaciers (=distance treshold for 'glacier lakes')

// import color pallets 
var palettes = require('users/gena/packages:palettes');

/// SRTM
// import SRTM, generate 'slope' from 'elevation'
var srtm = ee.Image('USGS/SRTMGL1_003');  // resolution 30md
var elevation = srtm.select('elevation'); // select 'elevation' bamd from SRTM
var slope = ee.Terrain.slope(elevation);  // calculate slope for each pixel based on 'elevation'

// Compute Mean from slope (3x3 pixels)
var slope_mean = slope.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: ee.Kernel.square({radius: 3, units: 'pixels', normalize: true})});



// add glims -> buffer 
var dataset = ee.FeatureCollection('GLIMS/current');
var glaciers = ee.Image().float().paint(dataset, 'db_area');

//buffer around glaciers 
var buffer1 = ee.Image(1).cumulativeCost({source: glaciers, maxDistance: 5000}).lt(5000); // 2000 - 2km buffer
//var palette5 = palettes.colorbrewer.YlGnBu[9].slice(5);
var params5 = {palette:palettes.colorbrewer.YlGnBu[9].slice(5), opacity: 0.5};

var buffer2 = ee.Image(1).cumulativeCost({source: glaciers, maxDistance: 2000}).lt(2000); // 500 - 0.5km buffer
//var palette2 = palettes.colorbrewer.YlGnBu[9].slice(6);
var params2 = {palette: palettes.colorbrewer.YlGnBu[9].slice(6), opacity: 0.5};

var buffer05 = ee.Image(1).cumulativeCost({source: glaciers, maxDistance: 500}).lt(500); // 500 - 0.5km buffer
//var palette05 = palettes.colorbrewer.YlGnBu[9].slice(7);
var params05 = {palette: palettes.colorbrewer.YlGnBu[9].slice(7), opacity: 0.5};

// add Glaciers (as raster) to the map
//var palette0 = palettes.colorbrewer.YlGnBu[9].slice(8);
var paramsgl = {palette: palettes.colorbrewer.YlGnBu[9].slice(8), opacity: 0.5};




//// ----------- Landsat 8 collection:

// Landsat 8 collection (surface reflectance)

var l8collection = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")      
//  .filter(ee.Filter.eq('WRS_PATH', 171))                  // filter path
//  .filter(ee.Filter.eq('WRS_ROW', 31))                    // filter row
  .filterBounds(caucasus)
  .filter(ee.Filter.lte('CLOUD_COVER', 15))                // filter cloud cover (should the threshold be different?)
//  .filter(ee.Filter.calendarRange(2013, 2018,'year'))      // filter years 
  .filter(ee.Filter.calendarRange(6, 9,'month'));          // filter month
//print('filtered collection info', l8collection); // print filtered collection info


// define functions using Landsat Bands 

// 1. Add 'NDWI_green_nir' as a new band to every img in a collection

var addNDWI1 = function(img) {
  var ndwi_band = img.normalizedDifference(['B3', 'B5']).rename('NDWI_1'); // vis green and NIR
  return img.addBands(ndwi_band); };
//  2. Add 'NDWI_nir_blue' as a new band to every img in a collection
var addNDWI2 = function(img) {
  // var ndwi_band2 = img.normalizedDifference(['B5', 'B2']).rename('NDWI_2'); // NIR and vis Blue
  var ndwi_band2 = img.normalizedDifference(['B3', 'B6']).rename('NDWI_2'); // vis green and SWIR - mNDWI
  return img.addBands(ndwi_band2); };  
/// 3. Green/Nir (band ratio 1)
var addRatio1 = function(img) {
  var ratio1 = img.select('B3').divide(img.select('B5')).rename('ratio1'); // vis Green/NIR
  return img.addBands(ratio1); };  
/// 4. NIR/SWIR [b5/b6]
var addRatio2 = function(img) {
  var ratio2 = img.select('B5').divide(img.select('B6')).rename('ratio2'); // NIR and SWIR
  return img.addBands(ratio2); };
/// 5. NDVI
var addNDVI = function(img) {
  var ndvi_band = img.normalizedDifference(['B5', 'B4']).rename('NDVI'); 
  return img.addBands(ndvi_band); };
  
  
/// create new collection with new bands 
var collection = l8collection.map(addNDWI1);  // 1.  NormDiff(vis green and NIR)
var collection = collection.map(addNDWI2);    // 2. NormDiff(vis green and SWIR)
var collection = collection.map(addRatio1);   // 3. vis Green/NIR
var collection = collection.map(addRatio2);   // 4. NIR / SWIR
var collection = collection.map(addNDVI);     // 5. NDVI

//  Add 'WATER' as a new band to every img in an updated collection

var maskNDWI = function(img) {
  var ndwi1 = img.select('NDWI_1');   // vis green and NIR (Norm Diff B3-B5)
  var ndwi2 = img.select('NDWI_2');   // NIR and vis Blue (Norm Diff B5-B2)
  var ndvi = img.select('NDVI');      // NIR and vis Red (Norm Diff B5-B4)
  var blue = img.select('B2');        // vis Blue (B2)
  var rat2 = img.select('ratio2');    // NIR and vis Blue (B5/B6)
  var mask = img.select('NDWI_1')
  .updateMask(ndwi1.gte(0.3)  
  .and(ndvi.lte(-0.2))  
  .and(blue.gte(300))  
  .and(slope_mean.lte(25)))  /// ??? 10 - [Gardelle, 2011]
  .rename('water1');
  return img.addBands(mask); };
var collection = collection.map(maskNDWI);  // Map over collection
print('L8 with WATER masked collection', collection); // print info 

/// set all WATER pixels to 1
//[water = 1], [everything else = 0]
var toBoolean = function(img) {
  var water = img.select('water1');
  var value = water.updateMask(water.neq(0)).multiply(0).add(1).rename('water');
  return img.addBands(value); };
var collection = collection.map(toBoolean);
//print(collection, 'binary code to water');


// pick collection only with NDWI_masked band 
var WaterCollection = collection.select('water'); 
// reduce collection -> summ all pixels 
var sum = WaterCollection.reduce(ee.Reducer.sum());


// ------------ LAYOUT and widgets for the app  

/*
 * Map setup
 */

var paletteLakes = palettes.colorbrewer.Spectral[11];
//Map.addLayer(sum, {palette:paletteLakes, min:1, max:25}, 'water occurence');


var vis = {min:1, max:25, palette:paletteLakes}; 
var composite = sum.visualize(vis);
var compositeLayer = ui.Map.Layer(composite).setName('water occurence');

var buffer11 = buffer1.mask(buffer1).visualize(params5);
var compositeBuffer5 = ui.Map.Layer(buffer11).setName('buffer - 5 km');

var buffer21 = buffer2.mask(buffer2).visualize(params2);
var compositeBuffer2 = ui.Map.Layer(buffer21).setName('buffer - 2 km');

var buffer051 = buffer05.mask(buffer05).visualize(params05);
var compositeBuffer05 = ui.Map.Layer(buffer051).setName('buffer - 0.5 km');

var glaciers1 = glaciers.visualize(paramsgl);
var compositeGlaciers = ui.Map.Layer(glaciers1).setName('Glaciers (source: GLIMS)');


// Create the main map and set layers.
var mapPanel = ui.Map();
var layers = mapPanel.layers();
layers.add(compositeBuffer5, 'buffer - 5 km');
layers.add(compositeBuffer2, 'buffer - 2 km');
layers.add(compositeBuffer05, 'buffer - 0.5 km');
layers.add(compositeGlaciers, 'Glaciers (source: GLIMS)');

layers.add(compositeLayer, 'water occurence');


/*
 * Panel setup
 */

// Create a panel to hold title, intro text, chart and legend components.
var toolPanel = ui.Panel({style: {width: '25%'}});

// Create an intro panel with labels.
var intro = ui.Panel([
  ui.Label({ 
    value: 'Glacier lakes area - Time Series Inspector',
    style: {fontSize: '18px', fontWeight: 'bold'}}),
  ui.Label({
    value: 'Create a polygon around the lake to see how water extent is changing in time',
    style: {fontSize: '16px'}}),
  ui.Label({
    value: 'Data source: Landsat 8',
    style: {fontSize: '16px'}}),
]);
toolPanel.add(intro);

ui.root.widgets().add(toolPanel); // keep that 

// Add placeholders for the chart and legend.
//toolPanel.add(ui.Label('[Chart]'));
//toolPanel.add(ui.Label('[Legend]'));

/*
 * Legend setup
 */

// Creates a color bar thumbnail image for use in legend from the given color
function makeColorBarParams(palette) {
  return {bbox: [0, 0, 25, 0.1], dimensions: '100x10', min: 1, max: 25, palette: palette};}
// Create the color bar for the legend.
var colorBar = ui.Thumbnail({
  image: ee.Image.pixelLonLat().select(0),
  params: makeColorBarParams(vis.palette),
  style: {stretch: 'horizontal', margin: '0px 8px', maxHeight: '24px'},});
// Create a panel with three numbers for the legend.
var legendLabels = ui.Panel({widgets: [ui.Label(vis.min, {margin: '4px 8px'}),
ui.Label((vis.max / 2), {margin: '4px 8px', textAlign: 'center', stretch: 'horizontal'}),
ui.Label(vis.max, {margin: '4px 8px'})],
layout: ui.Panel.Layout.flow('horizontal')});
var legendTitle = ui.Label({value: 'Map Legend: water occurence', style: {fontWeight: 'bold'}});

// Add the legendPanel to the map.

var legendPanel = ui.Panel([legendTitle, colorBar, legendLabels]);
toolPanel.widgets().set(3, legendPanel);



//// Chart setup for the polygon 

// Don't make imports that correspond to the drawn points.
mapPanel.drawingTools().setLinked(false);
// Limit the draw modes to points.
mapPanel.drawingTools().setDrawModes(['polygon']);
// Add an empty layer to hold the drawn points.
mapPanel.drawingTools().addLayer([]);
// Set the geometry type to be polygon.
mapPanel.drawingTools().setShape('polygon');
// Enter drawing mode.
mapPanel.drawingTools().draw();


// This function gets called when the geometry layer changes.
// Use debounce to call the function at most every 100 milliseconds.
var generateChart = ui.util.debounce(function() {
  var roi = mapPanel.drawingTools().layers().get(0).toGeometry();
  var timeseriesChart = ui.Chart.image.series(WaterCollection, roi, ee.Reducer.sum(), 30);

  timeseriesChart.setOptions({
    title: 'Water (area): time series',
    vAxis: {title: 'area, m sq.'},
    hAxis: {title: 'date', gridlines: {count: 7}},
    series: {
      0: {
        color: 'darkblue',
        lineWidth: 1.75,
        pointsVisible: true,
        pointSize: 2,
      },
    },
    legend: {position: 'right'},
  });
  // Add the chart at a fixed position, so that new charts overwrite older ones.
  toolPanel.widgets().set(2, timeseriesChart);
}, 100);

// Register a callback on the default map to be invoked when the map is clicked.
//mapPanel.drawingTools().onEdit(generateChart);
mapPanel.drawingTools().onDraw(generateChart);
mapPanel.drawingTools().onErase(generateChart);


/*
 * Map setup
 */

// Configure the map.
mapPanel.style().set('cursor', 'crosshair');

/*
 * Initialize the app
 */

// Replace the root with a SplitPanel that contains the inspector and map.
ui.root.clear();
ui.root.add(ui.SplitPanel(toolPanel, mapPanel));

mapPanel.setCenter(42.6941, 43.2185, 13);





