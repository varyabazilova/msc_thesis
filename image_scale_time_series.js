// https://code.earthengine.google.com/4582b516041f02bffb1990b3c73ee1de

/// find water on Landsat images, recode to 1/0, summ the water occurance (1s), make time-series plots
/// water detection is based on tresholds on (1) NDWI, (2) mNDWI, (3) NDVI, (4) band ratios (to exclude snow), 
// (5) vis Blue band (to exclude shadows) and (6) slope
/// GLIMS dataset is used to construct a buffer around glaciers (= distance treshold for 'glacier lakes')

/// Landsat 5 and Landsat 8 are imported separately (with 2 sets of thresholds on the same data for water detection) -> 2 "sum" layers are printed (L5 and L8) 
// -> 2 charts for chosen geometry (L5 and L8), corresponding for different time periods 

// high values in the beginning of June are considered snow 

// draw "geometry" polygon to inspect the area of interest 

// ------- 0. Geometries:

var caucasus = /* color: #d63000 */ee.Geometry.LineString(
        [[39.16712597811765, 43.96527833183177],
         [49.27454785311765, 41.1167473568263]]);

var roi = geometry


// -------- I. Classification ----------


// -------- 1. Importing "additional things" --------

// ---- import color pallets 
var palettes = require('users/gena/packages:palettes');
var paletteLakes = palettes.colorbrewer.Spectral[11];

// -------- 1.1 SRTM --------

// import SRTM, generate 'slope' from 'elevation'
var srtm = ee.Image('USGS/SRTMGL1_003');   // resolution 30md
var elevation = srtm.select('elevation');  // select 'elevation' bamd from SRTM
var slope = ee.Terrain.slope(elevation);   // calculate slope for each pixel based on 'elevation'
//Map.addLayer(slope, {min: 0, max: 90}, 'slope'); // print
//Map.addLayer(elevation, {min: 0, max: 90}, 'elevation', false);   // add to map

// compute Mean from slope (3x3 pixels)
var slope_mean = slope.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel:  ee.Kernel.square({radius: 3, units: 'pixels', normalize: true})});
//Map.addLayer(slope_mean, {min: 0, max: 90}, 'slope - mean', false);   // add to map


// -------- 1.2 GLACIERS --------

// add GLIMS dataset  
var dataset = ee.FeatureCollection('GLIMS/current');         // feature collection 
var glaciers = ee.Image().float().paint(dataset, 'db_area'); // convert to raster layer 

// add GLIMS/Glaciers (as raster) to the map
var paramsgl = {palette: palettes.colorbrewer.YlGnBu[9].slice(8), opacity: 0.5}; // visualization for glaciers 
//Map.addLayer(glaciers, paramsgl, 'GLIMS', false); 

// take only glaciers in Caucasus, described and submited to Glims by [Tielidze and Wheate, 2018] - source date = later then 3013
var glims2018 = dataset.filter(ee.Filter.gte('src_date', '2012-01-01'));//.clip(caucasus_pol);   
Map.addLayer(glims2018, null, 'GLIMS [Tielidze - 2018]', false); // vector layer to the map 


// buffer around GLIMS_feature collection

var buffer1km = function(feature) {
  var buff = feature.buffer(1000); 
  return ee.Feature(buff) };

var glac_buff_1km = glims2018.map(buffer1km);

//Map.addLayer(glac_buff_1km, null, 'buffer 1 km', false);



// ------------ 2. Import Landsat data ---------------- 

// ------------ 2.1 Landsat 5 ---------------- 

// Landsat 5 collection (surface reflectance)

var l5collection = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR")  // import data 
  .filterBounds(caucasus)                                        // filter region of interest 
  .filter(ee.Filter.calendarRange(6, 9,'month'))                 // filter month
  .filter(ee.Filter.lte('CLOUD_COVER', 20));                     // filter cloud cover 
  


// ------------ define functions using Landsat5 Bands 

// --------- define functions for L5 and L7

// 1. Add 'NDWI_green_nir' as a new band to every img in a collection
var addNDWI = function(img) {
  var ndwi_band = img.normalizedDifference(['B2', 'B4']).rename('NDWI'); // green and NIR - NDWI
  return img.addBands(ndwi_band); };
//  2. Add 'NDWI_nir_blue' as a new band to every img in a collection
var addNDSI = function(img) {
  var ndsi_band = img.normalizedDifference(['B2', 'B5']).rename('NDSI'); // green SWIR -  mNDWI/NDSI
  return img.addBands(ndsi_band); };  
// 4. Red/SWIR [b5/b6]
var addRatio2 = function(img) {
  var ratio2 = img.select('B3').divide(img.select('B5')).rename('ratio2'); // Red and SWIR
  return img.addBands(ratio2); };
// 5. NDVI
var addNDVI = function(img) {
  var ndvi_band = img.normalizedDifference(['B4', 'B3']).rename('NDVI');   // NIR and red - NDVI
  return img.addBands(ndvi_band); };
  

/// ------------ create new collection with new bands 
var collection5 = l5collection.map(addNDWI);   
var collection5 = collection5.map(addNDSI);    
var collection5 = collection5.map(addRatio2);  
var collection5 = collection5.map(addNDVI);    


// ----------- create a function to update water mask (define thresholds) 

var maskNDWI = function(img) {
  var ndwi = img.select('NDWI');      // vis green and NIR (Norm Diff)
  var ndsi = img.select('NDSI');      // NIR and vis Blue (Norm Diff)
  var ndvi = img.select('NDVI');      // NIR and vis Red (Norm Diff)
  var blue = img.select('B1');        // vis Blue 
  var rat2 = img.select('ratio2');    // NIR and vis Blue
  var mask = img.select('NDWI')
  .updateMask(
    (ndwi.gte(0.24).and(blue.gte(700)).and(blue.lte(2500))//.and(rat1.gte(2)).and(rat2.lte(16))
    .and(slope_mean.lte(25))))
  .rename('water');
  return img.addBands(mask); };
var collection5 = collection5.map(maskNDWI);  // Map over collection


// --------- recode to binary --------- 

// set all WATER pixels to 1
//[water = 1], [everything else = 0]
var toBoolean = function(img) {
  var water = img.select('water');
  var value = water.updateMask(water.neq(0)).multiply(0).add(1).rename('binary_water');
  return img.addBands(value); };
  
var collection5 = collection5.map(toBoolean);
print(collection5, 'Landsat5: binary code - water');


// ------------ 2.2 Landsat 7 ---------------- 

// ------------- Landsat 7 collection (SR)
var l7collection = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
  .filterBounds(caucasus)
  .filterDate('1987-01-01', '2003-05-31')
  .filter(ee.Filter.calendarRange(6, 9,'month'))           // filter month
  .filter(ee.Filter.lte('CLOUD_COVER', 20));               // filter cloud cover (should the threshold be different?)
//print('filtered L7 collection info', l7collection);        // print filtered collection info

/// ------------ create new collection with new bands (using same functions as Landsat5)
/// ------------ create new collection with new bands 
var collection7 = l7collection.map(addNDWI);   
var collection7 = collection7.map(addNDSI);    
var collection7 = collection7.map(addRatio2);  
var collection7 = collection7.map(addNDVI);    


// add 'water' band to the collection (same thresholds as Landsat 5)
var collection7 = collection7.map(maskNDWI);  


// --------- recode to binary 
var collection7 = collection7.map(toBoolean);




// ------------ 2.3 Landsat 8 ---------------- 

// ------------- Landsat 8 collection (SR)
var l8collection = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")      
  .filterBounds(caucasus)
  .filter(ee.Filter.lte('CLOUD_COVER', 20))                // filter cloud cover (should the threshold be different?)
  .filter(ee.Filter.calendarRange(6, 9,'month'));          // filter month
//print('filtered L8 collection info', l8collection); // print filtered collection info


// ------------- define functions using Landsat8 Bands 

// 1. mNDWI
var addMNDWI = function(img) {
  var mndwi_band = img.normalizedDifference(['B3', 'B5']).rename('mNDWI'); // vis green and NIR
  return img.addBands(mndwi_band); };
//  2. Add 'NDWI_nir_blue' as a new band to every img in a collection
var addNDWI = function(img) {
  var ndwi_band = img.normalizedDifference(['B5', 'B2']).rename('NDWI'); // NIR and vis Blue
  return img.addBands(ndwi_band); };  
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
var collection8 = l8collection.map(addMNDWI);   // 1. NormDiff(vis green and NIR)
var collection8 = collection8.map(addNDWI);     // 2. NormDiff(NIR and vis Blue)
var collection8 = collection8.map(addRatio1);   // 3. vis Green/NIR
var collection8 = collection8.map(addRatio2);   // 4. NIR / SWIR
var collection8 = collection8.map(addNDVI);     // 5. NDVI


// ------- define a function to update water mask

var maskNDWI = function(img) {
  var mndwi = img.select('mNDWI');    // vis green and NIR 
  var ndwi = img.select('NDWI');      // NIR and vis Blue
  var ndvi = img.select('NDVI');      // NIR and vis Red 
  var blue = img.select('B2');        // vis Blue 
  var rat1 = img.select('ratio1');    // vis Green/NIR
  var rat2 = img.select('ratio2');    // NIR/SWIR
  
  var mask = img.select('mNDWI')
  .updateMask(
    mndwi.gte(0.25).and(ndwi.lte(0)).and(blue.gte(300)).and(ndvi.lte(-0.1))               //.and(rat1.gte(1.4)).and(rat2.lte(30))
    .and(slope_mean.lte(25)))
  .rename('water');
  return img.addBands(mask); 
};
var collection8 = collection8.map(maskNDWI);  


print('L8 with WATER masked collection', collection8); 

// recode to binary 
/// set all WATER pixels to 1
//[water = 1], [everything else = 0]
var collection8 = collection8.map(toBoolean);
//print(collection8, 'binary code to water');


// ------- 3. Sum all WATER pixels to create "water occurance" layer -------

// make new collection only with 'binary_water' band
var WaterCollection5 = collection5.select('binary_water'); 
var WaterCollection7 = collection7.select('binary_water');
var WaterCollection8 = collection8.select('binary_water'); 

// merge Landsat 5 and Landsat7

var WaterCollection57 = WaterCollection5.merge(WaterCollection7);


// reduce collection -> summ all pixels 

var sum57 = WaterCollection57.reduce(ee.Reducer.sum());
var sum8 = WaterCollection8.reduce(ee.Reducer.sum());

// delete water_occurance = 1 (assume that thats error)

var sum57 = sum57.updateMask(sum57.neq(1));
var sum8  = sum8.updateMask(sum8.neq(1));


// --------- Add Layers to the Map 
//Map.addLayer(sum5, {palette: paletteLakes, min:1, max:100}, 'lake area - L5', false);  
//Map.addLayer(sum7, {palette: paletteLakes, min:1, max:20},  'lake area - L7', false);  
Map.addLayer(sum57, {palette: paletteLakes, min:1, max:120},  'lake area - L57', false);  
Map.addLayer(sum8,  {palette: paletteLakes, min:1, max:30},   'lake area - L8',  false);  




// ------------- 4. Time-Series charts for chosen geometry -------------

// Draw geometry of Interest: 

var roi = geometry

print(ui.Chart.image.series(WaterCollection57, roi, ee.Reducer.sum(), 30), 'Landsat 5-7'); 
print(ui.Chart.image.series(WaterCollection8,  roi, ee.Reducer.sum(), 30), 'Landsat 8');




