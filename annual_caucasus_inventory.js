// https://code.earthengine.google.com/34def5d9d56335ca29e3e5ed672cd4ba

// 1. define functions 
// 2. add bands to the collection
// 3. make an annual scale composite of SR bands and computed bands 
// 4. for map representation on the map chanhe (n) in the YEAR section 

// 5. define and apply thresholds to map water on an annual composites
// 6. make a 5-year average composite of water 



 
// -------- GLACIERS --------

// -------- 1.2 GLACIERS --------

// add GLIMS dataset  
var glims = ee.FeatureCollection('GLIMS/current');  // feature collection 
var glaciers = ee.Image().float().paint(glims);     // convert to raster layer 

// take only glaciers in Caucasus, described and submited to Glims by [Tielidze and Wheate, 2018] - source date = later then 3013
var glims2018 = glims.filter(ee.Filter.gte('src_date', '2012-01-01'));//.clip(caucasus_pol);   
//Map.addLayer(glims2018, null, 'GLIMS [Tielidze - 2018]', false); 

// functions for buffers around every feature 
var buffer05km = function(feature) {
  var buff = feature.buffer(500); // (meters)
  return ee.Feature(buff) };
  
var buffer1km = function(feature) {
  var buff = feature.buffer(1000); 
  return ee.Feature(buff) };

var glaciers05km = glims2018.map(buffer05km);
var glaciers1km = glims2018.map(buffer1km);





// ---------------- GLOBAL VARIABLES -------------------


// -------  names for landsat 5, 7 and landsat 8
var l57_bands = ['B1',   'B2',    'B3',  'B4',  'B5',   'B7'];
var l57_names = ['blue', 'green', 'red', 'nir', 'swir', 'swir2'];


var l8_bands = ['B2',   'B3',    'B4',  'B5',  'B6',   'B7'];
var l8_names = ['blue', 'green', 'red', 'nir', 'swir', 'swir2'];


var s2_bands = ['B2',   'B3',    'B4',  /*'B5',        'B6',        'B7',       */ 'B8',  'B11',  'B12'];
var s2_names = ['blue', 'green', 'red', /*'red_edge1', 'red_edge2', 'red_edge3',*/ 'nir', 'swir', 'swir2']; // swir = swir1

var reducer = ee.Reducer.minMax()
              .combine({
            //reducer2: ee.Reducer.percentile([10]),
              reducer2: ee.Reducer.median(),
              sharedInputs: true});

var yearListl5 = ee.List.sequence(1984, 2011);
var yearListl8 = ee.List.sequence(2013, 2021);


//// lists for everest area: 
//var yearListl5 = ee.List.sequence(1987, 2001).cat(ee.List.sequence(2003, 2011)); // 2002 is missing in Landsat 5
//var yearListl7 = ee.List.sequence(1999, 2002);




// ------------------ create a functions -------------  


// function to filter out failure on the edge of the 
var failure = function(img){
  var b2 = img.select('green') // B2
  var b3 = img.select('red')   // B3
  var b4 = img.select('nir')   // B4
  var b5 = img.select('swir')  // B5
  var mask = img.updateMask(b2.gte(0).and(b3.gte(0)).and(b4.gte(0)).and(b5.gte(0)))
  return(img.updateMask(mask))}

var twentyK = function(img){
  var b1 = img.select('blue') // blue
  var mask = b1.updateMask(b1.neq(20000))
  //var mask = b1.updateMask(b1.lte(19000))
  return img.updateMask(mask)}

// --------- define functions for L5 and L7

// Add 'NDWI_green_nir' as a new band to every img in a collection
var addNDWI = function(img) {
  var ndwi_band = img.normalizedDifference(['green', 'nir']).rename('NDWI'); 
  return img.addBands(ndwi_band); };
// Add 'NDWI_nir_blue' as a new band to every img in a collection
var addNDSI = function(img) {
  var ndsi_band = img.normalizedDifference(['green', 'swir']).rename('NDSI') ; 
  return img.addBands(ndsi_band); };  
// Red/SWIR [b5/b6]
var addRatio2 = function(img) {
  var ratio2 = img.select('red').divide(img.select('swir')).rename('ratio2'); 
  return img.addBands(ratio2); };
//  NDVI
var addNDVI = function(img) {
  var ndvi_band = img.normalizedDifference(['nir', 'red']).rename('NDVI');   
  return img.addBands(ndvi_band); };


// ratios:

//  water-ratio
var addNDWI_ratio = function(img) {
  var ndvi_ratio = img.select('green').divide(img.select('nir')).rename('water_ratio');   
  return img.addBands(ndvi_ratio); };

//  snow-ratio
var addNDSI_ratio = function(img) {
  var ndvi_ratio = img.select('green').divide(img.select('swir')).rename('snow_ratio');   
  return img.addBands(ndvi_ratio); };

//  vegetation-ratio
var addNDVI_ratio = function(img) {
  var ndvi_ratio = img.select('nir').divide(img.select('red')).rename('veg_ratio');   
  return img.addBands(ndvi_ratio); };

// product ratio: 
var product_ratio = function(img){
  var veg = img.select('veg_ratio').multiply(-1)
  var water = img.select('water_ratio')
  var snow = img.select('snow_ratio')
  var band = veg.multiply(water)/*.multiply(snow)*/.rename('ratio_product')
  return img.addBands(band);}

// product of ind
var product = function(img){
  var veg = img.select('NDVI').multiply(-1)
  var water = img.select('NDWI')
  var snow = img.select('NDSI')
  var prod = veg.multiply(water).multiply(snow).rename('ind_product')
  return img.addBands(prod)}

// AWEI(nsh)  
 var addAWEI1 = function(img) {
  var first = img.select('green').subtract(img.select('swir'))
  var first4 = first.multiply(ee.Image(2))
  var second = (img.select('nir').multiply(ee.Image(0.25))).add(img.select('swir2').multiply(ee.Image(2.75)))
  var band = first4.subtract(second).rename('AWEI_nsh');
  return img.addBands(band); }; 
  
// AWEI(sh)  
 var addAWEI2 = function(img) {
  var first = img.select('green').multiply(ee.Image(2.5))
  var second = (img.select('nir').add(img.select('swir'))).multiply(ee.Image(1.5))
  var third = img.select('swir2').multiply(ee.Image(0.25))
  var band = img.select('blue').add(first).add(second).add(third).rename('AWEI_sh')
  return img.addBands(band); }; 
 
  
  

  
// ------------------ IMPORT DATA ----------------- 

// import collection, filter, map functions 

var l5 = ee.ImageCollection("LANDSAT/LT05/C01/T1_SR").select(l57_bands, l57_names).filterBounds(caucasus)
var l7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR").select(l57_bands, l57_names).filterBounds(caucasus).filter(ee.Filter.calendarRange(1999, 2002,'year'))
var l8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR").select(l8_bands, l8_names)  .filterBounds(caucasus)


var l5collection = l5
               // errors:
               .map(failure)
               .map(twentyK)
               // NDI-s
               .map(addNDSI)     
               .map(addNDWI)     
               .map(addNDVI)
               .map(product)
               // ratios
               .map(addNDVI_ratio)
               .map(addNDWI_ratio)
               .map(addNDSI_ratio)
               .map(product_ratio)
               // ratio - ice snow
               .map(addRatio2)
               // AWEI-s
               .map(addAWEI1)
               .map(addAWEI2);

var l7collection = l7 //errors:
               .map(failure).map(twentyK)
               // NDI-s
               .map(addNDSI).map(addNDWI).map(addNDVI).map(product)
               // ratios
               .map(addNDVI_ratio).map(addNDWI_ratio).map(addNDSI_ratio).map(product_ratio)
               // ratio - ice snow
               .map(addRatio2)
               // AWEI-s
               .map(addAWEI1).map(addAWEI2);

var l8collection = l8
               // errors:
               .map(failure)
               .map(twentyK)
               // NDI-s
               .map(addNDSI)     
               .map(addNDWI)     
               .map(addNDVI)
               .map(product)
               // ratios
               .map(addNDVI_ratio)
               .map(addNDWI_ratio)
               .map(addNDSI_ratio)
               .map(product_ratio)
               // ratio - ice snow
               .map(addRatio2)
               // AWEI-s
               .map(addAWEI1)
               .map(addAWEI2);




// --------------------------------------------------------------------------------------------------------------

// -------------  YEARLY COMPOSITE ------------------

// -------------- LANDSAT 5 ---------------

// ------ 1st MIN ---------- 

// map over the list of years to generate a composite for each year.
var yearCompList_1stMin = yearListl5.map(function(year){
  
  var yearCol = l5collection.filter(ee.Filter.calendarRange(year, year, 'year'))
  var yearMin1 = yearCol.reduce(reducer);   // reducer
  
  var imgList = yearCol.aggregate_array('LANDSAT_ID');
  var n_img = imgList.size();
  var nBands = yearMin1.bandNames().size();
  return yearMin1.set({
    'year': year,
    'image_list': imgList,
    'n_bands': nBands,
    'n_img': n_img,
    'satellite': 'Landsat5',
    'method': '1st min',
    'system:time_start': year
  });
});
// create collection: 
// convert the annual composite image list to an ImageCollection and filter out years with no images.
var firstMinl5 = ee.ImageCollection.fromImages(yearCompList_1stMin).filter(ee.Filter.gt('n_bands', 0))
print('firstMinl5f', firstMinl5)





// ------ 2nd MIN ----------

var yearCompList_2ndMin = yearListl5.map(function(year){
  
  var yearCol = l5collection.filter(ee.Filter.calendarRange(year, year, 'year'))
  var yearMin1 = yearCol.reduce(ee.Reducer.min());
  var yearMax1 = yearCol.reduce(ee.Reducer.max());

  // var yearComp = yearCol.reduce(ee.Reducer.percentile([20]));
  var yearCol2 = yearCol.map(function (img) {
    var img_masked=img.updateMask(img.gt(yearMin1).and(img.lt(yearMax1)))
    return ee.Image(img_masked)});
  
  var yearMin2 = yearCol2.reduce(reducer);  // reducer
  
  var imgList = yearCol.aggregate_array('LANDSAT_ID');
  var n_img = imgList.size();
  var nBands = yearMin2.bandNames().size();
  
  return yearMin2.set({
    'year': year,
    'image_list': imgList,
    'n_bands': nBands,
    'n_img': n_img, 
    'satellite': 'Landsat5',
    'method': '2nd min',
    'system:time_start': year
  });
});
// create collection: 
var secondMinl5 = ee.ImageCollection.fromImages(yearCompList_2ndMin).filter(ee.Filter.gt('n_bands', 0))





// ----- 3 SIGMA empirical rule ------------
// 1. compute sigma 2. update mask
// rule: 95 of data fall in (MEAN-2sigma) range
// to perform .multiply(): make sure there are no 0-band images 


var year_l5_sigma = yearListl5.map(function(year){
  
  var yearCol = l5collection.filter(ee.Filter.calendarRange(year, year, 'year'))
  var yearMean = yearCol.reduce(ee.Reducer.mean());
  var yearSD = yearCol.reduce(ee.Reducer.stdDev());
  var year2SD = yearSD.multiply(ee.Image(2));
  var year_2sigma = yearMean.subtract(year2SD);
  var year_2SD_plus = yearMean.add(year2SD);
  var year_2SD_minus = yearMean.subtract(year2SD);
  
 // update mask ->[(mean-2sigma) < pixel < (mean+2sigma)]
 
  var yearCol2sigma = yearCol.map(function (img) {
    var img_masked = img.updateMask(img.gte(year_2SD_minus).and(img.lte(year_2SD_plus)))
    return ee.Image(img_masked)});
  
  var year2sigmal5 = yearCol2sigma.reduce(reducer);  // reducer

  var imgList = yearCol.aggregate_array('LANDSAT_ID');
  var n_img = imgList.size();
  var nBands = year2sigmal5.bandNames().size();
  
  return year2sigmal5.set({
    'year': year,
    'image_list': imgList,
    'n_bands': nBands,
    'n_img': n_img,
    'satellite': 'Landsat5',
    'method': 'sigma',
    'system:time_start': year
  });
});

// convert the annual composite image list to an ImageCollection and filter out years with no images.
var year2sigmal5 = ee.ImageCollection.fromImages(year_l5_sigma).filter(ee.Filter.gt('n_bands', 0))
print('year2sigma', year2sigmal5)



// --------------------------------------------------------------------------------------------------------------

// -------------- LANDSAT 8 ---------------

// ------ 1st MIN ---------- 

// map over the list of years to generate a composite for each year.
var year_l8_1stMin = yearListl8.map(function(year){
  
  var yearCol = l8collection.filter(ee.Filter.calendarRange(year, year, 'year'))
  var yearMin1 = yearCol.reduce(reducer);   // reducer
  
  var imgList = yearCol.aggregate_array('LANDSAT_ID');
  var n_img = imgList.size();
  var nBands = yearMin1.bandNames().size();
  return yearMin1.set({
    'year': year,
    'image_list': imgList,
    'n_bands': nBands,
    'n_img': n_img,
    'satellite': 'Landsat8',
    'method': '1st min',
    'system:time_start': year
  });
});
// create collection: 
// convert the annual composite image list to an ImageCollection and filter out years with no images.
var firstMinl8 = ee.ImageCollection.fromImages(year_l8_1stMin).filter(ee.Filter.gt('n_bands', 0))
print('firstMinl5f', firstMinl8)


// ------ 2nd MIN ----------

var year_l8_2ndMin = yearListl8.map(function(year){
  
  var yearCol = l8collection.filter(ee.Filter.calendarRange(year, year, 'year'))
  var yearMin1 = yearCol.reduce(ee.Reducer.min());
  var yearMax1 = yearCol.reduce(ee.Reducer.max());

  // var yearComp = yearCol.reduce(ee.Reducer.percentile([20]));
  var yearCol2 = yearCol.map(function (img) {
    var img_masked=img.updateMask(img.gt(yearMin1).and(img.lt(yearMax1)))
    return ee.Image(img_masked)});
  
  var yearMin2 = yearCol2.reduce(reducer);  // reducer
  
  var imgList = yearCol.aggregate_array('LANDSAT_ID');
  var n_img = imgList.size();
  var nBands = yearMin2.bandNames().size();
  
  return yearMin2.set({
    'year': year,
    'image_list': imgList,
    'n_bands': nBands,
    'n_img': n_img, 
    'satellite': 'Landsat8',
    'method': '2nd min',
    'system:time_start': year
  });
});
// create collection: 
var secondMinl8 = ee.ImageCollection.fromImages(year_l8_2ndMin).filter(ee.Filter.gt('n_bands', 0))




// ----- 3 SIGMA empirical rule ------------
// 1. compute sigma 2. update mask
// rule: 95 of data fall in (MEAN-2sigma) range
// to perform .multiply(): make sure there are no 0-band images 


var year_l8_sigma = yearListl8.map(function(year){
  
  var yearCol = l8collection.filter(ee.Filter.calendarRange(year, year, 'year'))
  var yearMean = yearCol.reduce(ee.Reducer.mean());
  var yearSD = yearCol.reduce(ee.Reducer.stdDev());
  var year2SD = yearSD.multiply(ee.Image(2));
  var year_2sigma = yearMean.subtract(year2SD);
  var year_2SD_plus = yearMean.add(year2SD);
  var year_2SD_minus = yearMean.subtract(year2SD);
  
 // update mask ->[(mean-2sigma) < pixel < (mean+2sigma)]
 
  var yearCol2sigma = yearCol.map(function (img) {
    var img_masked = img.updateMask(img.gte(year_2SD_minus).and(img.lte(year_2SD_plus)))
    return ee.Image(img_masked)});
  
  var year2sigmal8 = yearCol2sigma.reduce(reducer);  // reducer

  var imgList = yearCol.aggregate_array('LANDSAT_ID');
  var n_img = imgList.size();
  var nBands = year2sigmal8.bandNames().size();
  
  return year2sigmal8.set({
    'year': year,
    'image_list': imgList,
    'n_bands': nBands,
    'n_img': n_img,
    'satellite': 'Landsat8',
    'method': 'sigma',
    'system:time_start': year
  });
});


// convert the annual composite image list to an ImageCollection and filter out years with no images.
var year2sigmal8 = ee.ImageCollection.fromImages(year_l8_sigma).filter(ee.Filter.gt('n_bands', 0))

print('year2sigma', year2sigmal8)










// --------------------- MERGE YEARLY COMPOSITE -----------------

// 1st Min:
var firstMin_merged = firstMinl5.merge(firstMinl8).sort('year')
// 2nd Min:
var secondMin_merged = secondMinl5.merge(secondMinl8).sort('year')
// 2 sigma:
var sigma_merged = year2sigmal5.merge(year2sigmal8).sort('year')




// ------------------------------ YEAR -------------------------

//
// yearly collection to list:
var list = firstMin_merged.toList(firstMin_merged.size())
var list2nd = secondMin_merged.toList(secondMin_merged.size());
var list_sigma = sigma_merged.toList(sigma_merged.size());
print(list2nd)


// get image (1/year)

var visParams = {bands: ['red_min', 'green_min', 'blue_min'], min: 0, max: 3000, gamma: 1.4};
//var visParams = {bands: ['red_min', 'green_min', 'blue_min'], min: 0, max: 3000, gamma: 1.4};


// var year = ee.Image(list.get(15));
//var year = ee.Image(list2nd.get(18)); // 27 - bad one
 var year = ee.Image(list_sigma.get(21));
print(year)
print('image year:', year.get('year'))
print('image method:', year.get('method'))
print('years in total', list_sigma.size())

Map.addLayer(year, visParams, 'image (year)');


var ind = year.select('NDVI_min', 'NDVI_max', 'NDWI_min', 'NDWI_max', 'NDSI_min', 'NDSI_max', 
                      'ind_product_max', 'ind_product_min',
                      'water_ratio_max', 'water_ratio_min', 'snow_ratio_max', 'snow_ratio_min', 'veg_ratio_max', 'veg_ratio_min',
                      'ratio_product_max', 'ratio_product_min',
                      'blue_min', 'blue_max', 'green_min', 'red_min', 'nir_min', 'swir_min')


Map.addLayer(ind, {bands: ['NDVI_max', 'NDWI_max', 'NDSI_max'], min:[-1, -1, -1], max:[1, 1, 1]}, 'ind (year)', false);


// -------- FINAL WATER MASK Caucasus -------------
var water_final = ind1.expression('b("veg_ratio_min") < 0.41 && b("water_ratio_max") > 2.1') // caucasus thresholds

var water_final_mask = water_final.eq(1).selfMask().rename('waterMask')
print('water mask image', water_final_mask)
Map.addLayer(water_final_mask, {palette:['orange']}, '!! water_final_mask caucasus')


// ------------ CAREFULL!! -------
// ------------ THIS IS REGION DEPENDENT!!! ----------------
// CAUCASUS: 

// add water band to the collection 
// function to add water as a band
var addWater = function(img){
  var water = img.expression('b("veg_ratio_min") < 0.41 && b("water_ratio_max") > 2.1').rename('water') // caucasus thresholds
  var mask = water.updateMask(water.eq(1))
  return img.addBands(mask)
  //return img.addBands(water)
}


// --------- cut out everything outside the 1-km buffer: 
// convert buffer to raster 
var glaciers1kmRast = ee.Image().float().paint(glaciers1km, 'area')
                      .multiply(0).add(1)/*.unmask(0)*/.int(); // make binary
//Map.addLayer(glaciers1kmRast, {}, 'glaciers1kmRast')

var buff = function(img){
  return img.updateMask(glaciers1kmRast.eq(1))}

// Elbrus Blob - too bright pixels in the 8bit images 
var elbrus_blob_rast = ee.Image().float().paint(elbrus_blob)
                      .multiply(0).add(1).unmask(0); // make binary

var blob = function(img){
  return img.updateMask(elbrus_blob_rast.neq(1))}




// yearly collections:

// chose method:

//var col = firstMin_merged.map(addWater).map(buff).map(blob).select('water')
//var col = secondMin_merged.map(addWater).map(buff).map(blob).select('water')
var col = sigma_merged.map(addWater).map(buff).map(blob).select('water')
print('new collection with water', col);


// ------------ CAREFULL!! -------
// --------- THIS IS FOR CAUCASUS 
// hard-coded 5-year periods -> 1 image per 5 years 
// caucasus:

// collection to list
var col_list = col.toList(col.size())

// before 1989
var col89 = col.filter(ee.Filter.lte('year', 1989)).mean().float().set({'period': 1989})
// 1990-94
var col94 = col.filter(ee.Filter.eq('year', 1990))
     .merge(col.filter(ee.Filter.eq('year', 1991)))
     .merge(col.filter(ee.Filter.eq('year', 1992)))
     .merge(col.filter(ee.Filter.eq('year', 1993)))
     .merge(col.filter(ee.Filter.eq('year', 1994))).mean().float().set({'period': 1994})
//print(col94)

// 1995-99
var col99 = col.filter(ee.Filter.eq('year', 1995))
     .merge(col.filter(ee.Filter.eq('year', 1996)))
     .merge(col.filter(ee.Filter.eq('year', 1997)))
     .merge(col.filter(ee.Filter.eq('year', 1998)))
     .merge(col.filter(ee.Filter.eq('year', 1999))).mean().float().set({'period': 1999})

// 2000-04
var col04 = col.filter(ee.Filter.eq('year', 2000))
     .merge(col.filter(ee.Filter.eq('year', 2001)))
     .merge(col.filter(ee.Filter.eq('year', 2002)))
     .merge(col.filter(ee.Filter.eq('year', 2003)))
     .merge(col.filter(ee.Filter.eq('year', 2004))).mean().float().set({'period': 2004})
// 2005-09
var col09 = col.filter(ee.Filter.eq('year', 2005))
     .merge(col.filter(ee.Filter.eq('year', 2006)))
     .merge(col.filter(ee.Filter.eq('year', 2007)))
     .merge(col.filter(ee.Filter.eq('year', 2008)))
     .merge(col.filter(ee.Filter.eq('year', 2009))).mean().float().set({'period': 2009})
// 2010-14
var col14 = col.filter(ee.Filter.eq('year', 2010))
     .merge(col.filter(ee.Filter.eq('year', 2011)))
     .merge(col.filter(ee.Filter.eq('year', 2012)))
     .merge(col.filter(ee.Filter.eq('year', 2013)))
     .merge(col.filter(ee.Filter.eq('year', 2014))).mean().float().set({'period': 2014})
// 2015-19
var col19 = col.filter(ee.Filter.eq('year', 2015))
     .merge(col.filter(ee.Filter.eq('year', 2016)))
     .merge(col.filter(ee.Filter.eq('year', 2017)))
     .merge(col.filter(ee.Filter.eq('year', 2018)))
     .merge(col.filter(ee.Filter.eq('year', 2019))).mean().float().set({'period': 2019})
// 2020-on
var col20 = col.filter(ee.Filter.eq('year', 2020)).mean().float().set({'period': 2020})

// new collection with 5-year averages

var meanCol = ee.ImageCollection([col89, col94, col99, col04, col09, col14, col19, col20]).sort('period');
//print(meanCol)




Map.addLayer(col89, {palette:['red']},    'col89', false)
Map.addLayer(col94, {palette:['orange']}, 'col94', false)
Map.addLayer(col99, {palette:['yellow']}, 'col99', false)
Map.addLayer(col04, {palette:['green']},  'col04', false)
Map.addLayer(col09, {palette:['ADD8E6']}, 'col09', false)
Map.addLayer(col14, {palette:['blue']},   'col14', false)
Map.addLayer(col19, {palette:['purple']}, 'col19', false)
