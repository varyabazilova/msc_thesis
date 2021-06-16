// https://code.earthengine.google.com/2c81667726cec6646efd15391505ba95



// count landsat images available for the ROI, print histograms

var roi = geometry

var l5_coll = ee.ImageCollection('LANDSAT/LT5_L1T_TOA').filterBounds(roi)  
var l7_coll = ee.ImageCollection('LANDSAT/LE7_L1T_TOA')       
          .filter(ee.Filter.calendarRange(1999, 2002,'year'))
          .filterBounds(roi);  


var l8_coll = ee.ImageCollection('LANDSAT/LC08/C01/T1_TOA').filterBounds(roi)     




// LANDSAT 5 
var iC5 = l5_coll
  
iC5 = iC5.map(function(img){
  var year  = img.date().format("Y");   // get the acquisition year
  var month  = img.date().format("M")
  var CC = img.get('CLOUD_COVER')
  return img.set({
    'year': ee.Number.parse(year),
    'month': ee.Number.parse(month),
    'clouds': ee.Number.parse(CC),
    'satellite' : 'landsat 5'
  });
});


// LANDSAT 7 
var iC7 = l7_coll
  
iC7 = iC7.map(function(img){
  var year  = img.date().format("Y");   // get the acquisition year
  var month  = img.date().format("M")
  var CC = img.get('CLOUD_COVER')
  return img.set({
    'year': ee.Number.parse(year),
    'month': ee.Number.parse(month),
    'clouds': ee.Number.parse(CC),
    'satellite' : 'landsat 7'
  });
});

// LANDSAT 8 
var iC8 = l8_coll//.filter(ee.Filter.calendarRange(5, 9,'month'))
  
iC8 = iC8.map(function(img){
  var year  = img.date().format("Y");   // get the acquisition year
  var month  = img.date().format("M")
  var CC = img.get('CLOUD_COVER')
  return img.set({
    'year': ee.Number.parse(year),
    'month': ee.Number.parse(month),
    'clouds': ee.Number.parse(CC),
    'satellite' : 'landsat 8'
  });
});




// FEATURES

var iC_FC5 = ee.FeatureCollection(iC5);  
var iC_FC7 = ee.FeatureCollection(iC7);  
var iC_FC8 = ee.FeatureCollection(iC8);  



// Make histograms 


var iC_FC_size = iC_FC.size();

var options = {
    title: 'Landsat Mission GEE image availability',
    hAxis: {title: 'Year'},
    vAxis: {title: 'Image count'},
    colors: ['orange'],
  };

  
// Make the histogram, YEAR Landsat 5 
var histogram5 = ui.Chart.feature.histogram({
  features: iC_FC5,
  property: 'year',
  minBucketWidth: 1
}).setOptions(options);

print(histogram5)


// Make the histogram, YEAR Landsat 7
var histogram7 = ui.Chart.feature.histogram({
  features: iC_FC7,
  property: 'year',
  minBucketWidth: 1
}).setOptions(options);

print(histogram7)

// Make the histogram, YEAR Landsat 8
var histogram8 = ui.Chart.feature.histogram({
  features: iC_FC8,
  property: 'year',
  minBucketWidth: 1
}).setOptions(options);

print(histogram8)










