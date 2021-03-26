
var sentinel = ee.ImageCollection("COPERNICUS/S2_SR"),
    roi = 
    /* color: #98ff00 */
    /* shown: false */
    ee.Geometry.Point([72.43224767218874, 34.49236657665233]),
    sent = ee.ImageCollection("COPERNICUS/S1_GRD"),
    geometrys = 
    /* color: #98ff00 */
    /* shown: false */
    ee.Geometry.MultiPoint(
        [[72.16675473464745, 34.621388408588814],
         [72.17842770828027, 34.28109518135768],
         [72.7352972151162, 34.302085039608315],
         [72.72568417800683, 34.638903336034836]]),
    geometry = 
    /* color: #d63000 */
    /* shown: false */
    ee.Geometry.Point([72.42280913131358, 34.45189962454313]),
    country = 
    /* color: #98ff00 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[72.18899799162156, 34.639944827070565],
          [72.18899799162156, 34.302563524991655],
          [72.72664142423875, 34.302563524991655],
          [72.72664142423875, 34.639944827070565]]], null, false),
    region = 
    /* color: #d63000 */
    /* shown: false */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[72.26895252875929, 34.52851500996908],
          [72.26895252875929, 34.46669057005418],
          [72.40096012763624, 34.46669057005418],
          [72.40096012763624, 34.52851500996908]]], null, false);


var sentImage = ee.ImageCollection('COPERNICUS/S2_SR')
.filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 3))
.filter(ee.Filter.lt('NOT_VEGETATED_PERCENTAGE', 60))
.filterDate("2015-01-02", "2020-12-30")
.filterBounds(region);
//print(sentImage)
var sentImage = sentImage.median().clip(region);

//Select Sentinel - 2 Bands to use in the study
var sentbands = ['B2','B3','B4','B8','B11','B12'];
  //Sentinel - 2 L2A TCC
var trueColor = {
  bands: ["B4", "B3", "B2"],
  min: 0,
  max: 3000
};
Map.addLayer(sentImage, trueColor, "Sentinel 2 True Color");

var NDVI = sentImage.expression(
  "(NIR-RED) / (NIR+RED)",
  {
    RED: sentImage.select("B4"),
    NIR: sentImage.select("B8"),
    BLUE: sentImage.select("B2")
  })
  
  Map.addLayer(NDVI, {min: 0, max: 1}, "NDVI");
Export.image.toDrive({
  image: NDVI,
  description: 'NDVI',
  folder: "PlanB",
  scale: 15,
  region: region
});

// Nice visualization parameters for a vegetation index.
var vis = {min: 0, max: 1, palette: [
  'FFFFFF', 'CE7E45', 'FCD163', '66A000', '207401',
  '056201', '004C00', '023B01', '012E01', '011301']};
Map.addLayer(NDVI, vis, 'NDVI');
Map.addLayer(NDVI.gt(0.5), vis, 'NDVI Binarized');
Export.image.toDrive({
  image: NDVI.gt(0.5),
  description: 'NDVIgt',
  folder: "PlanB",
  scale: 15,
  region: region
});

//Ignore this, this is only used to avoid errors while concatinting eigen pairs
var eigenCollection = ee.Array([[0],[1],[2],[3],[4]]);
var tscale = 5;

// Decorrelation Stretching Main Function     
function decStr(bandsImage, location, scale){
  var bandNames = bandsImage.bandNames();
  // Naming the axis for intuition
  var dataAxis = 0;
  var bandsAxis = 1;
  // Calculate the mean for each band image
  var meansAll = bandsImage.reduceRegion(ee.Reducer.mean(), location, scale);
  // Generate an array (1D Matrix) of mean of each band
  var arrayOfmeans = ee.Image(meansAll.toArray());
  // Collapse the bands data such that each pixel is a matrix of pixel values of each band
  var pixelArrays = bandsImage.toArray();
  // Use the means array and the collapsed band data to center each pixel of each band by subtracting its corresponding mean from it
  var meanCent = pixelArrays.subtract(arrayOfmeans);
  // Calculate the Covariance matrix for the bands data
  var covar = meanCent.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: location,
    scale: scale
  });
  
  // Get the covariance in array format which shows the band-band covarince of the data
  var covarArray = ee.Array(covar.get('array'));
  // Perform eigen decomposition of the covariance matrix to obtain eigen values and eigen vector pairs
  var eigenPairs = covarArray.eigen();
  var eigenValues = eigenPairs.slice(bandsAxis, 0, 1); // slice(axis, start, end, step)
  var eigenVectors = eigenPairs.slice(bandsAxis, 1);
  // Rotate by the eigenvectors, scale to a variance of 30, and rotate back.
  //Store a diagonal matrix in i
  var i = ee.Array.identity(bandNames.length()); // i will be used to isolate each band data and scale its variance e.g i = [1,0,0,0,0] = isolate first band from 5 bands
  // Calculate variance from the eigenvalues ---> variance = 1/sqrt(eigenvalues)
  // matrixToDiag = Computes a square diagonal matrix from a single column matrix for multiplication purposes
  var variance = eigenValues.sqrt().matrixToDiag();
  //Multiply diagonal matrix i by 30 and divide by vaiance to obtain scaling variance matrix
  var scaled = i.multiply(30).divide(variance); //Changed from 30 -> 50, It was observed that changing variance scale increases contrast. Best contrast obtained for 30
  // Calculate a rotation matrix ---> rotationMatrix =  Eigenvect.Transpose * ScaledVariance * Eigenvect
  var rotation = eigenVectors.transpose()
    .matrixMultiply(scaled)
    .matrixMultiply(eigenVectors);
  // Convert 1-D nomalized array image data to 2-D and transpose it so it can be multiplied with rotation matrix
  var transposed = meanCent.arrayRepeat(bandsAxis, 1).arrayTranspose();
  // Multiply the transposed data with the rotation matrix
  return transposed.matrixMultiply(ee.Image(rotation))
    .arrayProject([bandsAxis])   //This drop unecessary axis from the transposed data and only retains 2 axis
    .arrayFlatten([bandNames])  //Flatten collections of collections
    .add(127).byte(); // Conver pixel values to 127 means so it can be visualized between 0 - 255 range.
    
    // .byte is used to force element wise operation
}

// Principal Component Analysis Main Function
function PCA(meanCent, scale, location){
  // Flatten the band image data in from 2D to a 1D array
  var arrays = meanCent.toArray();
  //print('PCA applying on', meanCent);
  // Calculate the covariance matrix for the bands data of the region
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: location,
    scale: scale,
    //tileScale: tscale,
    maxPixels: 1e9
  });
  // Get the band to band covariance of the region in 'array' format. Here .get('array') --> casts to an array
  var covarArray = ee.Array(covar.get('array'));
  // Perform an eigen analysis and slice apart the values and vectors.
  var eigenPairs = covarArray.eigen();
  // This is a P-length vector of Eigenvalues. Here P = number of PCs
  var eigenValues = eigenPairs.slice(1, 0, 1);
  // This is a PxP matrix with eigenvectors in rows.
  var eigenVectors = eigenPairs.slice(1, 1);
  //Print and store eigen pairs in eigenCollection variable and export to drive
 // print('eigen Values', eigenValues);
 // print('eigen Vector', eigenVectors);
    //Make feature collection out of eigenpairs so it can be exported to excel. From there we Convert it to a table using a python script
  eigenCollection = ee.Feature(null,{values:ee.Array.cat([eigenValues,eigenVectors],1)}); 
 // print('Eigen Collection Length',eigenCollection);
    // Export the FeatureCollection to excel sheet in drive
/*
  Export.table.toDrive({
  collection: ee.FeatureCollection([eigenCollection]),
  description: 'eigenAnalysis',
  fileFormat: 'CSV'
  });
*/
  // Convert the 1D image array back to 2D matrix for multiplication
  var imageMat = arrays.toArray(1);
  // To obtain PC = EigenVectors * 2D Image Matrix
  var PCs = ee.Image(eigenVectors).matrixMultiply(imageMat);
  // Turn the square roots of the Eigenvalues into a P-band image.
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);
  // Turn the PCs into a P-band image, normalized by SD.
  return PCs
    // Throw out an an unneeded dimension, [[]] -> [].
    .arrayProject([0])
    // Make the one band array image a multi-band image, [] -> image.
    .arrayFlatten([getNewBandNames('pc')])
    // Normalize the PCs by their SDs.
    .divide(sdImage);
}

  //Sentinel - 2 L2A TCC
var trueColor = {
  bands: ["B4", "B3", "B2"],
  min: 0,
  max: 3000
};
//Map.addLayer(sentImage, trueColor, "Sentinel 2 True Color");

  //Sentinel - 2 L2A FCC
var trueColor = {
  bands: ["B11", "B12", "B8"],
  min: 0,
  max: 3000
};
//Map.addLayer(sentImage, trueColor, "Sentinel 2 False Color");


              //APPLYING DS 
// Selecting bands to apply DS
var sentBandsImage = sentImage.select(sentbands);

//Obtain DS Results for All Satelites using dcs function
var DSsent = decStr(sentBandsImage, region, 1000);

//FCC of 3 bands of DS results for all satelites

var selectBands = [0,1,2,3,4,5]; 
//Map.addLayer(DSsent.select(selectBands), {}, 'DCS Sentinel Image');




            //PRINCIPAL COMPONENT ANALYSIS
   
  //Applying PCA on DS of Sentinel - 2
var region = sentImage.geometry();
var bands = [0,1,2,3,4,5];
var image =  DSsent.select(bands);
var scale = 30;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);
var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  });
};
var pcImageSDS = PCA(centered, scale, region);
/*Export.image.toDrive({
  image: pcImageDS,
  description: 'Sentinel2PCAofDS',
  folder: "GraniteExports",
  scale: 15,
  region: region
});*/

//Map.addLayer(pcImageSDS, {bands: ['pc2', 'pc3', 'pc6'], min: -2, max: 2}, 'Sentinel 2 L2C - PCA of DS '); //changed from PC1, PC2 and PC3

      
  //Sentinel PCA
var region = sentImage.geometry();
var image =  sentImage.select(sentbands);
var scale = 5;
var bandNames = image.bandNames();
var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
});
var means = ee.Image.constant(meanDict.values(bandNames));
var centered = image.subtract(means);


var pcImageS = PCA(centered, scale, region);

// Plot each PC as a new layerl
Map.addLayer(pcImageS, {bands: ['pc1', 'pc2', 'pc4'], min: -2, max: 2}, 'Sentinel 2 - PCA  used in paper');
Export.image.toDrive({
  image: pcImageS,
  description: 'Sentinel2PCA',
  folder: "PlanB",
  scale: 15,
  region: region
});

  // Plot each PC as a new layer
for (var i = 0; i < bandNames.length().getInfo(); i++) {
  var band = pcImageS.bandNames().get(i).getInfo();
  //Map.addLayer(pcImage.select([band]), {min: -2, max: 2}, band);
}



      //Sentinel - 2 Clustering
// Make the training dataset.
var training = sentImage.sample({
  region: region,
  scale: 20,
  numPixels: 8000
});

// Instantiate the clusterer and train it.
var clusterer = ee.Clusterer.wekaXMeans(2,15).train(training);
print('WekaXMeans Clusterer:', clusterer);
// Cluster the input using the trained clusterer.
var result = sentImage.cluster(clusterer);
print('Sentinel WekaXMeans', result)
// Display the clusters with random colors.
Map.addLayer(result.randomVisualizer(), {}, 'Sentinel WekaXMeans Results');

Export.image.toDrive({
  image: result,
  description: 'Sentinel2Xmeans',
  folder: "PlanB",
  scale: 15,
  region: region
});
Map.setCenter(72.30109838019655,34.528578508206145, 12);

// Load the data and select the bands of interest.
var image = ee.Image("UMD/hansen/global_forest_change_2019_v1_7");
var image = image.clip(region)
var lossImage = image.select(['loss']);
var gainImage = image.select(['gain']);

var areaImage = lossImage.multiply(ee.Image.pixelArea());
// Use the and() method to create the lossAndGain image.
var gainAndLoss = gainImage.and(lossImage);

// Show the loss and gain image.
//Map.addLayer(gainAndLoss.updateMask(gainAndLoss),{palette: 'FF00FF'}, 'Gain and Loss');
    
var treeCover = image.select(['treecover2000']);
Export.image.toDrive({
  image: treeCover,
  description: 'treeCover',
  folder: "PlanB",
  scale: 15,
  region: region
});
Map.addLayer(treeCover.updateMask(treeCover),{palette: ['green'], max: 100}, 'Forest Cover');
Export.image.toDrive({
  image: treeCover,
  description: 'treeCover',
  folder: "PlanB",
  scale: 15,
  region: region
});

// Add the loss layer in red.
Map.addLayer(lossImage.updateMask(lossImage),{palette: ['red']}, 'Loss');
Export.image.toDrive({
  image: lossImage,
  description: 'Loss',
  folder: "PlanB",
  scale: 15,
  region: region
});
// Add the gain layer in blue.
Map.addLayer(gainImage.updateMask(gainImage),{palette: ['blue']}, 'Gain');
Export.image.toDrive({
  image: gainImage,
  description: 'Gain',
  folder: "PlanB",
  scale: 15,
  region: region
});
// Sum the values of forest loss pixels in the region.
var stats_loss = lossImage.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: region,
  scale: 4
});
print('pixels representing loss: ', stats_loss.get('loss'), 'square meters');

var stats_gain = gainImage.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: region,
  scale: 4
});
print('pixels representing gain: ', stats_gain.get('gain'), 'square meters');

var treeLossVisParam = {bands: ['lossyear'],min: 0,  max: 5,palette: ['red','blue', 'black','yellow']};
Map.addLayer(image, treeLossVisParam, 'tree loss year');
Export.image.toDrive({
  image: treeLossVisParam,
  description: 'treelossvis',
  folder: "PlanB",
  scale: 15,
  region: region
});

/*

The code is for change detection using sentinel 1 
*/





function get_images(polarisations) {
  
    var img = ee.ImageCollection('COPERNICUS/S1_GRD')
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarisations[0]))
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', polarisations[1]))
            .filter(ee.Filter.eq('instrumentMode', 'IW'));
    return img;
}

var polarisation = 'VV and VH' ;//'VV'
// vh show it better though
var polarisations = ['VV', 'VH'];
var sentinel1_collection = get_images(polarisations);

// Different look angles
var ascending = sentinel1_collection.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
var descending = sentinel1_collection.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));

// Create a composite from means at different polarizations and look angles.

var before = ee.Filter.date('2018-04-01', '2019-05-30');
var after = ee.Filter.date('2019-08-01', '2020-09-15');

var composite_VV_VH_before = ee.Image.cat([
  ascending.filter(before).select('VH').mean(),
  ee.ImageCollection(ascending.filter(before).select('VV').merge(descending.filter(before).select('VV'))).mean(),
  descending.filter(before).select('VH').mean()
]).focal_median().clip(region);
Export.image.toDrive({
  image: composite_VV_VH_before,
  description: 'composite_VV_VH_before',
  folder: "PlanB",
  scale: 15,
  region: region
});
var composite_VV_VH_after = ee.Image.cat([
  ascending.filter(after).select('VH').mean(),
  ee.ImageCollection(ascending.filter(after).select('VV').merge(descending.filter(after).select('VV'))).mean(),
  descending.filter(after).select('VH').mean()
]).focal_median().clip(region);

Export.image.toDrive({
  image: composite_VV_VH_after,
  description: 'composite_VV_VH_after',
  folder: "PlanB",
  scale: 15,
  region: region
});
// Display as a composite of polarization and backscattering characteristics.
Map.addLayer(composite_VV_VH_before, {min: [-25, -20, -25], max: [0, 10, 0]}, 'Composite Before');
Map.addLayer(composite_VV_VH_after, {min: [-25, -20, -25], max: [0, 10, 0]}, 'Composite After');

Map.centerObject(region, 13);

// Make a composite
Map.addLayer(composite_VV_VH_before.addBands(composite_VV_VH_after), {min: -25, max: 0}, 'VV HV stack');


// Given that Sentinal-1 SAR data is already pre-processed, I would do a simple subtraction
var change = composite_VV_VH_after.subtract(composite_VV_VH_before).clip(region);
// Change
var params = {min: -2, max: 2};
Map.addLayer(change, params, 'Change');
// Here the whites - are the big changes
// Very pink fields - unchanged

Export.image.toDrive({
  image: change,
  description: 'change',
  folder: "PlanB",
  scale: 15,
  region: region
});
Map.addLayer(ee.Image().paint(region, 0, 2), {}, 'region');

var testAreaImage = ee.Image(1).clip(region);
var scaleforTestArea = 5;
var reducer = testAreaImage.reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: region,
  crs: 'EPSG:32645', // WGS Zone N 45
  scale: scaleforTestArea,
  maxPixels: 1E13
});
// km square
var area = ee.Number(reducer.get('constant')).multiply(scaleforTestArea).multiply(scaleforTestArea).divide(1000000);
// gives an area  km2
print('area of region using pixel count method: ', area.getInfo() + ' km2');
