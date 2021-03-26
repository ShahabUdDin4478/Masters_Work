var sentinel = ee.ImageCollection("COPERNICUS/S2_SR"),
    roi = 
    /* color: #98ff00 */
    /* shown: false */
    ee.Geometry.Point([72.43224767218874, 34.49236657665233]),
    sent = ee.ImageCollection("COPERNICUS/S1_GRD"),
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
var collection = sentImage.median().clip(region);

//Select Sentinel - 2 Bands to use in the study
var sentbands = ['B2','B3','B4','B8','B11','B12'];


function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = (1 << 3);
  var cloudsBitMask = (1 << 5);
  // Get the pixel QA band.
  var qa = image.select('pixel_qa');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0));
  return image.updateMask(mask); 
}

// Choose whether you would like to use the first image from the collection or create a median value image from the collection
var chosenimage = ee.Image(collection)//.first()); //Use ".first();" or ".median();"

// ****Input the known resolution of you dataset (in meters)***
var resolution = ee.Number(30);


/*--------------------------------------------------------*/
// Section 1.2 Spectral Subset Option
/*
On which band range(s) would you like to perform the transformation?
Define them, pairwise (inclusive), in the array below.
Ideally, all bands from an image would be used. However, if you're performing the analysis on imagery with many bands, 
(e.g. hyperspectral data) consider subsetting the data to lessen the computational intensity of the algorithm.
E.g.    [[1,10]] would select bands 1-10;
		
        [[2,5],
        [8,12],
        [16,20]] would select bands 2-5,8-12, and 16-20;
        
        [[1,1],
        [2,2],
        [3,3]] would select bands 1, 2, and 3;

*/

var brarray =  [[1,7]];

/*--------------------------------------------------------*/
// 1.3. Band Selection for deriving the Shift Difference area
/*
Below is a separate band subset selection array used to subset bands for use in the shift difference area derivation.
Ideally, the shift difference would be performed using all of the selected bands from above, but given the computational intensity
of the calculation within GEE the script times-out if it is run on more than a few bands. The current suggestion is to select 
the true color bands. (Most analysts performing an MNF transformation would select a shift difference area by visually inspecting 
the image with these bands; as such, the true color bands for LandSat 8 appear below as the default values.)

If you would like to make changes, define them in a pairwise fashion (inclusive) in the array below.
*/
  
var noisebrarray = [[2,2],
                    [3,3],
                    [4,4]];




/*-------------------------------------------------------------------*/
// Section 2. Formatting Reserves
/*-------------------------------------------------------------------*/
// Preferred color schemes can be designed here for later use

var paletteRMS = ['3366FF','000033','ffffff','3366CC','000000'];

var palette_blue2purp = ['ff00ff','006eff','00ffd0','459619'];

///////////////////////NO USER INPUT AFTER THIS POINT///////////////////////

/*---------------------------------------------------*/
// Section 3: Image Preparation
/*---------------------------------------------------*/

// This function performs the band range selection using the arrays defined in the section above.
var bandrangecat = function(inputarray){
  
  var openlist = [];
  
  var numofranges = inputarray.length-1;
  
  for (var i = 0; i <= numofranges; i++){
    
    var start = inputarray[i][0];
    var end = inputarray[i][1];
    
    for (var j = start; j <= end; j++){
      openlist.push(j); 
    }
}
return openlist};

// Use the helper function from above to form a concatenated list of image bands.
var catbands = bandrangecat(brarray);

// Format the concatenated list of image bands as an array.
var catbandsarray = ee.Array(catbands).subtract(1).toList();

// Use the helper function from above to form a concatenated list of noise bands.
var catnoisebands = bandrangecat(noisebrarray);

// Format the concatenated list of noise bands as an array.
var catnoisebandsarray = ee.Array(catnoisebands).subtract(1).toList();


/*--------------------------------------------------------*/
// Format input data


// Spectrally subset input data
var originalImage = chosenimage.select(catbandsarray).clip(region);

// Get band names of bands in use
var bands = originalImage.bandNames();

// Compute the total number of image bands and store this number
var numberofbands = bands.length();

// Subset the image for the noise shift difference calculations
var SubsetNoiseImage = chosenimage.select(catnoisebandsarray).clip(region);

// Retrieve a list of noise band names
var noisebands = SubsetNoiseImage.bandNames();

// Compute the total number of noise bands and store this number
var numberofnoisebands = noisebands.length();




/*----------------------------------------------------*/
// Section 4: Homogenous Region Calculation and Selection
/*----------------------------------------------------*/


// Format a kernel to use when finding the homogenous area
var square_kernel = ee.Kernel.square(resolution.multiply(3), "meters");

// Find standard deviation for each pixel and its determined neighborhood
var stdDev = SubsetNoiseImage.reduceNeighborhood(ee.Reducer.stdDev(), square_kernel);
// Map.addLayer(stdDev, {min:0, max:0.1}, 'Std Dev');
// Map.addLayer(stdDev.select([1]), {min:0, max:1000,palette:['ff00ff','006eff', '00ffd0', '459619']}, 'Std Dev');

// Compute the quadratic mean (root mean square / RMS) of the neighborhood standard deviation through all bands for each pixel
// and then sum these values for each pixel
var RMS = stdDev.multiply(stdDev).reduce(ee.Reducer.sum()).divide(SubsetNoiseImage.bandNames().length()).sqrt();
// Map.addLayer(RMS, {min:0, max:1}, 'RMS');

// Find and store the minimum and maximum variance values (and their range) within the area of interest
var RMSDict = RMS
  .reduceRegion({
  reducer: ee.Reducer.minMax(),
  geometry: region, 
  scale: resolution
});
var dictrmsmax = ee.Number(RMSDict.get('sum_max'));
var dictrmsmin = ee.Number(RMSDict.get('sum_min'));
var dictrange = dictrmsmax.subtract(dictrmsmin);


/*--------------------------------------------------------*/
// Find the area with the lowest variance

// Display the quadratic mean layer using the computed minimum and maximum
var RMS_vis = RMS.visualize({min:dictrmsmin, max:dictrmsmax, palette:palette_blue2purp});
// Map.addLayer(RMS_vis,{},'RMS Full Range', false);

// Define the threshold of how much variance is acceptable for consideration for a shift difference area
var desiredpercentage = 0.05;
var percentofrange = dictrange.multiply(desiredpercentage);
var bottompercent = dictrmsmin.add(percentofrange);
var threshold = RMS.select(['sum']).lt(bottompercent);

// Display pixels that have variance below the defined threshold
var lowRMS = threshold.mask(threshold);
// Map.addLayer(lowRMS, {min:0, max:1, palette:palette_blue2purp}, 'Pixels with low RMS', false);


/*--------------------------------------------------------*/
// Make each area an island

// Define the kernel used in the connectedComponents() call
var kernel_fixed = ee.Kernel.fixed(3, 3,[[1,1,1],
                                        [1,0,1],
                                        [1,1,1]]);

// Connect all pixels that meet the threshold defined above
var connected = lowRMS.connectedComponents(kernel_fixed,256);
// Map.addLayer(connected,{}, 'Connected Low RMS Islands');

// Determine the minimum number of pixels in the areas of interest
var pixelmin = ee.Algorithms.If(numberofbands.lte(50),ee.Number(50),numberofbands);

// Compute the number of connected pixels in each island
var count = connected.connectedPixelCount(ee.Array(pixelmin).multiply(4), true);
// Map.addLayer(count, null, 'Connected Count');

// Reproject the layer and disregard all areas below the minimum pixel threshold
var precount = count.reproject(
  chosenimage.projection().crs(),
  null,
  resolution).gt(ee.Number(pixelmin));
// Map.addLayer(precount, null, 'Connected Count Reprojected');

// Run a mask to leave only the islands of interest then restore the bands
var countFV = precount.mask(precount);
var RTV = countFV.addBands(countFV);

// Turn each island into a vector feature and filter out islands that do not
// meet the minimum area requirement
var islandsprep = RTV.reduceToVectors({
  reducer: ee.Reducer.first(),
  geometry: region, 
  geometryType: 'polygon',
  scale: resolution, 
  bestEffort: true 
  });
// Map.addLayer(islandsprep, {color:'e51bc7'}, 'Islands', false);

var islandsWithArea = islandsprep.map(function(f) {
return f.set('area', f.geometry().area(5))});
// Map.addLayer(withArea);

var islands = islandsWithArea.filterMetadata('area', 'greater_than', resolution.multiply(resolution).multiply(numberofbands));
// Map.addLayer(islands, {color: '900000'}, 'Big Islands');


/*--------------------------------------------------------*/
// Sum each region's variance then sort the regions

//Find the total variance within each box and choose the box with the least variance
var buffvariance = RMS.reduceRegions(islands, ee.Reducer.sum(), resolution);
// Map.addLayer(buffvariance,{color:'00ff00'},'Buffered Sums of RMS St.d.');

// Select the buffered area with the lowest variance to get shift difference area
var buffvarImage = buffvariance.limit(1, 'sum', true);
// Map.addLayer(buffvarImage, {color:'ff0000'}, 'Final Buffered Area');




/*-------------------------------------------------------------------*/
// Section 5: Perform the Shift Difference Calculation
/*-------------------------------------------------------------------*/

// Clip the image to the chosen shift difference area
var kernelarea = originalImage.clip(region);
// Map.addLayer(kernelarea,{},'Noise Kernel Area',false);


// Define kernels that link a pixel to its neighbors immediately above and to the left of it
var kernel_left = ee.Kernel.fixed(3, 3, [[0,0,0],
                                        [1,0,0],
                                        [0,0,0]]);

var kernel_up   = ee.Kernel.fixed(3, 3, [[0,1,0],
                                        [0,0,0],
                                        [0,0,0]]);

// Create a layer stack of niehgboring pixel values in order to perform math between pixel values
var kernelimage_left = kernelarea.neighborhoodToBands(kernel_left);
var kernelimage_up = kernelarea.neighborhoodToBands(kernel_up);

var diff_left = kernelimage_left.subtract(kernelarea);
var diff_up = kernelimage_up.subtract(kernelarea);

// Find average difference between pixels for the whole shift difference area
var diff = diff_left.add(diff_up).divide(2).clip(kernelarea.geometry());
// Map.addLayer(diff, {min:-100, max:100}, 'diff');




/*-------------------------------------------------------------------*/
// Section 6: Finalize the MNF Transformation
/*-------------------------------------------------------------------*/

// Find the covariance matrix of the finalized shift difference area
var covardict = diff.toArray().reduceRegion(ee.Reducer.covariance(), null, resolution, null, null, false, 800000000);

// Convert the covariance matrix into an array
var noisecovariancematrix = ee.Array(covardict.get('array'));
// print("Noise Coavariance Matrix", noisecovariancematrix);

// Decompose the matrix into eigenvalues and eigenvectors
var eigendecomplist = noisecovariancematrix.eigen();


/*----------------------------------------*/
// MNF Matrix

// Use the results of the decomposition to formulate the required matrices for the subsequent mathematics
var eigenvalues = eigendecomplist.slice(1, 0, 1);
var eigenvectors = eigendecomplist.slice(1, 1);

var matrixr = eigenvalues.sqrt().pow(-1).matrixToDiag();
var matrixcmnf = eigenvalues.pow(-1).matrixToDiag();


/*--------------------------------------------------------*/
// Noise-whiten the dataset

// Convert the image to an array
var arrayimage = originalImage.toArray();
// Map.addLayer(arrayimage);

// Find the mean value of the bands in the whole image
var meanimage = originalImage.reduceRegion('mean', region, resolution, null, null, false, 800000000);

// Make an array from the image’s band means for each band
var meanarray = ee.Array(meanimage.values(bands));

// Mean correct the image
var meancenteredimage = arrayimage.subtract(meanarray);
// Map.addLayer(meancenteredimage, {}, "Mean Centered Image");

// Multiply the mean centered image by the noise eigenvectors then scale the data by the noise standard deviation values
var nwarrayimage = meancenteredimage.arrayRepeat(1,1).arrayTranspose()
                    .matrixMultiply(eigenvectors.transpose())
                    .matrixMultiply(matrixr);
// Map.addLayer(nwarrayimage, {min:0, max:10}, "Noise Whitened Array Image");

// Covariate the noise-whitened dataset
var nwimage = nwarrayimage.arrayFlatten([['Noise-Whitened'],bands]);
// Map.addLayer(nwimage);
var nwcovardict = nwarrayimage.arrayProject([1]).reduceRegion(ee.Reducer.covariance(), region, resolution, null, null, false, 800000000);

// Eigendecompose the covariance matrix
var nwcovariancematrix = ee.Array(nwcovardict.get('array'));
var nweigendecomplist = nwcovariancematrix.eigen();

// Retrieve the eigenvalues and eigenvectors for each MNF transformed band
var nweigenvectors = nweigendecomplist.slice(1, 1);
var nweigenvalues = nweigendecomplist.slice(1, 0, 1);

// Finalize the MNF Transformation by multiplying the second eigenvector matrix by the noise-whitened data
var mnfdata = nwarrayimage.matrixMultiply(nweigenvectors.transpose());
// Map.addLayer(mnfdata);

// Use a map function to retrieve the band names for the image
var bl = ee.List.sequence(1, numberofbands, 1);
var fbands = bl.map(function(n) {
  return ee.String('Band ').cat(ee.Number(n).int());
});

// Flatten the array image back into a normal image and add results to map
var mnfimage = mnfdata.arrayFlatten([['MNF Transformed'],fbands]," ");
Map.addLayer(mnfimage, {min:-3, max:8}, "MNF Transformed Data");

// **End of the MNF Transformation**

var mnfcovardict = mnfdata.arrayProject([1]).reduceRegion(ee.Reducer.covariance(), null, resolution, null, null, false, 3000000000);
var mnfcovariance = ee.Array(mnfcovardict.get('array'));
// print(mnfcovariance);
var mnfeigendecomp = mnfcovariance.eigen();
var mnfeigenvalues = mnfeigendecomp.slice(1, 0, 1);




/*-------------------------------------------------------------------*/
// Section 7: Display MNF bands and chart the eigenvalues
/*-------------------------------------------------------------------*/
 
var eigenValueArray = ee.Array(mnfeigenvalues).repeat(0,1);

var charty = Chart.array.values(eigenValueArray, 0, bl).setSeriesNames(["Eigen Values"]);
charty = charty.setOptions({
  title: 'Eigenvalues For MNF Bands',
  hAxis: {
    title: 'MNF Bands'
  },
  vAxis: {
    title: 'Eigenvalue'
  },
  lineWidth: 1,
  pointSize: 4,
  series: {
    0: {color: 'darkgreen'}
  }
});
print(charty);
