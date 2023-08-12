# S.-alterniflora-removal-detection
var S2 = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
var L8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
//L8  = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")

var bandIn_S2 = ['B3','B4','B8','B11','cloud']
var bandOut = ['green','red','nir','swir1','cloud']
//bandIn_L8 = ['SR_B3','SR_B4','SR_B5','SR_B6','QA_PIXEL']
var bandIn_L8 = ['B3','B4','B5','B6','pixel_qa']
