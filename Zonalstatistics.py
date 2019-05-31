import gdal, osr
from skimage.graph import route_through_array
import numpy as np
import ogr, gdal, osr, os
import numpy as np
import itertools
from math import sqrt,ceil
import time

def main(CostSurfacefn, NewCostSurfacefn):

	costSurfaceArray, noDataValue = raster2array(CostSurfacefn) # creates array from cost surface raster
	costSurfaceArray = zonalstatistics(costSurfaceArray, noDataValue)
	array2raster (CostSurfacefn, NewCostSurfacefn, costSurfaceArray)


def raster2array(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray()
    print ('array made')
    print (band.GetNoDataValue())
    return array, band.GetNoDataValue()

def zonalstatistics (array, noDataValue):
    print ('doing statistics')
    print (len(array))
    for i in range(len(array)-1):
        print (i)        
        for j in range (len(array[i])-1):
            if int(array[i][j]) == int(noDataValue):
                if array [i-1][j-1] != int(noDataValue) or array [i-1][j+1] != int(noDataValue) or array [i+1][j-1] != int(noDataValue) or array [i+1][j+1] != int(noDataValue):
                    values = []
                    if array [i-1][j-1] != noDataValue:
                        values.append(array [i-1][j-1])
                    if array [i-1][j] != noDataValue:
                        values.append(array [i-1][j])
                    if array [i-1][j+1] != noDataValue:
                        values.append(array [i-1][j+1])
                    if array [i][j-1] != noDataValue:
                        values.append(array [i][j-1])
                    if array [i][j+1] != noDataValue:
                        values.append(array [i][j+1])
                    if array [i+1][j-1] != noDataValue:
                        values.append(array [i+1][j-1])
                    if array [i+1][j] != noDataValue:
                        values.append(array[i+1][j])
                    if array [i+1][j+1] != noDataValue:
                        values.append(array [i+1][j+1])
                    if len(values) > 4:
                        values.sort()
                        array[i][j] = values[int(len(values)/2)]
            
            if int(array[i][j]) == int(999):

                    #check really quickly if it's not a 'field' of 999
                    if array [i-1][j-1] != int(999) or array [i-1][j+1] != int(999) or array [i+1][j-1] != int(999) or array [i+1][j+1] != int(999):
                        values = []
                        values.append(array [i-1][j-1])
                        values.append(array [i-1][j])
                        values.append(array [i-1][j+1])
                        values.append(array [i][j-1])
                        values.append(array [i][j+1])
                        values.append(array [i+1][j-1])
                        values.append(array [i+1][j])
                        values.append(array [i+1][j+1])
                        if values.count(4) >= 4:
                            array[i][j] = int(4)

    return array

def array2raster(rasterfn,newRasterfn,array):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    cols = raster.RasterXSize
    rows = raster.RasterYSize

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromWkt(raster.GetProjectionRef())
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()


if __name__ == "__main__":
	start_time = time.time()
	CostSurfacefn = '/Users/daveyoldenburg 1/Desktop//Test for LCA/CostRaster_size1.tif'
	NewCostSurfacefn = '/Users/daveyoldenburg 1/Desktop/Test for LCA/CostRaster_size1_filled.tif'
	main(CostSurfacefn, NewCostSurfacefn)
	print("--- %s seconds ---" % (time.time() - start_time))