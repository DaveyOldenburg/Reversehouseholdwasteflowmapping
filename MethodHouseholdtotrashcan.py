from shapely.geometry import Point, LineString, Polygon
import shapely.wkt
from skimage.graph import route_through_array 
import numpy as np
import ogr, gdal, osr, os
import itertools
from math import sqrt,ceil, atan, tan, degrees, cos, sin
import rasterio
from rasterio import features
import time
import geojson
import geopandas as gpd
 
#ox.config(log_console=False, use_cache=False)

def main(CostSurfacefn,startCoordShp,stopCoordShp,pixelValue, searchdist, WKT_file, closest_trash_dist, newRasterfn, linebuff, Hotspot, rasteradjustment):

    #G = ox.graph_from_place('Willemsparkbuurt, Netherlands', network_type='walk')
    costSurfaceArray, noDataValue = raster2array(CostSurfacefn) # creates array from cost surface raster
    #get coordinates out of shp HOUSEHOLDS
    file_house = ogr.Open(startCoordShp)
    shape_house = file_house.GetLayer(0)
    #get coordinates out of shp TrashCan)
    file_trash = ogr.Open(stopCoordShp)
    shape_trash = file_trash.GetLayer(0)
    #open hotspot file and put into list
    file_hotspot = ogr.Open(Hotspot)
    shape_hotspot = file_hotspot.GetLayer(0)
    
    hotspotlist = []
    for i in range(shape_hotspot.GetFeatureCount()):
        feature_hotspot = shape_hotspot.GetFeature(i)
        feature_json_hotspot = feature_hotspot.ExportToJson()
        feature_json_hotspot = eval(feature_json_hotspot)
        hotcoordX = feature_json_hotspot['geometry']['coordinates'][0]
        hotcoordY = feature_json_hotspot['geometry']['coordinates'][1]
        hotCoord = (hotcoordX,hotcoordY)
        hotspotlist.append(hotCoord)

    #errorfix
    null = None

    
    with open(WKT_file, 'w') as f:
        #write first line of WKT file
        f.write(str('Geometry') + '\t' +'Identifica_house' + '\t' + 'length' + '\t' + 'LOCATIENR'+'\t' + 'geom_pap'+'\t' + 'LOCATIE_pap'+'\t' +'geom_glas'+'\t' + 'LOCATIE_glas'+'\t' +'geom_plas'+'\t' + 'LOCATIE_plas'+'\n')

 
        already_made = {}
        #every household
        #for i in range(shape_house.GetFeatureCount()):
        for i in range(1):
            feature_house = shape_house.GetFeature(i)
            feature_json_house = feature_house.ExportToJson()
            feature_json_house = eval(feature_json_house)
            startcoordX = feature_json_house['geometry']['coordinates'][0]
            startcoordY = feature_json_house['geometry']['coordinates'][1]
            startCoord = (startcoordX,startcoordY)
            #to make sure that the costsurface is reset every time
            costSurfaceArray2 = costSurfaceArray
            """find the closest hotspot area, for now, its this"""
            """with in certain distances make a list of closest hotspots"""
            hotspotdistancedict = {}
            small = 999999
            for hotspot in hotspotlist:
                distance = distanceR2(startcoordX, startcoordY, hotspot[0], hotspot[1])
                hotspotdistancedict[distance] = hotspot
                if distance < small:
                    small = distance

            hotspotxy = hotspotdistancedict[small]
            housexy = (startcoordX, startcoordY)
            #shortestpath_wkt = shortestpathtohotspot (G, housexy, hotspotxy, node_buff, edge_buff)
            direction_wkt = directiontohotspot (housexy, hotspotxy, linebuff)
            direction_geojson = WKT2GEOJSON (str(direction_wkt))
            hotspotinraster = Bresenham_with_rasterio(CostSurfacefn, direction_geojson, rasteradjustment)

            """This is where the 'preference' is added ot the regular ] raster, this is done for every household"""
            costandhotspotarray =  costSurfaceArray2 * hotspotinraster
            #array2raster(CostSurfacefn,newRasterfn,costandhotspotarray)

            considerable_trashcans = [] 
            #every trashcan
            # if point in polygon which is already calculated, use that one.
            if feature_json_house['properties']['id_pand'] in already_made:
                print ('==========duplicate found=========')
                lists = already_made[feature_json_house['properties']['id_pand']]
                f.write(str(lists[1])+ '\t'+ str(lists[2]) +'\t'+ str(lists[0]) +'\t' + str(lists[3])  +'\t' + str(lists[4])+'\t' + str(lists[5])+'\t' + str(lists[6])+'\t' + str(lists[7])+'\t' + str(lists[8])+'\t' + str(lists[9])+'\n')
            

            else:
                
                for i in range(shape_trash.GetFeatureCount()):
                    feature_trash = shape_trash.GetFeature(i)
                    feature_json_trash = feature_trash.ExportToJson()
                    feature_json_trash = eval(feature_json_trash)
                    if 'rest' in feature_json_trash['properties']['FRACTIE1']:
                        #weird reason it's multipoint so all the same point, but more of them
                        stopcoordX = feature_json_trash['geometry']['coordinates'][0][0]
                        stopcoordY = feature_json_trash['geometry']['coordinates'][0][1]

                        #define search radius, only possible if CRS is in meters.
                        if (float(stopcoordX) > startcoordX - searchdist and float(stopcoordX) < startcoordX + searchdist) and (float(stopcoordY) > startcoordY - searchdist and float(stopcoordY) < startcoordY + searchdist):
                            #if distance is super close then just do this:
                            if distanceR2 (startcoordX, startcoordY, stopcoordX,stopcoordY) < closest_trash_dist:
                                print ('super close trash found!')
                                stopCoord = (stopcoordX,stopcoordY)
                                #for test underlaying 'costandhotspot' is changed back into the previous array'
                                pathArray, costvalue = createPath(CostSurfacefn,costandhotspotarray,startCoord,stopCoord)
                                Multiline = array2multiline(pathArray,CostSurfacefn,pixelValue)
                                geom = ogr.CreateGeometryFromWkt(str(Multiline))
                                length = geom.Length()
                                trashinfo = []
                                trashinfo.append (dict(feature_json_house))
                                trashinfo.append (float(length))
                                trashinfo.append (str(Multiline))
                                trashinfo.append (str(feature_json_trash['properties']['LOCATIENR']))
                                trashinfo.append (float(costvalue))
                                trashinfo.append (str(feature_json_trash['properties']['FRACTIE1']))
                                considerable_trashcans.append (trashinfo)
                                break
                            else:
                                stopCoord = (stopcoordX,stopcoordY)
                                pathArray, costvalue = createPath(CostSurfacefn,costandhotspotarray,startCoord,stopCoord)
                                Multiline = array2multiline(pathArray,CostSurfacefn,pixelValue)
                                geom = ogr.CreateGeometryFromWkt(str(Multiline))
                                length = geom.Length()
                                trashinfo = []
                                trashinfo.append (dict(feature_json_house))
                                trashinfo.append (float(length))
                                trashinfo.append (str(Multiline))
                                trashinfo.append (str(feature_json_trash['properties']['LOCATIENR']))
                                trashinfo.append (float(costvalue))
                                trashinfo.append (str(feature_json_trash['properties']['FRACTIE1']))
                                considerable_trashcans.append (trashinfo)



                #find closest in list & check if empty, make BIGGER search (double?) radius and pic the closest (THIS IS SLOW, BUT NEVER NECESARRY)
                if considerable_trashcans == []:
                    print ('==========Nothing in radius')
                    further_trashcans = {}
                    for i in range(shape_trash.GetFeatureCount()):
                        feature_trash = shape_trash.GetFeature(i)
                        feature_json_trash = feature_trash.ExportToJson()
                        feature_json_trash = eval(feature_json_trash)
                        #weird reason it's multipoint so all the same point, but more of them
                        stopcoordX = feature_json_trash['geometry']['coordinates'][0][0]
                        stopcoordY = feature_json_trash['geometry']['coordinates'][0][1]
                        biggersearch = 200
                        if (float(stopcoordX) > startcoordX - biggersearch and float(stopcoordX) < startcoordX + biggersearch) and (float(stopcoordY) > startcoordY - biggersearch and float(stopcoordY) < startcoordY + biggersearch):
                            stopCoord = (stopcoordX,stopcoordY)
                            distance = float(distanceR2 (startcoordX,startcoordY,stopcoordX,stopcoordY))
                            further_trashcans[distance] = (stopCoord, str(feature_json_trash['properties']['LOCATIENR']))
                    closest_trash = min(further_trashcans)
                    stopCoordNew = further_trashcans[closest_trash][0]
                    LOCATIENR = further_trashcans[closest_trash][1]
                    #for test underlaying 'costandhotspot' is changed back into the previous array'
                    pathArray, costvalue = createPath(CostSurfacefn,costandhotspotarray,startCoord,stopCoordNew)
                    Multiline = array2multiline(pathArray,CostSurfacefn,pixelValue)
                    geom = ogr.CreateGeometryFromWkt(str(Multiline))
                    length = geom.Length()
                    trashinfo = []
                    trashinfo.append (dict(feature_json_house))
                    trashinfo.append (float(length))
                    trashinfo.append (str(Multiline))
                    trashinfo.append (str(LOCATIENR))
                    trashinfo.append (str(feature_json_trash['properties']['FRACTIE1']))
                    f.write(str(trashinfo[2]) + '\t' + str(trashinfo[0]["properties"]["Identifica"]) + '\t' + str(float(trashinfo[1])) + '\t' + str(trashinfo[3]))

                    papier = []
                    glas = []
                    plastic = []
                    for i in range(shape_trash.GetFeatureCount()):
                        feature_trash = shape_trash.GetFeature(i)
                        feature_json_trash = feature_trash.ExportToJson()
                        feature_json_trash = eval(feature_json_trash)
                        stopcoordX = feature_json_trash['geometry']['coordinates'][0][0]
                        stopcoordY = feature_json_trash['geometry']['coordinates'][0][1]
                        biggersearch = 250
                        if (float(stopcoordX) > startcoordX - biggersearch and float(stopcoordX) < startcoordX + biggersearch) and (float(stopcoordY) > startcoordY - biggersearch and float(stopcoordY) < startcoordY + biggersearch):
                            if 'papier' not in trashinfo[4]:
                                if 'papier' in feature_json_trash['properties']['FRACTIE1']:
                                    distance = int(distanceR2(startcoordX,startcoordY,stopcoordX,stopcoordY))
                                    papier.append((feature_json_trash['properties']['LOCATIENR'], distance, (stopcoordX,stopcoordY)))
                            if 'glass' not in trashinfo[4]:
                                if 'glass' in feature_json_trash['properties']['FRACTIE1']:
                                    distance = int(distanceR2(startcoordX,startcoordY,stopcoordX,stopcoordY))
                                    glas.append((feature_json_trash['properties']['LOCATIENR'], distance, (stopcoordX,stopcoordY)))
                            if 'kunststof' not in trashinfo[4]:
                                if 'kunststof' in feature_json_trash['properties']['FRACTIE1']:
                                    distance = int(distanceR2(startcoordX,startcoordY,stopcoordX,stopcoordY))
                                    plastic.append((feature_json_trash['properties']['LOCATIENR'], distance, (stopcoordX,stopcoordY)))
                    
                    if papier ==[]:
                        f.write('\t' + '' + '\t'+ '')
                        Multiline_papier = ''
                        closest_papier = '',''
                    else:
                        closest_papier = min(papier, key = lambda t: t[1])
                        stopCoord= closest_papier[2]
                        pathArray_papier, costvalue = createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord)
                        Multiline_papier = array2multiline(pathArray_papier,CostSurfacefn,pixelValue)
                        f.write ('\t' + str(Multiline_papier) + '\t'+ str(closest_papier[0]))
                    
                    if glas ==[]:
                        f.write('\t' + '' + '\t'+ '')
                        Multiline_glas = ''
                        closest_glas = '',''
                    else:
                        closest_glas = min(glas, key = lambda t: t[1])
                        stopCoord= closest_glas[2]
                        pathArray_glas, costvalue = createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord)
                        Multiline_glas = array2multiline(pathArray_glas,CostSurfacefn,pixelValue)
                        f.write ('\t' + str(Multiline_glas) + '\t'+ str(closest_glas[0]))
                    
                    if plastic ==[]:
                        f.write('\t' + '' + '\t'+ '')
                        Multiline_plastic = ''
                        closest_plastic = '', ''
                    else:
                        closest_plastic = min(plastic, key = lambda t: t[1])
                        stopCoord= closest_plastic[2]
                        pathArray_plastic, costvalue = createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord)
                        Multiline_plastic = array2multiline(pathArray_plastic,CostSurfacefn,pixelValue)
                        f.write ('\t' + str(Multiline_plastic) + '\t'+ str(closest_plastic[0]))
                    f.write('\n')
                    already_made[trashinfo[0]['properties']['id_pand']] = (float(trashinfo[1]),str(trashinfo[2]),str(trashinfo[0]["properties"]["Identifica"]),int(trashinfo[3]), Multiline_papier,closest_papier[0],Multiline_glas, closest_glas[0], Multiline_plastic, closest_plastic[0] )
                    print ('==========exception found')
                else:
                    closest_trash = min(considerable_trashcans, key = lambda t: t[4])
                    f.write(str(closest_trash[2]) + '\t' + str(closest_trash[0]["properties"]["Identifica"]) + '\t' + str(float(closest_trash[1])) + '\t' + str(closest_trash[3]))
                    #here we're checking which additional trashcans can be added, based on the fact that we need all 4 of the kinds (glass, paper, plastic)
                    papier = []
                    glas = []
                    plastic = []
                    for i in range(shape_trash.GetFeatureCount()):
                        feature_trash = shape_trash.GetFeature(i)
                        feature_json_trash = feature_trash.ExportToJson()
                        feature_json_trash = eval(feature_json_trash)
                        stopcoordX = feature_json_trash['geometry']['coordinates'][0][0]
                        stopcoordY = feature_json_trash['geometry']['coordinates'][0][1]
                        biggersearch = 250
                        if (float(stopcoordX) > startcoordX - biggersearch and float(stopcoordX) < startcoordX + biggersearch) and (float(stopcoordY) > startcoordY - biggersearch and float(stopcoordY) < startcoordY + biggersearch):
                            if 'papier' not in closest_trash[5]:
                                if 'papier' in feature_json_trash['properties']['FRACTIE1']:
                                    distance = int(distanceR2(startcoordX,startcoordY,stopcoordX,stopcoordY))
                                    papier.append((feature_json_trash['properties']['LOCATIENR'], distance, (stopcoordX,stopcoordY)))
                            if 'glas' not in closest_trash[5]:
                                if 'glas' in feature_json_trash['properties']['FRACTIE1']:
                                    distance = int(distanceR2(startcoordX,startcoordY,stopcoordX,stopcoordY))
                                    glas.append((feature_json_trash['properties']['LOCATIENR'], distance, (stopcoordX,stopcoordY)))
                            if 'kunststof' not in closest_trash[5]:
                                if 'kunststof' in feature_json_trash['properties']['FRACTIE1']:
                                    distance = int(distanceR2(startcoordX,startcoordY,stopcoordX,stopcoordY))
                                    plastic.append((feature_json_trash['properties']['LOCATIENR'], distance, (stopcoordX,stopcoordY)))
                    
                    if papier ==[]:
                        f.write('\t' + '' + '\t'+ '')
                        Multiline_papier = ''
                        closest_papier = '',''
                    else:
                        closest_papier = min(papier, key = lambda t: t[1])
                        stopCoord= closest_papier[2]
                        pathArray_papier, costvalue = createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord)
                        Multiline_papier = array2multiline(pathArray_papier,CostSurfacefn,pixelValue)
                        f.write ('\t' + str(Multiline_papier) + '\t'+ str(closest_papier[0]))
                    
                    if glas ==[]:
                        f.write('\t' + '' + '\t'+ '')
                        Multiline_glas = ''
                        closest_glas = '',''
                    else:
                        closest_glas = min(glas, key = lambda t: t[1])
                        stopCoord= closest_glas[2]
                        pathArray_glas, costvalue = createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord)
                        Multiline_glas = array2multiline(pathArray_glas,CostSurfacefn,pixelValue)
                        f.write ('\t' + str(Multiline_glas) + '\t'+ str(closest_glas[0]))
                    
                    if plastic ==[]:
                        f.write('\t' + '' + '\t'+ '')
                        Multiline_plastic = ''
                        closest_plastic = '', ''
                    else:
                        closest_plastic = min(plastic, key = lambda t: t[1])
                        stopCoord= closest_plastic[2]
                        pathArray_plastic, costvalue = createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord)
                        Multiline_plastic = array2multiline(pathArray_plastic,CostSurfacefn,pixelValue)
                        f.write ('\t' + str(Multiline_plastic) + '\t'+ str(closest_plastic[0]))
                    f.write('\n')

                    already_made[closest_trash[0]['properties']['id_pand']] = (float(closest_trash[1]),str(closest_trash[2]),str(closest_trash[0]["properties"]["Identifica"]),int(closest_trash[3]), Multiline_papier,closest_papier[0],Multiline_glas, closest_glas[0], Multiline_plastic, closest_plastic[0] )
                    print ('==========best found=========') 

    #for now not important yet
    array2raster(CostSurfacefn,newRasterfn,costandhotspotarray)

def raster2array(rasterfn):
    raster = gdal.Open(rasterfn)
    band = raster.GetRasterBand(1)
    array = band.ReadAsArray()
    print ('array made')
    return array, band.GetNoDataValue()

def coord2pixelOffset(rasterfn,x,y):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    xOffset = int((x - originX)/pixelWidth)
    yOffset = int((y - originY)/pixelHeight)
    return xOffset,yOffset

def createPath(CostSurfacefn,costSurfaceArray,startCoord,stopCoord):

    # coordinates to array index
    startCoordX = startCoord[0]
    startCoordY = startCoord[1]
    startIndexX,startIndexY = coord2pixelOffset(CostSurfacefn,startCoordX,startCoordY)

    stopCoordX = stopCoord[0]
    stopCoordY = stopCoord[1]
    stopIndexX,stopIndexY = coord2pixelOffset(CostSurfacefn,stopCoordX,stopCoordY)

    # create path
    indices, weight = route_through_array(costSurfaceArray, (startIndexY,startIndexX), (stopIndexY,stopIndexX),geometric=True,fully_connected=True)
    indices = np.array(indices).T
    path = np.zeros_like(costSurfaceArray)
    path[indices[0], indices[1]] = 1
    print ('path found')
    return path, weight

def pixelOffset2coord(rasterfn,xOffset,yOffset):
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    originX = geotransform[0]
    originY = geotransform[3]
    pixelWidth = geotransform[1]
    pixelHeight = geotransform[5]
    coordX = originX+pixelWidth*xOffset
    coordY = originY+pixelHeight*yOffset
    return coordX, coordY

def array2multiline(array,rasterfn,pixelValue):

    # max distance between points
    raster = gdal.Open(rasterfn)
    geotransform = raster.GetGeoTransform()
    pixelWidth = geotransform[1]
    maxDistance = ceil(sqrt(2*pixelWidth*pixelWidth))

    # array2dict
    count = 0
    roadList = np.where(array == pixelValue)
    multipoint = ogr.Geometry(ogr.wkbMultiLineString)
    pointDict = {}
    for indexY in roadList[0]:
        indexX = roadList[1][count]
        Xcoord, Ycoord = pixelOffset2coord(rasterfn,indexX,indexY)
        pointDict[count] = (Xcoord, Ycoord)
        count += 1

    multiline = ogr.Geometry(ogr.wkbMultiLineString)
    for i in itertools.combinations(pointDict.values(), 2):
        point1 = ogr.Geometry(ogr.wkbPoint)
        point1.AddPoint(i[0][0],i[0][1])
        point2 = ogr.Geometry(ogr.wkbPoint)
        point2.AddPoint(i[1][0],i[1][1])

        distance = point1.Distance(point2)
        if distance < maxDistance:
            line = ogr.Geometry(ogr.wkbLineString)
            line.AddPoint(i[0][0],i[0][1])
            line.AddPoint(i[1][0],i[1][1])
            multiline.AddGeometry(line)
    print ('multiline made')
    return (multiline)

def Bresenham_with_rasterio(CostSurfacefn, shortestpath_geojson, cellvalue):
    """
    Returns rasterised line as an input-raster-sized array with zeroes,
    and ones for the line pixels.
    Input is a list with (x,y) coordinates, already given in CRS of TIF.
    """
    d = rasterio.open(CostSurfacefn)
    #if you multiply this value, I make sure the outcome is simular so. The non-touched values are 1! not 0!. So to get 0.8 is -0.2
    cellvalue = cellvalue - 1

    shapes = [(shortestpath_geojson, cellvalue)]
    re = features.rasterize(shapes, out_shape=d.shape, all_touched=False,transform=d.transform)
    re = re + np.ones_like(re)
    """Return is an array!!"""
    return re

def distanceR2(x1, y1, x2, y2):
    """
    Distance in R^2.
    """
    return sqrt( (x2 - x1)**2 + (y2 - y1)**2 )

def array2raster(rasterfn,newRasterfn,array):
    print ('made raster')
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

def directiontohotspot (housexy, hotspotxy, buff):

    #if hotspot is left of house
    if housexy[0] > hotspotxy[0]:
        angle_rad = atan((housexy[1]-hotspotxy[1])/(housexy[0]-hotspotxy[0]))
        angle_deg = degrees(angle_rad)
        distance = distanceR2(housexy[0],housexy[1], hotspotxy[0], hotspotxy[1])
        new_distance = distance + 5
        change_in_y = sin(angle_rad)*new_distance
        change_in_x = cos(angle_rad)*new_distance
        housexy = (hotspotxy[0] + change_in_x, hotspotxy[1]+change_in_y)

    #if hotspot is right of house:
    if housexy[0] < hotspotxy[0]:
        angle_rad = atan((hotspotxy[1]-housexy[1])/(hotspotxy[0]-housexy[0]))
        angle_deg = degrees(angle_rad)
        distance = distanceR2(housexy[0],housexy[1], hotspotxy[0], hotspotxy[1])
        new_distance = distance + 5
        change_in_y = sin(angle_rad)*new_distance
        change_in_x = cos(angle_rad)*new_distance
        housexy = (hotspotxy[0] - change_in_x, hotspotxy[1]-change_in_y)

    line = LineString([Point(housexy[0], housexy[1]), Point(hotspotxy[0], hotspotxy[1])])
    
    polygon = line.buffer(buff, cap_style=2)
    print (polygon)
    return polygon

def shortestpathtohotspot (G, housexy, hotspotxy, node_buff, edge_buff):
    source_point = Point(housexy[0], housexy[1])
    target_point = Point(hotspotxy[0], hotspotxy[1])
    source_pointLL = RD2ll (source_point)
    source_pointLL = shapely.wkt.loads(source_pointLL)
    target_pointLL = RD2ll (target_point)
    target_pointLL = shapely.wkt.loads(target_pointLL)
    startnode = find_nearest_node (G, source_pointLL.x, source_pointLL.y)
    targetnode = find_nearest_node (G, target_pointLL.x, target_pointLL.y)

    shortest_path = nx.dijkstra_path(G, source=startnode, target=targetnode, weight=((data['distance']) for data in G.edges(data=True, keys=True)))
    shortest_path_edge = zip(shortest_path,shortest_path[1:])
    subgraph = G.subgraph (shortest_path)
    nodegraph = G.subgraph(startnode)
    node_points = [Point((data['x'], data['y'])) for node, data in subgraph.nodes(data=True)]
    node_points_RD = []
    """For loop under is specially to transform from LL to RD..
    because the nodes are to RD the edges go automatically """
    for node_point in node_points:
        node_point_RD = ll2RD(Point(node_point))
        node_points_RD.append(shapely.wkt.loads(node_point_RD))

    nodes_gdf = gpd.GeoDataFrame({'id': subgraph.nodes()}, geometry=node_points_RD)
    nodes_gdf = nodes_gdf.set_index('id')
    edge_lines = []
    for n_fr, n_to in subgraph.edges():
        f = nodes_gdf.loc[n_fr].geometry
        t = nodes_gdf.loc[n_to].geometry
        edge_lines.append(LineString([f,t]))
    #can we add the line 'closest node' -> housepoint here?
    close_node = [Point((data['x'], data['y'])) for node, data in nodegraph.nodes(data=True)]
    close_node = close_node[0]
    close_node = ll2RD(Point(node_point))
    close_node = shapely.wkt.loads(close_node)
    edge_lines.append(LineString([source_point,close_node]))
    n = nodes_gdf.buffer(node_buff).geometry
    e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
    all_gs = list(n) + list(e)
    new_path = gpd.GeoSeries(all_gs).unary_union
    return (new_path)

def WKT2GEOJSON (WKT):
    g1 = shapely.wkt.loads(WKT)
    g2 = geojson.Feature(geometry=g1, properties={})
    return g2.geometry

def ll2RD (geometry):
    #4326 is crs from OSM
    source = osr.SpatialReference()
    source.ImportFromEPSG(4326)

    target = osr.SpatialReference()
    target.ImportFromEPSG(28992)

    transform = osr.CoordinateTransformation(source, target)

    gdalgeom = ogr.CreateGeometryFromWkt(str(geometry))
    gdalgeom.Transform(transform)

    return (gdalgeom.ExportToWkt())

def RD2ll (geometry):
    #4326 is crs from OSM
    source = osr.SpatialReference()
    source.ImportFromEPSG(28992)

    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(source, target)

    gdalgeom = ogr.CreateGeometryFromWkt(str(geometry))
    gdalgeom.Transform(transform)

    return (gdalgeom.ExportToWkt())

def find_nearest_node (G, x, y):
    node = ox.get_nearest_node(G, (y, x))
    return node

if __name__ == "__main__":
    start_time = time.time()
    CostSurfacefn = '/Users/daveyoldenburg 1/Desktop/Test for LCA/CostRaster_size1_filled.tif'
    startCoordShp = '/Users/daveyoldenburg 1/Desktop/Test for LCA/Verblijfsobject_ams_test.shp'
    newRasterfn = '/Users/daveyoldenburg 1/Desktop/Test for LCA/Costandhotspotarrayrast.tif'
    stopCoordShp = '/Users/daveyoldenburg 1/Desktop/Thesis/Gemeente Containers/Containerlocations_GM.shp'
    WKT_file = '/Users/daveyoldenburg 1/Desktop/Test for LCA/Method1_size1_0.1.tsv'
    Hotspot = '/Users/daveyoldenburg 1/Desktop/Test for LCA/Density Clusters Middlepoint/Middlepoint of clusters.shp'
    pixelValue = 1
    searchdist = 100
    closest_trash_dist = 20
    linebuff = 200
    #these are precentages now!
    rasteradjustment = 0.1
    main(CostSurfacefn,startCoordShp,stopCoordShp,pixelValue, searchdist, WKT_file, closest_trash_dist, newRasterfn, linebuff, Hotspot, rasteradjustment)
    print("--- %s seconds ---" % (time.time() - start_time))



