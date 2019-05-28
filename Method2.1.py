import osmnx as ox, networkx as nx, geopandas as gpd
from shapely.geometry import Point, LineString, Polygon
from shapely.ops import cascaded_union
import shapely.wkt
from descartes import PolygonPatch
import ogr, osr
import time
import math

ox.config(log_console=False, use_cache=False)

def main(place, speed_dict, houseshp, trip_time, actorsshp, WKT_file,Actors_household_file,Hotspot):

	#get houses
	file_house = ogr.Open(houseshp)
	shape_house = file_house.GetLayer(0)
	#error fix

	file_actor = ogr.Open(actorsshp)
	shape_actor = file_actor.GetLayer(0)

	null = None
	graph_walk = get_network (place, 'bike')
	print ('Downloaded network')

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

    #biking network for now!

	with open(WKT_file, 'w') as f:
		with open(Actors_household_file, 'w') as w:
			"""Write the first line of the WKT File"""
			f.write(str('Geometry') + '\t' +'Identifica' + '\t' + 'id_supplier_1' + '\t' + 'id_supplier_2'+ '\t' + 'id_supplier_3'+'\n')
			w.write(str('house_id') + '\t' +'store_id' + '\t' + 'influence'+ '\t' + 'id_pand'+'\n')

			already_made = {}
			already_made2 = {}
			#every household
			#for i in range(shape_house.GetFeatureCount()):
			for i in range (5):
				feature_house = shape_house.GetFeature(i)
				feature_json_house = feature_house.ExportToJson()
				feature_json_house = eval(feature_json_house)
				startcoordX = feature_json_house['geometry']['coordinates'][0]
				startcoordY = feature_json_house['geometry']['coordinates'][1]

				"""find_nearest_node takes LL so transform"""
				point = Point(startcoordX, startcoordY)
				point = RD2ll(point)
				point_ll = shapely.wkt.loads(point)

				#finding closest hotspot
				hotspotdistancedict = {}
				small = 999999
				for hotspot in hotspotlist:
					distance = distanceR2(startcoordX, startcoordY, hotspot[0], hotspot[1])
					hotspotdistancedict[distance] = hotspot
					if distance < small:
				   		small = distance
				hotspotxy = hotspotdistancedict[small]

				#check if polygon is already made
				if feature_json_house['properties']['id_pand'] in already_made:
					print ('duplicate found')
					f.write(str(already_made[feature_json_house['properties']['id_pand']][0])+ '\t' + str(feature_json_house['properties']["Identifica"])  +'\t'+ ",".join(str(i[0]) for i in already_made[feature_json_house['properties']['id_pand']][1])+'\t'+ ",".join(str(i[0]) for i in already_made[feature_json_house['properties']['id_pand']][2])+'\t'+ ",".join(str(i[0]) for i in already_made[feature_json_house['properties']['id_pand']][3])+'\n')

					importance_actors_dict = already_made2[feature_json_house['properties']['id_pand']]
					for actor in importance_actors_dict:
						w.write(str(feature_json_house['properties']["Identifica"])+ '\t' + str(actor) + '\t' + str(importance_actors_dict[actor][0])+ '\t' + str(feature_json_house['properties']['id_pand'])+'\n')


				else:
					prep_walk, house_node = find_nearest_node(graph_walk, point_ll.x, point_ll.y)
					print (house_node)
					G_time_walk = add_attribute_time (prep_walk, speed_dict['walk'])

					isocrone1 = make_iso_poly(G_time_walk, house_node, 0.33*trip_time, edge_buff=35, node_buff=50, infill=True)
					isocrone2 = make_iso_poly(G_time_walk, house_node, 0.66*trip_time, edge_buff=35, node_buff=50, infill=True)
					isocrone3 = make_iso_poly(G_time_walk, house_node, trip_time, edge_buff=35, node_buff=50, infill=True)
					shortestpath_wkt = shortestpathtohotspot (G_time_walk, (startcoordX,startcoordY), hotspotxy, 50, 35)
					

					isocrone1_geom = shapely.wkt.loads(str(isocrone1))
					isocrone2_geom = shapely.wkt.loads(str(isocrone2))
					isocrone3_geom = shapely.wkt.loads(str(isocrone3))
					shortestpath_geom = shapely.wkt.loads(str(shortestpath_wkt))
					isocrone3_geom = cascaded_union([isocrone3_geom, shortestpath_geom])

					isocrone2_geom = isocrone2_geom.difference(isocrone1_geom)
					isocrone3_geom = (isocrone3_geom.difference(isocrone2_geom)).difference(isocrone1_geom)

					#under this takes super long
					actors_in_isocrone1, actors_in_isocrone2, actors_in_isocrone3 = get_actors_in_isocrone(isocrone1_geom, isocrone2_geom, isocrone3_geom, shape_actor)
					
					f.write(str(isocrone1) + '\t' + str(feature_json_house['properties']['Identifica']) + '\t' + ",".join(str(i[0]) for i in actors_in_isocrone1) +'\t'+ ",".join(str(i[0]) for i in actors_in_isocrone2) +'\t'+ ",".join(str(i[0]) for i in actors_in_isocrone3) +'\n')
					already_made[feature_json_house['properties']['id_pand']] = str(isocrone1), actors_in_isocrone1, actors_in_isocrone2, actors_in_isocrone3
					print ('=========== actors found! ============')

					importance_dict = {}
					if actors_in_isocrone1 != []:
						
						for actor in actors_in_isocrone1:
							feature_actor = shape_actor.GetFeature(actor[1])
							feature_json_actor = feature_actor.ExportToJson()
							#bugfix
							null = None
							feature_json_actor = eval(feature_json_actor)
							
							importance_dict[feature_json_actor['properties']['OBJECTID'], 1] = feature_json_actor['properties']['Oppervlakt']
							#w.write(str(feature_json_house['properties']["Identifica"])+ '\t'+ str(feature_json_actor['properties']['OBJECTID']) +'\t'+ str(importance) +'\t'+ str(feature_json_house['properties']['id_pand'])+ '\n')
					# if the first one is empty, put the 50% importance in the second isocrone (smart?, NOPE)
					
					if actors_in_isocrone2 != []:
						
						for actor in actors_in_isocrone2:
							feature_actor = shape_actor.GetFeature(actor[1])
							feature_json_actor = feature_actor.ExportToJson()
							#bugfix
							null = None
							feature_json_actor = eval(feature_json_actor)
							
							importance_dict[feature_json_actor['properties']['OBJECTID'], 2] = feature_json_actor['properties']['Oppervlakt']
							#w.write(str(feature_json_house['properties']["Identifica"])+ '\t'+ str(feature_json_actor['properties']['OBJECTID']) +'\t'+ str(importance) +'\t'+ str(feature_json_house['properties']['id_pand'])+ '\n')

					if actors_in_isocrone3 != []:
						
						for actor in actors_in_isocrone3:
							feature_actor = shape_actor.GetFeature(actor[1])
							feature_json_actor = feature_actor.ExportToJson()
							#bugfix
							null = None
							feature_json_actor = eval(feature_json_actor)
							
							importance_dict[feature_json_actor['properties']['OBJECTID'], 3] = feature_json_actor['properties']['Oppervlakt']
							#w.write(str(feature_json_house['properties']["Identifica"])+ '\t'+ str(feature_json_actor['properties']['OBJECTID']) +'\t'+ str(importance) +'\t'+ str(feature_json_house['properties']['id_pand'])+ '\n')

					total = 0
					for actor in importance_dict:
						if actor[1] == 1:
							total += 3 * importance_dict[actor]
						elif actor[1] == 2:
							total += 2 * importance_dict[actor]
						elif actor[1] == 3:
							total += importance_dict[actor]

					importance_actors_dict = {}

					for actor in importance_dict:
						if actor[1] == 1:
							importance_normal = 3 * importance_dict[actor]
						elif actor[1] == 2:
							importance_normal = 2 * importance_dict[actor]
						elif actor[1] == 3:
							importance_normal = 1 * importance_dict[actor]

						importance = (importance_normal/total) *100
						importance_actors_dict[actor[0]]= importance, actor[1]
						w.write(str(feature_json_house['properties']["Identifica"]) + '\t' + str(actor[0]) + '\t' + str(importance)+ '\t' + str(feature_json_house['properties']['id_pand'])+'\n')

					already_made2[feature_json_house['properties']['id_pand']] = importance_actors_dict



def get_actors_in_isocrone(isocrone1_geom, isocrone2_geom, isocrone3_geom, shape_actor):
	"""takes isocrone in shapely format, actors as shape, 
	returns list of actors"""
	actors_in_isocrone1 = []
	actors_in_isocrone2 = []
	actors_in_isocrone3 = []

	for i in range(shape_actor.GetFeatureCount()):
		feature_actor = shape_actor.GetFeature(i)
		feature_json_actor = feature_actor.ExportToJson()
		#bugfix
		null = None
		feature_json_actor = eval(feature_json_actor)
		#below statement is not necesarry after filtering of shapefile, but for the idea:
		if feature_json_actor['properties']['Lvl1'] == 'G WHOLESALE AND RETAIL TRADE; REPAIR OF MOTOR VEHICLES AND MOTORCYCLES':
			actor_x = feature_json_actor['geometry']['coordinates'][0]
			actor_y = feature_json_actor['geometry']['coordinates'][1]
			actor_point = Point(actor_x, actor_y)
			#we're gonna append the 'i' here too, in this case we can find the feature easily back. Quicker later
			if isocrone1_geom.contains(actor_point):
				actors_in_isocrone1.append((feature_json_actor['properties']['OBJECTID'], i))
			elif isocrone2_geom.contains(actor_point):
				actors_in_isocrone2.append((feature_json_actor['properties']['OBJECTID'], i))
			elif isocrone3_geom.contains(actor_point):
				actors_in_isocrone3.append((feature_json_actor['properties']['OBJECTID'], i))
	return actors_in_isocrone1, actors_in_isocrone2, actors_in_isocrone3

def get_network (place, mode):
	"""Returns network in LL"""
	G = ox.graph_from_place(place, network_type=mode)
	#makes UTM of LAT LONG
	#G = ox.project_graph(G)
	#ox.plot_graph(G)
	return G

def find_nearest_node (G, x, y):
	node = ox.get_nearest_node(G, (y, x))

	return G, node

def add_attribute_time (G, mode_speed):
	meters_per_minute = mode_speed * 1000/60
	for u, v, k, data in G.edges(data=True, keys=True):
		data['time'] = data['length'] / meters_per_minute
	return G

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

def make_iso_poly(G, house_node, trip_time, edge_buff, node_buff, infill):
	"""Returns a polygon in wkt """
	subgraph = nx.ego_graph(G, house_node, radius=trip_time, distance='time')
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
	n = nodes_gdf.buffer(node_buff).geometry
	e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
	all_gs = list(n) + list(e)
	new_iso = gpd.GeoSeries(all_gs).unary_union
    
    # try to fill in surrounded areas so shapes will appear solid and blocks without white space inside them
	if infill:
		new_iso = Polygon(new_iso.exterior)
	return '{}'.format(new_iso)

	
def shortestpathtohotspot (G, housexy, hotspotxy, node_buff, edge_buff):
	source_point = Point(housexy[0], housexy[1])
	target_point = Point(hotspotxy[0], hotspotxy[1])
	source_pointLL = RD2ll (source_point)
	source_pointLL = shapely.wkt.loads(source_pointLL)
	target_pointLL = RD2ll (target_point)
	target_pointLL = shapely.wkt.loads(target_pointLL)
	startnode = find_nearest_node (G, source_pointLL.x, source_pointLL.y)
	targetnode = find_nearest_node (G, target_pointLL.x, target_pointLL.y)

	shortest_path = nx.dijkstra_path(G, source=startnode[1], target=targetnode[1], weight=((data['distance']) for data in G.edges(data=True, keys=True)))
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
	n = nodes_gdf.buffer(node_buff).geometry
	e = gpd.GeoSeries(edge_lines).buffer(edge_buff).geometry
	all_gs = list(n) + list(e)
	new_path = gpd.GeoSeries(all_gs).unary_union
	return '{}'.format(new_path)

def distanceR2(x1, y1, x2, y2):
	"""
	Distance in R^2.
	"""
	return math.sqrt( (x2 - x1)**2 + (y2 - y1)**2 )


if __name__ == "__main__":
	start_time = time.time()
	place = 'Amsterdam , Netherlands'
	speed_dict = {'walk': 4.5, 'bike': 13.5, 'drive': 15.0}
	trip_time = 15
	houseshp = '/Users/daveyoldenburg 1/Desktop/Test for LCA/Verblijfsobject_ams_test.shp'
	actorsshp =  '/Users/daveyoldenburg 1/Desktop/Test for Method2/T3.2_AMA_Companies_food_waste_chain_located_RD_new_FILTERED_opp.shp'
	WKT_file = '/Users/daveyoldenburg 1/Desktop/Test for Method2/Method2.tsv'
	Actors_household_file = '/Users/daveyoldenburg 1/Desktop/Test for Method2/Influence_actors_try.tsv' 
	Hotspot = '/Users/daveyoldenburg 1/Desktop/Test for LCA/Density Clusters Middlepoint/Middlepoint of clusters.shp'
	main(place, speed_dict, houseshp, trip_time, actorsshp, WKT_file,Actors_household_file, Hotspot)
	print("--- %s seconds ---" % (time.time() - start_time))

