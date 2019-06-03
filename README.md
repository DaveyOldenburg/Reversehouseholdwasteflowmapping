# Reversehouseholdwasteflowmapping
Geomatics MSc Thesis: Reverse Household Waste Flow Mapping

Our linear way of producing has to come to an end. The world population is estimated to grow to 9 billion people in 2050. All these people have needs and to answer to these needs, we require raw materials. However, almost everything we consume is finite. These limitations result in fluctuations in prices and scarcity. More and more companies and government institutions see the value of circularity of raw materials. 

The direction of this thesis is created together with the REPAiR program by the AMS Institute. This thesis is building up upon research and data provided by the H2020 project “REPAiR - Resource Management in Peri-Urban Areas. Going Beyond Urban Metabolism” led by prof. Arjan van Timmeren. In this research, REPAiR, among other things, visualised a few steps of the waste chain and mapped them for Amsterdam from the level of the neighbourhood to the waste treatment plant. 

This thesis maps the waste flows in the opposite direction. The starting point is the waste storage in public space (this was also the start point in the REPAiR research) and follows the waste back to the supplier. The households are found that use the starting point (waste storage) (Method 1, household to trashcan), whereafter the suppliers of the food from the households are found (Method 2, household to supplier).

Since mapping on household level is never done before, potential behaviour data is missing. To be able to map on this level of detail, the relation between different nodes from the food waste chain are found. This relation is therefor mapped based on environmental factors, like where to walk, and where to go. 

Main question: How can reverse waste flows for food waste be mapped based on environmental factors?

Method 1:
To connect the household to the supplier, the following input is used: the households, the trashcans, the hotspot analysis and the least cost raster. The last two are based on environmental factors, so is the hotspot analysis created by density based matching based on food suppliers and the least cost raster created by rasterising the BGT giving walkability values.  

Method 2:
To connect the household to the food suppliers, the following input is used: the households, the food suppliers and an OSM road network. By creating an isochrone based on a certain bike-time, we can find the food suppliers that are potentially used by the household. Knowing that bigger and closer stores are more likely to be used, possibility values can be given to every potential food supplier.

Validation (Correlation.R):
To validate method 1, the correlation is found between the amount of trash in a certain trashcan and the number of households that use that trashcan. The manual parameter used in method 1 is adjusted until this correlation is the highest. Method 2 is considered a demonstration, since an accurate validation method is not available, due to missing data.


Algorithms:

MethodHouseholdtotrashcan.py - Creating the spatial connection between households and their used trashcans
MethodHouseholdtosupplier.py - Creating the spatial connection between households and their food suppliers

Additional Algorithms:

Zonalstatistics.py - Filling gaps in the created least cost raster and fixing minor errors

Correlations.R - Plots the correlation between the number of households assigned to a trashcan and the mass per year of this trashcan.

