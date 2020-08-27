# raster_extraction

This script is used to read multiple rasters from tif files, 
using a loop to get them all in the same crs and extent so they can be packed into a rasterStack
From the rasterStacks, estimates are extracted within IU polygons
The rasterStacks are then extracted to 5kmx5km pixels to produce a covariates matrix which is an efficient storage medium
