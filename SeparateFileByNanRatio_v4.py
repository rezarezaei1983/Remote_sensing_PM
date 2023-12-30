

# Ref : https://stackoverflow.com/questions/70241007/mask-netcdf-using-shapefile-and-calculate-average-and-anomaly-for-all-polygons-w
import numpy as np
import regionmask
import xarray as xr
import geopandas as gpd


# load polygons of US states
countris = regionmask.defined_regions.natural_earth_v5_0_0.countries_110

print(countris)
print(countris["TR"])

TRshape = countris["TR"]
