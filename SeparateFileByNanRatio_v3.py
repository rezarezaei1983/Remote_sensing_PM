

# Ref : https://www.earthdatascience.org/courses/use-data-open-source-python/intro-raster-data-python/raster-data-processing/crop-raster-data-with-shapefile-in-python/

# Alternative method : 
#       https://mygeoblog.com/2019/06/25/mask-netcdf-using-shp-file/  
          

# Shapefile and NetCDF files are in the same coordinate system (epsg:4326)

import os
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import mapping
#import rioxarray as rxr
import xarray as xr
import geopandas as gpd
from netCDF4 import Dataset

data = Dataset("C:/Users/Reza/Desktop/Deniz_data/chirps-v2.0.2022.days_p05.nc", "r")

#print("\n>>> Variables : \n", data.variables.keys())
#print("\n>>> Data : \n", data)

precip = data.variables["precip"][0,:,:]
lat = data.variables["latitude"][:]
lon = data.variables["longitude"][:]
times = data.variables["time"]

#lats, lons = np.meshgrid(lat, lon)




crop_extent = gpd.read_file("C:/Users/Reza/Desktop/Deniz_data/shapefile/TUR_adm0.shp")

country_name = list(crop_extent["NAME_ISO"])[0]

### draw shapefile
fig, ax = plt.subplots(figsize=(16,10))
crop_extent.plot(ax=ax, column = "NAME_ISO") # get "NAME_ISO" from print(countries)
ax.set_title(f"{country_name} Shapefile", fontsize=16)






#print('crop extent crs: ', crop_extent.crs)


###### how to clip/crop netcdf file using shapefile


