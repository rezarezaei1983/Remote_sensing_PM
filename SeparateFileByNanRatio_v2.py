
# Ref :  https://www.youtube.com/watch?v=RdAFsMhdMOY

#       https://www.guillaumedueymes.com/post/shapefiles_country/
#       https://code.mpimet.mpg.de/boards/2/topics/13006?r=13008

import os, sys
from netCDF4 import Dataset
import pandas as pd
import datetime
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore")



data = Dataset("C:/Users/Reza/Desktop/Deniz_data/chirps-v2.0.2022.days_p05.nc", "r")


#print("\n>>> Variables : \n", data.variables.keys())
#print("\n>>> Data : \n", data)

precip = data.variables["precip"][0,:,:]
lat = data.variables["latitude"][:]
lon = data.variables["longitude"][:]
times = data.variables["time"]

lats, lons = np.meshgrid(lat, lon)




import xarray as xr
import numpy as np
import regionmask
import geopandas as gpd
import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

 
countries = gpd.read_file("C:/Users/Reza/Desktop/Deniz_data/shapefile/TUR_adm0.shp")
#print(countries)

### draw shapefile
fig, ax = plt.subplots(figsize=(16,10))
countries.plot(ax=ax, column = "NAME_ISO") # get "NAME_ISO" from print(countries)


### If we had more then one country/state in the shape file, wew can select the
#    desired sountry/state running the following lines:
my_list = list(countries["ISO"])
print("\n>>> my_list : ")
print(my_list)
my_list_unique = set(my_list)
print("\n>>> my_list_unique : ")
print(my_list_unique)

indexes = [my_list.index(x) for x in my_list_unique]
print("\n>>> indexes : ")
print(indexes)


countries_mask_poly = regionmask.Regions_cls(name="NAME_ISO", numbers=indexes, 
                                             names=countries.NAME_ISO[indexes], 
                                             abbrevs=countries.NAME_ISO[indexes], 
                                             outlines=list(countries.geometry.values[i] for i in range(0,countries.shape[0])))

print("\n>>> countries_mask_poly : ")
print(countries_mask_poly)



# Ref : https://www.youtube.com/watch?v=RdAFsMhdMOY


#             Up to 28':30"