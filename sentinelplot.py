import netCDF4
import glob
import os
import xarray as xr
import rioxarray
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

eur_extent={'name': 'eur_extent',
    'lonmin':-10.0,
            'lonmax':30.0,
            'latmin':35.0,
            'latmax':60.0}

tr_extent={'name': 'tr_extent',
    'lonmin':25.0,
            'lonmax':45.0,
            'latmin':20.0,
            'latmax':50.0}

# Italy (Po Valley, Milano)
# povalley_extent={'name': 'povalley_extent',
#     'lonmin':7.0,
#             'lonmax':13.0,
#             'latmin':44.0,
#             'latmax':47.0}



# China (Hubei, Wuhan)
# hubei_extent={'name': 'hubei_extent',
#     'lonmin':108.3,
#             'lonmax':116.1,
#             'latmin':29.1,
#             'latmax':33.3}


region_list =[eur_extent,
             tr_extent]


fileList = glob.glob('/media/elviradeniz/thesis/data/L2_data/S5P_OFFL_L2__AER_AI_20180628T102407_20180628T120537_03661_01_010002_20180704T095226.nc')
fileList_sort = sorted(fileList)

# Create the bins based on a 0.036 x 0.036 grid
lat_bins = np.arange(-90,90+(0.036/2),0.036)
lon_bins = np.arange(-180,180+(0.036/2), 0.036)

# define a label for each bin corresponding to the central latitude
lat_center = np.arange(-90+(0.036/2),90,0.036)
lon_center = np.arange(-180+(0.036/2),180,0.036)

unit = 'mol/cm2'
longname= 'Aerosol Index'

for i in range(0,len(fileList_sort)):
        print(i)
        # Open the COG file with rasterio
        tmp = xr.open_rasterio(fileList_sort[i])
        tmp_name = tmp.split('_')[8]
        print(tmp_name)
        tmp= tmp.rename({'x': 'lon', 'y':'lat'})
        # Flag out negative values
        tmp_flag = tmp.where(tmp>0,np.nan)
        
        # Bring NO2 values onto a regular latitude grid and create the average of multiple values per cell
        tmp_regrid_lat = tmp_flag.groupby_bins('lat', lat_bins, labels=lat_center).mean()
        # Bring NO2 values onto a regular longitude grid
        tmp_regrid = tmp_regrid_lat.groupby_bins('lon', lon_bins, labels=lon_center).mean()
        
        # Create a xarray.DataArray with the regridded array and the define latituden and longitude bins
        data_array = xr.DataArray(
                tmp_regrid.isel(band=0).values,
                dims=['lat','lon'],
                coords={
                    'time': pd.to_datetime(tmp_name),
                    'lat':(['lat'],tmp_regrid.isel(band=0).lat_bins),
                    'lon':(['lon'],tmp_regrid.isel(band=0).lon_bins)
                },
                attrs={'long_name': longname, 'units': unit},
                name='AI'
            )
        
        # Save the created xarray.DatArray as netCDF file
        data_array.to_netcdf('./S5P_OFFL_L2__AER_AI'+tmp_name+'.nc', 'w')  


aai= xr.open_mfdataset('/media/elviradeniz/thesis/data/L2_data/S5P_*', concat_dim='time', combine='nested')
aai=ndjfmam.aai
aai



aai

aai_monthly= aai.resample(time='M',skipna=False).mean()


conversion_factor=6.02214*1e+19

aai = aai_monthly*conversion_factor



region = region_list[0] # Choose first region in region list

# Rename lat/lon to latitude and longitude
aai_monthly = aai_monthly.rename({'lon': 'longitude', 'lat': 'latitude'})


eur = generate_geographical_subset(xarray=aai_monthly,
                                   latmin=region['latmin'],
                                   latmax=region['latmax'],
                                   lonmin=region['lonmin'],
                                   lonmax=region['lonmax'])



eur

month=5 #Month April
visualize_pcolormesh(data_array=eur.isel(time=month)*1e-15, 
                     longitude=eur.longitude, 
                     latitude=eur.latitude, 
                     projection=ccrs.PlateCarree(), 
                     color_scale='viridis', 
                     unit='*1e-15 molec./cm2', 
                     long_name=aai.long_name + ' April 2020', 
                     vmin=0, 
                     vmax=8, 
                     lonmin=region['lonmin'],
                     lonmax=region['lonmax'],
                     latmin=region['latmin'],
                     latmax=region['latmax'],
                     set_global=False)

month = 5 # April
anomaly_1920 = aai.isel(time=month) - aai_monthly.isel(time=month)
anomaly_1920

region = region_list[0]
eur_anomaly = generate_geographical_subset(xarray=anomaly_1920,
                                           latmin=region['latmin'],
                                           latmax=region['latmax'],
                                           lonmin=region['lonmin'],
                                           lonmax=region['lonmax'])


visualize_pcolormesh(data_array=eur_anomaly*1e-15, 
                     longitude=eur_anomaly.longitude, 
                     latitude=eur_anomaly.latitude, 
                     projection=ccrs.PlateCarree(), 
                     color_scale='RdBu_r', 
                     unit='*1e-15 mol/cm2', 
                     long_name=no2_tropo_2019.long_name + 'Anomaly April 2020 to April 2019', 
                     vmin=-3, 
                     vmax=3, 
                     lonmin=region['lonmin'],
                     lonmax=region['lonmax'],
                     latmin=region['latmin'],
                     latmax=region['latmax'],
                     set_global=False)