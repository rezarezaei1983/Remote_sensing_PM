
import string
import pandas as pd
from datetime import datetime
import glob
import netCDF4 as nc
import os
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
from matplotlib.axes import Axes
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
import warnings
warnings.simplefilter(action = "ignore", category = UserWarning)

def visualize_pcolormesh(data_array, longitude, latitude, projection, color_scale, 
                         vmin, vmax, var_name, set_global=True, lonmin=-180, 
                         lonmax=180, latmin=-90, latmax=90):

    fig=plt.figure(figsize=(20, 10))

    ax = plt.axes(projection=projection)
   
    img = plt.pcolormesh(longitude, latitude, data_array, 
                        cmap=plt.get_cmap(color_scale), transform=ccrs.PlateCarree(),
                        vmin=vmin,
                        vmax=vmax,
                        shading='auto')

    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=1)

    if (projection==ccrs.PlateCarree()):
        ax.set_extent([lonmin, lonmax, latmin, latmax], projection)
        gl = ax.gridlines(draw_labels=True, linestyle='--')
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter=LONGITUDE_FORMATTER
        gl.yformatter=LATITUDE_FORMATTER
        gl.xlabel_style={'size':14}
        gl.ylabel_style={'size':14}

    if(set_global):
        ax.set_global()
        ax.gridlines()

    cbar = fig.colorbar(img, ax=ax, orientation='horizontal', fraction=0.04, pad=0.1)
    barlbl = var_name + " "
    cbar.set_label(barlbl, fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    
    return fig

      

                  
dt_path = glob.glob("C:/Users/Reza/Desktop/Deniz_data/l3/S5P_OFFL_L3__AER_AI_*.nc")

ds = xr.open_mfdataset(dt_path, combine="nested", compat='override')
#print(ds["absorbing_aerosol_index"].values.shape)

"""
da = ds.to_array('absorbing_aerosol_index')

aai = da.to_dataframe(dim_order=None)

aai['absorbing_aerosol_index'] = pd.to_datetime(aai['absorbing_aerosol_index'])
annual_aai = aai.resample('W', on='absorbing_aerosol_index').mean()
#aai_2019 = aai.resample('M', on='absorbing_aerosol_index').mean()
aai_final = annual_aai.to_xarray()
aai_final.to_netcdf('yearly_combined.nc')

latmin=30.
latmax=60.
lonmin=35.
lonmax=55.

var_name = 'absorbing_aerosol_index'
file_name = 'yearly_combined.nc'
ds.to_array(dim='absorbing_aerosol_index', name=None)
#new_file = xr.DataArray('yearly_combined.nc')
visualize_pcolormesh(data_array=ds, 
                         longitude=ds.longitude, 
                         latitude=ds.latitude,
                         projection=ccrs.PlateCarree(), color_scale='YlGn',
                         vmin=-5, vmax=5, var_name='absorbing_aerosol_index',
                         lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax, 
                         set_global=False)
out_file = "{0}_{1}.png".format(str(var_name).replace(" ", "_"), file_name)
plt.savefig(out_file)      
plt.close()
"""