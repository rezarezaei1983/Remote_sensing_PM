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

%run ../functions.ipynb

def generate_geographical_subset(xarray, latmin, latmax, lonmin, lonmax, reassign=False):
    """ 
    Generates a geographical subset of a xarray.DataArray and if kwarg reassign=True, shifts the longitude grid 
    from a 0-360 to a -180 to 180 deg grid.
    
    Parameters:
        xarray(xarray.DataArray): a xarray DataArray with latitude and longitude coordinates
        latmin, latmax, lonmin, lonmax(int): lat/lon boundaries of the geographical subset
        reassign(boolean): default is False
        
    Returns:
        Geographical subset of a xarray.DataArray.
    """   
    if(reassign):
        xarray = xarray.assign_coords(longitude=(((xarray.longitude + 180) % 360) - 180))
    return xarray.where((xarray.latitude < latmax) & (xarray.latitude > latmin) & (xarray.longitude < lonmax) & (xarray.longitude > lonmin),drop=True)


def visualize_pcolormesh(data_array, longitude, latitude, projection, color_scale, vmin, vmax, 
                         set_global=True, lonmin=-180, lonmax=180, latmin=-90, latmax=90):
    """ 
    Visualizes a xarray.DataArray with matplotlib's pcolormesh function.
    
    Parameters:
        data_array(xarray.DataArray): xarray.DataArray holding the data values
        longitude(xarray.DataArray): xarray.DataArray holding the longitude values
        latitude(xarray.DataArray): xarray.DataArray holding the latitude values
        projection(str): a projection provided by the cartopy library, e.g. ccrs.PlateCarree()
        color_scale(str): string taken from matplotlib's color ramp reference
        unit(str): the unit of the parameter, taken from the NetCDF file if possible
        long_name(str): long name of the parameter, taken from the NetCDF file if possible
        vmin(int): minimum number on visualisation legend
        vmax(int): maximum number on visualisation legend
        set_global(boolean): optional kwarg, default is True
        lonmin,lonmax,latmin,latmax(float): optional kwarg, set geographic extent is set_global kwarg is set to 
                                            False

    """
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
    cbar.set_label(unit, fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    

 #   plt.show()
    return fig

file = xr.open_mfdataset('/media/elviradeniz/thesis/data/L3_data/S5P_OFFL_L3__AER_AI_20180628T102407_20180628T120537_03661_01_010002_20180704T095226.nc')
file

s5p_aai = file['absorbing_aerosol_index']
s5p_aai

s5p_aai_time = s5p_aai[0,:,:]
s5p_aai_time


latitude = s5p_aai_time.latitude
print(latitude)
longitude = s5p_aai_time.longitude
print(longitude)

units = s5p_aai_time.units
print(units)
#long_name = s5p_aai_time.long_name
#print(long_name)

multiplication_factor=s5p_aai_time.multiplication_factor_to_convert_to_molecules_percm2
print(multiplication_factor)

latmin=20.
latmax=50.
lonmin=25.
lonmax=45.

s5p_aai_subset = generate_geographical_subset(xarray=s5p_aai_time,
                                             latmin=latmin,
                                             latmax=latmax,
                                             lonmin=lonmin,
                                             lonmax=lonmax)

s5p_aai_subset

visualize_pcolormesh(data_array=s5p_aai_subset,
                     longitude=s5p_aai_subset.longitude,
                     latitude=s5p_aai_subset.latitude,
                     projection=ccrs.PlateCarree(),
                     color_scale='viridis',
                     vmin=0, 
                     vmax=0.2,
                     lonmin=lonmin,
                     lonmax=lonmax,
                     latmin=latmin,
                     latmax=latmax,
                     set_global=False)

multiplication_factor = s5p_aai_subset.multiplication_factor_to_convert_to_molecules_percm2
print(multiplication_factor)

s5p_aai_subset_converted = s5p_aai_subset*multiplication_factor
print(s5p_aai_subset_converted)


visualize_pcolormesh(data_array=s5p_aai,
                     longitude=longitude,
                     latitude=latitude,
                     projection=ccrs.PlateCarree(),
                     color_scale='viridis',
                     vmin=0, 
                     vmax=10,
                     lonmin=lonmin,
                     lonmax=lonmax,
                     latmin=latmin,
                     latmax=latmax)