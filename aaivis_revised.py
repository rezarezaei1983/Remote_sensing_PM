
import string
from datetime import datetime
import glob

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

#%run ../functions.ipynb

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


def visualize_pcolormesh(data_array, longitude, latitude, projection, color_scale, 
                         vmin, vmax, unit, var_name, set_global=True, lonmin=-180, 
                         lonmax=180, latmin=-90, latmax=90):
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
        lonmin,lonmax,latmin,latmax(float): optional kwarg, set geographic extent 
        is set_global kwarg is set to False
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
    barlbl = var_name + " " + unit
    cbar.set_label(barlbl, fontsize=14)
    cbar.ax.tick_params(labelsize=14)
    
    return fig

        
      

nc_path = "C:\\Users\\Reza\\Desktop\\Deniz_data\\l3\\"
output_path = "C:\\Users\\Reza\\Desktop\\Deniz_data\\l3\\"

nc_files = glob.glob(nc_path + "*.nc")

for count, item in enumerate(nc_files):
    file = xr.open_dataset(item)
    name = str(item).split("_")[7].replace("T", "_") 
    print("Name: ", name)
    file_name = name[:11] + "-" + name[11:13]+ "-" + name[13:]
    print(file_name)

    s5p_aai = file['absorbing_aerosol_index']
    s5p_aai_time = s5p_aai[0,:,:]
    unit = s5p_aai_time.units
    long_name = s5p_aai.description
    var_name = string.capwords(long_name)
    
    latmin=20.
    latmax=50.
    lonmin=25.
    lonmax=45.
    
    s5p_aai_subset = generate_geographical_subset(xarray=s5p_aai_time, 
                                                  latmin=latmin, latmax=latmax, 
                                                  lonmin=lonmin, lonmax=lonmax)

    visualize_pcolormesh(data_array=s5p_aai_subset, 
                         longitude=s5p_aai_subset.longitude, 
                         latitude=s5p_aai_subset.latitude,
                         projection=ccrs.PlateCarree(), color_scale='YlGn',
                         vmin=0, vmax=0.2, unit=unit, var_name=var_name, 
                         lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax, 
                         set_global=False)
    
    ############################## set title 
    date = str(item).split("_")[7].replace("T", " ")
    date = date[2:]
    date = datetime.strptime(date, '%y%m%d %H%M%S')
    plt.title("{0} ({1})".format(var_name, date), fontsize=16)
    ##############################
    
    out_file = "{0}_{1}.png".format(str(long_name).replace(" ", "_"), file_name)
    plt.savefig(output_path + out_file)      
    plt.close()