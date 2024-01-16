"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                             <<< SOME NOTES >>>                              #
#                                                                             #
#>>> This script separates Sentinel-5P files based on the given maximum Nan   #
#    percentage ('max_nans_percent').                                         #
#>>> L2 and L3 files must be in the same directory.                           #
#>>> L2 files are global, while L3 files are cropped data (regional).         #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@author : Reza Rezaei
email   : rezarezaei2008@gmail.com
version : 0.7
year    : 2022
"""

import os
import sys
import glob
import netCDF4
import numpy as np
from datetime import datetime
from osgeo import gdal,osr,ogr
import matplotlib.pyplot as plt
import pprint

# I got this error:   ERROR 1: PROJ: proj_as_wkt: Cannot find proj.db
# To solve the error I added the below line:
os.environ['PROJ_LIB'] = r"C:\\Users\\Reza\\anaconda3\\envs\\env39\\Library\\share\\proj"


class GetMaskedImage():
    def __init__(self, data_dir, shapefile_dir, max_nans_percent):
        self.datadir = data_dir
        self.shapedir = shapefile_dir
        self.max_nans = max_nans_percent
        self.accepted = "accepted_data"
        self.rejected = "rejected_data"
    
    def open_nc(self, nc_type):
        if os.path.exists(self.datadir):
            os.chdir(self.datadir)
            files = glob.glob(nc_type + "*.nc")     
        elif not os.path.exists(self.datadir):
            print("\n\033[1;31;47mError :", f"\033[30m'{self.datadir}' not found.\n")
            sys.exit()
            
        return files
        
    def ulc_coor(self, lon, lat):
        """ Get the Upper Left Corner coordinate. """
        mesh_lons, mesh_lats = np.meshgrid(lon, lat)
        ulc_lon = mesh_lons[-1,0]
        ulc_lat = mesh_lats[-1,0]
        
        return ulc_lon, ulc_lat 
    
    def makeMask(self, lons,lats,lon_size,lat_size,res):
        ulc_lon, ulc_lat = self.ulc_coor(lons, lats)
        source_ds = gdal.OpenEx(self.shapedir)
        src_band = source_ds.GetRasterBand(1)
        
        source_layer = source_ds.GetLayer(0)   
        number_of_bands = 1
        driver = gdal.GetDriverByName("MEM").Create(                           
            '', lon_size, lat_size, number_of_bands, gdal.GDT_Byte)
            
        top_left_x = ulc_lon
        top_left_y = ulc_lat
        x_rotation, y_rotation = 0, 0
        cell_width, cell_height = res, -res
        geotransform = [top_left_x, cell_width, x_rotation, top_left_y, y_rotation, cell_height ]
        driver.SetGeoTransform(geotransform)  
        gdal.RasterizeLayer(driver, [1], source_layer, burn_values=[1], options=["ALL_TOUCHED=TRUE"])        
        band = driver.GetRasterBand(1)
        array = band.ReadAsArray()
        array = np.flipud(array)    
        driver = None
        band = None
           
        return array
    
    def CalculationOnMaskedArr(self, file, var_name):
        cur_date = file.split("_")[6]
        cur_date = cur_date[:4] + "." + cur_date[4:6] + "." + cur_date[6:]
        cur_date = cur_date[:-4] + ":" + cur_date[-4:-2] + ":" + cur_date[-2:]
        cur_date = datetime.strptime(cur_date, '%Y.%m.%d %H:%M:%S')
            
        nc = netCDF4.Dataset(file,'r')
        abs_aer_ind = nc.variables[var_name][0,:,:]
        lons = nc.variables['longitude'][:]
        lats = nc.variables['latitude'][:]

        plt.imshow(abs_aer_ind, origin='lower')
        plt.title(cur_date)
        plt.show()
            
        lon_size = lons.size
        lat_size = lats.size
        cellsize = lons[:][1] - lons[:][0]
        
        mask = self.makeMask(lons,lats,lon_size,lat_size,cellsize) 
        mskd_abs_aer_ind = np.ma.masked_where(mask==0,abs_aer_ind)
        non_mskd_grids = mskd_abs_aer_ind.count()
        non_mskd_nans = np.count_nonzero(np.isnan(mskd_abs_aer_ind[~mskd_abs_aer_ind.mask]))            
        nan_pernentage = (100 * non_mskd_nans) / non_mskd_grids
        print("\nnan_pernentage : ", nan_pernentage)
            
        if nan_pernentage > self.max_nans:
            print("\n>>> Should be removed")
        elif nan_pernentage < self.max_nans:
            print("\n>>> Should be used.")
            
        plt.imshow(mskd_abs_aer_ind, origin='lower')
        plt.title(cur_date)
        plt.show()
            
    def RunCalculation(self):
        L3_files = self.open_nc("S5P_OFFL_L3__AER_AI_")
        
        for file in L3_files:
            L3 = self.CalculationOnMaskedArr(file=file, var_name="absorbing_aerosol_index")
        




#    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      Run      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ins = GetMaskedImage(data_dir="C:/Users/Reza/Desktop/data/",
                     shapefile_dir="C:/Users/Reza/Desktop/shapefile/TUR_adm0.shp",
                     max_nans_percent = 15)        
ins.RunCalculation()
