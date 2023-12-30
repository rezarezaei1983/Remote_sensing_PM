"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                             <<< SOME NOTES >>>                              #
#                                                                             #
#>>> This script determines the NaNs percentage of the data within the Turkey # 
#    shapefile boundaries and separates the data based on the NaN percentage. #
#                                                                             #
#>>> L2 and L3 files must be in the same directory.                           #
#                                                                             #
#>>> After the installation of 'osgeo', edit 'PROJ_LIB' path in line 35.      #
#                                                                             #
#>>> Edit paths in lines 191-192, and set the maximum acceptable NaNs         #
#    percentage {max_nans_percent}.                                           #
#                                                                             #
#>>> NaNs% will be written in a CSV file named 'nan_percentage.csv' in the    #
#    'data_dir/' directory.                                                   #
#                                                                             #
#>>> The masked input data, which are used to determine NaNs% will be plotted #
#    and saved in the 'shape_plots/' directory, to check the results          #
#    visoually. Comment lines 133-138 to deactivate data plotting.            #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
"""

import os
import sys
import glob
import shutil
import netCDF4
import numpy as np
import pandas as pd
from osgeo import gdal 
import matplotlib.pyplot as plt

#os.environ['PROJ_LIB'] = r"C:\\Users\\Reza\\anaconda3\\envs\\env39\\Library\\share\\proj"

class GetMaskedImage():
    def __init__(self, data_dir, shapefile_dir, max_nans_percent):
        self.datadir = data_dir
        self.shapedir = shapefile_dir
        self.max_nans = max_nans_percent
        self.plotdir = "shape_plots"               
        self.accepted = "accepted_data"
        self.rejected = "rejected_data"
    
    def open_nc(self, nc_type):
        if os.path.exists(self.datadir):
            os.chdir(self.datadir)
        elif not os.path.exists(self.datadir):
            print("\n\033[1;31;47mError :", 
                  f"\033[30m'{self.datadir}' not found.\n")
            sys.exit()
        
        if os.path.exists(self.accepted):
            files = glob.glob(self.datadir + self.accepted + "/*.nc")
            for file in files:
                os.replace(file, self.datadir + "/" + os.path.basename(file))
            shutil.rmtree(self.datadir + self.accepted )
        
        if os.path.exists(self.rejected):
            files = glob.glob(self.datadir + self.rejected + "/*.nc")
            for file in files:
                os.replace(file, self.datadir + "/" + os.path.basename(file))
            shutil.rmtree(self.datadir + self.rejected )
        
        if os.path.exists(self.plotdir):      
            shutil.rmtree(self.plotdir)
        
        files = glob.glob(nc_type + "*.nc") 
                       
        return files
        
    def ulc_coor(self, lon, lat):
        mesh_lons, mesh_lats = np.meshgrid(lon, lat)
        ulc_lon = mesh_lons[-1,0]
        ulc_lat = mesh_lats[-1,0]
        
        return ulc_lon, ulc_lat 
    
    def makeMask(self, lons,lats,lon_size,lat_size,res):
        ulc_lon, ulc_lat = self.ulc_coor(lons, lats)
        source_ds = gdal.OpenEx(self.shapedir)
        
        source_layer = source_ds.GetLayer(0)   
        number_of_bands = 1
        driver = gdal.GetDriverByName("MEM").Create(                           
            '', lon_size, lat_size, number_of_bands, gdal.GDT_Byte)
            
        top_left_x = ulc_lon
        top_left_y = ulc_lat
        x_rotation, y_rotation = 0, 0
        cell_width, cell_height = res, -res
        geotransform = [top_left_x, cell_width, x_rotation, 
                        top_left_y, y_rotation, cell_height ]
        driver.SetGeoTransform(geotransform)  
             
        gdal.RasterizeLayer(driver, [1], source_layer, burn_values=[1], 
                            options=["ALL_TOUCHED=TRUE"])
        band = driver.GetRasterBand(1)
        array = band.ReadAsArray()
        array = np.flipud(array)    
        
        driver = None
        band = None
           
        return array
    
    def CalculationOnMaskedArr(self, file, var_name):
        cur_date = file.split("_")[6].replace("T", " ")
        cur_date = cur_date[:4] + "." + cur_date[4:6] + "." + cur_date[6:]
        cur_date = cur_date[:-4] + ":" + cur_date[-4:-2] + ":" + cur_date[-2:]     
         
        nc = netCDF4.Dataset(file,'r')
        abs_aer_ind = nc.variables[var_name][0,:,:]
        lons = nc.variables['longitude'][:]
        lats = nc.variables['latitude'][:]

        lon_size = lons.size
        lat_size = lats.size
        cellsize = lons[:][1] - lons[:][0]
            
        mask = self.makeMask(lons,lats,lon_size,lat_size,cellsize)
        mskd_abs_aer_ind = np.ma.masked_where(mask==0,abs_aer_ind)
        non_mskd_grids = mskd_abs_aer_ind.count()
        non_mskd_nans = np.count_nonzero(np.isnan(
            mskd_abs_aer_ind[~mskd_abs_aer_ind.mask]))            
            
        nan_pernentage = (100 * non_mskd_nans) / non_mskd_grids
        nan_pernentage = round(nan_pernentage, 2)
        nc.close()    
        
        # comment this block to deactivate the data plotting
        plt.imshow(mskd_abs_aer_ind, origin='lower')
        plt.title(cur_date)   
        figname = file[:-3]
        outfile = os.path.join(self.datadir, self.plotdir, figname)  
        plt.savefig(outfile + ".png") 
        plt.close()
        
        return cur_date, nan_pernentage
            
    def MakeCSV(self):
        L3_files = self.open_nc("S5P_OFFL_L3__AER_AI_")
        L2_files = self.open_nc("S5P_OFFL_L2__AER_AI_")
        os.makedirs(self.plotdir)                 
        os.makedirs(self.accepted)
        os.makedirs(self.rejected)
        
        msg = "NOTE : NaNs percentages were calculated by using the grid "\
              "values placed within the Turkey shapefile boundaries."
        print("\n", msg)
        
        date_arr = []
        time_arr = []
        perc_arr = []
        for file in L3_files:
            cur_date, nan_pernentage = self.CalculationOnMaskedArr(file=file, 
                                            var_name="absorbing_aerosol_index")
            date_arr.append(cur_date.split(" ")[0])
            time_arr.append(cur_date.split(" ")[1])
            perc_arr.append(nan_pernentage)
            
            cur_L2_file = "S5P_OFFL_L2__AER_AI_" + file.split("_")[6]
            L2_file = next((L2 for L2 in L2_files if cur_L2_file in L2), None)

            if nan_pernentage > self.max_nans:
                os.replace(os.getcwd() + "/" + file, self.rejected + "/" + file) 
                os.replace(os.getcwd() + "/" + L2_file, self.rejected + "/" + L2_file)
                
            else:
                os.replace(os.getcwd() + "/" + file, self.accepted + "/" + file) 
                os.replace(os.getcwd() + "/" + L2_file, self.accepted + "/" + L2_file)
        
        data = {"Date" : date_arr, "Time" : time_arr, "NaNs%" : perc_arr}
        df = pd.DataFrame(data)
        outfile = os.path.join(self.datadir, "nan_percentage.csv")
        
        with open(outfile, "w", newline="") as f:
            f.write("\n" + msg + "\n\n")
            df.to_csv(f, index=False)
            f.close()
        
        print(f"\n The data were filtered based on max NaNs = {self.max_nans}")
        
        print("\n The results are saved in : ", outfile)
        

        
#####################################  Run  ###################################

ins = GetMaskedImage(data_dir="C:/Users/Reza/Desktop/Deniz_data/Deniz_data_2/",
                     shapefile_dir="C:/Users/Reza/Desktop/Deniz_data/shapefile/TUR_adm0.shp",
                     max_nans_percent=99)        

ins.MakeCSV()