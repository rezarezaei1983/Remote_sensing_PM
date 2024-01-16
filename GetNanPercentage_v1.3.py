"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                             <<< SOME NOTES >>>                              #
#                                                                             #
#>>> This script determines the NaNs percentage of all tile grids or grids    # 
#    within the Turkey shapefile boundaries, and separates the data based on  #
#    the NaN percentage.                                                      #
#    Also, the script calculates some statistical parameters such as max, min,#
#    quantiles, and Cumulative Probability.                                   #
#                                                                             #
#>>> L2 and L3 files must be in the same directory.                           #
#                                                                             #
#>>> After the installation of 'osgeo', edit 'PROJ_LIB' path in line 40, if   #
#    necessary.                                                               #
#                                                                             #
#>>> Edit paths in lines 247-249, and set the maximum acceptable NaNs         #
#    percentage {max_nans_percent}.                                           #
#                                                                             #
#>>> The outputs will be written in a CSV file named 'nan_percentage.csv' in  #
#    the 'data_dir/' directory.                                               #
#                                                                             #
#>>> Set 'use_shapefile_subset=True' if the calculations must be done based on# 
#    the Turkey shapefile. Otherwise, set 'use_shapefile_subset=False'.       #
#                                                                             #
#>>> Input data (tile data or masked based on Turkey shapefile), which are    #
#    used to determine NaNs%, will be plotted if 'create_map=True'. The plots #
#    will be saved in the 'shape_plots/' directory. Set 'create_map=False' to #
#    deactivate the data plotting.                                            #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@author : Reza Rezaei
email   : rezarezaei2008@gmail.com
version : 1.3
year    : 2022
"""

import os
import sys
import glob
import shutil
import netCDF4
import numpy as np
import pandas as pd
from osgeo import gdal 
from datetime import datetime
import matplotlib.pyplot as plt


#os.environ['PROJ_LIB'] = r"C:\\Users\\Reza\\anaconda3\\envs\\env39\\Library\\share\\proj"

class GetMaskedImage():
    def __init__(self, data_dir, use_shapefile_subset, shapefile_dir, 
                 max_nans_percent, create_map):
        self.datadir = data_dir
        self.subset = use_shapefile_subset
        self.shapedir = shapefile_dir
        self.max_nans = max_nans_percent
        self.map = create_map
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
        os.makedirs(self.accepted)
        
        if os.path.exists(self.rejected):
            files = glob.glob(self.datadir + self.rejected + "/*.nc")
            for file in files:
                os.replace(file, self.datadir + "/" + os.path.basename(file))
            shutil.rmtree(self.datadir + self.rejected )
        os.makedirs(self.rejected)
        
        if os.path.exists(self.plotdir):     
            shutil.rmtree(self.plotdir)    
        os.makedirs(self.plotdir) 
        
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
    
    def CreateMap(self, data, title, file_name):
        plt.imshow(data, origin='lower')
        plt.title(title)   
        figname = file_name[:-3]
        outfile = os.path.join(self.datadir, self.plotdir, figname)  
        plt.savefig(outfile + ".png") 
        plt.close()
    
    def CalculationOnMaskedArr(self, file, var_name):
        cur_date = file.split("_")[6].replace("T", " ")
        cur_date = cur_date[:4] + "." + cur_date[4:6] + "." + cur_date[6:]
        cur_date = cur_date[:-4] + ":" + cur_date[-4:-2] + ":" + cur_date[-2:]     
        cur_date = datetime.strptime(cur_date, '%Y.%m.%d %H:%M:%S')
                
        nc = netCDF4.Dataset(file,'r')
        abs_aer_ind = nc.variables[var_name][0,:,:]
        lons = nc.variables['longitude'][:]
        lats = nc.variables['latitude'][:]
        lon_size = lons.size
        lat_size = lats.size
        cellsize = lons[:][1] - lons[:][0]
        
        if self.subset == True:
            mask = self.makeMask(lons,lats,lon_size,lat_size,cellsize)
            mskd_abs_aer_ind = np.ma.masked_where(mask==0,abs_aer_ind)
            non_mskd_grids = mskd_abs_aer_ind.count()
            non_mskd_nans = np.count_nonzero(np.isnan(
                mskd_abs_aer_ind[~mskd_abs_aer_ind.mask]))            
            nan_pernentage = (100 * non_mskd_nans) / non_mskd_grids
            nan_pernentage = round(nan_pernentage, 2)
            data = mskd_abs_aer_ind[mskd_abs_aer_ind.mask == False]
        elif self.subset == False:
            grids = abs_aer_ind.count()
            non_nans = np.count_nonzero(~np.isnan(abs_aer_ind))
            nans = grids - non_nans
            nan_pernentage = (100 * nans) / grids
            nan_pernentage = round(nan_pernentage, 2)
            data = abs_aer_ind
                               
        data = np.copy(data)
        data = data[~np.isnan(data)]
        
        if len(data) > 0:        
            max_val = np.max(data)            
            min_val = np.min(data)            
            q_50 = np.quantile(data, 0.50)            
            q_75 = np.quantile(data, 0.75)            
            q_90 = np.quantile(data, 0.90)
            q_99 = np.quantile(data, 0.99)
            q_999 = np.quantile(data, 0.999)
        else:
            max_val = "-"            
            min_val = "-"          
            q_50 = "-"           
            q_75 = "-"          
            q_90 = "-"
            q_99 = "-"
            q_999 = "-"
        
        nc.close()    
        
        if self.map == True:
            if self.subset == True:
                self.CreateMap(data=mskd_abs_aer_ind, title=cur_date, file_name=file)
            elif self.subset == False:
                self.CreateMap(data=abs_aer_ind, title=cur_date, file_name=file)
                
        return cur_date, nan_pernentage, max_val, min_val, q_50, q_75, q_90, q_99, q_999
            
    def MakeCSV(self):
        L3_files = self.open_nc("S5P_OFFL_L3__AER_AI_")
        L2_files = self.open_nc("S5P_OFFL_L2__AER_AI_")          
        var_name= "absorbing_aerosol_index"
        
        if self.subset == True:
            msg1 = "NOTE : NaNs percentages were calculated by using the grid "\
                   "values placed within the Turkey shapefile boundaries."
        elif self.subset == False:
            msg1 = "NOTE : NaNs percentages were calculated using all grids "\
                   "within the tile."
            
        msg2 = "NOTE : NaNs are ignored in calculating the quantiles."
        print("\n", msg1, "\n", msg2)
        
        col_names = ["Year", "Month", "Day", "Hour", "Minute", "Second", "NaNs%", 
                     "Maximum", "Minimum", "0.5th quantile", "0.75th quantile",
                     "0.90th quantile", "0.99th quantile", "0.999th quantile"]
        analysis_results = pd.DataFrame(columns=col_names)
        
        for idx, file in enumerate(L3_files):
            cur_date, nan_pernentage, max_val, min_val, q_50, q_75, q_90, q_99, q_999 = \
                                        self.CalculationOnMaskedArr(file=file, 
                                        var_name=var_name)
            
            analysis_results.loc[idx, "Year"] = cur_date.year
            analysis_results.loc[idx, "Month"] = cur_date.month
            analysis_results.loc[idx, "Day"] = cur_date.day
            analysis_results.loc[idx, "Hour"] = cur_date.hour
            analysis_results.loc[idx, "Minute"] = cur_date.minute
            analysis_results.loc[idx, "Second"] = cur_date.second
            analysis_results.loc[idx, "NaNs%"] = nan_pernentage
            analysis_results.loc[idx, "Maximum"] = max_val
            analysis_results.loc[idx, "Minimum"] = min_val
            analysis_results.loc[idx, "0.5th quantile"] = q_50
            analysis_results.loc[idx, "0.75th quantile"] = q_75
            analysis_results.loc[idx, "0.90th quantile"] = q_90
            analysis_results.loc[idx, "0.99th quantile"] = q_99
            analysis_results.loc[idx, "0.999th quantile"] = q_999
            
            cur_L2_file = "S5P_OFFL_L2__AER_AI_" + file.split("_")[6]
            L2_file = next((L2 for L2 in L2_files if cur_L2_file in L2), None)

            if nan_pernentage > self.max_nans:
                os.replace(os.getcwd() + "/" + file, self.rejected + "/" + file) 
                os.replace(os.getcwd() + "/" + L2_file, self.rejected + "/" + L2_file)
                
            else:
                os.replace(os.getcwd() + "/" + file, self.accepted + "/" + file) 
                os.replace(os.getcwd() + "/" + L2_file, self.accepted + "/" + L2_file)
        
        analysis_results = analysis_results.sort_values(by=["NaNs%"])
        
        rank = np.arange(1, (len(analysis_results["Year"]) + 1), 1)
        
        cum_prob = []
        for i in rank:
            cp = (i - 0.5)/np.max(rank)
            cum_prob.append(cp)

        analysis_results["Rank"] = rank
        analysis_results["Cumulative Probability"] = cum_prob
        
        outfile = os.path.join(self.datadir, "nan_percentage.csv")
        with open(outfile, "w", newline="") as f:
            f.write("\n" + "Variable : " + var_name + "\n" + msg1 + "\n" + msg2 + "\n\n")
            analysis_results.to_csv(f, index=False)
            f.close()
        
        print(f"\n The data were filtered based on max NaNs = {self.max_nans}")
        print("\n The results are saved in : ", outfile)
        

        
#####################################  Run  ###################################

ins = GetMaskedImage(data_dir="C:/Users/Reza/Desktop/data/",
                     use_shapefile_subset=False,
                     shapefile_dir="C:/Users/Reza/Desktop/data/shapefile/TUR_adm0.shp",
                     max_nans_percent=100,
                     create_map=True)

ins.MakeCSV()
