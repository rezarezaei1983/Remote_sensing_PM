"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              <<< SOME NOTES >>>                               #
#                                                                               #
#>>> This script prepares the spatial and temporal subsets of the S5P OFFL (L3) # 
#    data and visualizes the mean, maximum, or minimum of the temporal subsets. #
#                                                                               #
#>>> To define the subset domain, set the latitude and longitude boundaries     #
#    (lat_bound, lon_bound).                                                    #
#                                                                               #
#>>> To eliminate the input files with NaNs percentages higher than a specific  #
#    value, set the acceptable maximum percentage of NaNs in the                #
#    'max_nans_percent' argument. Otherwise, set 'max_nans_percent'=100 to use  #
#    all files without considering the NaNs percentage.                         #
#                                                                               #
#>>> The script is only applicable for daily("d"), monthly("m"), seasonal("s"), #
#    and yearly("y") visualization.                                             # 
#                                                                               #
#>>> In seasonal subsetting, every year's winter season includes December from  #
#    that year as well as January and February from the following year.         #
#                                                                               #
#>>> The script can visualize the mean of the period (day, month, season, and   #
#    year), or the maximum or minimum of the period. Select the desired         #
#    statistical measure option ("mean", "max", "min") in the 'stat_measure'    #
#    argument.                                                                  #
#                                                                               #
#>>> Set the 'threshold' argument (give a number) to visualize the grids        #
#    containing values higher than a specific value (threshold), only.          #
#    Otherwise set the 'threshold'=None to visualize all grids (do not apply a  #
#    threshold value). The threshold value is applicable to mean, maximum, and  #
#    minimum measures.                                                          #
#                                                                               #
#>>> The outputs will be saved in "{input_dir}/Plots". Also, a file named       #
#    "day_number_of_periods.csv" could be found in the same path, which shows   #
#    the number of files that were used to create each period measure ("mean",  # 
#    "max", "min").                                                             #
#                                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@author : Reza Rezaei
email   : rezarezaei2008@gmail.com
version : 1.2
year    : 2022
"""

import os
import sys
import glob
import math
import pprint
import string
import warnings
import numpy as np
import xarray as xr
import pandas as pd
import cartopy.crs as ccrs
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import cartopy.feature as cfeature
from cartopy.mpl.geoaxes import GeoAxes
GeoAxes._pcolormesh_patched = Axes.pcolormesh
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

warnings.simplefilter(action = "ignore", category = UserWarning)
warnings.filterwarnings('ignore')


class handleS5Pl3(): 
    def __init__(self, input_dir, lat_bound, lon_bound, max_nans_percent, 
                 period=["d","m","s","y"], stat_measure=["mean","max","min"], 
                 threshold=None):
        self.indir = input_dir
        self.var = "absorbing_aerosol_index"
        self.lat_bnd = lat_bound
        self.lon_bnd = lon_bound
        self.max_nans = max_nans_percent
        self.period = period
        self.measure = stat_measure
        self.threshold = threshold
        
    def makeSpatialSubset(self):
        """ This function make the geographical subset of the main data using 
            the given boundary coordinates. """
        if os.path.exists(self.indir):
            files = os.path.join(self.indir, "S5P_OFFL_L3__AER_AI_*.nc")
            infiles = glob.glob(files)
        else:
            print("Error: Input files directory not found.")
            sys.exit()
        
        df_dict = {"values": [], "times": [], "lat": None, "lon": None, "NaNs%": []}
        for count, file in enumerate(infiles):
            cur_file = xr.open_dataset(file)
            lats = cur_file["latitude"][:].values
            lons = cur_file["longitude"][:].values
            idx_lat=(lats>=self.lat_bnd[0])*(lats<=self.lat_bnd[1])
            idx_lon=(lons>=self.lon_bnd[0])*(lons<=self.lon_bnd[1])
            idxlat=np.nonzero(idx_lat)[0]
            idxlon=np.nonzero(idx_lon)[0]             
            s5p_aai = cur_file[self.var][0, idxlat[0]:idxlat[-1], idxlon[0]:idxlon[-1]]
            aai_vals = s5p_aai.values
            total_grids = aai_vals.size
            nans = np.count_nonzero(np.isnan(s5p_aai))
            nan_pernentage = (100 * nans) / total_grids
            
            if count == 0:
                lat_arr = cur_file["latitude"][idxlat[0]:idxlat[-1]].values
                lon_arr = cur_file["longitude"][idxlon[0]:idxlon[-1]].values
                df_dict["lat"] = lat_arr
                df_dict["lon"] = lon_arr
                self.name = string.capwords(s5p_aai.description)    
                
            if nan_pernentage <= self.max_nans:
                df_dict["NaNs%"].append(nan_pernentage)
                df_dict["values"].append(aai_vals)
                date = str(file).split("_")[7].replace("T", " ") 
                datetime_str = date[:11] + ":" + date[11:13]+ ":" + date[13:]
                datetime_obj = datetime.strptime(datetime_str, '%Y%m%d %H:%M:%S')
                df_dict["times"].append(datetime_obj)
            cur_file.close()
        
        return df_dict
            
    def setSubsetDict(self):
        """ Prepares dictionary to store daily, monthly, or yearly data in the 
            next function 'append2SubsetDict()'. """
        df_dict = self.makeSpatialSubset()
        years = list(np.unique([time.year for time in df_dict["times"]]))
        
        subset_dict = {}
        for year in years:
            subset_dict[year] = {"data":{}, "lat":None, "lon":None, 
                                 "date":[], "file no":[]}
                
        if self.period.casefold() == "D".casefold():               
            for time in df_dict["times"]:
                subset_dict[time.year]["data"][str(time.date())]=[]
        elif self.period.casefold() == "M".casefold():
            for time in df_dict["times"]:
                subset_dict[time.year]["data"][time.strftime("%b")]=[]
        elif self.period.casefold() == "Y".casefold():
            for time in df_dict["times"]:
                subset_dict[time.year]["data"][time.year]=[]
        elif self.period.casefold() == "S".casefold():
            seasons = ["Spring", "Summer", "Autumn", "Winter"]
            for year in years:
                for season in seasons:
                    subset_dict[year]["data"][season]=[]                    
        else:
            print("\nError: Only 'D' (Daily), 'M' (Monthly), 'S' (seasonal), " \
                  "and 'Y' (Yearly) values are acceptable as inputs for the " \
                  "variable named 'period'.")
            sys.exit()  
        
        return subset_dict
    
    def append2SubsetDict(self):    
        """ This function appends (stores) the data to the relevent day, month,  
            or year dictionary, without averaging the data. """
        df_dict = self.makeSpatialSubset()
        temporal_subset = self.setSubsetDict()                             
        
        if self.period.casefold() == "D".casefold() or \
            self.period.casefold() == "M".casefold() or \
                self.period.casefold() == "Y".casefold():
            for count, dt in enumerate(df_dict["times"]): 
                cur_year = dt.year
                if self.period.casefold() == "D".casefold():
                    cur_date = str(dt.date())
                elif self.period.casefold() == "M".casefold():
                    cur_date = dt.strftime("%b")
                elif self.period.casefold() == "Y".casefold():
                    cur_date = dt.year
                temporal_subset[cur_year]["data"][cur_date].append(df_dict["values"][count])
                temporal_subset[cur_year]["lat"] = df_dict["lat"]
                temporal_subset[cur_year]["lon"] = df_dict["lon"]
                if dt.date() not in temporal_subset[cur_year]["date"]:
                    temporal_subset[cur_year]["date"].append(dt.date())
        
        elif self.period.casefold() == "S".casefold():   
            season_dict = {"Spring":[3,4,5], "Summer":[6,7,8],
                           "Autumn":[9,10,11], "Winter":[12,1,2]}
            for count, dt in enumerate(df_dict["times"]):
                cur_year = dt.year
                if self.period.casefold() == "S".casefold():
                    if dt.month in season_dict["Spring"]:
                        cur_date = "Spring"
                        temporal_subset[cur_year]["data"][cur_date].append(df_dict["values"][count])
                        temporal_subset[cur_year]["lat"] = df_dict["lat"]
                        temporal_subset[cur_year]["lon"] = df_dict["lon"]
                        if dt.date() not in temporal_subset[cur_year]["date"]:
                            temporal_subset[cur_year]["date"].append(dt.date())                             
                    elif dt.month in season_dict["Summer"]:
                        cur_date = "Summer"
                        temporal_subset[cur_year]["data"][cur_date].append(df_dict["values"][count])
                        temporal_subset[cur_year]["lat"] = df_dict["lat"]
                        temporal_subset[cur_year]["lon"] = df_dict["lon"]
                        if dt.date() not in temporal_subset[cur_year]["date"]:
                            temporal_subset[cur_year]["date"].append(dt.date())     
                    elif dt.month in season_dict["Autumn"]:
                        cur_date = "Autumn"
                        temporal_subset[cur_year]["data"][cur_date].append(df_dict["values"][count])
                        temporal_subset[cur_year]["lat"] = df_dict["lat"]
                        temporal_subset[cur_year]["lon"] = df_dict["lon"]
                        if dt.date() not in temporal_subset[cur_year]["date"]:
                            temporal_subset[cur_year]["date"].append(dt.date())
                    elif dt.month in season_dict["Winter"]:
                        cur_date = "Winter"
                        if dt.month == 12:
                            temporal_subset[cur_year]["data"][cur_date].append(df_dict["values"][count])
                            temporal_subset[cur_year]["lat"] = df_dict["lat"]
                            temporal_subset[cur_year]["lon"] = df_dict["lon"]
                            if dt.date() not in temporal_subset[cur_year]["date"]:
                                temporal_subset[cur_year]["date"].append(dt.date())
                        elif dt.month == 1 or dt.month == 2:
                            if cur_year - 1 in temporal_subset.keys():
                                temporal_subset[cur_year - 1]["data"][cur_date].append(df_dict["values"][count])
                                temporal_subset[cur_year - 1]["lat"] = df_dict["lat"]
                                temporal_subset[cur_year - 1]["lon"] = df_dict["lon"]
                                if dt.date() not in temporal_subset[cur_year - 1]["date"]:
                                    temporal_subset[cur_year - 1]["date"].append(dt.date())

        for year in temporal_subset.keys():
            for date in temporal_subset[year]["data"].keys():
                file_no = len(temporal_subset[year]["data"][date])
                temporal_subset[year]["file no"].append(file_no)
        
        del df_dict
        return temporal_subset
    
    def getVminVmax(self):  
        """ This function finds the max and in value in all time-steps of the 
            data and adds a margine in both limits. The max and min values will
            be used as the upper and lower limits of the colorbar of the map. """
        temporal_subset = self.append2SubsetDict()
        
        mins = []
        maxs = []  
        for year in temporal_subset.keys():
            for date, value in temporal_subset[year]["data"].items():
                value = np.array(value)
                if len(value) > 0:
                    min_list = np.nanmin(value, axis=(1,2))
                    max_list = np.nanmax(value, axis=(1,2))
                    for i in range(len(value)):
                        mins.append(min_list[i])
                        maxs.append(max_list[i])
        vmin = math.floor(np.nanmin(np.array(mins)))   
        vmax = math.ceil(np.nanmax(np.array(maxs)))       

        return vmin, vmax
        
    def calculatePeriodicStatMeasure(self):    
        """ This function calculates the periodic (daily, monthly, seasonal, or   
            yearly) statistical measure (mean, max, min) of the values. """   
        vmin, vmax = self.getVminVmax()                    
        periodic_val = self.append2SubsetDict()                     
                
        for year in periodic_val.keys():
            for date, value in periodic_val[year]["data"].items():
                if len(value) > 0:
                    arr = np.array(value)                                     
                    if self.measure.casefold() == "mean".casefold():
                        val = np.nanmean(arr, axis=0)    
                        if self.threshold != None:
                            val = np.where(val >= self.threshold, val, vmin)
                    elif self.measure.casefold() == "max".casefold():
                        val = np.nanmax(arr, axis=0) 
                        if self.threshold != None:
                            val = np.where(val >= self.threshold, val, vmin)
                    elif self.measure.casefold() == "min".casefold():
                        val = np.nanmin(arr, axis=0)   
                        if self.threshold != None:
                            val = np.where(val >= self.threshold, val, vmin)
                    periodic_val[year]["data"][date] = val
        
        return periodic_val
           
    def visualizeOnMap(self, arr, lon, lat, vmin, vmax, measure, period, date, out_dir):      
        fig=plt.figure(figsize=(20, 15))
        projection = ccrs.PlateCarree()
        ax = plt.axes(projection=projection)
   
        img = plt.pcolormesh(lon, lat, arr, 
                            cmap=plt.get_cmap("jet"),
                            transform=ccrs.PlateCarree(),
                            vmin=vmin,
                            vmax=vmax,
                            shading="auto")
        
        ax.add_feature(cfeature.BORDERS, edgecolor="black", linewidth=1)
        ax.add_feature(cfeature.COASTLINE, edgecolor="black", linewidth=1)
        ax.set_extent([self.lon_bnd[0], self.lon_bnd[1], self.lat_bnd[0], 
                       self.lat_bnd[1]], ccrs.PlateCarree())
        
        gl = ax.gridlines(draw_labels=True, linewidth=1, color="black",
                          alpha=0.8, linestyle='--')
        gl.top_labels=False
        gl.right_labels=False
        gl.xformatter=LONGITUDE_FORMATTER
        gl.yformatter=LATITUDE_FORMATTER
        gl.xlabel_style={"size":14}
        gl.ylabel_style={"size":14}
        plot_title = f"{period} {measure} {self.name}\n({date})"
        plt.title(plot_title, fontsize=20)
        #ax.set_global()
        #ax.gridlines()
        
        cbar = fig.colorbar(img, ax=ax, orientation='horizontal', fraction=0.04, pad=0.1)
        cbar.set_label(self.name, fontsize=14)
        cbar.ax.tick_params(labelsize=14)
        
        out_file_name = f"{out_dir}/{period}_{measure}_{date}.png"
        plt.savefig(out_file_name, dpi=100, bbox_inches="tight") 
        plt.close()       
        
    def runVisualization(self):
        vmin, vmax = self.getVminVmax()
        periodic_val = self.calculatePeriodicStatMeasure()
        
        path = os.path.join(self.indir, "Plots")
        if not os.path.exists(path):
            os.makedirs(path)
        
        measures = {"mean": "Mean", "max": "Maximum", "min": "Minimum"}
        periods = {"d": "Daily", "m": "Monthly", "s": "Seasonal", "y": "Yearly"}
        report = {"Period":[], "Number of used files":[]}
        
        for year in periodic_val.keys():
            counter=0
            for period, value in periodic_val[year]["data"].items():
                if len(value) > 1:
                    lats = periodic_val[year]["lat"]
                    lons = periodic_val[year]["lon"]
                    msr = measures[self.measure.lower()]
                    prd = periods[self.period.lower()]
                    if self.period.casefold() == "m".casefold() or \
                        self.period.casefold() == "s".casefold():
                        cur_date = f"{period} {year}"
                    elif self.period.casefold() == "d".casefold():
                        cur_date = period
                    elif self.period.casefold() == "y".casefold():
                        cur_date = year
                    
                    report["Period"].append(cur_date)
                    file_no = periodic_val[year]["file no"][counter]
                    report["Number of used files"].append(file_no)
                    
                    self.visualizeOnMap(arr=value, lon=lons, lat=lats, vmin=vmin, 
                                        vmax=vmax, measure=msr, period=prd, 
                                        date=cur_date, out_dir=path)
                counter += 1
                
        df = pd.DataFrame.from_dict(report)
        out_file_name = f"{path}/day_number_of_periods.csv"
        df.to_csv(out_file_name , index=False)
        del periodic_val
        
        print(f"\nPlots directory:\n    {path}\n")
        
        msg =   """
                       + + + + + + + + + + + + + + + + + + +
                       +                                   +
                       +         SUCCESSFULLY DONE         +
                       +                                   +
                       + + + + + + + + + + + + + + + + + + +
                """
        print(msg)



#================================== Instance ==================================
df_dir = "C:/Users/Reza/Desktop/data/l3/"               
                 
ins = handleS5Pl3(input_dir=df_dir,
          lat_bound=[35, 43], 
          lon_bound=[25, 46],    
          max_nans_percent=25,
          period="m",              # options: d (daily), m (monthly), s (seasonal), y (yearly)
          stat_measure="mean",      # options: mean, max, min
          threshold=None)          # options: None or an integer

ins.runVisualization()
