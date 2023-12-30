"""
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                              <<< SOME NOTES >>>                             #
#                                                                             #
# If action = "get_inappropriate_file_name" :     The names of inappropriate  #
#            file will be saved in a file named "inappropriate_files_list.csv"#
#                                                                             #
# If action = "remove_inappropriate_file" :       All of nappropriate files   #  
#             will be deleted.                                                #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
"""

import os
import glob
import math
from netCDF4 import Dataset
import numpy as np
import csv

class data_separation():
    def __init__(self, files_path, llc, ulc, urc, lrc, 
                 action=["get_inappropriate_file_name", "remove_inappropriate_file"]):
        self.path = files_path
        self.llc = llc
        self.ulc = ulc
        self.urc = urc
        self.lrc = lrc
        self.action = action
    
    def get_problematic_data(self):
        inappropriate_files_list = []
        
        nc_files = glob.glob(self.path + "*.nc")  

        for count, item in enumerate(nc_files):
            ds = Dataset(item, "r")
            llc_abs_ind =  ds.variables["absorbing_aerosol_index"][0,self.llc[0],self.llc[1]]
            ulc_abs_ind =  ds.variables["absorbing_aerosol_index"][0,self.ulc[0],self.ulc[1]]
            urc_abs_ind =  ds.variables["absorbing_aerosol_index"][0,self.urc[0],self.urc[1]]
            lrc_abs_ind =  ds.variables["absorbing_aerosol_index"][0,self.lrc[0],self.lrc[1]]
            """
            print("llc_abs_ind : ", llc_abs_ind)
            print("ulc_abs_ind : ", ulc_abs_ind)
            print("urc_abs_ind : ", urc_abs_ind)
            print("lrc_abs_ind : ", lrc_abs_ind)
            print("\n--------------------------------------------------------")
            """
            if math.isnan(llc_abs_ind) or math.isnan(ulc_abs_ind) or \
                math.isnan(urc_abs_ind) or math.isnan(lrc_abs_ind):
                    inappropriate_files_list.append(str(item))
                    ds.close()
        
        return inappropriate_files_list
        
    def take_action(self):
        inappropriate_files_list = self.get_problematic_data()
        
        if self.action == "get_inappropriate_file_name":
            resultFile = open("inappropriate_files_list.csv", 'a')
            wr = csv.writer(resultFile, dialect='excel')
            i = 0
            while i!=len(inappropriate_files_list): 
                wr.writerow(inappropriate_files_list[i:i+1])   
                i+=1
            resultFile.close()
            
        elif self.action == "remove_inappropriate_file":
            for file in inappropriate_files_list:
                os.remove(file)

        
# The index of non NaN values (in corners)                         
llc_ind = [0, 0]         #lat_lon = [25.005, 20.005] 
ulc_ind = [1999, 0]      #lat_lon = [44.995, 20.005]
urc_ind = [1999, 2012]   #lat_lon = [44.995, 41.025]
lrc_ind = [0, 2334]      #lat_lon = [25.005, 43.345]


ins = data_separation(files_path="C:/Users/Reza/Desktop/Deniz_data/l3/", 
                        llc=llc_ind, ulc=ulc_ind, urc=urc_ind, lrc=lrc_ind,
                        action="remove_inappropriate_file")

ins.take_action()
