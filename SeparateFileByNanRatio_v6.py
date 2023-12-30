
import os
import sys
import glob
import netCDF4
import numpy as np
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
            print(print("\n\033[1;31;47mError :", f"\033[30m'{self.datadir}' not found.\n"))
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
        
        ### get the layer
        source_layer = source_ds.GetLayer(0)   
        """
        ### get layer count
        layer_count = source_ds.GetLayerCount()
        print("\nlayer_count : \n", layer_count)
        
        ### get the layer name
        layer = source_ds.GetLayerByIndex(0)
        print('The layer is named: {n}\n'.format(n=layer.GetName()))
           
        ### get feature number
        numFeatures = source_layer.GetFeatureCount()
        print("\nnumFeatures : ", numFeatures)
        
        ### get a feature
        feature = source_layer.GetFeature(0)
        print("\nfeature : \n", feature)
        
        ### get attriute
        attribute = feature.GetFieldAsString('ISO')
        print("\nattribute : \n", attribute)
        
        ### get reference system of the shapefile
        spatialRef = source_layer.GetSpatialRef()  
        print("\nReference System : \n", spatialRef)
        """
        number_of_bands = 1
        driver = gdal.GetDriverByName("MEM").Create(                           
            '', lon_size, lat_size, number_of_bands, gdal.GDT_Byte)
            
        top_left_x = ulc_lon
        top_left_y = ulc_lat
        x_rotation, y_rotation = 0, 0
        #y_rotation = 0 
        cell_width, cell_height = res, -res
        #cell_height = -res
        geotransform = [top_left_x, cell_width, x_rotation, top_left_y, y_rotation, cell_height ]
        driver.SetGeoTransform(geotransform)  
             
        # Rasterize shapefile to grid
        gdal.RasterizeLayer(driver, [1], source_layer, burn_values=[1], options=["ALL_TOUCHED=TRUE"])
            
        ### ****************************
        #target_dsSRS = osr.SpatialReference()
        #target_dsSRS.ImportFromEPSG(4326)
        #driver.SetProjection(target_dsSRS.ExportToWkt())
        ### ****************************
        
        band = driver.GetRasterBand(1)
        
        # Get rasterized shapefile as numpy array
        array = band.ReadAsArray()
        array = np.flipud(array)    
        """
        ####################### Write Shapefile as TIFF #######################
        proj = driver.GetProjection()
        
        #print("\nproj : \n", proj)
        
        driver2 = gdal.GetDriverByName("GTiff")  
        
        driver2.Register()
        
        outds = driver2.Create(
            "C:/Users/Reza/Desktop/Deniz_data/testraster222.tif", 
            lon_size, lat_size, number_of_bands, gdal.GDT_Byte) 
        
        #outds.SetGeoTransform(geotransform)  
        #outds.SetProjection(proj)
        outband = outds.GetRasterBand(1).WriteArray(array) 
        #######################################################################
        """
        # Flush memory file
        driver = None
        band = None
           
        return array
    
    def CalculationOnMaskedArr(self):
        L3_files = self.open_nc("S5P_OFFL_L3__AER_AI_")
    
        for file in L3_files:
            cur_date = file.split("_")[6]
            print("\n\n---------------------------\n", file)
            print(cur_date)
            
            nc = netCDF4.Dataset(file,'r')
            abs_aer_ind = nc.variables["absorbing_aerosol_index"][0,:,:]
            lons = nc.variables['longitude'][:]
            lats = nc.variables['latitude'][:]

            ### show the abs_aer_ind
            plt.imshow(abs_aer_ind, origin='lower')  # , origin='lower'   #################################
            plt.title(cur_date)
            plt.show()
            
            lon_size = lons.size
            lat_size = lats.size
            cellsize = lons[:][1] - lons[:][0]
            
            ### create the mask
            mask = self.makeMask(lons,lats,lon_size,lat_size,cellsize) 
            
            ### show the mask array
            #plt.imshow(mask, origin='lower')
            #plt.show()
            
            ### mask the abs_aer_ind data
            mskd_abs_aer_ind = np.ma.masked_where(mask==0,abs_aer_ind)
            #print("\nmskd_abs_aer_ind : \n", mskd_abs_aer_ind)
            
            ### count the total number of grids inside the non-masked region of array (including NaNs and non-NaNs)
            non_mskd_grids = mskd_abs_aer_ind.count()
        
            ### count the number on nans inside the non-masked region of array
            non_mskd_nans = np.count_nonzero(np.isnan(mskd_abs_aer_ind[~mskd_abs_aer_ind.mask]))            
            
            nan_pernentage = (100 * non_mskd_nans) / non_mskd_grids
            print("\nnan_pernentage : ", nan_pernentage)
            
            if nan_pernentage > self.max_nans:
                print("\n>>> Should be removed")
            elif nan_pernentage < self.max_nans:
                print("\n>>> Should be used.")
            
            plt.imshow(mskd_abs_aer_ind, origin='lower')  # , extent =[0, lon_size, 0, lat_size]
            plt.title(cur_date)
            plt.show()




#        %%%%%%%%%%%%%%%%%%%%%      Instance      %%%%%%%%%%%%%%%%%%%%%


ins = GetMaskedImage(data_dir="C:/Users/Reza/Desktop/Deniz_data/Deniz_data_2/",
                     shapefile_dir="C:/Users/Reza/Desktop/Deniz_data/shapefile/TUR_adm0.shp",
                     max_nans_percent = 15)        

ins.CalculationOnMaskedArr()