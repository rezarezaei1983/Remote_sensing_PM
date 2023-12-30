


# Ref : https://mygeoblog.com/2019/06/25/mask-netcdf-using-shp-file/
#       https://gis.stackexchange.com/questions/286634/gdal-raster-is-rotated-flipped-incorrectly

#       chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides1.pdf
#       https://notebook.community/ritviksahajpal/open-geo-tutorial/Python/chapters/chapter_4_vector

import os
import netCDF4
import numpy as np
from osgeo import gdal,osr,ogr
import matplotlib.pyplot as plt
import pprint



# I got this error:   ERROR 1: PROJ: proj_as_wkt: Cannot find proj.db
# To solve the error I added the below line:
os.environ['PROJ_LIB'] = r"C:\\Users\\Reza\\anaconda3\\envs\\env39\\Library\\share\\proj"


def ulc_ccor(lon,lat):
    mesh_lons, mesh_lats = np.meshgrid(lon, lat)
    
    ulc_lon = mesh_lons[-1,0]   
    ulc_lat = mesh_lats[-1,0]   

    return ulc_lon, ulc_lat 


def makeMask(lons,lats,lon_size,lat_size,res):
    ulc_lon, ulc_lat = ulc_ccor(lons, lats)
        
    #source_ds = ogr.Open(shapefile)                 ###############        
    source_ds = gdal.OpenEx(shapefile)
    src_band = source_ds.GetRasterBand(1)#.ReadAsArray()
    
    ### get the layer
    source_layer = source_ds.GetLayer(0)   #GetLayer("TUR_adm0")
    #print("\nsource layer : \n", source_layer)
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
    x_rotation = 0 
    y_rotation = 0 
    cell_width = res
    cell_height = -res
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
    
    array = np.flipud(array)    ##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    #### https://gis.stackexchange.com/questions/227406/reading-a-shapefile-as-an-array-using-python
   
    ######################### Write Shapefile as TIFF #########################
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
    ###########################################################################
    """
    driver.Register()
    outband = driver.GetRasterBand(1).WriteArray(array) 
    """
    
    #print("\n\n----------------- GetGeoTransform -----------------")
    #print(driver.GetGeoTransform())
     
    # Flush memory file
    driver = None
    band = None
        
    return array
 
# set the data directories
datadir = "C:/Users/Reza/Desktop/Deniz_data/Deniz_data_2/"
shapefile = "C:/Users/Reza/Desktop/Deniz_data/shapefile/TUR_adm0.shp"
infile = "S5P_OFFL_L3__AER_AI_20180628T102407_20180628T120537_03661_01_010002_20180704T095226.nc"
ncs = datadir + infile



nc = netCDF4.Dataset(ncs,'r')
abs_aer_ind = nc.variables["absorbing_aerosol_index"][0,:,:]

### show the abs_aer_ind
plt.imshow(abs_aer_ind, origin='lower')  # , origin='lower'   #################################
plt.show()

lons = nc.variables['longitude'][:]
lats = nc.variables['latitude'][:]

lon_size = lons.size
lat_size = lats.size

### calculate the cellsize
cellsize = lons[:][1] - lons[:][0]
#y_cellsize = lats[:][1] - lats[:][0]

#print("\ncellsize : \n", cellsize)
#print("\ny_cellsize : \n", y_cellsize)


### create the mask
mask = makeMask(lons,lats,lon_size,lat_size,cellsize) 

### show the mask
plt.imshow(mask, origin='lower')
plt.show()

### mask the abs_aer_ind data
mskd_abs_aer_ind = np.ma.masked_where(mask==0,abs_aer_ind)

print("\nmskd_abs_aer_ind : \n", mskd_abs_aer_ind)


plt.imshow(mskd_abs_aer_ind, origin='lower')  # , extent =[0, lon_size, 0, lat_size]
plt.show()


""" 
# print some stats
print(np.min(mskd_abs_aer_ind), np.mean(mskd_abs_aer_ind), np.max(mskd_abs_aer_ind))
"""