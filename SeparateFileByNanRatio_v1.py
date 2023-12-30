

from osgeo import ogr
import geopandas as gpd


class SeparateSuitableFiles():
    def __init__(self, shape_path):
        self.shp = shape_path
        
    
    def read_shape(self):        
        #shp_file = ogr.Open(self.shp)
        #shape = shp_file.GetLayer(0)
        shapefile = gpd.read_file(self.shp)
        
        
        
        
        
    


ins = SeparateSuitableFiles(
    shape_path="C:/Users/Reza/Desktop/Deniz_data/shapefile/TUR_adm0.shp")

ins.read_shape()