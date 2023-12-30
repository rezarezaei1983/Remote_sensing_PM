
import calendar
import numpy as np
import pprint
import warnings
import re

warnings.filterwarnings('ignore')
"""
data = {
        "a":{"a1":10, "a2":11, "a3":12},
        "b":{"b1":20, "b2":21, "b3":22}
        }

#print(data["a"].keys())


if "a1" in data["a"].keys():
    print("ok")
else:
    print("no")
"""    

"""
data2 = {
         "a":[10,11,12,13,14,15,16,17,18,19],
         "b":[20,21,22,23,24,25,26,27,28,29]
         }


print(data2["a"][3])
"""

"""
year_dict = {"2018": {'D': {'date': [], 'data': {"a":[1,2,3], "b":[10,11,12]}},
                      'M': {'Jan': [],
                            'Feb': [],
                            'Mar': [],
                            'Apr': [],
                            'May': [],
                            'Jun': [],
                            'Jul': [],
                            'Aug': [],
                            'Sep': [],
                            'Oct': [],
                            'Nov': [],
                            'Dec': []}},
            "2019": {'D': {'date': [], 'data': {}},
                     'M': {'Jan': [],
                           'Feb': [],
                           'Mar': [],
                           'Apr': [],
                           'May': [],
                           'Jun': [],
                           'Jul': [],
                           'Aug': [],
                           'Sep': [],
                           'Oct': [],
                           'Nov': [],
                           'Dec': []}}}

print(year_dict["2018"]["D"]["data"].keys())
"""

"""
year_list = [2019, 2020, 2021]
temporal_subset = {}

for year in year_list:
    temporal_subset[year] = {"data":{}}
    
pprint.pprint(temporal_subset)
"""
"""
a = [[[10,11,12],
      [13,14,15],
      [16,17,np.nan]],
     
     [[20,21,22],
      [23,24,25],
      [26,27,np.nan]],
     
     [[30,31,32],
      [33,34,35],
      [36,37,38]],
     
     [[40,41,42],
      [43,44,45],
      [46,47,48]]]

arr = np.array(a)
#print(arr.shape)

#print("\n")          #  RuntimeWarning

#mean = arr.mean(axis=(0))
#print(mean)

#print("\n")

mean = np.nanmean(arr, axis=0)
#print("mean:\n", mean)

print("\n")

maxim = np.nanmax(arr, axis=(1,2))
#print("maxim:\n", maxim)

print("\n")

minim = np.nanmin(arr, axis=(1,2))
#print("minim:\n", minim)

print("\n")

quantile = np.nanquantile(arr, 0.5, axis=0)
#print("quantile:\n", quantile)
"""
"""
b = [[10,11,12],
     [13,14,15],
     [16,17,np.nan]]

arr = np.array(b)
maxim = np.nanmax(arr, axis=(0,1))
print("maxim:\n", maxim)

print("\n")

minim = np.nanmin(arr, axis=(0,1))
print("minim:\n", minim)
"""
"""
inp_str = "percentile_0.5"
#num = re.findall(r'\d+', inp_str) 
num = float(inp_str.split("_")[1])

print(type(num))
"""

"""
seasons = [month%12 // 3 + 1 for month in range(1, 13)]
#print(seasons)


seasons = [calendar.month_abbr[month%12 // 3 + 1] for month in range(1, 13)]
#print(seasons)

month = 3

season = month%12 // 3 + 1
#print(season)


season_dict = {"Spring":[3,4,5], "Summer":[6,7,8],
               "Autumn":[9,10,11], "Winter":[12,1,2]}

a = [1, 5, 11, 8]

for i in a:
    if i in season_dict["Winter"]:
        print(i)
"""

"""
a = [[[10,11,12],
      [13,14,15],
      [16,17,np.nan]],
     
     [[20,21,22],
      [23,24,25],
      [26,27,np.nan]],
     
     [[30,31,32],
      [33,34,35],
      [36,37,38]],
     
     [[40,41,42],
      [43,44,45],
      [46,47,48]]]


a = np.array(a)
a = np.where(a >= 30, a, 0)
    
print(a)
"""

import math

x = np.nan
print(math.isnan(x))
