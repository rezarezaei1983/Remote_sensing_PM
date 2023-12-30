import numpy as np
import numpy.ma as ma

# Create an array with int elements using the numpy.array() method
arr = np.array([[65, 68, np.nan, np.nan, np.nan], 
                [93, 33, 39, 63, np.nan], 
                [73, 88, 59, 71, 51], 
                [np.nan, 45, 39, 41, 67],
                [np.nan, 41, 47, 34, 50]])

print("Array...\n", arr)

# Create a masked array and mask some of them as invalid
maskArr = ma.masked_array(arr, mask =[[0, 0, 0, 0, 0], 
                                      [0, 0, 0, 0, 0],
                                      [1, 0, 0, 0, 0], 
                                      [1, 1, 0, 0, 0], 
                                      [1, 1, 1, 1, 0]])
print("\nOur Masked Array\n", maskArr)


##############

##############


# Get the number of elements of the Masked Array
print("\nElements in the Masked Array...\n",maskArr.size)


# To count the non-masked elements of the masked array, use the ma.MaskedArray.count() method
print("\nNon-masked elements...\n", maskArr.count())

 
print("\nNumber of masked elements : ", ma.count_masked(maskArr))


print("\nNumber of NaNs : ", np.count_nonzero(np.isnan(maskArr)))  # ~

### get the non-masked part as a 1D array
#non_mskd_data = maskArr[~maskArr.mask]
print("\nNumber of NaNs inside the non-masked area : \n", 
      np.count_nonzero(np.isnan(maskArr[~maskArr.mask])))

