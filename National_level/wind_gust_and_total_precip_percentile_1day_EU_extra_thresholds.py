
#################
#
#
# This script will load in the whole of ERA5 and find a set of useful percentiles of 10m wind gust and total precipitation over the European domain.
# The 90th, 95th and 98th 99th 99.5th percentile
#
##################


# import libraries
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt # not used except for on the fly map tests.


# predefine an array.
octmar_days = 31+30+31+31+29+31 # include leap day.

# use ones instead of zeros due to precip having a lot of zeros anyway so harder
# to remove
array_of_totals_gust = np.ones([41,octmar_days,145,201])
array_of_totals_precip = np.ones([41,octmar_days,145,201])

data_dir =  '/path_to_data/'

for yr in range(1980,2021):
    print(yr)
    tp = []
    i10 = []
    for mon in [10,11,12,1,2,3]:

        # load in data
        fname = 'ERA5_EU_1hr_i10mg_tp_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
        dataset = Dataset(data_dir + fname ,mode='r')
        lat = dataset.variables['latitude'][:]
        lon = dataset.variables['longitude'][:]
        gust = dataset.variables['i10fg'][:] # instantaneous
        precip = dataset.variables['tp'][:]*1000 # convert to mm
        time = dataset.variables['time'][:]
        dataset.close()
        
        # get daily max
        gust_reshape = np.reshape(gust,(int(len(gust)/24),24,len(lat),len(lon)))
        gust_max = np.max(gust_reshape,axis=1)

        # get daily total
        precip_reshape = np.reshape(precip,(int(len(precip)/24),24,len(lat),len(lon)))
        precip_max = np.sum(precip_reshape,axis=1)

        i10.append(gust_max)
        tp.append(precip_max)


    i10_for_array = np.vstack(i10)
    tp_for_array = np.vstack(tp)

    array_of_totals_gust[yr-1980,0:len(i10_for_array),:,:] = i10_for_array
    array_of_totals_precip[yr-1980,0:len(i10_for_array),:,:] = tp_for_array



array_of_totals_gust[array_of_totals_gust == 1.] = np.nan
# commenting this out for now as a lot of the time it's not raining so i'm not diluting the sample too much here...think if there is a better way!
array_of_totals_precip[array_of_totals_precip == 1.] = np.nan


percentiles_precip = np.zeros([5,len(lat),len(lon)])
percentiles_gust = np.zeros([5,len(lat),len(lon)])


# calculate percentiles at each gridpoint.
for i in range(0,len(lat)):
    print(i)
    for j in range(0,len(lon)):

        # an earlier version has P99 in (so dont be too confused if cant find file source)
        percentiles_precip[:,i,j] = np.nanpercentile(array_of_totals_precip[:,:,i,j],(90,95,98,99,99.5))
        percentiles_gust[:,i,j] = np.nanpercentile(array_of_totals_gust[:,:,i,j],(90,95,98,99,99.5))
 

# save the outputs

#np.save('total_precip_p90_95_98_99_995_1day_EU',percentiles_precip)
#np.save('inst_10m_gust_p90_95_98_99_995_1day_EU',percentiles_gust)



