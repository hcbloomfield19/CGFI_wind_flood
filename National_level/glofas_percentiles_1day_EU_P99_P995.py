#################
#
#
# This script will load in the whole of glofas and find a set of useful percentiles of river runoff 
# over the European domain.
# The 90th, 95th 98th and 99th and 99.5th percentile
#
# Also saving all the percentiles for ease.
##################


# import libraries
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt # not used except for on the fly map tests.


# predefine an array.
days = 366 # include leap day.

# use ones instead of zeros due to precip having a lot of zeros anyway so harder
# to remove

array_of_totals_flow = np.ones([40,days,351,552])


# note you need to download this before you start using the script in the /Data/ folder.

data_dir =  '/path_to_glofas_data/'

for yr in range(1980,2020):
    print(yr)
    tflow = []
        
    for mon in range(1,13): 

        fname = 'glofas_EU_daily_river_discharge_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
        dataset = Dataset(data_dir + fname ,mode='r')
        lat = dataset.variables['lat'][:]
        lon = dataset.variables['lon'][:]
        discharge = dataset.variables['dis24'][:] # m^3 s^-1
        time = dataset.variables['time'][:]
        dataset.close()
        
        
        tflow.append(discharge)


    tflow_for_array = np.vstack(tflow)
    print(np.shape(tflow_for_array))

    array_of_totals_flow[yr-1980,0:len(tflow_for_array),:,:] = tflow_for_array



# set all ones to NaNs as there is no flow here. 
array_of_totals_flow[array_of_totals_flow == 1.] = np.nan



# now calculate the percentiles 

percentiles_flow = np.zeros([5,len(lat),len(lon)])

for i in range(0,len(lat)):
    print(i)
    for j in range(0,len(lon)):

        percentiles_flow[:,i,j] = np.nanpercentile(array_of_totals_flow[:,:,i,j],(90,95,98,99,99.5))



# save the data 

#np.save('glofas_EU_flows_p90_95_98_99_995_1day_annual',percentiles_flow)

