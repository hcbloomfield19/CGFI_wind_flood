#
# This script calculates the damage metric from Priestley et al., (2018)
#
#



# import libraries
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt # not used except for on the fly map tests.
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.cm as cm
import cartopy.io.shapereader as shpreader
import shapely.geometry
import shapefile
import scipy.interpolate 
import geopandas as gpd

# firstly load in the percentiles data. and setup the masks

perc_data = np.load('/home/fz21704/code_folder/CGFI/data/glofas_EU_flows_p90_95_98_99_995_1day_annual.npy')
# shape = [ month, threshold, lats, lons]
perc_data = perc_data[4,:,:] # 2 for 98, 3 for 99, 4 for 99.5
# get rid of any odd datapoints
perc_data[perc_data<0.] = 0.





# load in the latitudes and longitudes for plotting the maps
data_dir = '/scratch/hydro1/fz21704/EU_glofas/'
fname = 'glofas_EU_daily_river_discharge_2015_12.nc'
dataset = Dataset(data_dir + fname ,mode='r')
lat = dataset.variables['lat'][:] # not a regular grid so take their meshgrid
lon = dataset.variables['lon'][:] -360 # not a regular grid so take their meshgrid
#print(dataset.variables.keys())
data = dataset.variables['dis24'][:] # m^3s^-1
dataset.close()

LONS,LATS = np.meshgrid(lon,lat)
x,y = LONS.flatten(), LATS.flatten()
points = np.vstack((x,y)).T


print('loading population data')


# load in population data
# 5km data from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11/data-download
dataset = Dataset('/home/fz21704/Data/population_data/gpw_v4_population_density_adjusted_rev11_2pt5_min.nc',mode='r')
pop_lons = dataset.variables['longitude'][4000:5500] # isolate over Europe (-49.5E:50.25E)
pop_lats = dataset.variables['latitude'][400:1500] # isolate 4.75N:80.75N (similar to the ERA5 gridded data.
pop_totals = dataset.variables['UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][4,400:1500,4000:5500] 
dataset.close()

pops = np.ma.getdata(pop_totals)
pops[pops <=0.] = 0

# now lets interpolate this onto the same grid as ERA5.
POP_LONS,POP_LATS = np.meshgrid(pop_lons,pop_lats)

print('doing interpolation')
points_array = np.array([POP_LONS.flatten(),POP_LATS.flatten()])
points_to_interp_array = np.array([LONS.flatten(),LATS.flatten()])

# do the interpolation
interp_nonstandard = scipy.interpolate.griddata(points_array.transpose(),pops.flatten(),points_to_interp_array.transpose())

Population_totals_interp = np.reshape(interp_nonstandard,(np.shape(LATS)))


# and now the masks for NUTS1 regions

dataset = Dataset('NUTS_1_masks_GLOFAS.nc',mode='r')
NUTS1_keys = dataset.variables['NUTS keys'][:]
NUTS1_masks = dataset.variables['NUTS zones'][:]
dataset.close


COUNTRY_LIST = NUTS1_keys


for COUNTRY in  COUNTRY_LIST:

    print(COUNTRY)

    count = np.where(NUTS1_keys == COUNTRY)
    counter = count[0][0]
 
    MASK_MATRIX_RESHAPE = NUTS1_masks[counter,:,:]

    if np.sum(MASK_MATRIX_RESHAPE) == 0:
        array_of_totals_flow = np.zeros([41,366])
        array_of_totals_FSI_pops = np.zeros([41,366])

    else:

        #fig = plt.figure(figsize=(8,6))
        #ax = plt.subplot(111,projection=ccrs.PlateCarree())
        #ax.coastlines(resolution='50m')
        #ax.set_extent([-9,20,40,60])
        #s = plt.pcolormesh(LONS,LATS,MASK_MATRIX_RESHAPE,cmap='Blues')
        #cb =fig.colorbar(s)
        #cb.set_label('Mask value',fontsize=14)
        #ax.set_title(COUNTRY,fontsize=14)
        #ax.set_aspect('auto',adjustable=None)
        #plt.show()

        pop_mask = Population_totals_interp*MASK_MATRIX_RESHAPE


        #fig = plt.figure(figsize=(8,6))
        #ax = plt.subplot(111,projection=ccrs.PlateCarree())
        #ax.coastlines(resolution='50m')
        #ax.set_extent([-9,20,40,60])
        #s = plt.pcolormesh(LONS,LATS,pop_mask,cmap='Blues',vmin=0,vmax=10)
        #cb =fig.colorbar(s)
        #cb.set_label('Flow (m$^{3}$s$^{-1}$)',fontsize=14)
        #ax.set_title('Glofas example flow',fontsize=14)
        #ax.set_aspect('auto',adjustable=None)
        #plt.show()


        
        array_of_totals_flow = np.zeros([41,366])
        array_of_totals_FSI_pops = np.zeros([41,366])

        for yr in range(1980,2021): 
            
            flow_list = []
            FSI_list_pops = []
            print(yr)

            for mon in range(1,13): # change back to 13
            
                perc_thresh = perc_data

                fname = 'glofas_EU_daily_river_discharge_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
                dataset = Dataset(data_dir + fname ,mode='r')
                data = dataset.variables['dis24'][:] # subset area for ease of computation
                dataset.close()
        
                masked_variable = np.zeros(np.shape(data))
            
                diff_of_p9X = np.zeros(np.shape(data))
                FSI_pops = np.zeros(np.shape(data))
            
                for i_day in range(0,len(data)):
                
                    masked_variable[i_day,:,:] = data[i_day,:,:]*MASK_MATRIX_RESHAPE
                    # sort out the masking and remove very large negative numbers in the masked region
                    masked_variable[masked_variable <0.] =0.
                    diff_of_p9X[i_day,:,:] = ((masked_variable[i_day,:,:] / perc_thresh) -1)
                    # set infinity to zero on mountainous region in spain where there are issues with the riverflows over mountain (3 gridpoints, Northern Spain)
                    diff_of_p9X[diff_of_p9X == np.inf] = 0.

         

                    # if a certin proportion of the land area is covered... (3140 GB land points in glofas) so 0.01$ of land is 3 gridpoints.... so need a value > 0 (ie exceeding percentile) for 3 gridpoints.
                    if len(np.where(diff_of_p9X[i_day,:,:] > 0.)[0]) > 0.: # for now have this as >0. as all countries are different sizes!
                        for i_lat in range(0,len(lat)):
                            for i_lon in range(0,len(lon)):
                                if diff_of_p9X[i_day,i_lat,i_lon] > 0.:
                                    # calculate the difference and normalise by the mean flow in that gridbox
                                    # if the variable in question < P9X of the variable:
                                    # calculate the difference and normalise by the mean flow in that gridbox
                                    FSI_pops[i_day,i_lat,i_lon] = (diff_of_p9X[i_day, i_lat, i_lon] )*pop_mask[i_lat,i_lon]
                        # if not enough is covered FSI is zero.
                    else:
                        FSI_pops[i_day,:,:] = np.zeros(np.shape(diff_of_p9X[i_day,:,:]))
     


                # sum over all the lats and lons
                flow_month = np.sum(np.nansum(masked_variable,axis=2),axis=1)
                flow_list.append(flow_month)

                FSI_month_pops= np.sum(np.nansum(FSI_pops,axis=2),axis=1)
                FSI_list_pops.append(FSI_month_pops)


            flow_for_array = np.concatenate(flow_list)
            array_of_totals_flow[yr-1980,0:len(flow_for_array)] = flow_for_array

            FSI_for_array_pops = np.concatenate(FSI_list_pops)
            array_of_totals_FSI_pops[yr-1980,0:len(FSI_for_array_pops)] = FSI_for_array_pops


    np.save('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/NUTS1_data/' + COUNTRY + '_daily_sum_glofas_riverflow_ALL_MON',array_of_totals_flow)
    np.save('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/NUTS1_data/' + COUNTRY + '_daily_sum_glofas_FSI_995_ALL_MON_gridthr_popweight',array_of_totals_FSI_pops)


