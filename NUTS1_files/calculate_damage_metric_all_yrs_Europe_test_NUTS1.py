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
from scipy.interpolate import interp2d
import geopandas as gpd

# firstly load in the percentiles data. and setup the masks

#perc_data = np.load('/home/fz21704/code_folder/CGFI/data/inst_10m_gust_p90_95_98_1day_EU.npy')
perc_data = np.load('/home/fz21704/code_folder/CGFI/data/inst_10m_gust_p90_95_98_99_995_1day_EU.npy')

perc_98 = perc_data[4,:,:] # 3= 99, 4=99.5 






# load in a land mask for the appropriate data:

# load in the latitudes and longitudes for plotting the maps
data_dir = '/scratch/hydro1/fz21704/ERA5/'
fname = 'ERA5_EU_1hr_i10mg_tp_2015_12.nc'
dataset = Dataset(data_dir + fname ,mode='r')
lat = dataset.variables['latitude'][:]
lon = dataset.variables['longitude'][:]
data = dataset.variables['i10fg'][:]
dataset.close()

# now lets turn it into a grid of 1s and 0s for use with the data.
LONS,LATS = np.meshgrid(lon,lat)
x,y = LONS.flatten(), LATS.flatten()

points = np.vstack((x,y)).T
MASK_MATRIX = np.zeros((len(x),1))


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
# create interpolation function.
interp_function = interp2d(pop_lons,pop_lats,pops)

# we need the flip here so latitudes are increasing for the function.
Population_totals_interp = np.flip(interp_function(lon,lat),axis=0)




oct_mar_days = 31+30+31+31+29+31 # include leap day.


# now load in the NUTS1 data
dataset = Dataset('NUTS_1_masks.nc',mode='r')
NUTS1_keys = dataset.variables['NUTS keys'][:]
NUTS1_masks = dataset.variables['NUTS zones'][:]
dataset.close()






####################################################
#
#
#
#
######################################################

COUNTRY_LIST = NUTS1_keys


for COUNTRY in COUNTRY_LIST:

    count = np.where(NUTS1_keys == COUNTRY)
    counter = count[0][0]
    print(COUNTRY)
    MASK_MATRIX_RESHAPE = NUTS1_masks[counter,:,:]

    if np.sum(MASK_MATRIX_RESHAPE) == 0:
        array_of_totals_SSI_pop = np.zeros([41,oct_mar_days])
        array_of_totals_gust = np.zeros([41,oct_mar_days])
        array_of_totals_tp = np.zeros([41,oct_mar_days])
 

    else:

        '''
        fig = plt.figure(figsize=(8,6))
        ax = plt.subplot(111,projection=ccrs.PlateCarree())
        ax.coastlines(resolution='50m')
        ax.set_extent([-9,20,40,60])
        s = plt.pcolormesh(LONS,LATS,MASK_MATRIX_RESHAPE,cmap='Blues')
        cb =fig.colorbar(s)
        cb.set_label('Mask value',fontsize=14)
        ax.set_title('Land Mask',fontsize=14)
        ax.set_aspect('auto',adjustable=None)
        plt.show()


        fig = plt.figure(figsize=(8,6))
        ax = plt.subplot(111,projection=ccrs.PlateCarree())
        ax.coastlines(resolution='50m')
        ax.set_extent([-9,20,40,60])
        s = plt.pcolormesh(LONS,LATS,data[0,:,:],cmap='Blues',vmin=0,vmax=10)
        cb =fig.colorbar(s)
        cb.set_label('Flow (m$^{3}$s$^{-1}$)',fontsize=14)
        ax.set_title('Glofas example flow',fontsize=14)
        ax.set_aspect('auto',adjustable=None)
        plt.show()
    


        fig = plt.figure(figsize=(8,6))
        ax = plt.subplot(111,projection=ccrs.PlateCarree())
        ax.coastlines(resolution='50m')
        ax.set_extent([-9,20,40,60])
        s = plt.pcolormesh(LONS,LATS,data[0,:,:]*MASK_MATRIX_RESHAPE,cmap='Blues',vmin=0,vmax=10)
        cb =fig.colorbar(s)
        cb.set_label('Flow (m$^{3}$s$^{-1}$)',fontsize=14)
        ax.set_title('Glofas example flow',fontsize=14)
        ax.set_aspect('auto',adjustable=None)
        plt.show()
        '''


        pop_mask = Population_totals_interp*MASK_MATRIX_RESHAPE

        #fig = plt.figure(figsize=(8,6))
        #ax = plt.subplot(111,projection=ccrs.PlateCarree())
        #ax.coastlines(resolution='50m')
        ##ax.set_extent([-9,20,40,60])
        #s = plt.pcolormesh(LONS,LATS,pop_mask,cmap='Blues',vmin=0,vmax=10)
        #cb =fig.colorbar(s)
        #cb.set_label('Flow (m$^{3}$s$^{-1}$)',fontsize=14)
        #ax.set_title('Glofas example flow',fontsize=14)
        #ax.set_aspect('auto',adjustable=None)
        #plt.show()


        array_of_totals_SSI_pop = np.zeros([41,oct_mar_days])
        array_of_totals_gust = np.zeros([41,oct_mar_days])
        array_of_totals_tp = np.zeros([41,oct_mar_days])

        for yr in range(1980,2021):

            print(yr)
            SSI_list_pop = []
            gust_list = []
            tp_list = []

            for mon in [10,11,12,1,2,3]:
            
                #print(mon)
                fname = 'ERA5_EU_1hr_i10mg_tp_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
                dataset = Dataset(data_dir + fname ,mode='r')
                gust = dataset.variables['i10fg'][:] # instantaneous
                tp = dataset.variables['tp'][:]*1000 # get into mm
                dataset.close()
                gust_reshape = np.reshape(gust,(int(len(gust)/24),24,len(lat),len(lon)))
                tp_reshape = np.reshape(tp,(int(len(tp)/24),24,len(lat),len(lon)))

                data = np.max(gust_reshape,axis=1)
                data_tp = np.sum(tp_reshape,axis=1)
 
                masked_variable = np.zeros(np.shape(data))
                masked_variable_tp = np.zeros(np.shape(data_tp))
                ratio_of_p98 = np.zeros(np.shape(data))
                SSI_pop = np.zeros(np.shape(data))
        

                #print('Masking variables')
                for i in range(0,len(data_tp)):
                    masked_variable[i,:,:] = data[i,:,:]*MASK_MATRIX_RESHAPE
                    masked_variable_tp[i,:,:] = data_tp[i,:,:]*MASK_MATRIX_RESHAPE
                    ratio_of_p98[i,:,:] = (masked_variable[i,:,:]/perc_98) - 1

                
                #print('calculating SSI')
                for i_day in range(0,len(data)):
                    #print(i_day)
                    for i_lat in range(0,len(lat)):
                        for i_lon in range(0,len(lon)):
                            # if the variable in question < P98 of the variable:
                            if ratio_of_p98[i_day,i_lat,i_lon] > 0.:
                                cubic_data = ratio_of_p98[i_day, i_lat, i_lon]**3
                                SSI_pop[i_day, i_lat , i_lon] = cubic_data*pop_mask[i_lat,i_lon]
                        
                

                SSI_month_pop = np.sum(np.sum(SSI_pop,axis=2),axis=1)
                gust_day = np.sum(np.sum(masked_variable,axis=2),axis=1) / np.sum(MASK_MATRIX_RESHAPE)
                tp_day = np.sum(np.sum(masked_variable_tp,axis=2),axis=1) / np.sum(MASK_MATRIX_RESHAPE)


                gust_list.append(gust_day)
                tp_list.append(tp_day)
                SSI_list_pop.append(SSI_month_pop)

            SSI_for_array_pop = np.concatenate(SSI_list_pop)
            gust_for_array = np.concatenate(gust_list)
            tp_for_array = np.concatenate(tp_list)

            array_of_totals_SSI_pop[yr-1980,0:len(gust_for_array)] = SSI_for_array_pop
            array_of_totals_gust[yr-1980,0:len(gust_for_array)] = gust_for_array
            array_of_totals_tp[yr-1980,0:len(tp_for_array)] = tp_for_array



    np.save('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/NUTS1_data/' + COUNTRY +'_daily_sum_tp.npy',array_of_totals_tp)
    np.save('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/NUTS1_data/' + COUNTRY +'_daily_max_gust.npy',array_of_totals_gust)
    np.save('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/NUTS1_data_sens_test/' + COUNTRY + '_daily_SSI995_pop_weighted.npy',array_of_totals_SSI_pop)
    
