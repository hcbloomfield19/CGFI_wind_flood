#
# This script calculates the National aggregate total precipitation from ERA5 for creating corrleation metrics.
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


# load in the latitudes and longitudes of the gridded ERA5 data for calculating masks.
# note you need to download this before you start using the script in the /Data/ folder.

data_dir = '/path_to_glofas_data/'
fname = 'ERA5_EU_1hr_i10mg_tp_2015_12.nc'
dataset = Dataset(data_dir + fname ,mode='r')
lat = dataset.variables['latitude'][:]
lon = dataset.variables['longitude'][:]
dataset.close()


# this will sort all countries except GB and Ireland where we need some special shapefiles (in Data folder)
# These are inbuilt files in python from natural_earth.

countries_shp = shpreader.natural_earth(resolution='10m',category='cultural',name='admin_0_countries')
country_shapely=[]

# now lets turn it into a grid of 1s and 0s for use with the data.
LONS,LATS = np.meshgrid(lon,lat)
x,y = LONS.flatten(), LATS.flatten()

points = np.vstack((x,y)).T
MASK_MATRIX = np.zeros((len(x),1))


# load in population data
# 5km data from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11/data-download
# Not included in the repository due to size, but very easy to download. 
dataset = Dataset('/home/fz21704/Data/population_data/gpw_v4_population_density_adjusted_rev11_2pt5_min.nc',mode='r')
pop_lons = dataset.variables['longitude'][4000:4800] # isolate over Europe (-49.5E:50.25E)
pop_lats = dataset.variables['latitude'][600:1300] # isolate 4.75N:80.75N (similar to the ERA5 gridded data.
pop_totals = dataset.variables['UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][4,600:1300,4000:4800] 
dataset.close()

pops = np.ma.getdata(pop_totals)
pops[pops <=0.] = 0

# now lets interpolate this onto the same grid as ERA5.
POP_LONS,POP_LATS = np.meshgrid(pop_lons,pop_lats)

interp_function = interp2d(pop_lons,pop_lats,pops)

# we need the flip here (have tested in plots that it's working - we need it)
# some quirk of the interp-function
Population_totals_interp = np.flip(interp_function(lon,lat),axis=0)


# set up number of days to work over

oct_mar_days = 31+30+31+31+29+31 # include leap day.


####################################################
#
# NOW LETS WORK THROUGH ALL COUNTRIES AND CALCULATE NATIONAL AGGREGATE PRECIPITATION
#
#
######################################################

COUNTRY_LIST = ['Great Britain','All Ireland','France','Germany','Spain','Portugal','Italy','Belgium','Luxembourg','Netherlands','Norway','Sweden','Denmark','Poland','Czechia','Austria','Hungary','Switzerland','Latvia','Lithuania','Belarus','Ukraine','Slovakia','Slovenia','Romania','Bulgaria','Croatia','Greece','Montenegro','Moldova','Finland', 'Bosnia and Her','Estonia','Poland','Serbia','North Macedoni','Albania'] # remember names need spaces in


for COUNTRY in COUNTRY_LIST:

    # initialise the mask matrix 
    MASK_MATRIX = np.zeros((len(x),1))
    # create an empty array to store country shapefiles
    country_shapely=[]
    print(COUNTRY)

    # set up special conditions for the GB and All Ireland data (we've done this as collections of countries across Islands to make more sense with the river flow data later)

    if COUNTRY == 'Great Britain':
        GB_data = gpd.read_file("/home/fz21704/Data/NUTS/NUTS_1/Great_Britain.shp")
        country_shapely.append(GB_data.geometry)
    elif COUNTRY == 'All Ireland':
        IE_data = gpd.read_file("/home/fz21704/Data/NUTS/NUTS_1/All_Ireland.shp")
        country_shapely.append(IE_data.geometry)


   # standard natural_earth shapefiles for the rest of Europe

    else:
        for country in shpreader.Reader(countries_shp).records():
            #print(country.attributes["NAME"])
            if country.attributes["NAME"][0:len(COUNTRY)] == COUNTRY:
                print('Found Country')
                country_shapely.append(country.geometry)

    for i in range(0,len(x)):
        my_point = shapely.geometry.Point(x[i],y[i])

        if COUNTRY in ['Great Britain','All Ireland']:
            if country_shapely[0][0].contains(my_point) == True:
                MASK_MATRIX[i,0] = 1.0
        else:
            if country_shapely[0].contains(my_point) == True:
                MASK_MATRIX[i,0] = 1.0


    MASK_MATRIX_RESHAPE = np.reshape(MASK_MATRIX,(len(lat),len(lon)))

 
 
    pop_mask = Population_totals_interp*MASK_MATRIX_RESHAPE


    # setup arrays to fill up
    array_of_totals_tp = np.zeros([41,oct_mar_days])
    array_of_totals_tp_pop = np.zeros([41,oct_mar_days])

    for yr in range(1980,2021):


        # initialise lists to store interim files

        print(yr)
        tp_list = []
        tp_list_pop = []

        for mon in [10,11,12,1,2,3]:
        
            #print(mon)
            fname = 'ERA5_EU_1hr_i10mg_tp_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
            dataset = Dataset(data_dir + fname ,mode='r')
            tp = dataset.variables['tp'][:]*1000 # convert to mm
            dataset.close()
            tp_reshape = np.reshape(tp,(int(len(tp)/24),24,len(lat),len(lon)))
            data = np.sum(tp_reshape,axis=1)
 
            masked_variable = np.zeros(np.shape(data))
            pop_masked_variable = np.zeros(np.shape(data))

            # mask the data
            for i in range(0,len(data)):
                masked_variable[i,:,:] = data[i,:,:]*MASK_MATRIX_RESHAPE
                pop_masked_variable[i,:,:] = data[i,:,:]*pop_mask

            tp_day = np.sum(np.sum(masked_variable,axis=2),axis=1) / np.sum(MASK_MATRIX_RESHAPE)
            tp_day_pop = np.nansum(np.nansum(pop_masked_variable,axis=2),axis=1) / np.nansum(pop_mask)

            tp_list.append(tp_day)
            tp_list_pop.append(tp_day_pop)

        tp_for_array = np.concatenate(tp_list)
        tp_for_array_pop = np.concatenate(tp_list_pop)

        array_of_totals_tp[yr-1980,0:len(tp_for_array)] = tp_for_array
        array_of_totals_tp_pop[yr-1980,0:len(tp_for_array)] = tp_for_array_pop


    # uncomment plot to check all is working well
    
    #fig = plt.figure()
    #ax = fig.add_subplot(1,2,1)
   
    #x_plot =plt.pcolormesh(LONS,LATS,masked_variable[0,:,:],vmin=0,vmax=18,cmap='Blues')

    #ax = fig.add_subplot(1,2,2)

    #x_plot = plt.pcolormesh(LONS,LATS,data[0,:,:],vmin=0,vmax=18,cmap='Blues')
    #plt.show()
    
    # save the datafiles as numpy files for easy processing later. This can be changed to your preferred file format. 

    #np.save(str(COUNTRY) + '_daily_sum_tp.npy',array_of_totals_tp)
    #np.save(str(COUNTRY) + '_daily_sum_tp_pop_weighted.npy',array_of_totals_tp_pop)


