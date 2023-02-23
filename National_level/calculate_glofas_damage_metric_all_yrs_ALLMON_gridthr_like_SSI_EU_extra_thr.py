#
# This script calculates the Flood Severity Index metric from Bloomfield et al., (2023) designed to be similar to 
 # the Storm Severity Index metric from Priestley et al., (2018)
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

# firstly load in the percentiles data. and setup the masks.
# multiple percentile thresholds are available here if you wish to do sensitivity tests but it should be 
# the 99.5th percentile for the metric

perc_data = np.load('/home/code_for_github_reposutory_CLEANED/Data/glofas_EU_flows_p90_95_98_99_995_1day_annual.npy')

# shape = [ month, threshold, lats, lons]
perc_data = perc_data[4,:,:] # 2 for 98, 3 for 99, 4 for 99.5

# get rid of any odd datapoints
perc_data[perc_data<0.] = 0.

# load in the latitudes and longitudes of the gridded glofas data for calculating masks.
# note you need to download this before you start using the script in the /Data/ folder.
data_dir = '/path_to_glofas_data/'
fname = 'glofas_EU_daily_river_discharge_2015_12.nc'
dataset = Dataset(data_dir + fname ,mode='r')
lat = dataset.variables['lat'][:] # not a regular grid so take their meshgrid
lon = dataset.variables['lon'][:] -360 # not a regular grid so take their meshgrid
dataset.close()

# now lets turn it into a grid for making the country masks
LONS,LATS = np.meshgrid(lon,lat)
x,y = LONS.flatten(), LATS.flatten()
points = np.vstack((x,y)).T

# this will sort all countries except GB and Ireland where we need some special shapefiles (in Data folder)

# These are inbuilt files in python from natural_earth.

countries_shp = shpreader.natural_earth(resolution='10m',category='cultural',name='admin_0_countries')


# load in population data
# 5km data from https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-adjusted-to-2015-unwpp-country-totals-rev11/data-download
# Not included in the repository due to size, but very easy to download. 
dataset = Dataset('/home/fz21704/Data/population_data/gpw_v4_population_density_adjusted_rev11_2pt5_min.nc',mode='r')
pop_lons = dataset.variables['longitude'][4000:5500] # isolate over Europe (-49.5E:50.25E)
pop_lats = dataset.variables['latitude'][400:1500] # isolate 4.75N:80.75N (similar to the ERA5 gridded data.
pop_totals = dataset.variables['UN WPP-Adjusted Population Density, v4.11 (2000, 2005, 2010, 2015, 2020): 2.5 arc-minutes'][4,400:1500,4000:5500] 
dataset.close()

# set any negative values to zero as they mess up interpoltation.
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


####################################################
#
#
# NOW LETS WORK THROUGH ALL COUNTRIES AND CREATE THE FSI METRIC (this could be sped up a lot using a library like xarray, cfpython or iris for a bulk data load)
#
######################################################


COUNTRY_LIST = ['Great Britain','All Ireland','France','Germany','Spain','Portugal','Italy','Belgium','Luxembourg','Netherlands','Norway','Sweden','Denmark','Poland','Czechia','Austria','Hungary','Switzerland','Latvia','Lithuania','Belarus','Ukraine','Slovakia','Slovenia','Romania','Bulgaria','Croatia','Greece','Montenegro','Moldova','Finland', 'Bosnia and Her','Estonia','Poland','Serbia','North Macedoni','Albania'] 



for COUNTRY in COUNTRY_LIST:


    # create an empty array to store country shapefiles
    country_shapely=[]



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

    # initialise the mask matrix 
    MASK_MATRIX = np.zeros((len(x),1))

    # there are a lot of points so this bit may take some time!
    for i in range(0,len(x)):
        my_point = shapely.geometry.Point(x[i],y[i])

        if COUNTRY in ['Great Britain','All Ireland']:
            if country_shapely[0][0].contains(my_point) == True:
                MASK_MATRIX[i,0] = 1.0
        else:
            if country_shapely[0].contains(my_point) == True:
                MASK_MATRIX[i,0] = 1.0

 
    MASK_MATRIX_RESHAPE = np.reshape(MASK_MATRIX,(np.shape(LATS)[0],np.shape(LATS)[1]))


   
    # uncomment plots if you want to check everything is working!  
    
    #fig = plt.figure(figsize=(8,6))
    #ax = plt.subplot(111,projection=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m')
    #ax.set_extent([-9,20,40,60])
    #s = plt.pcolormesh(LONS,LATS,MASK_MATRIX_RESHAPE,cmap='Blues')
    #cb =fig.colorbar(s)
    #cb.set_label('Mask value',fontsize=14)
    #ax.set_title('Land Mask',fontsize=14)
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


    # initalise arrays to fill. Note our percentiles are only calculated from data from Jan-Dec 1980-2020. So we generate FSI for the full year.
    # This was to match some data provided to us from CEH which used annual percentiles of river flows (results are very similar if you restrict
    # to Oct-Mar percentiles like the SSI metric)

    array_of_totals_flow = np.zeros([41,366])
    array_of_totals_FSI = np.zeros([41,366])
    array_of_totals_FSI_pops = np.zeros([41,366])

    for yr in range(1980,2021): 


        # initialise lists to store interim files

        flow_list = []
        FSI_list = []
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
            FSI = np.zeros(np.shape(data))
            FSI_pops = np.zeros(np.shape(data))
            


            # Mask the data
            for i_day in range(0,len(data)):
                
                masked_variable[i_day,:,:] = data[i_day,:,:]*MASK_MATRIX_RESHAPE
                # sort out the masking and remove very large negative numbers in the masked region
                masked_variable[masked_variable <0.] =0.
                diff_of_p9X[i_day,:,:] = ((masked_variable[i_day,:,:] / perc_thresh) -1)
                # set infinity to zero on mountainous region in spain where there are issues with the riverflows over mountain (3 gridpoints, Northern Spain)
                diff_of_p9X[diff_of_p9X == np.inf] = 0.

         
                # Calculate the FSI

                # extra condition here about needeing 0.01% of land to be 'flooded' as this is a relatively fine grid compared to ERA5. 

                # if a certin proportion of the land area is covered... (3140 GB land points in glofas) so 0.01% of land is 3 gridpoints.... so need a value > 0 (ie exceeding percentile) for 3 gridpoints.
                if len(np.where(diff_of_p9X[i_day,:,:] > 0.)[0]) > 0.: # for now have this as >0. as all countries are different sizes!
                    for i_lat in range(0,len(lat)):
                        for i_lon in range(0,len(lon)):
                            if diff_of_p9X[i_day,i_lat,i_lon] > 0.:
                            # calculate the difference and normalise by the mean flow in that gridbox
                            # if the variable in question < P9X of the variable:
                            # calculate the difference and normalise by the mean flow in that gridbox
                                FSI[i_day,i_lat,i_lon] = (diff_of_p9X[i_day, i_lat, i_lon] )
                                FSI_pops[i_day,i_lat,i_lon] = (diff_of_p9X[i_day, i_lat, i_lon] )*pop_mask[i_lat,i_lon]
                # if not enough is covered FSI is zero.
                else:
                    FSI[i_day,:,:] = np.zeros(np.shape(diff_of_p9X[i_day,:,:]))
                    FSI_pops[i_day,:,:] = np.zeros(np.shape(diff_of_p9X[i_day,:,:]))
     


            # sum over all the lats and lons
            flow_month = np.sum(np.nansum(masked_variable,axis=2),axis=1)
            flow_list.append(flow_month)

            FSI_month = np.sum(np.nansum(FSI,axis=2),axis=1)
            FSI_list.append(FSI_month)

            FSI_month_pops= np.sum(np.nansum(FSI_pops,axis=2),axis=1)
            FSI_list_pops.append(FSI_month_pops)


        flow_for_array = np.concatenate(flow_list)
        #print(np.shape(flow_for_array))
        array_of_totals_flow[yr-1980,0:len(flow_for_array)] = flow_for_array

        FSI_for_array = np.concatenate(FSI_list)
        #print(np.shape(FSI_for_array))
        array_of_totals_FSI[yr-1980,0:len(FSI_for_array)] = FSI_for_array

        FSI_for_array_pops = np.concatenate(FSI_list_pops)
        #print(np.shape(FSI_for_array_pops))
        array_of_totals_FSI_pops[yr-1980,0:len(FSI_for_array_pops)] = FSI_for_array_pops



    # save the datafiles as numpy files for easy processing later. This can be changed to your preferred file format. 

    #np.save(COUNTRY + '_daily_sum_glofas_riverflow_ALL_MON',array_of_totals_flow)
    #np.save(COUNTRY + '_daily_sum_glofas_FSI_995_ALL_MON_gridthr',array_of_totals_FSI)
    #np.save(COUNTRY + '_NEW_daily_sum_glofas_FSI_99_ALL_MON_gridthr_popweight',array_of_totals_FSI_pops)



