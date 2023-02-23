#
# This script calculates the Storm Severity Index metric from Priestley et al., (2018)
# Using ERA5 data - this will need to be downloaded from the climate data store to run the code below in advance. 
#



# import libraries
import numpy as np
import os
from netCDF4 import Dataset
import matplotlib.pyplot as plt 
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.cm as cm
import cartopy.io.shapereader as shpreader
import shapely.geometry
import shapefile
from scipy.interpolate import interp2d
import geopandas as gpd

# firstly load in the percentiles data. and setup the masks.
# multiple percentile thresholds are available here if you wish to do sensitivity tests but it should be 
# the 98th percentile for the metric

perc_data = np.load('/home/code_for_github_reposutory_CLEANED/Data/inst_10m_gust_p90_95_98_99_995_1day_EU.npy')
perc_98 = perc_data[2,:,:] # 3= 99, 4=99.5 



# load in the latitudes and longitudes of the gridded ERA5 data for calculating masks.
# note you need to download this before you start using the script in the /Data/ folder.

data_dir = '/path_to_ERA5_data/'
fname = 'ERA5_EU_1hr_i10mg_tp_2015_12.nc'
dataset = Dataset(data_dir + fname ,mode='r')
lat = dataset.variables['latitude'][:]
lon = dataset.variables['longitude'][:]
dataset.close()

# now lets turn it into a grid for making the country masks
LONS,LATS = np.meshgrid(lon,lat)
x,y = LONS.flatten(), LATS.flatten()

points = np.vstack((x,y)).T
MASK_MATRIX = np.zeros((len(x),1))

# this will sort all countries except GB and Ireland where we need some special shapefiles (in Data folder)

# These are inbuilt files in python from natural_earth.
countries_shp = shpreader.natural_earth(resolution='10m',category='cultural',name='admin_0_countries')
country_shapely=[]

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
# create interpolation function.
interp_function = interp2d(pop_lons,pop_lats,pops)

# we need the flip here so latitudes are increasing for the function.
Population_totals_interp = np.flip(interp_function(lon,lat),axis=0)

oct_mar_days = 31+30+31+31+29+31 # include leap day.





####################################################
#
#
# NOW LETS WORK THROUGH ALL COUNTRIES AND CREATE THE SSI METRIC (this could be sped up a lot using a library like xarray, cfpython or iris for a bulk data load)
#
######################################################

COUNTRY_LIST = ['Great Britain','All Ireland','France','Germany','Spain','Portugal','Italy','Belgium','Luxembourg','Netherlands','Norway','Sweden','Denmark','Poland','Czechia','Austria','Hungary','Switzerland','Latvia','Lithuania','Belarus','Ukraine','Slovakia','Slovenia','Romania','Bulgaria','Croatia','Greece','Montenegro','Moldova','Finland', 'Bosnia and Her','Estonia','Poland','Serbia','North Macedoni','Albania']



for COUNTRY in COUNTRY_LIST:


    # initialise the mask matrix
    MASK_MATRIX = np.zeros((len(x),1))
    # create an empty array to store country shapefiles
    country_shapely=[]

    print(COUNTRY)


    # set up special conditions for the GB and All Ireland data (we've done this as collections of countries across Islands to make more sense with the river flow data later)
    if COUNTRY == 'Great Britain':
        GB_data = gpd.read_file("/home/code_for_github_reposutory_CLEANED/Data/NUTS/NUTS_1/Great_Britain.shp")
        country_shapely.append(GB_data.geometry)
    elif COUNTRY == 'All Ireland':
        IE_data = gpd.read_file("/home/code_for_github_reposutory_CLEANED/Data/NUTS/NUTS_1/All_Ireland.shp")
        country_shapely.append(IE_data.geometry)

    # standard natural_earth shapefiles for the rest of Europe
    else:
        for country in shpreader.Reader(countries_shp).records():
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



    # you can uncomment the plots below if you want to test everything is working ok with your ERA5 data! 

    #plt.pcolormesh(MASK_MATRIX_RESHAPE)
    #plt.show()

    
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


    #fig = plt.figure(figsize=(8,6))
    #ax = plt.subplot(111,projection=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m')
    #ax.set_extent([-9,20,40,60])
    #s = plt.pcolormesh(LONS,LATS,data[0,:,:],cmap='Blues',vmin=0,vmax=10)
    #cb =fig.colorbar(s)
    #cb.set_label('Flow (m$^{3}$s$^{-1}$)',fontsize=14)
    #ax.set_title('Glofas example flow',fontsize=14)
    #ax.set_aspect('auto',adjustable=None)
    #plt.show()
    


    #fig = plt.figure(figsize=(8,6))
    #ax = plt.subplot(111,projection=ccrs.PlateCarree())
    #ax.coastlines(resolution='50m')
    #ax.set_extent([-9,20,40,60])
    #s = plt.pcolormesh(LONS,LATS,data[0,:,:]*MASK_MATRIX_RESHAPE,cmap='Blues',vmin=0,vmax=10)
    #cb =fig.colorbar(s)
    #cb.set_label('Flow (m$^{3}$s$^{-1}$)',fontsize=14)
    #ax.set_title('Glofas example flow',fontsize=14)
    #ax.set_aspect('auto',adjustable=None)
    #plt.show()
    

    # create population mask
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


    # initalise arrays to fill. Note our percentiles are only calculated from data from Oct-Mar 1980-2020. So we only generate SSI for Oct-Mar.

    array_of_totals_SSI = np.zeros([41,oct_mar_days])
    array_of_totals_SSI_pop = np.zeros([41,oct_mar_days])
    array_of_totals_gust = np.zeros([41,oct_mar_days])
    array_of_totals_gust_pop = np.zeros([41,oct_mar_days])

    for yr in range(1980,2021):

        # initialise lists to store interim files
        print(yr)
        SSI_list = []
        SSI_list_pop = []
        gust_list = []
        gust_list_pop = []


        # Oct- Mar
        for mon in [10,11,12,1,2,3]:
            
            # load in data
            fname = 'ERA5_EU_1hr_i10mg_tp_' + str(yr) + '_' + str(mon).zfill(2) + '.nc'
            dataset = Dataset(data_dir + fname ,mode='r')
            gust = dataset.variables['i10fg'][:] # instantaneous
            dataset.close()
            gust_reshape = np.reshape(gust,(int(len(gust)/24),24,len(lat),len(lon)))
            data = np.max(gust_reshape,axis=1)
 
            # initialise arrays to fill for that month
            masked_variable = np.zeros(np.shape(data))
            pop_masked_variable = np.zeros(np.shape(data))
            ratio_of_p98 = np.zeros(np.shape(data))
            SSI = np.zeros(np.shape(data))
            SSI_pop = np.zeros(np.shape(data))
            


            # mask the data
            for i in range(0,len(data)):
                masked_variable[i,:,:] = data[i,:,:]*MASK_MATRIX_RESHAPE
                pop_masked_variable[i,:,:] = data[i,:,:]*pop_mask

                ratio_of_p98[i,:,:] = (masked_variable[i,:,:]/perc_98) - 1


            # calculate the Storm Severity Index
            for i_day in range(0,len(data)):
                #print(i_day)
                for i_lat in range(0,len(lat)):
                    for i_lon in range(0,len(lon)):
                        # if the variable in question < P98 of the variable:
                        if ratio_of_p98[i_day,i_lat,i_lon] > 0.:
                            cubic_data = ratio_of_p98[i_day, i_lat, i_lon]**3
                            SSI_pop[i_day, i_lat , i_lon] = cubic_data*pop_mask[i_lat,i_lon]
                            SSI[i_day,i_lat,i_lon] = cubic_data*MASK_MATRIX_RESHAPE[i_lat,i_lon]
                        

            # aggregate ovcer gridpoints and weight by country area.
            SSI_month = np.sum(np.sum(SSI,axis=2),axis=1)
            SSI_month_pop = np.sum(np.sum(SSI_pop,axis=2),axis=1)
            gust_day = np.sum(np.sum(masked_variable,axis=2),axis=1) / np.sum(MASK_MATRIX_RESHAPE)
            gust_day_pop = np.sum(np.sum(pop_masked_variable,axis=2),axis=1) / np.sum(pop_mask)



            # store to apprpriate lists
            gust_list.append(gust_day)
            gust_list_pop.append(gust_day_pop)
            SSI_list.append(SSI_month)
            SSI_list_pop.append(SSI_month_pop)

        SSI_for_array = np.concatenate(SSI_list)
        SSI_for_array_pop = np.concatenate(SSI_list_pop)
        gust_for_array = np.concatenate(gust_list)
        gust_for_array_pop = np.concatenate(gust_list_pop)

        array_of_totals_SSI[yr-1980,0:len(gust_for_array)] = SSI_for_array
        array_of_totals_SSI_pop[yr-1980,0:len(gust_for_array)] = SSI_for_array_pop
        array_of_totals_gust[yr-1980,0:len(gust_for_array)] = gust_for_array
        array_of_totals_gust_pop[yr-1980,0:len(gust_for_array)] = gust_for_array_pop



    # save the datafiles as numpy files for easy processing later. This can be changed to your preferred file format. 

    #np.save(COUNTRY +'_daily_max_gust.npy',array_of_totals_gust)
    #np.save(COUNTRY +'_daily_max_gust_pop_weighted.npy',array_of_totals_gust_pop)
    #np.save(COUNTRY +'_daily_SSI.npy',array_of_totals_SSI)
    #np.save(COUNTRY + '_daily_SSI95_pop_weighted.npy',array_of_totals_SSI_pop)

