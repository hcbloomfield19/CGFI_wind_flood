######
#
#
# THis script does the correlation analysis shown in Bloomfield et al., (2023) but also a very similar version (slightly different uncertainty sampling) 
 # to that on https://the-iea.github.io/cgfi-wind-flood/
#
######

# import libraries
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scipy.stats as sp
import matplotlib as mpl
import pandas as pd
import cartopy.crs as ccrs # this is on jasmin2
import cartopy.io.shapereader as shpreader
import geopandas as gpd

octmar_days = 31+30+31+31+29+31 # include leap day.
yrs = 40 # 1980-2020 data included



# countries processed
COUNTRY_LIST = ['Great Britain','All Ireland','France','Germany','Spain','Portugal','Italy','Belgium','Luxembourg','Netherlands','Norway','Sweden','Denmark','Poland','Czechia','Austria','Hungary','Switzerland','Latvia','Lithuania','Belarus','Ukraine','Slovakia','Slovenia','Romania','Bulgaria','Croatia','Greece','Montenegro','Moldova','Finland', 'Bosnia and Her','Estonia','Poland','Serbia','North Macedoni','Albania']



#setup empty lists for storing correlations
day_corr = []
week_corr = []
mon_corr = []
seas_corr = []

# identify key timescales of interest to calcualte correlations for - we could do more or less of these if required. You can pick any of interest to you
timesteps_xday = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30,40,50,60,90,120,150,180] 


for COUNTRY in COUNTRY_LIST:

    print(COUNTRY)
    # load in the precipiation, wind gusts and river flow (country average)
    #Oct-Mar 
    precip = np.load('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/' + COUNTRY + '_daily_sum_tp.npy')
    gust = np.load('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/' + COUNTRY + '_daily_max_gust.npy')

    # Load in the storm severity index (SSI) data  Oct- Mar
    gust_damage = np.load('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/' + COUNTRY + '_daily_SSI_pop_weighted.npy')

    # Load in glofas and FSI data/ Jan-Dec year, 366 days (last one blank if not a leap year)
    glofas = np.load('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/' + COUNTRY + '_daily_sum_glofas_riverflow_ALL_MON.npy')
    glofas_damage = np.load('/home/fz21704/code_folder/CGFI/find_percentiles_of_weather_data/' + COUNTRY + '_daily_sum_glofas_FSI_995_ALL_MON_gridthr_popweight.npy')



    # reshpe these into  [ years * days ] note the glofas data is a different shape - (Jan-Dec) not (Oct-Mar).
    # You can simplify this a lot if you just calculate data over the same time periods!)

    reordered_precip = np.zeros([yrs,octmar_days])
    reordered_gust = np.zeros([yrs,octmar_days])
    reordered_glofas= np.zeros([yrs,octmar_days])
    reordered_FSI = np.zeros([yrs,octmar_days])
    reordered_SSI = np.zeros([yrs,octmar_days])


    # loop through years and reorder
    leapyrs = [0,4,8,12,16,20,24,28,32,36,40]

    # get the data into Oct-Mar winter seasons.
    for qyr in range(0,yrs):

        reordered_precip[qyr,0:31+30+31] = precip[qyr,0:31+30+31]
        reordered_precip[qyr,31+30+31:octmar_days] = precip[qyr+1,31+30+31:octmar_days]
        reordered_gust[qyr,0:31+30+31] = gust[qyr,0:31+30+31]
        reordered_gust[qyr,31+30+31:octmar_days] = gust[qyr+1,31+30+31:octmar_days]
        reordered_SSI[qyr,0:31+30+31] = gust_damage[qyr,0:31+30+31]
        reordered_SSI[qyr,31+30+31:octmar_days] = gust_damage[qyr+1,31+30+31:octmar_days]
        if qyr in leapyrs:
            reordered_FSI[qyr,0:31+30+31] = glofas_damage[qyr,274:366]
            reordered_FSI[qyr,31+30+31:octmar_days] = glofas_damage[qyr+1,0:91]
            reordered_glofas[qyr,0:31+30+31] = glofas[qyr,274:366]
            reordered_glofas[qyr,31+30+31:octmar_days] = glofas[qyr+1,0:91]
        else:
            reordered_FSI[qyr,0:31+30+30] = glofas_damage[qyr,274:365]
            reordered_FSI[qyr,31+30+30:octmar_days] = glofas_damage[qyr+1,0:92]
            reordered_glofas[qyr,0:31+30+30] = glofas[qyr,274:365]
            reordered_glofas[qyr,31+30+30:octmar_days] = glofas[qyr+1,0:92]
 
 
 


    # define arrays for storing bootstrapped results
    bootstrapped_correlations1 = np.zeros([1000,len(timesteps_xday)])
    bootstrapped_correlations2 = np.zeros([1000,len(timesteps_xday)])
    bootstrapped_correlations3 = np.zeros([1000,len(timesteps_xday)])

    # more lists to store the values for the different variable lines
    corr_xday = []
    corr_xday2 = []
    corr_xday3 = []
    pval_xday = []
    pval_xday2 = []
    pval_xday3 = []

    
    counter = 0
    
    for day_scales in timesteps_xday:


        # divide all of our data into chunks based on how many years we have.
        # don't re-use chunks.
        ds = int(day_scales)
        divisor = int(180/ds)
        total_points_ERA5 = 40*divisor
        data_to_use1 = np.reshape(reordered_gust[:,0:divisor*ds],(40,divisor,ds))
        xday_data1 = np.mean(data_to_use1,axis=2)
        data_to_use2 = np.reshape(reordered_precip[:,0:divisor*ds],(40,divisor,ds))
        xday_data2 = np.sum(data_to_use2,axis=2)

        data_to_use3 = np.reshape(reordered_glofas[:,0:divisor*ds],(40,divisor,ds))
        xday_data3 = np.sum(data_to_use3,axis=2)

        data_to_use4 = np.reshape(reordered_FSI[:,0:divisor*ds],(40,divisor,ds))
        xday_data4 = np.sum(data_to_use4,axis=2)

        data_to_use5 = np.reshape(reordered_SSI[:,0:divisor*ds],(40,divisor,ds))
        xday_data5 = np.sum(data_to_use5,axis=2)


        # calculate spearmans rank corrleation
        a = sp.spearmanr(xday_data1.flatten(),xday_data2.flatten())

        corr_xday.append(a[0])
        pval_xday.append(a[1])

        b = sp.spearmanr(xday_data1.flatten(),xday_data3.flatten())

        corr_xday2.append(b[0])
        pval_xday2.append(b[1])


        c = sp.spearmanr(xday_data4.flatten(),xday_data5.flatten())

        corr_xday3.append(c[0])
        pval_xday3.append(c[1])


        # boostrap the data (96% confidence interval) to see how robust results are to sample size.
        for samp_i in range(0,1000):
            rand_samp = np.random.randint(low = 0, high=total_points_ERA5,size=int(total_points_ERA5*0.95))
            d1 =xday_data1.flatten()
            d2 = xday_data2.flatten()
            d3 = xday_data3.flatten()
            d4 = xday_data4.flatten()
            d5 = xday_data5.flatten()
 

            samp1 = d1[rand_samp]
            samp2 = d2[rand_samp]
            samp3 = d3[rand_samp]
            samp4 = d4[rand_samp]
            samp5 = d5[rand_samp]
 

            sampa = sp.spearmanr(samp1,samp2)
            bootstrapped_correlations1[samp_i,counter] = sampa[0]

            sampb = sp.spearmanr(samp1,samp3)
            bootstrapped_correlations2[samp_i,counter] = sampb[0]

            sampc = sp.spearmanr(samp4,samp5)
            bootstrapped_correlations3[samp_i,counter] = sampc[0]


        counter += 1

     
    std_bs = np.std(bootstrapped_correlations1,axis=0)
    std_bs2 = np.std(bootstrapped_correlations2,axis=0)
    std_bs3 = np.std(bootstrapped_correlations3,axis=0)


    # plot it all out as a pretty line graph
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(1,1,1)
    plt.plot(np.array(timesteps_xday),np.array(corr_xday),color='purple',label='Wind gust vs. precip')
    plt.fill_between(np.array(timesteps_xday), np.array(corr_xday) + std_bs, np.array(corr_xday) - std_bs, alpha=0.5, color='purple')

    plt.plot(np.array(timesteps_xday),np.array(corr_xday2),color='orange',label='Wind gust vs. riverflow')
    plt.fill_between(np.array(timesteps_xday), np.array(corr_xday2) + std_bs2, np.array(corr_xday2) - std_bs2, alpha=0.5, color='orange')

    plt.plot(np.array(timesteps_xday),np.array(corr_xday3),color='lightblue',label='SSI vs. FSI')
    plt.fill_between(np.array(timesteps_xday), np.array(corr_xday3) + std_bs3, np.array(corr_xday3) - std_bs3, alpha=0.5, color='lightblue')

    plt.ylim([-0.2,1])
    plt.xlim([1,180])

    plt.xlabel('Timescale of interest (days)',fontsize=14)
    plt.ylabel('Spearmans Rank correlation',fontsize=14)
    plt.title('ERA5 ' + COUNTRY ,fontsize=16)

    plt.legend(frameon=False,loc='lower left')
    plt.savefig(COUNTRY + '_correlations.png')
    plt.close()


    # append the correlations for maps later (you could change the numbers here if you care about different timescales!)
    day_corr.append(corr_xday2[0])
    week_corr.append(corr_xday2[7])
    mon_corr.append(corr_xday2[17])
    seas_corr.append(corr_xday2[24])



# create a dataframe of the daily/weekly/monthly/seasonal correlations

d_clusters_EUROPE = {'Country Name': COUNTRY_LIST ,'Days': day_corr ,'Weeks': week_corr,'Months': mon_corr,'Seasons':seas_corr}
df_clusters_EUROPE = pd.DataFrame(data=d_clusters_EUROPE)




# setup colorbar
bounds = np.linspace(0,1.0,11)
nbins = np.float(len(bounds)) -2 # must not be an integer
bins = bounds
cmap = plt.cm.YlGnBu
norm = mpl.colors.BoundaryNorm(bounds,cmap.N)


values_country = df_clusters_EUROPE['Country Name']
values_cluster_EUROPE = df_clusters_EUROPE['Days']
df_clusters_EUROPE['bin'] = np.digitize(values_cluster_EUROPE, bins,right=True)


# plot out the correlations as a map

fig = plt.figure(figsize=(7,7))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-15,40,32,72])
ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
countries_shp = shpreader.natural_earth(resolution='10m',category='cultural', name='admin_0_countries')


    
for i in range(0,len(values_country)): # loop through each country

    color = cmap(df_clusters_EUROPE['bin'][i]/nbins)

    if i==0:
        COUNTRY == 'Great Britain'
        GB_data = gpd.read_file("/home/Data/NUTS/NUTS_1/Great_Britain.shp")
        ax.add_geometries(GB_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
    elif i ==1:
        COUNTRY == 'All Ireland'
        IE_data = gpd.read_file("/home/Data/NUTS/NUTS_1/All_Ireland.shp")
        ax.add_geometries(IE_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')

    else:
        for country in shpreader.Reader(countries_shp).records():
            if country.attributes['NAME'][0:len(values_country[i])] == values_country[i]:
                print(values_country[i])

                if country.geometry.geom_type=='MultiPolygon':
                    ax.add_geometries(country.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
                elif country.geometry.geom_type=='Polygon': 
                    # this is a single geometry so need to put it in a list to plot
                    ax.add_geometries([country.geometry], ccrs.PlateCarree(),facecolor=color,edgecolor='k')

                 
ax1=fig.add_axes([0.2,0.1,0.6,0.05]) # [left, bottom, width, height !!]
cb1 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,norm=norm,extend='min',orientation='horizontal',ticks = [-0.2,0,0.2,0.4,0.6,0.8,1])
cb1.set_label('Daily Correlation ',fontsize=20)
cb1.ax.tick_params(labelsize=20)
plt.savefig('Europe_daily_correlations.png')
plt.close()



values_country = df_clusters_EUROPE['Country Name']
values_cluster_EUROPE = df_clusters_EUROPE['Weeks']
df_clusters_EUROPE['bin'] = np.digitize(values_cluster_EUROPE, bins,right=True)




fig = plt.figure(figsize=(7,7))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-15,40,32,72])
ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
countries_shp = shpreader.natural_earth(resolution='10m',category='cultural', name='admin_0_countries')
   

for i in range(0,len(values_country)): # loop through each country

    color = cmap(df_clusters_EUROPE['bin'][i]/nbins)

    if i==0:
        COUNTRY == 'Great Britain'
        GB_data = gpd.read_file("/home/Data/NUTS/NUTS_1/Great_Britain.shp")
        ax.add_geometries(GB_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
    elif i ==1:
        COUNTRY == 'All Ireland'
        IE_data = gpd.read_file("/home/Data/NUTS/NUTS_1/All_Ireland.shp")
        ax.add_geometries(IE_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')

    else:
        for country in shpreader.Reader(countries_shp).records():
            if country.attributes['NAME'][0:len(values_country[i])] == values_country[i]:
                print(values_country[i])

                if country.geometry.geom_type=='MultiPolygon':
                    ax.add_geometries(country.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
                elif country.geometry.geom_type=='Polygon': 
                    # this is a single geometry so need to put it in a list to plot
                    ax.add_geometries([country.geometry], ccrs.PlateCarree(),facecolor=color,edgecolor='k')

 
                 
ax1=fig.add_axes([0.2,0.1,0.6,0.05]) # [left, bottom, width, height !!]
cb1 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,norm=norm,extend='min',orientation='horizontal',ticks = [-0.2,0,0.2,0.4,0.6,0.8,1])
cb1.set_label('Weekly Correlation ',fontsize=20)
cb1.ax.tick_params(labelsize=20)
plt.savefig('Europe_weekly_correlations.png')
plt.close()


values_country = df_clusters_EUROPE['Country Name']
values_cluster_EUROPE = df_clusters_EUROPE['Months']
df_clusters_EUROPE['bin'] = np.digitize(values_cluster_EUROPE, bins,right=True)


fig = plt.figure(figsize=(7,7))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-15,40,32,72])
ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
countries_shp = shpreader.natural_earth(resolution='10m',category='cultural', name='admin_0_countries')
    
for i in range(0,len(values_country)): # loop through each country

    color = cmap(df_clusters_EUROPE['bin'][i]/nbins)

    if i==0:
        COUNTRY == 'Great Britain'
        GB_data = gpd.read_file("/home/Data/NUTS/NUTS_1/Great_Britain.shp")
        ax.add_geometries(GB_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
    elif i ==1:
        COUNTRY == 'All Ireland'
        IE_data = gpd.read_file("/home/Data/NUTS/NUTS_1/All_Ireland.shp")
        ax.add_geometries(IE_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')

    else:
        for country in shpreader.Reader(countries_shp).records():
            if country.attributes['NAME'][0:len(values_country[i])] == values_country[i]:
                print(values_country[i])

                if country.geometry.geom_type=='MultiPolygon':
                    ax.add_geometries(country.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
                elif country.geometry.geom_type=='Polygon': 
                    # this is a single geometry so need to put it in a list to plot
                    ax.add_geometries([country.geometry], ccrs.PlateCarree(),facecolor=color,edgecolor='k')

 
                 
ax1=fig.add_axes([0.2,0.1,0.6,0.05]) # [left, bottom, width, height !!]
cb1 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,norm=norm,extend='min',orientation='horizontal',ticks = [-0.2,0,0.2,0.4,0.6,0.8,1])
cb1.set_label('Monthly Correlation ',fontsize=20)
cb1.ax.tick_params(labelsize=20)
plt.savefig('Europe_monthly_correlations.png')
plt.close()


values_country = df_clusters_EUROPE['Country Name']
values_cluster_EUROPE = df_clusters_EUROPE['Seasons']
df_clusters_EUROPE['bin'] = np.digitize(values_cluster_EUROPE, bins,right=True)


fig = plt.figure(figsize=(7,7))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-15,40,32,72])

ax.background_patch.set_visible(False)
ax.outline_patch.set_visible(False)
countries_shp = shpreader.natural_earth(resolution='10m',category='cultural', name='admin_0_countries')
    
for i in range(0,len(values_country)): # loop through each country

    color = cmap(df_clusters_EUROPE['bin'][i]/nbins)

    if i==0:
        COUNTRY == 'Great Britain'
        GB_data = gpd.read_file("/home/Data/NUTS/NUTS_1/Great_Britain.shp")
        ax.add_geometries(GB_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
    elif i ==1:
        COUNTRY == 'All Ireland'
        IE_data = gpd.read_file("/home/Data/NUTS/NUTS_1/All_Ireland.shp")
        ax.add_geometries(IE_data.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')

    else:
        for country in shpreader.Reader(countries_shp).records():
            if country.attributes['NAME'][0:len(values_country[i])] == values_country[i]:
                print(values_country[i])

                if country.geometry.geom_type=='MultiPolygon':
                    ax.add_geometries(country.geometry, ccrs.PlateCarree(),facecolor=color,edgecolor='k')
                elif country.geometry.geom_type=='Polygon': 
                    # this is a single geometry so need to put it in a list to plot
                    ax.add_geometries([country.geometry], ccrs.PlateCarree(),facecolor=color,edgecolor='k')

 
                 
ax1=fig.add_axes([0.2,0.1,0.6,0.05]) # [left, bottom, width, height !!]
cb1 = mpl.colorbar.ColorbarBase(ax1,cmap=cmap,norm=norm,extend='min',orientation='horizontal',ticks = [-0.2,0,0.2,0.4,0.6,0.8,1])
cb1.set_label('Seasonal Correlation ',fontsize=20)
cb1.ax.tick_params(labelsize=20)
plt.savefig('Europe_seasonal_correlations.png')
plt.close()


