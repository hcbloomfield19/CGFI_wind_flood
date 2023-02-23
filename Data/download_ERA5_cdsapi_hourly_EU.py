#
#
# script to download the ERA5 gusts and precip
#
#
#
# before you do this you need to do make sure you've installed the python cdsapi module


import cdsapi



# loop through all months and years to download the data. 

for YEAR in range(1980,2022):
    for MONTH in [10,11,12,1,2,3]:

        # sort out the different days in a month
        m = str(MONTH).zfill(2) # make sure it is 01, 02 etc
        if m in ['04','06','09','11']:
            days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30']
        elif m in ['01','03','05','07','08','10','12']:
            days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']
        else:

            # special conditions for leap years and february.
            if YEAR in [1962,1966,1972,1976,1980,1984,1988,1992,1996,2000,2004,2008,2012,2016,2020]:
                days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']
            else:
                days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28']

        y = str(YEAR)

        c = cdsapi.Client()
       	c.retrieve('reanalysis-era5-single-levels',
                       {'variable':['instantaneous_10m_wind_gust','total_precipitation'],
                        'product_type':'reanalysis',
                        'year':[y],
                        'month':[m],
                        'day':days,
                        'time':['00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00','08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00','17:00','18:00','19:00','20:00','21:00','22:00','23:00'],
                        'area' : [70,-15,34,35],# N/W/S/E 
                        'format':'netcdf'},'ERA5_EU_1hr_i10mg_tp_' +str(y) + '_' + str(m) + '.nc')


