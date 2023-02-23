#
#
# script to download the ERA5 10m winds,T2m and MSLP
#
#
#
# before you do this you need to have the python cdsapi module installed 

import cdsapi

# loop through all months and years
for YEAR in range(1980,2021):
    for MONTH in [1,2,3,4,5,6,7,8,9,10,11,12]:
        m = str(MONTH).zfill(2) # make sure it is 01, 02 etc


        # special conditions for days per month
        if m in ['04','06','09','11']:
            days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30']
        elif m in ['01','03','05','07','08','10','12']:
            days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31']
        else:

            # special conditions for feburary leap years.
            if YEAR in [1980,1984,1988,1992,1996,2000,2004,2008,2012,2016,2020]:
                days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29']
            else:
                days = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28']


        if m == '01':
            MON = 'january'
        elif m == '02':
            MON = 'february'
        elif m == '03':
            MON = 'march'
        elif m == '04':
            MON = 'april'
        elif m == '05':
            MON = 'may'
        elif m == '06':
            MON = 'june'
        elif m == '07':
            MON = 'july'
        elif m == '08':
            MON = 'august'
        elif m == '09':
            MON = 'september'
        elif m == '10':
            MON = 'october'
        elif m == '11':
            MON = 'november'
        else:
            MON = 'december'

        y = str(YEAR)
        c = cdsapi.Client()
        r = c.retrieve('cems-glofas-historical',
                       
                       {'variable':['river_discharge_in_the_last_24_hours'],
                        'system_version':'version_3_1',
                        'hydrological_model':'lisflood',
                        'product_type':'consolidated',
                        'hyear':[y],
                        'hmonth':[MON],
                        'hday':days,
                        'area':[70,-15,35,40], #N/W/S/E
                        'format':'grib'}) # same area as seasonal forecasts
        r.download('glofas_EU_daily_river_discharge_' +str(y) + '_' + str(m) + '.grib')


