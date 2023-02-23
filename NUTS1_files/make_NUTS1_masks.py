import numpy as np 
import geopandas as gpd
import xarray as xr
import regionmask
from netCDF4 import Dataset
import matplotlib.pyplot as plt

NUTLEV =1
path_to_shapefile = '/home/Data/NUTS/NUTS_' +str(NUTLEV) + '/NUTS_RG_01M_2021_4326_LEVL_' + str(NUTLEV) + '.shp'
nuts = gpd.read_file(path_to_shapefile)
print(nuts.head())

len_field = len(nuts)

folder_with_netcdf ='/path_to_ERA5/ERA5_EU_1hr_i10mg_tp_2015_12.nc'

d = xr.open_mfdataset(folder_with_netcdf)
data = d['i10fg']
test =data.longitude.data


# all 125 NUTS1 regions
codes_to_do = list(range(0,125))

nuts_mask_poly = regionmask.Regions(name='nuts_mask',numbers =codes_to_do,names = list(nuts.NUTS_ID) ,abbrevs = list(nuts.NUTS_ID),outlines = list(nuts.geometry.values[i] for i in codes_to_do))

print(nuts_mask_poly)


mask = nuts_mask_poly.mask(data.isel(time=0),lat_name='latitude',lon_name='longitude')

plt.figure(figsize=(12,8))
ax=plt.axes()
mask.plot(ax=ax)
nuts.plot(ax=ax,alpha=0.8,facecolor='none',lw=1)
plt.title('NUTS '+ str(NUTLEV),fontsize=16)
plt.show()

country_strings=[]

# create final masks
final_mask_array = np.zeros([len_field,len(mask.latitude.values),len(mask.longitude.values)])
LATS = mask.latitude.values
LONS = mask.longitude.values
 
for i in range(0,len_field):
    print(i)
    sel_mask = mask.where(mask==i).values
    # change it to [1,0] rather than [number,nan]
    sel_mask[sel_mask == i] = 1.
    sel_mask[np.isnan(sel_mask)] = 0.
    # save into an array we can load from later.
    final_mask_array[i,:,:] = sel_mask
    country_strings.append(str(nuts.NUTS_ID[i]))

    
plt.imshow(final_mask_array[106,:,:],cmap='Reds')
plt.colorbar()
plt.title(country_strings[106])
plt.show()


# save as a netcdf
dataset_1 = Dataset('NUTS_' + str(NUTLEV) + '_masks.nc','w',format='NETCDF4')
lat_dim = dataset_1.createDimension('lat',len(LATS))
lon_dim = dataset_1.createDimension('lon',len(LONS))
NUTS_dim = dataset_1.createDimension('NUTS',len_field)

lats = dataset_1.createVariable('lat',np.float32,('lat',))   # (var_name. type of data, Dimension)
lons = dataset_1.createVariable('lon',np.float32,('lon',))
if NUTLEV == 0:
    NUTS_keys = dataset_1.createVariable('NUTS keys','S2',('NUTS',))
if NUTLEV == 1:
    NUTS_keys = dataset_1.createVariable('NUTS keys','S3',('NUTS',))
if NUTLEV == 2:
    NUTS_keys = dataset_1.createVariable('NUTS keys','S4',('NUTS',))


NUTS_array = dataset_1.createVariable('NUTS zones',np.float32,('NUTS','lat','lon'))

NUTS_array[:] = final_mask_array
lats[:] = LATS
lons[:] = LONS
for i in range(0,len_field):
    NUTS_keys[i] = country_strings[i]

dataset_1.close()



