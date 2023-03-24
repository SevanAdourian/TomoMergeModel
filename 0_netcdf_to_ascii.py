import numpy as np
import pyshtools
import sys
import matplotlib.pyplot as plt
import yaml
import pandas as pd
import netCDF4 as nc
import xarray as xr
import pdb
import cartopy.crs as ccrs
from Model1D import Model1D
from scipy.interpolate import RectSphereBivariateSpline, RectBivariateSpline


def interpolate_model(ds, target_lats, target_lons, target_depths):
     # Create a new dataset for interpolated values
     new_ds = xr.Dataset()
     
     # Decreasing latitude for increasing colatitudes
     ds = ds.sortby('latitude', ascending=False)
     
     # Interpolate over depth using builtin interp function
     tmp_ds = ds.interp(depth=target_depths, method='linear')
     
     # Get the latitudes, longitudes, and depths of the variable
     lats = tmp_ds.coords['latitude'].values
     lons = tmp_ds.coords['longitude'].values
     depths = tmp_ds.coords['depth'].values
     
     # Convert to colatitude and radians
     colats = np.pi/2 - np.deg2rad(lats)
     lons = np.deg2rad(lons)
     
     target_colats = np.pi/2 - np.deg2rad(np.flipud(target_lats))
     target_lons = np.deg2rad(target_lons)
     
     nlat = len(target_lats)
     nlon = len(target_lons)
     ndepths = len(target_depths)
     
     # Create mask to interpolate on, rest will be filled with NaNs
     min_colat = np.min(colats)
     max_colat = np.max(colats)
     min_lon = np.min(lons)
     max_lon = np.max(lons)
     lat_mask = np.logical_and(target_colats >= min_colat, target_colats <= max_colat)
     lon_mask = np.logical_and(target_lons >= min_lon, target_lons <= max_lon)
     lat_indices, lon_indices = np.where(lat_mask)[0], np.where(lon_mask)[0]
     
     # Loop over each variable in the dataset
     for varname in tmp_ds.data_vars:
          # Interpolate the variable at the target latitudes and longitudes
          # interpolated_data = np.zeros((nlat, nlon, ndepths))
          # interpolated_data = np.empty((nlat, nlon, ndepths)).fill(np.nan)
          interpolated_data = np.full([nlat, nlon, ndepths], np.nan)
          
          # Loop over depths
          for i_d, depth in enumerate(tmp_ds.depth):
               # Get the data array for the variable
               var_data = tmp_ds[varname].isel(depth=i_d).values
               
               # Create an interpolation function for the variable at given depth
               interp_func = RectSphereBivariateSpline(colats, lons, var_data)
               
               # Inteprolate on new grid for given depth
               # interpolated_data[lat_mask,lon_mask,i_d] = interp_func(target_colats[lat_mask],\
               #                                                        target_lons[lon_mask])
               raw_interp = interp_func(target_colats[lat_mask],\
                                        target_lons[lon_mask])
               # breakpoint()
               interpolated_slice = interpolated_data[:,:,i_d]
               interpolated_slice[np.ix_(lat_indices, lon_indices)] = raw_interp
               interpolated_data[:,:,i_d] = interpolated_slice
                                                  

               # breakpoint()
          # Update the data array for the variable in the new dataset
          new_ds[varname]= xr.DataArray(interpolated_data, \
                                        coords=[('latitude', np.flipud(target_lats)), \
                                                ('longitude', np.rad2deg(target_lons)), \
                                                ('depth', target_depths)]) 
          
          new_ds = new_ds.sortby('latitude', ascending=True)
          
     return new_ds


# sys.exit()
# Reading the yaml file to 
for conf_file in sys.argv[1:]:
    print('using configuration loaded from: %s' % (conf_file))
    with open(conf_file) as f:
        conf = yaml.safe_load(f)

regional_model_file = conf['path_to_regional_model']
global_model_file   = conf['path_to_global_model']
spline_file         = conf['depth_file']
ref_model           = conf['path_to_ref_model']

# Load netcdf from file
regional_model = xr.open_dataset(regional_model_file)
global_model   = xr.open_dataset(global_model_file)

# ax = plt.axes(projection=ccrs.EckertIII())
# ax.set_global()
# global_model.vsv.interp(depth=50,method='linear').plot(ax=ax, transform=ccrs.PlateCarree(), \
#                                                       cmap='plasma')
# ax.coastlines()
# plt.show()

# Load spline file
spline_knots_radius =  np.flipud(np.loadtxt(spline_file, skiprows=1))
spline_knots = regional_model.radius_in_meters/1000.0 - spline_knots_radius

# Homogenize units (depth in km, velocitites in km/s)
# Longitude from -180/180 to 0/360 (add 360 to negative values)
# Regional is in m/s and depth in meters
if conf['units_regional_depth'] == 'm':
    print('Converting regional model depth unit to km.')
    regional_model = regional_model.assign_coords(depth=regional_model.depth/1000.0)
elif conf['units_regional_depth'] == 'km':
    print('The regional model is in correct depth unit.')
else:
    sys.exit("Incorrect unit for the depth of the regional model.")
    
# Setting longitude to [0-360]
# lon_0_360 = []
# for lon in regional_model.longitude.to_numpy():
#     if lon < 0:
#         lon_r = lon+360.0
#     else:
#         lon_r = lon
#     lon_0_360.append(lon_r)
# regional_model = regional_model.assign_coords(longitude=lon_0_360)

# lon_0_360 = []
# for lon in global_model.longitude.to_numpy():
#     if lon < 0:
#         lon_g = lon+360.0
#         print(lon, lon_g)
#     else:
#         lon_g = lon
#         print(lon, lon_g)
#     lon_0_360.append(lon_g)

# global_model = global_model.assign_coords(longitude=lon_0_360)

# Dealing with velocity units we want km/s
regional_model = regional_model.assign(VSV_km=lambda x: regional_model.VSV / 1000.0)
regional_model = regional_model.assign(VSH_km=lambda x: regional_model.VSH / 1000.0) 

# Defining target latitude and longitude for equally spaced grid
dlon_reg = regional_model.longitude[1]-regional_model.longitude[0]
dlat_reg = regional_model.latitude[1]-regional_model.latitude[0]

dlon_glo = global_model.longitude[1]-global_model.longitude[0]
dlat_glo = global_model.latitude[1]-global_model.latitude[0]

target_reg_increment = int(np.floor(360./max(dlon_reg,dlat_reg)))
lon_reg = np.linspace(-180,180.,target_reg_increment)
lat_reg = np.linspace(-90.,90.,target_reg_increment)

# print(lon_reg)
# print(lat_reg)
target_glo_increment = max(dlon_glo,dlat_glo)
lon_glo = np.arange(0.,360.+target_glo_increment,target_glo_increment)
lat_glo = np.arange(90.,-90.-target_glo_increment,-target_glo_increment)

# Find relevant splines for regional model
spline_knots_reg = spline_knots[np.where(np.logical_and \
                                         (spline_knots > min(regional_model.depth.to_numpy()), \
                                          spline_knots < max(regional_model.depth.to_numpy())))]
spline_knots_reg_radius = regional_model.radius_in_meters/1000.0 - spline_knots_reg

# Interpolate models on equally spaced grid and chosen depths
# regional_int = regional_model.interp(latitude=lat_reg, longitude=lon_reg, \
#                                      depth=spline_knots_reg, method='linear')
regional_int = interpolate_model(regional_model, lat_reg, lon_reg, spline_knots_reg)
# global_int   = global_model.interp(latitude=lat_glo, longitude = lon_glo, \
#                                    depth=spline_knots, method='linear')
global_int   = global_model.interp(depth=spline_knots, method='linear')

# Extract variables to numpy arrays
vsv_reg = regional_int.VSV_km.to_numpy()
vsh_reg = regional_int.VSH_km.to_numpy()
vsv_glo = global_int.vsv.to_numpy()
vsh_glo = global_int.vsh.to_numpy()

# Read 1D model
ref = Model1D(ref_model)
ref.load_from_file()
vsv = ref.get_values(1000 * spline_knots_radius, parameter='vsv')
vsh = ref.get_values(1000 * spline_knots_radius, parameter='vsh')
vs0 = 0.001 * np.sqrt((2 * vsv ** 2 + vsh ** 2) / 3)

# Compute Voigt average for Vs
vs0_reg_abs = np.sqrt((2 * vsh_reg ** 2 + vsv_reg ** 2) / 3)
vs0_glo_abs = np.sqrt((2 * vsh_glo ** 2 + vsv_glo ** 2) / 3)

# Compute perturbations
sizes_glo = global_int.vsv.sizes
dims_glo  = global_int.vsv.dims
sizes_reg = regional_int.VSV_km.sizes
dims_reg  = regional_int.VSV_km.dims

vs0_glo = np.repeat(vs0,((global_int.latitude.size*global_int.longitude.size))).\
    reshape(sizes_glo[dims_glo[0]],sizes_glo[dims_glo[1]],sizes_glo[dims_glo[2]])
vs0_reg = np.repeat(vs0[:-13:],((regional_int.latitude.size*regional_int.longitude.size))).\
    reshape(sizes_reg[dims_reg[0]],sizes_reg[dims_reg[1]],sizes_reg[dims_reg[2]])

# pdb.set_trace()
vs0_reg = (vs0_reg_abs - vs0_reg) / vs0_reg
vs0_glo = (vs0_glo_abs - vs0_glo) / vs0_glo

# 
regional_int = regional_int.assign(vs0=(('latitude','longitude','depth'), vs0_reg))
global_int = global_int.assign(vs0=(('depth','latitude','longitude'), vs0_glo))
global_int = global_int.sortby('longitude')

# regional_int.sel(depth=30).plot(cmap='plasma')
# ax = plt.axes(projection=ccrs.EckertIII())
# ax.set_global()
# global_int.vsv.interp(depth=250,method='linear').plot(ax=ax, transform=ccrs.PlateCarree(), \
#                                                       cmap='plasma')
# ax.coastlines()
# plt.show()

# ax = plt.axes(projection=ccrs.EckertIII())
# ax.set_global()
# global_int.vs0.interp(depth=250,method='linear').plot(ax=ax, transform=ccrs.PlateCarree(), \
#                                                       cmap='plasma')
# ax.coastlines()
# plt.show()

# pdb.set_trace()

# Write ascii files
base_global_ascii_files = conf['path_to_ascii_files']+'/global_'
base_reg_ascii_files    = conf['path_to_ascii_files']+'/regional_'

# for i,dep in enumerate(spline_knots):
#     global_ascii_filename   = base_global_ascii_files+str(format(int(dep),'04d'))+".xyz"
#     global_int.sel(depth=dep).to_dataframe().swaplevel(0).to_csv(global_ascii_filename, \
#                                                                  sep='\t', na_rep=0, \
#                                                                  float_format='%8.3f', \
#                                                                  columns=['vs0'],header=False)

for i,dep in enumerate(spline_knots_reg):
    regional_ascii_filename = base_reg_ascii_files+str(format(int(dep),'04d'))+".xyz"
    regional_int.sel(depth=dep).to_dataframe().swaplevel(0).to_csv(regional_ascii_filename, \
                                                                   sep='\t', na_rep=0, \
                                                                   float_format='%8.3f', \
                                                                   columns=['vs0'],header=False)
    
