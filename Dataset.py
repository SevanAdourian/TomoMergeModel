# imports here
import os, yaml
import numpy as np
import pyshtools
import sys
import matplotlib.pyplot as plt
import pandas as pd
import netCDF4 as nc
import xarray as xr
import pdb
import cartopy.crs as ccrs
from Model1D import Model1D
from scipy.interpolate import RectSphereBivariateSpline, RectBivariateSpline

class Dataset():
    REGIONAL = 0
    GLOBAL = 1
    # constructors
    def __init__(self, filePath: str, modelType: int, splineFilePath = None, globalModel = None):

        # Input validation for file path
        if filePath is None:
            raise ValueError("Need to pass in file path")
        if type(filePath) is not str:
            raise ValueError("File path must be a string")
        if not os.path.exists(filePath):
            raise ValueError("Incorrect file path: file could not be found")
        
        # Input validation for model type
        if modelType not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        
        if modelType == Dataset.REGIONAL:
            if splineFilePath == None or globalModel == None or globalModel.getModelType() == Dataset.REGIONAL:
                raise ValueError("Global model needed for initalization of regional models")
            
        
        # assign the file path and model type, and parse file for other attributes
        self.filePath = filePath
        self.modelType = modelType

        self.dataset = xr.open_dataset(filePath)

        if self.modelType == Dataset.REGIONAL:
            target_lats, target_lons, target_depths = getInterpolationParameters(globalModel, splineFilePath)

            self.dataset = Dataset.interpolate_model(self.dataset, target_lats, target_lons, target_depths)

    # Initialize a Dataset from the conf.yaml file
    @classmethod
    def initFromConf(modelType: int):

        # Input validation for model type
        if modelType not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        
        # load the conf.yaml file
        try:
            with open("conf.yaml", "r") as file:
                conf = yaml.safe_load(file)
        except OSError:
            raise ValueError("conf.yaml not found or not in working directory")
        # return to the user, a Dataset instance depending on the model type and file path in conf.yaml
        if modelType == Dataset.REGIONAL:
            return Dataset(conf["path_to_regional_model"], modelType)
        if modelType == Dataset.GLOBAL:
            return Dataset(conf["path_to_global_model"], modelType)


    # getters and setters here
    def getFilePath(self) -> str:
        return self.filePath

    def setFilePath(self, path: str):
        if path is None:
            raise ValueError("Need to pass in file path")
        if not isinstance(path, str):
            raise ValueError("File path must be a string")
        if not os.path.exists(path):
            raise ValueError("Incorrect file path: file could not be found")
        self.filePath = path

    def getModelType(self) -> int:
        return self.modelType

    def setModelType(self, mt: int):
        if mt not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        self.modelType = mt
    
    @staticmethod
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


    # dataset specific functions here
    pass