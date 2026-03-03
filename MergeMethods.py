# imports here
from Dataset import Dataset
import multiprocessing
import os
import concurrent.futures
import time
import sys
import yaml
import pyshtools
import numpy as np
import pandas as pd
import xarray as xr
import pdb
import matplotlib.pyplot as plt

class MergeMethods():
    # constructors
    def __init__(self, modelOne: Dataset, modelTwo: Dataset, confFilePath, regionalVariable, globalVariable):

        # input validation
        if modelOne is None or modelTwo is None:
            raise ValueError("Null model(s) provided")
        
        # checking to see that there is one regional and global model, and assigning them
        regional_mod = None
        global_mod = None
        if modelOne.getModelType() == Dataset.REGIONAL:
            regional_mod = modelOne
        if modelOne.getModelType() == Dataset.GLOBAL:
            global_mod = modelOne
        if modelTwo.getModelType() == Dataset.REGIONAL:
            regional_mod = modelTwo
        if modelTwo.getModelType() == Dataset.GLOBAL:
            global_mod = modelTwo

        # raise an error for two global or two regional models
        if regional_mod is None or global_mod is None:
            raise ValueError("Constructor requires one REGIONAL and one GLOBAL Dataset")

        # store the models
        self.regional_model = regional_mod
        self.global_model = global_mod
        self.merge_model = None

        # store config for merging variables
        if confFilePath is None:
            raise ValueError("Must provide a configuration file containing merging variables")

        try:
            self.conf = self.parseConfFile(confFilePath=confFilePath)
        except ValueError as e:
            raise ValueError(e)
        except:
            raise ValueError("Unable to parse config file")

        if regionalVariable not in list(self.regional_model.getDataset().data_vars):
            raise ValueError(f"Could not find regional variable {regionalVariable} in regional model")
        if globalVariable not in list(self.global_model.getDataset().data_vars):
            raise ValueError(f"Could not find global variable {globalVariable} in global model")
        self.regionalVariable = regionalVariable
        self.globalVariable = globalVariable

    def containsAllParams(self, conf):
        params = ['path_to_regional_model', 
                  'path_to_global_model',
                  'depth_file',
                  'max_depth_regional_model',
                  'units_regional_depth',
                  'reg_lmax',
                  'win_lmax',
                  'win_eff_lmax',
                  'lon_min_mask',
                  'lon_max_mask',
                  'lat_min_mask',
                  'lat_max_mask',
                  'radius_in_meters']
        for param in params:
            if param not in conf:
                raise ValueError(f"Param not found: {param}")

    def parseConfFile(self, confFilePath):
        with open(confFilePath) as f:
            conf = yaml.safe_load(f)
        self.containsAllParams(conf)
        depth_file = conf['depth_file'] 
        #  - Mask parameters
        # Part of the regional model that will be preserved in the merging.
        # Read in the depth file
        conf['depth_knots_radius'] =  np.array(self.regional_model.depths)
        conf['depth_knots'] = conf['radius_in_meters']/1000.0 - conf['depth_knots_radius']
        conf['base_global_ascii_files'] = conf['path_to_ascii_files']+'/global_'
        conf['base_regional_ascii_files']    = conf['path_to_ascii_files']+'/regional_'
        conf['base_merged_ascii_files'] = conf['path_to_ascii_files']+'/merged_'

        return conf


    # getters and setters here
    def getRegionalModel(self) -> Dataset:
        return self.regional_model

    def getGlobalModel(self) -> Dataset:
        return self.global_model

    def getConf(self):
        return self.conf

    def setRegionalModel(self, model: Dataset) -> None:
        if model is None:
            raise ValueError("Null model provided")
        if model.getModelType() != Dataset.REGIONAL:
            raise ValueError("setRegionalModel requires a REGIONAL Dataset")
        self.regional_model = model

    def setGlobalModel(self, model: Dataset) -> None:
        if model is None:
            raise ValueError("Null model provided")
        if model.getModelType() != Dataset.GLOBAL:
            raise ValueError("setGlobalModel requires a GLOBAL Dataset")
        self.global_model = model
    
    def setConf(self, confFilePath):
        if confFilePath is None:
            raise ValueError("Must provide a configuration file for merging variables")
        
        try:
            self.conf = self.parseConfFile(confFilePath=confFilePath)
        except:
            raise ValueError("Unable to parse config file")

    # merge function and helper functions here
    def convert_to_spherical_harmonics(self, zmesh, reg_lmax):
        #   - Convert to spherical harmonics with pySHtools
        mask = np.isnan(zmesh)
        zmesh[mask] = np.nanmean(np.nanmean(zmesh))
        grid = pyshtools.SHGrid.from_array(zmesh, grid= 'DH')
        clm = pyshtools.SHGrid.expand(grid)
        clm = clm.pad(reg_lmax)  #Pad to match global clm

        # Trim coefficients in SPH to match number used in regional model
        grid = pyshtools.SHCoeffs.expand(clm,lmax=reg_lmax)

        return grid, clm

    def create_window(self, reg_field, global_clm):
        '''
        Create a mask of spherical harmonic windows using the range given above
        Sharpness of the edges of the mask are controlled by reg_nwin
        Normalise mask between 0 and 1 and then apply to the regional data
        '''
        
        #-----------
        xvar = 'longitude'
        xlen = (len(reg_field[xvar]))
        yvar = 'latitude'
        ylen = (len(reg_field[yvar]))

        lon_left=self.conf['lon_min_mask']
        lon_right=self.conf['lon_max_mask']
        lat_bottom = self.conf['lat_min_mask'] # possible range: -90, 90 deg
        lat_top    =  self.conf['lat_max_mask'] # possible range: -90, 90 deg

        #  - Regional mask - set 1s inside region and 0s outside
        reg_zmask= np.where( ( (reg_field[xvar] >= lon_left) & (reg_field[xvar] <= lon_right) )&\
                            ( (reg_field[yvar] >= lat_bottom) & (reg_field[yvar] <= lat_top) ),1,0)
        # 
        reg_zmesh_mask=(np.reshape(reg_zmask,(ylen,xlen)))

        print(f"Inside create_window")
        pdb.set_trace()
        #-----------
        if self.conf['win_type'] == 'spherical': # spherical or rectangular'
            #   - Construct spherical harmonic window function from mask

            reg_win=pyshtools.SHWindow.from_mask(reg_zmesh_mask,lwin=self.conf['win_lmax'])
            reg_win_clm=pyshtools.SHWindow.to_shcoeffs(reg_win,0)
            reg_win_clm_pad=reg_win_clm.pad(global_clm.lmax)  #Pad to match global clm
            
            reg_win_energy = (reg_win.to_shgrid(0).to_array())**2
            for i in range(1, self.conf['win_eff_lmax']):
                reg_win_energy += (reg_win.to_shgrid(i).to_array())**2
            reg_win_energy_grid = pyshtools.SHGrid.from_array(reg_win_energy)
                
            #  - Normalise window (0 to 1) so that we can mask outside and inside
            valmax=(np.amax(reg_win_energy_grid.data))
            reg_win_energy_grid=reg_win_energy_grid/valmax
        
            reg_win_energy_clm = pyshtools.SHGrid.expand(reg_win_energy_grid)
            reg_win_energy_clm=reg_win_energy_clm.pad(self.conf['reg_lmax'])  #Pad to match global clm
            reg_win_energy_grid = pyshtools.SHCoeffs.expand(reg_win_energy_clm)  #Grid of mask

        elif self.conf['win_type'] == 'rectangular': # spherical or rectangular
            # Construct rectangular window from mask (no smoothing at edges, sharp window)
            reg_win_grid=pyshtools.SHGrid.from_array(reg_zmesh_mask)
            reg_win_clm = pyshtools.SHGrid.expand(reg_win_grid)
            reg_win_clm_pad=reg_win_clm.pad(global_clm.lmax)  #Pad to match reg clm
            reg_win_energy_grid=pyshtools.SHCoeffs.expand(reg_win_clm_pad)

            
        return reg_win_energy_grid
    
    def apply_window(self, global_grid, global_clm, reg_grid, reg_win_energy_grid):

        # - Multiply grid by mask
        reg_grid_masked=reg_grid*reg_win_energy_grid  
        reg_clm_masked=pyshtools.SHGrid.expand(reg_grid_masked)

        # Plot map and spectra for reg_clm_masked
        self.plot_map_and_spectra(reg_grid_masked, reg_clm_masked, "plots_178_lmax/region.png")

        # Sum regional spectrum inside box (masked) with global spectra outside box (unmasked)
        summed_clm=reg_clm_masked+global_clm

        # Plot map and spectra for summed_clm

        summed_grid=pyshtools.SHCoeffs.expand(summed_clm)
        self.plot_map_and_spectra(summed_grid, summed_clm, "plots_178_lmax/merged.png")
        return summed_grid
    
    def process_slice(self, depth):
    # Reading in global tomography model
        zmesh_global = self.reshape_field(self.global_model.getDataset(), depth, self.globalVariable)
        
        global_grid, global_clm = self.convert_to_spherical_harmonics(zmesh_global.values, self.conf["reg_lmax"])
        # Plot map and spectra for global_clm

        if depth == 240:
            self.plot_map_and_spectra(global_grid, global_clm, "plots_178_lmax/global_240.png")
        
        if depth < self.conf['max_depth_regional_model']:
            # Above where the regional model is defined in depth, actual merging
            # Reading in regional tomography model
            zmesh_regional = self.reshape_field(self.regional_model.getDataset(), depth, self.regionalVariable)
            reg_grid, reg_clm = self.convert_to_spherical_harmonics(zmesh_regional.values, self.conf["reg_lmax"])
            # Doing mask windowing
            pdb.set_trace()
            # reg_win_energy_grid = self.create_window(self.regional_model.getDataset().sel(depth=depth), global_clm)
            reg_win_energy_grid = self.create_window(zmesh_regional, global_clm)
            summed_grid = self.apply_window(global_grid, global_clm, reg_grid, reg_win_energy_grid)
            pdb.set_trace()
        else:
            # Below where the regional model is defined, we just write the global model
            summed_grid = global_grid

        return self.write_model(depth=depth, grid=summed_grid)
    
    def write_model(self, depth, grid):
        '''
        Write out masked, merged model as 3 columns: lon, lat, dv 
        '''
        m_dv = grid.data
        lats = grid.lats()
        lons = grid.lons()

        # Adding a 3rd dimension 
        m_dv = np.expand_dims(m_dv, axis=-1)
        da = xr.DataArray(
            data=m_dv,
            dims=("latitude", "longitude", "depth"),
            coords={"latitude": lats, "longitude": lons, "depth": [depth]},
            name="dv"
        )

        return da

    def plot_map_and_spectra(self, grid_object, clm_object, file_name):
        fig, (col1, col2) = plt.subplots(2, 1)
        grid_object.plot(ax=col1, colorbar='right', cb_label='Power', show=False)
        clm_object.plot_spectrum(ax=col2)
        fig.legend(loc = 'upper right')
        fig.savefig(file_name, dpi=400)
        plt.close(fig)
        return

    def reshape_field(self, lon_lat_field, depth, varname):
        zmesh = lon_lat_field[varname].sel(depth=float(depth))
        zmesh = zmesh.sortby('longitude')
        zmesh.assign_coords(latitude=zmesh.latitude[::-1])
        # zmesh_values = np.flipud(zmesh.values)
        return zmesh

    def merge(self):
        depths = list(self.conf['depth_knots'])

        # On Windows, process-based multiprocessing uses "spawn" which re-imports modules.
        # Using a process Pool with bound methods/self can also hang due to pickling.
        # A thread pool avoids both issues and is reliable here.
        merged_arrays = []
        for depth in depths:
            depth_val = float(depth)
            print("depth - ", depth_val)
            merged_arrays.append(self.process_slice(depth_val))

        # Actually concatneation returns a DataArray, not a DataSet
        merged_all_arrays = xr.concat(merged_arrays, dim="depth")
        self.merge_model = merged_all_arrays.to_dataset(name="merged")

        # Needed for multiprocessing
        self.merge_model = self.merge_model.sortby("depth")
        # self.merge_model = self.merge_model.assign_coords(latitude=self.merge_model.latitude[::-1]).sortby("latitude")
        # self.merge_model = self.merge_model.assign_coords(longitude=((self.merge_model.longitude + 180) % 360)).sortby("longitude")
        self.merge_model = Dataset("", Dataset.GLOBAL, xrDataset=self.merge_model, depthUnits='km')

        return self.merge_model
