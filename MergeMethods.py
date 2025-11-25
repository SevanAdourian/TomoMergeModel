# imports here
from Dataset import Dataset
import multiprocessing
import time
import sys
import yaml
import pyshtools
import numpy as np
import pandas as pd
import xarray as xr

class MergeMethods():
    # constructors
    def __init__(self, modelOne: Dataset, modelTwo: Dataset, confFilePath):

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
        conf['depth_knots_radius'] =  np.flipud(np.loadtxt(depth_file, skiprows=1))
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
        xvar = np.unique(reg_field[:,0])
        xlen = (len(xvar))
        yvar = np.unique(reg_field[:,1])
        ylen = (len(yvar))

        lon_left=self.conf['lon_min_mask']
        lon_right=self.conf['lon_max_mask']
        lat_bottom = self.conf['lat_min_mask'] # possible range: -90, 90 deg
        lat_top    =  self.conf['lat_max_mask'] # possible range: -90, 90 deg

        #  - Regional mask - set 1s inside region and 0s outside
        reg_zmask= np.where( ( (reg_field[:,0] >= lon_left) & (reg_field[:,0] <= lon_right) )&\
                            ( (reg_field[:,1] >= lat_bottom) & (reg_field[:,1] <= lat_top) ),1,0)
        # 
        reg_zmesh_mask=(np.reshape(reg_zmask,(ylen,xlen)))
        reg_mask_xyz=np.column_stack((reg_field[:,0],reg_field[:,1],reg_zmask))
        
        #-----------
        if self.conf['win_type'] == 'spherical': # spherical or rectangular'
            #   - Construct spherical harmonic window function from mask
            reg_win=pyshtools.SHWindow.from_mask(reg_zmesh_mask,lwin=self.conf['win_lmax'])
            reg_win_clm=pyshtools.SHWindow.to_shcoeffs(reg_win,0)
            reg_win_clm_pad=reg_win_clm.pad(global_clm.lmax)  #Pad to match global clm
            
            reg_win_energy = (reg_win.to_shgrid(0).to_array())**2
            for i in range(1,self.conf['win_eff_lmax']):
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
        # pdb.set_trace()
        reg_grid_masked=reg_grid*reg_win_energy_grid  
        reg_clm_masked=pyshtools.SHGrid.expand(reg_grid_masked)
        
        #   - Apply windows to global grid
        # pdb.set_trace()
        global_grid_masked=global_grid*reg_win_energy_grid
        global_clm_masked=pyshtools.SHGrid.expand(global_grid_masked)
        
        # Sum regional and global spectra inside box
        # pdb.set_trace()
        summed_clm_masked=reg_clm_masked+global_clm_masked

        # Sum regional spectrum inside box (masked) with global spectra outside box (unmasked)
        # pdb.set_trace()
        summed_clm=reg_clm_masked+global_clm
        summed_grid=pyshtools.SHCoeffs.expand(summed_clm)

        return summed_grid
    
    def process_slice(self, depth):
    # Reading in global tomography model
        zmesh_global = self.reshape_field(self.global_model.getDataset())

        global_grid, global_clm = self.convert_to_spherical_harmonics(zmesh_global, self.conf["reg_lmax"])
        
        if depth < self.conf['max_depth_regional_model']:
            # Above where the regional model is defined in depth, actual merging
            # Reading in regional tomography model
            zmesh_regional = self.reshape_field(self.regional_model.getDataset())
            
            reg_grid, reg_clm = self.convert_to_spherical_harmonics(zmesh_regional, self.conf["reg_lmax"])
            # pdb.set_trace()
            # Doing mask windowing
            reg_win_energy_grid = self.create_window(self.regional_model.getDataset(), global_clm)
            summed_grid = self.apply_window(global_grid, global_clm, reg_grid, reg_win_energy_grid)
            # ForkedPdb().set_trace()
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

        da = xr.DataArray(
            data=m_dv,
            dims=("lat", "lon"),
            coords={"lat": lats, "lon": lons, "depth": [depth]},
            name="dv"
        )

        return da

    def reshape_field(self, lon_lat_field):
        #  - Create mesh of z
        xvar = np.unique(lon_lat_field[:,0])
        xlen = (len(xvar))
        yvar = np.unique(lon_lat_field[:,1])
        ylen = (len(yvar))
        zmesh = (np.reshape(lon_lat_field[:,2],(ylen,xlen)))

        return zmesh

    def merge(self):
        # Initiate computing slices on multiple processes
        multiprocessing.set_start_method('spawn')

        # Loop over depths with multiprocessing:
        depths = self.conf['depth_knots']  # get depth values from xarray dataset
        n_depths = len(depths)  # number of times to run process_slice
        
        with multiprocessing.Pool() as pool:
            # Use map to apply process_slice to each depth and collect results directly
            merged_arrays = pool.map(self.process_slice, depths)

        self.merge_model = xr.concat(merged_arrays, dim="depth")

        self.merge_model = self.merge_model.sortby("depth")

        self.merge_model = Dataset("", Dataset.GLOBAL, xrDataset=self.merge_model)

        return self.merge_model