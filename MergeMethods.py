# imports here
from Dataset import Dataset
from ConfigParams import ConfigParams
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


class MergeMethods:
    # constructors
    def __init__(
        self,
        modelOne: Dataset,
        modelTwo: Dataset,
        configParams: ConfigParams,
        regionalVariable: str,
        globalVariable: str,
    ):
        """Initialize with one regional and one global dataset, config parameters, and variable names."""

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
        if configParams is None:
            raise ValueError(
                "Must provide a ConfigParmas class containing merging variables"
            )

        self.conf = configParams

        if regionalVariable not in list(self.regional_model.getDataset().data_vars):
            raise ValueError(
                f"Could not find regional variable {regionalVariable} in regional model"
            )
        if globalVariable not in list(self.global_model.getDataset().data_vars):
            raise ValueError(
                f"Could not find global variable {globalVariable} in global model"
            )
        self.regionalVariable = regionalVariable
        self.globalVariable = globalVariable

    # getters and setters here
    def getRegionalModel(self) -> Dataset:
        """Return the regional model dataset."""

        return self.regional_model

    def getGlobalModel(self) -> Dataset:
        """Return the global model dataset."""

        return self.global_model

    def getConf(self):
        """Return the configuration parameters."""

        return self.conf

    def setRegionalModel(self, model: Dataset) -> None:
        """Set the regional model after validation."""

        if model is None:
            raise ValueError("Null model provided")
        if model.getModelType() != Dataset.REGIONAL:
            raise ValueError("setRegionalModel requires a REGIONAL Dataset")
        self.regional_model = model

    def setGlobalModel(self, model: Dataset) -> None:
        """Set the global model after validation."""

        if model is None:
            raise ValueError("Null model provided")
        if model.getModelType() != Dataset.GLOBAL:
            raise ValueError("setGlobalModel requires a GLOBAL Dataset")
        self.global_model = model

    def setConf(self, configParams):
        """Set the configuration parameters after validation."""

        if configParams is None:
            raise ValueError("Must provide a configuration file for merging variables")

        self.conf = configParams

    # merge function and helper functions here
    def convert_to_spherical_harmonics(self, zmesh, reg_lmax):
        """Convert a 2D grid to spherical harmonics coefficients up to specified lmax."""

        # Set nan values to the mean of the region
        mask = np.isnan(zmesh)
        zmesh[mask] = np.nanmean(np.nanmean(zmesh))

        # Create the grid and coefficients
        grid = pyshtools.SHGrid.from_array(zmesh, grid="DH")
        clm = pyshtools.SHGrid.expand(grid)
        clm = clm.pad(reg_lmax)  # Pad to match global clm

        # Trim coefficients in SPH to match number used in regional model
        grid = pyshtools.SHCoeffs.expand(clm, lmax=reg_lmax)

        return grid, clm

    def create_window(self, reg_field, global_clm):
        """
        Create a mask of spherical harmonic windows using the range given above
        Sharpness of the edges of the mask are controlled by reg_nwin
        Normalise mask between 0 and 1 and then apply to the regional data
        """

        reg_field = reg_field.sortby("longitude")

        xvar = "longitude"
        yvar = "latitude"

        # get the window bounds
        lon_left = self.conf.lon_min_mask
        lon_right = self.conf.lon_max_mask
        lat_bottom = self.conf.lat_min_mask
        lat_top = self.conf.lat_max_mask

        #  - Regional mask - set 1s inside region and 0s outside
        lon = reg_field[xvar].values
        lat = reg_field[yvar].values

        lon2d, lat2d = np.meshgrid(lon, lat)

        reg_zmesh_mask = np.where(
            (lon2d >= lon_left)
            & (lon2d <= lon_right)
            & (lat2d >= lat_bottom)
            & (lat2d <= lat_top),
            1,
            0,
        )

        if self.conf.win_type == "spherical":  # spherical or rectangular'
            #   - Construct spherical harmonic window function from mask

            reg_win = pyshtools.SHWindow.from_mask(
                reg_zmesh_mask, lwin=self.conf.win_lmax
            )
            reg_win_clm = pyshtools.SHWindow.to_shcoeffs(reg_win, 0)
            reg_win_clm_pad = reg_win_clm.pad(
                global_clm.lmax
            )  # Pad to match global clm

            reg_win_energy = (reg_win.to_shgrid(0).to_array()) ** 2
            for i in range(1, self.conf.win_eff_lmax):
                reg_win_energy += (reg_win.to_shgrid(i).to_array()) ** 2
            reg_win_energy_grid = pyshtools.SHGrid.from_array(reg_win_energy)

            #  - Normalise window (0 to 1) so that we can mask outside and inside
            valmax = np.amax(reg_win_energy_grid.data)
            reg_win_energy_grid = reg_win_energy_grid / valmax

            reg_win_energy_clm = pyshtools.SHGrid.expand(reg_win_energy_grid)
            reg_win_energy_clm = reg_win_energy_clm.pad(
                self.conf.reg_lmax
            )  # Pad to match global clm
            reg_win_energy_grid = pyshtools.SHCoeffs.expand(
                reg_win_energy_clm
            )  # Grid of mask

        elif self.conf.win_type == "rectangular":  # spherical or rectangular
            # Construct rectangular window from mask (no smoothing at edges, sharp window)
            reg_win_grid = pyshtools.SHGrid.from_array(reg_zmesh_mask)
            reg_win_clm = pyshtools.SHGrid.expand(reg_win_grid)
            reg_win_clm_pad = reg_win_clm.pad(global_clm.lmax)  # Pad to match reg clm
            reg_win_energy_grid = pyshtools.SHCoeffs.expand(reg_win_clm_pad)

        reg_win_energy_grid.data = np.flipud(reg_win_energy_grid.data)

        return reg_win_energy_grid

    def apply_window(self, global_grid, global_clm, reg_grid, reg_win_energy_grid):
        """Apply windowing mask to regional grid and merge with global model."""

        # - Multiply grid by mask
        reg_grid_masked = reg_grid * reg_win_energy_grid
        reg_clm_masked = pyshtools.SHGrid.expand(reg_grid_masked)

        # Sum regional spectrum inside box (masked) with global spectra outside box (unmasked)
        summed_clm = reg_clm_masked + global_clm

        # Plot map and spectra for summed_clm

        summed_grid = pyshtools.SHCoeffs.expand(summed_clm)
        return summed_grid

    def process_slice(self, depth):
        """Process and merge regional/global models at a single depth slice."""

        # Reading in global tomography model
        zmesh_global = self.reshape_field(
            self.global_model.getDataset(), depth, self.globalVariable
        )

        global_grid, global_clm = self.convert_to_spherical_harmonics(
            zmesh_global, self.conf.reg_lmax
        )

        if depth <= max(list(self.regional_model.depths)):
            # Above where the regional model is defined in depth, actual merging
            # Reading in regional tomography model
            zmesh_regional = self.reshape_field(
                self.regional_model.getDataset(), depth, self.regionalVariable
            )
            reg_grid, reg_clm = self.convert_to_spherical_harmonics(
                zmesh_regional, self.conf.reg_lmax
            )
            # Doing mask windowing
            # reg_win_energy_grid = self.create_window(self.regional_model.getDataset().sel(depth=depth), global_clm)
            reg_win_energy_grid = self.create_window(
                self.regional_model.getDataset().sel(depth=depth), global_clm
            )
            summed_grid = self.apply_window(
                global_grid, global_clm, reg_grid, reg_win_energy_grid
            )
        else:
            # Below where the regional model is defined, we just write the global model
            summed_grid = global_grid

        return self.write_model(depth=depth, grid=summed_grid)

    def write_model(self, depth, grid):
        """
        Write out masked, merged model as 3 columns: lon, lat, dv
        """
        m_dv = grid.data
        lats = grid.lats()
        lons = grid.lons()

        # Adding a 3rd dimension
        m_dv = np.expand_dims(m_dv, axis=-1)
        da = xr.DataArray(
            data=m_dv,
            dims=("latitude", "longitude", "depth"),
            coords={"latitude": lats, "longitude": lons, "depth": [depth]},
            name="dv",
        )

        return da

    def plot_map_and_spectra(self, grid_object, clm_object, file_name):
        """Plot spatial grid and power spectrum side by side and save to file."""

        fig, (col1, col2) = plt.subplots(2, 1)
        grid_object.plot(ax=col1, colorbar="right", cb_label="Power", show=False)
        clm_object.plot_spectrum(ax=col2)
        fig.legend(loc="upper right")
        fig.savefig(file_name, dpi=400)
        plt.close(fig)
        return

    def reshape_field(self, lon_lat_field, depth, varname):
        """Extract and reshape a 2D field at a given depth for spherical harmonic expansion."""

        zmesh = lon_lat_field[varname].sel(depth=float(depth))
        zmesh = zmesh.sortby("longitude")
        # zmesh.assign_coords(latitude=zmesh.latitude[::-1])
        zmesh = np.flipud(zmesh.values)
        return zmesh

    def merge(self):
        """Merge regional and global models across all depth slices and return merged Dataset."""

        depths = self.regional_model.depths

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
        self.merge_model = Dataset(
            "", Dataset.GLOBAL, xrDataset=self.merge_model, depthUnits="km"
        )

        return self.merge_model
