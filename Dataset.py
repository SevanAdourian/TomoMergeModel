from __future__ import annotations

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


class Dataset:
    REGIONAL = 0  # Regional Dataset object
    GLOBAL = 1  # Global Dataset object
    RADIUS_IN_METERS = 6371000

    # constructors
    def __init__(
        self,
        filePath: str,
        modelType: int,
        depthUnits: str = None,
        xrDataset: xr.Dataset = None,
        globalModel: Dataset = None,
        depths: list[int] = None,
    ):

        # Input validation for file path
        if filePath is None:
            raise ValueError("Need to pass in file path")

        # Input validation for model type
        if modelType not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")

        # If model is regional, and a spline file and global model are given, depth units are needed
        if modelType == Dataset.REGIONAL:
            if depths is not None and globalModel is not None and depthUnits is None:
                raise ValueError(
                    "If interpolation on depths is needed, depth units must be given"
                )

        # assign the file path and model type, and parse file for other attributes
        self.filePath: str = filePath
        self.modelType: int = modelType
        self.depths: list[int] = depths
        if depths == None and modelType == Dataset.REGIONAL:
            self.depths = globalModel.getDataset().depth.values

        # open the file, reassign variable names, homogenize units, and convert longitude values
        if filePath == "":
            if xrDataset is None:
                raise ValueError("If empty filePath is given, xrDataset must be given")
            self.dataset: xr.Dataset = xrDataset
        else:
            self.dataset: xr.Dataset = xr.open_dataset(filePath)

        # assign names to varialbes and convert the units
        self.assign_names()
        self.convert_units(depthUnits)

        if self.modelType == Dataset.REGIONAL:  # interpolating the regional model
            target_lats, target_lons, target_depths = self.getInterpolationParameters(
                depths, globalModel=globalModel
            )
            self.dataset: xr.Dataset = Dataset.interpolate_model(
                self.dataset, target_lats, target_lons, target_depths
            )

        # convert the longitude
        self.convert_longitude()

    # Initialize a Dataset from the conf.yaml file
    @classmethod
    def initFromConf(modelType: int, globalModel: Dataset = None) -> Dataset:

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
            return Dataset(
                conf["path_to_regional_model"],
                modelType,
                conf["depth_file"],
                None,
                conf["path_to_global_model"],
                conf["units_regional_depth"],
                globalModel=globalModel,
            )
        if modelType == Dataset.GLOBAL:
            return Dataset(conf["path_to_global_model"], modelType)

    # get the interpolation parameters based on the regional and global model
    def getInterpolationParameters(
        self, depths: list[int], globalModel: Dataset
    ) -> tuple:
        # get target latitude and longitude
        glo: xr.Dataset = globalModel.getDataset()
        dlon_reg: xr.DataArray = glo.longitude[1] - glo.longitude[0]
        dlat_reg: xr.DataArray = glo.latitude[1] - glo.latitude[0]

        # get the regional increment
        target_reg_increment: int = int(np.floor(360.0 / max(dlon_reg, dlat_reg)))
        lon_reg: np.ndarray = np.linspace(-180, 180.0, target_reg_increment)
        lat_reg: np.ndarray = np.linspace(-90.0, 90.0, target_reg_increment)

        spline_knots_reg = None
        if depths is not None:
            # get spline knots from the spline file
            spline_knots_radius: np.ndarray = np.array(depths)
            spline_knots: np.ndarray = (
                Dataset.RADIUS_IN_METERS / 1000.0 - spline_knots_radius
            )

            # Get the relevant splines for the regional model
            spline_knots_reg: np.ndarray = spline_knots[
                np.where(
                    np.logical_and(
                        spline_knots > min(self.dataset.depth.to_numpy()),
                        spline_knots < max(self.dataset.depth.to_numpy()),
                    )
                )
            ]
        return (lat_reg, lon_reg, spline_knots_reg)

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

    def getDataset(self) -> xr.Dataset:
        return self.dataset

    def getModelType(self) -> int:
        return self.modelType

    def setModelType(self, modelType: int):
        if modelType not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        self.modelType = modelType

    # interpolate a depth slice
    def interpolate_slice(
        colats,
        lons,
        var_data,
        target_colats,
        lat_mask,
        target_lons,
        lon_mask,
        interpolated_data,
        i_d,
        lat_indices,
        lon_indices,
    ):
        interp_func = RectSphereBivariateSpline(colats, lons, var_data)

        # Inteprolate on new grid for given depth
        raw_interp = interp_func(target_colats[lat_mask], target_lons[lon_mask])
        if i_d is not None:
            interpolated_slice = interpolated_data[:, :, i_d]
        else:
            interpolated_slice = interpolated_data[:, :]

        interpolated_slice[np.ix_(lat_indices, lon_indices)] = raw_interp
        return interpolated_slice

    @staticmethod
    # Interpolate all depth slices
    def interpolate_model(
        ds: xr.Dataset,
        target_lats: np.ndarray,
        target_lons: np.ndarray,
        target_depths: np.ndarray = None,
    ):
        # Create a new dataset for interpolated values
        new_ds: xr.Dataset = xr.Dataset()

        # Decreasing latitude for increasing colatitudes
        ds: xr.Dataset = ds.sortby("latitude", ascending=False)

        # Interpolate over depth using builtin interp function
        if target_depths is not None:
            tmp_ds: xr.Dataset = ds.interp(depth=target_depths, method="linear")
        else:
            tmp_ds: xr.Dataset = ds.copy()

        # Get the latitudes, longitudes, and depths of the variable
        lats: np.ndarray = tmp_ds.coords["latitude"].values
        lons: np.ndarray = tmp_ds.coords["longitude"].values

        # Convert to colatitude and radians
        colats: np.ndarray = np.pi / 2 - np.deg2rad(lats)
        lons: np.ndarray = np.deg2rad(lons)

        target_colats: np.ndarray = np.pi / 2 - np.deg2rad(np.flipud(target_lats))
        target_lons: np.ndarray = np.deg2rad(target_lons)

        nlat: int = len(target_lats)
        nlon: int = len(target_lons)
        ndepths = None
        if target_depths is not None:
            ndepths: int = len(target_depths)

        # Create mask to interpolate on, rest will be filled with NaNs
        min_colat = np.min(colats)
        max_colat = np.max(colats)
        min_lon = np.min(lons)
        max_lon = np.max(lons)
        lat_mask = np.logical_and(
            target_colats >= min_colat, target_colats <= max_colat
        )
        lon_mask = np.logical_and(target_lons >= min_lon, target_lons <= max_lon)
        lat_indices, lon_indices = np.where(lat_mask)[0], np.where(lon_mask)[0]

        # Loop over each variable in the dataset
        for varname in tmp_ds.data_vars:
            # Interpolate the variable at the target latitudes and longitudes
            shape = [nlat, nlon, ndepths] if ndepths is not None else [nlat, nlon]
            interpolated_data = np.full(shape, np.nan)

            # Loop over depths if there are depths, otherwise just do it for the one layer
            if ndepths is not None:
                for i_d, depth in enumerate(tmp_ds.depth):
                    # Get the data array for the variable
                    var_data = tmp_ds[varname].isel(depth=i_d).values

                    # Create an interpolation function for the variable at given depth
                    interpolated_data[:, :, i_d] = Dataset.interpolate_slice(
                        colats,
                        lons,
                        var_data,
                        target_colats,
                        lat_mask,
                        target_lons,
                        lon_mask,
                        interpolated_data,
                        i_d,
                        lat_indices,
                        lon_indices,
                    )
            else:
                var_data = tmp_ds[varname].values
                interpolated_data[:, :] = Dataset.interpolate_slice(
                    colats,
                    lons,
                    var_data,
                    target_colats,
                    lat_mask,
                    target_lons,
                    lon_mask,
                    interpolated_data,
                    None,
                    lat_indices,
                    lon_indices,
                )

            # Update the data array for the variable in the new dataset
            new_ds[varname] = xr.DataArray(
                interpolated_data,
                coords=[
                    ("latitude", np.flipud(target_lats)),
                    ("longitude", np.rad2deg(target_lons)),
                    ("depth", target_depths),
                ],
            )

            new_ds = new_ds.sortby("latitude", ascending=True)

        return new_ds

    def convert_units(self, depthUnit: str = "m"):

        if depthUnit == "m":
            self.dataset = self.dataset.assign_coords(depth=self.dataset.depth / 1000.0)
            if "VSV" in self.dataset:
                self.dataset = self.dataset.assign(
                    VSV=lambda x: self.dataset.VSV / 1000.0
                )

            if "VSH" in self.dataset:
                self.dataset = self.dataset.assign(
                    VSH=lambda x: self.dataset.VSH / 1000.0
                )

        elif depthUnit != "km":
            raise ValueError("Depth Units must be m or km")

    def convert_longitude(self):
        lon_0_360 = []
        for lon in self.dataset.longitude.to_numpy():
            if lon < 0:
                lon_r = lon + 360.0
            else:
                lon_r = lon
            lon_0_360.append(lon_r)

        self.dataset = self.dataset.assign_coords(longitude=lon_0_360)

    def assign_names(self):
        name_map = {
            "latitude": ["lat", "latitude", "Latitude", "LAT"],
            "longitude": ["lon", "longitude", "Longitude", "LON"],
            "depth": ["depth", "DEPTH", "z", "level"],
        }

        for key in name_map.keys():
            found = False
            for value in name_map[key]:
                if value in self.dataset.variables:
                    self.dataset = self.dataset.rename({value: key})
                    found = True
                    break

            if not found:
                raise ValueError(f"Could not detect variable name for {key}")

    def save_model(self, filePath: str):
        if filePath is None:
            raise ValueError("File path to save model must be given")
        self.dataset.to_netcdf(filePath)

    def plot_variable(
        self,
        varname: str,
        depth=None,
        ax=None,
        cmap="viridis",
        vmin=None,
        vmax=None,
        title: str = None,
        show_colorbar=True,
        figsize=(10, 5),
    ):
        """
        Plot one variable from the dataset on a PlateCarree map.
        - varname: name of a data variable in the xr.Dataset
        - depth: either an index (int) or a value (float) to pick nearest depth for 3D variables.
                 If None and variable has a depth dimension, the first level is used.
        Returns (fig, ax).
        """
        if varname not in self.dataset:
            raise ValueError(f"{varname} not found in dataset")

        da = self.dataset[varname]

        # handle depth dimension if present
        if "depth" in da.dims:
            if depth is None:
                da_sel = da.isel(depth=0)
                depth_label = str(da.depth.values[0])
            elif isinstance(depth, (int, np.integer)):
                da_sel = da.isel(depth=int(depth))
                depth_label = str(da.depth.values[int(depth)])
            else:
                # assume numeric value -> find nearest depth
                idx = int(np.argmin(np.abs(da.depth.values - float(depth))))
                da_sel = da.isel(depth=idx)
                depth_label = str(da.depth.values[idx])
        else:
            da_sel = da
            depth_label = None

        lats = da_sel["latitude"].values
        lons = da_sel["longitude"].values
        data = da_sel.values

        # prepare figure / axis
        created_fig = False
        if ax is None:
            fig, ax = plt.subplots(
                figsize=figsize, subplot_kw={"projection": ccrs.PlateCarree()}
            )
            created_fig = True
        else:
            fig = ax.figure

        # pcolormesh expects 2D lon/lat grids or 1D coords with matching shapes
        # Use transform to PlateCarree for proper georeferencing
        pcm = ax.pcolormesh(
            lons,
            lats,
            data,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        ax.coastlines()
        try:
            ax.set_global()
        except Exception:
            pass

        if title is None:
            title = varname
            if depth_label is not None:
                title = f"{title} (depth={depth_label})"
        ax.set_title(title)

        if show_colorbar:
            cbar = fig.colorbar(
                pcm, ax=ax, orientation="vertical", pad=0.02, fraction=0.04
            )
            cbar.set_label(varname)

        if created_fig:
            plt.show()

        return fig, ax

    def plot_all_variables(
        self,
        depth=None,
        save_dir: str = None,
        cmap="viridis",
        vmin=None,
        vmax=None,
        figsize=(10, 5),
        show=True,
    ):
        """
        Iterate through all data variables in the dataset and plot each.
        - depth: forwarded to plot_variable for 3D variables
        - save_dir: if provided, PNG files are written to this directory (created if needed)
        - show: whether to call plt.show() after each plot (useful to turn off when saving many files)
        Returns list of (varname, filepath_or_None).
        """
        results = []
        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)

        for var in self.dataset.data_vars:
            fig, ax = self.plot_variable(
                var,
                depth=depth,
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
                figsize=figsize,
                show_colorbar=True,
            )
            filepath = None
            if save_dir is not None:
                safe_name = var.replace("/", "_").replace(" ", "_")
                filepath = os.path.join(save_dir, f"{safe_name}_{str(depth)}.png")
                fig.savefig(filepath, bbox_inches="tight", dpi=150)
            if show:
                plt.show()
            else:
                plt.close(fig)
            results.append((var, filepath))
        return results
