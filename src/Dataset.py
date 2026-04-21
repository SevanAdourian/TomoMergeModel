from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
from scipy.interpolate import RectSphereBivariateSpline


class Dataset:
    REGIONAL = 0  # Regional Dataset object
    GLOBAL = 1  # Global Dataset object
    RADIUS_IN_METERS = 6371000

    def __init__(
        self,
        file_path: str,
        model_type: int,
        depth_units: str = None,
        xr_dataset: xr.Dataset = None,
        global_model: Dataset = None,
        depths: list[int] = None,
    ):
        """Initialize a dataset and normalize its coordinates and units."""

        # Input validation for file path
        if file_path is None:
            raise ValueError("Need to pass in file path")

        # Input validation for model type
        if model_type not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")

        # If model is regional, and a spline file and global model are given, depth units are needed
        if model_type == Dataset.REGIONAL:
            if depths is not None and global_model is not None and depth_units is None:
                raise ValueError(
                    "If interpolation on depths is needed, depth units must be given"
                )

        # assign the file path and model type, and parse file for other attributes
        self.file_path: str = file_path
        self.model_type: int = model_type
        self.depths: list[int] = depths
        if depths is None and model_type == Dataset.REGIONAL:
            if global_model is None:
                raise ValueError(
                    "globalModel must be provided for REGIONAL models when depths is None"
                )
            self.depths = global_model.getDataset().depth.values

        # Load from an existing xr.Dataset or open from file
        if file_path == "":
            if xr_dataset is None:
                raise ValueError("If empty filePath is given, xrDataset must be given")
            self.dataset: xr.Dataset = xr_dataset
        else:
            self.dataset: xr.Dataset = xr.open_dataset(file_path)

        # Normalize coordinate names, units, and longitude convention
        self.assign_names()
        if depth_units is not None:
            self.convert_units(depth_units)

        if self.model_type == Dataset.REGIONAL:
            target_lats, target_lons, target_depths = self.getInterpolationParameters(
                self.depths, global_model=global_model
            )
            self.dataset: xr.Dataset = Dataset.interpolate_model(
                self.dataset, target_lats, target_lons, target_depths
            )

        # convert the longitude
        self.convert_longitude()

    def getInterpolationParameters(
        self, depths: list[int], global_model: Dataset
    ) -> tuple:
        """Return target latitude/longitude grids and optional depth knots."""

        # get target latitude and longitude
        glo: xr.Dataset = global_model.getDataset()
        dlon_reg: xr.DataArray = glo.longitude[1] - glo.longitude[0]
        dlat_reg: xr.DataArray = glo.latitude[1] - glo.latitude[0]

        # get the regional increment
        target_reg_increment: int = int(np.floor(360.0 / max(dlon_reg, dlat_reg)))
        lon_reg: np.ndarray = np.linspace(-180, 180.0, target_reg_increment)
        lat_reg: np.ndarray = np.linspace(-90.0, 90.0, target_reg_increment)

        spline_knots_reg = None
        if depths is not None:
            spline_knots: np.ndarray = np.array(depths)

            # Filter to depths that fall within the regional model's depth range
            spline_knots_reg: np.ndarray = spline_knots[
                np.where(
                    np.logical_and(
                        spline_knots > min(self.dataset.depth.to_numpy()),
                        spline_knots < max(self.dataset.depth.to_numpy()),
                    )
                )
            ]
        return (lat_reg, lon_reg, spline_knots_reg)

    def getFilePath(self) -> str:
        """Return the source file path for this dataset."""

        return self.file_path

    def setFilePath(self, path: str):
        """Set the dataset file path after validating it exists."""

        if path is None:
            raise ValueError("Need to pass in file path")
        if not isinstance(path, str):
            raise ValueError("File path must be a string")
        if not os.path.exists(path):
            raise ValueError("Incorrect file path: file could not be found")
        self.file_path = path

    def getDataset(self) -> xr.Dataset:
        """Return the underlying xarray dataset."""

        return self.dataset

    def getModelType(self) -> int:
        """Return the model type constant for this dataset."""

        return self.model_type

    def setModelType(self, model_type: int):
        """Set the model type after validating allowed values."""

        if model_type not in [Dataset.REGIONAL, Dataset.GLOBAL]:
            raise ValueError("Model type must be Dataset.REGIONAL or Dataset.GLOBAL")
        self.model_type = model_type

    @staticmethod
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
        """Interpolate one 2D depth slice onto the target spherical grid."""

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
    def interpolate_model(
        ds: xr.Dataset,
        target_lats: np.ndarray,
        target_lons: np.ndarray,
        target_depths: np.ndarray = None,
    ):
        """Interpolate all variables to target latitude/longitude and optional depth grids."""

        new_ds: xr.Dataset = xr.Dataset()

        # Sort descending in latitude so colatitudes are ascending
        ds: xr.Dataset = ds.sortby("latitude", ascending=False)

        # Interpolate to target depths (if provided) before spatial interpolation
        if target_depths is not None:
            tmp_ds: xr.Dataset = ds.interp(depth=target_depths, method="linear")
        else:
            tmp_ds: xr.Dataset = ds.copy()

        lats: np.ndarray = tmp_ds.coords["latitude"].values
        lons: np.ndarray = tmp_ds.coords["longitude"].values

        # Convert lat/lon to colatitude (radians) required by RectSphereBivariateSpline
        colats: np.ndarray = np.pi / 2 - np.deg2rad(lats)
        lons: np.ndarray = np.deg2rad(lons)

        target_colats: np.ndarray = np.pi / 2 - np.deg2rad(np.flipud(target_lats))
        target_lons: np.ndarray = np.deg2rad(target_lons)

        nlat: int = len(target_lats)
        nlon: int = len(target_lons)
        ndepths = None
        if target_depths is not None:
            ndepths: int = len(target_depths)

        # Build lat/lon masks so points outside the regional extent stay NaN
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
            shape = [nlat, nlon, ndepths] if ndepths is not None else [nlat, nlon]
            interpolated_data = np.full(shape, np.nan)

            if ndepths is not None:
                for i_d, depth in enumerate(tmp_ds.depth):
                    var_data = tmp_ds[varname].isel(depth=i_d).values
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

            if target_depths is not None:
                new_ds[varname] = xr.DataArray(
                    interpolated_data,
                    coords=[
                        ("latitude", np.flipud(target_lats)),
                        ("longitude", np.rad2deg(target_lons)),
                        ("depth", target_depths),
                    ],
                )
            else:
                new_ds[varname] = xr.DataArray(
                    interpolated_data,
                    coords=[
                        ("latitude", np.flipud(target_lats)),
                        ("longitude", np.rad2deg(target_lons)),
                    ],
                )

        new_ds = new_ds.sortby("latitude", ascending=True)

        return new_ds

    def convert_units(self, depth_unit: str = "m"):
        """Convert depth coordinates and velocity fields from meters to kilometers.

        Args:
            depth_unit: The unit of the source data. Use ``'m'`` to trigger
                conversion to km; ``'km'`` is a no-op. Any other value raises.
        """
        if depth_unit == "m":
            self.dataset = self.dataset.assign_coords(depth=self.dataset.depth / 1000.0)
            if "VSV" in self.dataset:
                self.dataset = self.dataset.assign(
                    VSV=lambda x: self.dataset.VSV / 1000.0
                )
            if "VSH" in self.dataset:
                self.dataset = self.dataset.assign(
                    VSH=lambda x: self.dataset.VSH / 1000.0
                )
        elif depth_unit != "km":
            raise ValueError("depth_unit must be 'm' or 'km'")

    def convert_longitude(self):
        """Remap longitude coordinates from [-180, 180] to [0, 360] and sort."""
        lon_0_360 = [
            lon + 360.0 if lon < 0 else lon
            for lon in self.dataset.longitude.to_numpy()
        ]
        self.dataset = self.dataset.assign_coords(longitude=lon_0_360)
        self.dataset = self.dataset.sortby("longitude", ascending=True)

    def assign_names(self):
        """Rename coordinate variables to standard latitude, longitude, and depth names."""

        name_map = {
            "latitude": ["lat", "latitude", "y"],
            "longitude": ["lon", "longitude", "x"],
            "depth": ["depth", "z", "level"],
        }

        for key in name_map.keys():
            found = False
            for value in name_map[key]:
                if value.lower() in self.dataset.variables:
                    self.dataset = self.dataset.rename({value: key})
                    found = True
                    break

            if not found:
                raise ValueError(f"Could not detect variable name for {key}")

    def save_model(self, file_path: str):
        """Save the current dataset to a NetCDF file path."""

        if file_path is None:
            raise ValueError("File path to save model must be given")
        self.dataset.to_netcdf(file_path)

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
        """Plot a single data variable on a global PlateCarree map.

        Args:
            varname: Name of the variable in the dataset to plot.
            depth: Depth level to visualize for 3-D variables. Accepts an
                integer index, a float value (nearest depth is selected), or
                ``None`` to use the first depth level.
            ax: Existing matplotlib axes to draw on. A new figure is created
                when ``None``.
            cmap: Matplotlib colormap name.
            vmin: Lower bound for the color scale.
            vmax: Upper bound for the color scale.
            title: Axes title. Defaults to ``"<varname> (depth=<depth>)"``.
            show_colorbar: Whether to attach a colorbar to the plot.
            figsize: Figure size passed to ``plt.subplots``.

        Returns:
            Tuple of ``(fig, ax)``.
        """
        if varname not in self.dataset:
            raise ValueError(f"{varname} not found in dataset")

        da = self.dataset[varname]

        if "depth" in da.dims:
            if depth is None:
                da_sel = da.isel(depth=0)
                depth_label = str(da.depth.values[0])
            elif isinstance(depth, (int, np.integer)):
                da_sel = da.isel(depth=int(depth))
                depth_label = str(da.depth.values[int(depth)])
            else:
                # Nearest-depth selection for float values
                idx = int(np.argmin(np.abs(da.depth.values - float(depth))))
                da_sel = da.isel(depth=idx)
                depth_label = str(da.depth.values[idx])
        else:
            da_sel = da
            depth_label = None

        da_sel = da_sel.sortby("latitude")

        lats = da_sel["latitude"].values[::-1]
        lons = da_sel["longitude"].values
        data = da_sel.values

        created_fig = False
        if ax is None:
            fig, ax = plt.subplots(
                figsize=figsize, subplot_kw={"projection": ccrs.PlateCarree()}
            )
            created_fig = True
        else:
            fig = ax.figure

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
        """Plot every data variable in the dataset at a given depth.

        Args:
            depth: Depth level forwarded to :meth:`plot_variable`.
            save_dir: Directory to write PNG files into. Created automatically
                if it does not exist. When ``None``, figures are not saved.
            cmap: Matplotlib colormap name applied to all plots.
            vmin: Shared lower bound for the color scale.
            vmax: Shared upper bound for the color scale.
            figsize: Figure size passed to each :meth:`plot_variable` call.
            show: Call ``plt.show()`` after each plot. Set to ``False`` when
                batch-saving to avoid blocking.

        Returns:
            List of ``(varname, filepath_or_None)`` tuples.
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
