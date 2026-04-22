from .Dataset import Dataset
from .ConfigParams import ConfigParams

import concurrent.futures
import numpy as np
import pyshtools
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


class MergeMethods:
    def __init__(
        self,
        model_one: Dataset,
        model_two: Dataset,
        config_params: ConfigParams,
        regional_variable: str,
        global_variable: str,
    ):
        """Initialize with one regional and one global dataset, config parameters, and variable names."""

        # input validation
        if model_one is None or model_two is None:
            raise ValueError("Null model(s) provided")

        # checking to see that there is one regional and one global model, and assigning them
        model_one_type = model_one.getModelType()
        model_two_type = model_two.getModelType()

        regional_mod = None
        global_mod = None
        if model_one_type == Dataset.REGIONAL and model_two_type == Dataset.GLOBAL:
            regional_mod = model_one
            global_mod = model_two
        elif model_one_type == Dataset.GLOBAL and model_two_type == Dataset.REGIONAL:
            regional_mod = model_two
            global_mod = model_one

        # raise an error for two global or two regional models
        if regional_mod is None or global_mod is None:
            raise ValueError("Constructor requires one REGIONAL and one GLOBAL Dataset")

        self.regional_model = regional_mod
        self.global_model = global_mod
        self.merge_model = None

        if config_params is None:
            raise ValueError(
                "Must provide a ConfigParmas class containing merging variables"
            )

        self.conf = config_params

        if regional_variable not in list(self.regional_model.getDataset().data_vars):
            raise ValueError(
                f"Could not find regional variable {regional_variable} in regional model"
            )
        if global_variable not in list(self.global_model.getDataset().data_vars):
            raise ValueError(
                f"Could not find global variable {global_variable} in global model"
            )
        self.regional_variable = regional_variable
        self.global_variable = global_variable

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

    def setConf(self, config_params):
        """Set the configuration parameters after validation."""

        if config_params is None:
            raise ValueError("Must provide a configuration file for merging variables")

        self.conf = config_params

    # merge function and helper functions here
    def convert_to_spherical_harmonics(self, zmesh, reg_lmax):
        """Convert a 2D grid to spherical harmonics coefficients up to specified lmax."""

        # Work on a copy to avoid mutating caller input
        zmesh_clean = np.array(zmesh, copy=True)

        # Set nan values to the mean of the region
        mask = np.isnan(zmesh_clean)
        if np.any(mask):
            if np.all(mask):
                raise ValueError("Cannot convert all-NaN grid to spherical harmonics")
            zmesh_clean[mask] = np.nanmean(zmesh_clean)

        # Create the grid and coefficients
        grid = pyshtools.SHGrid.from_array(zmesh_clean, grid="DH")
        clm = pyshtools.SHGrid.expand(grid)
        clm = clm.pad(reg_lmax)  # Pad to match global clm

        # Trim coefficients in SPH to match number used in regional model
        grid = pyshtools.SHCoeffs.expand(clm, lmax=reg_lmax)

        return grid, clm

    def create_window(self, reg_field, global_clm):
        """Build a normalized spherical harmonic window function over the regional area.

        Constructs a binary geographic mask, then lifts it into either a
        multi-taper spherical harmonic window (``win_type='spherical'``) or a
        sharp rectangular window (``win_type='rectangular'``). The result is a
        normalised grid (values in [0, 1]) suitable for spatially blending the
        regional and global models.

        Args:
            reg_field: Regional dataset slice (xr.Dataset) used to define
                the spatial extent of the window.
            global_clm: Global SHCoeffs object whose ``lmax`` is used to pad
                the window coefficients.

        Returns:
            pyshtools SHGrid representing the normalized window.
        """

        reg_zmesh_mask = self.build_binary_mask(
            reg_field=reg_field.copy(),
            mask_mode=getattr(self.conf, "mask_mode", "bounds"),
        )
        # pdb.set_trace()

        if self.conf.win_type == "spherical":
            # Multi-taper spherical harmonic window — smooth, tapered edges
            reg_win = pyshtools.SHWindow.from_mask(
                reg_zmesh_mask, lwin=self.conf.win_lmax
            )
            reg_win_clm = pyshtools.SHWindow.to_shcoeffs(reg_win, 0)
            reg_win_clm.pad(global_clm.lmax)

            reg_win_energy = (reg_win.to_shgrid(0).to_array()) ** 2
            for i in range(1, self.conf.win_lmax):
                reg_win_energy += (reg_win.to_shgrid(i).to_array()) ** 2
            reg_win_energy_grid = pyshtools.SHGrid.from_array(reg_win_energy)

            # Normalise to [0, 1] for use as a blend weight
            valmax = np.amax(reg_win_energy_grid.data)
            reg_win_energy_grid = reg_win_energy_grid / valmax

            reg_win_energy_clm = pyshtools.SHGrid.expand(reg_win_energy_grid)
            reg_win_energy_clm = reg_win_energy_clm.pad(self.conf.reg_lmax)
            reg_win_energy_grid = pyshtools.SHCoeffs.expand(reg_win_energy_clm)

        elif self.conf.win_type == "rectangular":
            # Sharp binary window — no edge smoothing
            reg_win_grid = pyshtools.SHGrid.from_array(reg_zmesh_mask)
            reg_win_clm = pyshtools.SHGrid.expand(reg_win_grid)
            reg_win_clm_pad = reg_win_clm.pad(global_clm.lmax)
            reg_win_energy_grid = pyshtools.SHCoeffs.expand(reg_win_clm_pad)
        else:
            raise ValueError(
                f"Unsupported win_type '{self.conf.win_type}'. Expected 'spherical' or 'rectangular'."
            )

        # reg_win_energy_grid.data = np.flipud(reg_win_energy_grid.data)

        return reg_win_energy_grid

    def build_binary_mask(self, reg_field, mask_mode="bounds"):
        """Build a geographic binary mask for the provided regional field."""

        reg_field = reg_field.sortby("longitude")
        lon = reg_field["longitude"].values
        lat = reg_field["latitude"].values
        lon2d, lat2d = np.meshgrid(lon, lat)

        mode = (mask_mode or "bounds").lower()
        if mode == "bounds":
            mask = np.where(
                (lon2d >= self.conf.lon_min_mask)
                & (lon2d <= self.conf.lon_max_mask)
                & (lat2d >= self.conf.lat_min_mask)
                & (lat2d <= self.conf.lat_max_mask),
                1,
                0,
            )
            return mask.astype(np.uint8)

        if mode == "continent":
            from cartopy.io import shapereader
            from shapely.geometry import Point
            from shapely.ops import unary_union
            from shapely.prepared import prep

            resolution = getattr(self.conf, "mask_resolution", "110m")
            target = getattr(self.conf, "mask_target", "land").lower()
            continents = getattr(self.conf, "mask_continents", []) or []

            shp = shapereader.natural_earth(
                resolution=resolution,
                category="cultural",
                name="admin_0_countries",
            )
            reader = shapereader.Reader(shp)
            selected = set(c.lower() for c in continents) if continents else None
            geoms = []
            for rec in reader.records():
                cont = (rec.attributes.get("CONTINENT") or "").lower()
                if selected is None or cont in selected:
                    geoms.append(rec.geometry)

            if not geoms:
                raise ValueError("No geometries found for requested continents")

            geom = prep(unary_union(geoms))
            lon180 = ((lon2d + 180.0) % 360.0) - 180.0
            inside = np.fromiter(
                (
                    geom.contains(Point(float(x), float(y)))
                    for x, y in zip(lon180.ravel(), lat2d.ravel())
                ),
                dtype=bool,
                count=lon180.size,
            ).reshape(lon2d.shape)

            if target == "ocean":
                inside = ~inside

            return inside.astype(np.uint8)

        raise ValueError(
            f"Unsupported mask_mode '{mode}'. Expected 'bounds' or 'continent'."
        )

    # def apply_window(self, global_grid, global_clm, reg_grid, reg_win_energy_grid):
    #     """Apply windowing mask to regional grid and merge with global model."""

    #     # - Multiply grid by mask
    #     reg_grid_masked = reg_grid * reg_win_energy_grid
    #     reg_clm_masked = pyshtools.SHGrid.expand(reg_grid_masked)

    #     # Global grid 

    #     # Sum regional spectrum inside box (masked) with global spectra outside box (unmasked)
    #     summed_clm = reg_clm_masked + global_clm

    #     # Plot map and spectra for summed_clm

    #     summed_grid = pyshtools.SHCoeffs.expand(summed_clm)
    #     return summed_grid

    def apply_window(self, global_grid, global_clm, reg_grid, reg_win_grid,
                     lcut=60, delta=5):
        """Merge global and regional models via spectral blending and geographic windowing.

        Performs a two-stage merge:
        1. **Spectral blending** — regional and global SH coefficients are
           mixed per degree using a logistic weight centred at ``lcut``,
           favouring global structure at long wavelengths and regional detail
           at short wavelengths.
        2. **Spatial blending** — the spectrally blended grid is composited
           with the original global grid using the normalised window
           ``reg_win_grid`` as a per-pixel blend weight.

        Args:
            global_grid: Global model as a pyshtools SHGrid.
            global_clm: Global model SH coefficients (pyshtools SHCoeffs).
            reg_grid: Regional model as a pyshtools SHGrid.
            reg_win_grid: Normalised window grid (values in [0, 1]) defining
                the regional footprint; output of :meth:`create_window`.
            lcut: Spherical harmonic degree at which the logistic blend
                transitions from global-dominated to regional-dominated.
            delta: Steepness of the logistic transition (lower = sharper).

        Returns:
            Merged model as a pyshtools SHGrid.
        """
        
        lmax = global_clm.lmax

        # Expand regional grid to SH coefficients and match lmax to global
        reg_clm = reg_grid.expand().pad(lmax)

        G = global_clm.to_array()
        R = reg_clm.to_array()

        # Logistic weight: favours global at low degrees, regional at high degrees
        degrees = np.arange(lmax + 1)
        w = 1 / (1 + np.exp(-(degrees - lcut) / delta))

        # Blend SH coefficients per degree
        C = np.zeros_like(G)
        for l in range(lmax + 1):
            C[:, l, :l+1] = (1 - w[l]) * G[:, l, :l+1] + w[l] * R[:, l, :l+1]

        blended_clm = pyshtools.SHCoeffs.from_array(C)
        blended_grid = blended_clm.expand()

        # Composite: global outside the window, blended inside
        merged_grid = global_grid * (1 - reg_win_grid) + blended_grid * reg_win_grid
        merged_clm = merged_grid.expand().pad(lmax)

        return merged_grid, reg_clm, blended_clm, merged_clm

    def process_slice(self, depth):
        """Process and merge regional and global models at a single depth slice.

        When the requested depth falls within the regional model's depth range,
        the regional field is windowed and blended into the global model.
        Outside that range, the global model is returned unchanged.

        Args:
            depth: Depth value (km) to process.

        Returns:
            xr.DataArray with dimensions ``(latitude, longitude, depth)``.
        """

        # Reading in global tomography model
        zmesh_global = self.reshape_field(
            self.global_model.getDataset(), depth, self.global_variable
        )

        global_grid, global_clm = self.convert_to_spherical_harmonics(
            zmesh_global, self.conf.reg_lmax
        )

        regional_depths = np.asarray(self.regional_model.depths, dtype=float)
        if regional_depths.size == 0:
            raise ValueError("Regional model depths are empty")

        regional_min_depth = float(np.min(regional_depths))
        regional_max_depth = float(np.max(regional_depths))

        if regional_min_depth <= float(depth) <= regional_max_depth:
            # Depth is within regional coverage — perform the merge
            zmesh_regional = self.reshape_field(
                self.regional_model.getDataset(), depth, self.regional_variable
            )
            reg_grid, reg_clm = self.convert_to_spherical_harmonics(
                zmesh_regional, self.conf.reg_lmax
            )
            reg_win_energy_grid = self.create_window(
                self.regional_model.getDataset().sel(
                    depth=float(depth)
                ),
                global_clm,
            )
            summed_grid, reg_clm, blended_clm, merged_clm = self.apply_window(
                global_grid, global_clm, reg_grid, reg_win_energy_grid
            )
            self.plot_combined_spectra(
                spectra_series=[
                    ("regional", reg_clm),
                    ("global", global_clm),
                    ("blended", blended_clm),
                    ("merged", merged_clm),
                ],
                file_name=f"spectra_{float(depth):g}km.png",
                title=f"Power Spectra Comparison ({float(depth):g} km)",
            )
        else:
            # Depth is below regional coverage — pass through global model unchanged
            summed_grid = global_grid

        return self.write_model(depth=depth, grid=summed_grid)

    def write_model(self, depth, grid):
        """
        Write out merged model as an xarray DataArray with latitude, longitude, depth dimensions.
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

        fig, ax = plt.subplots(
            1, 1, subplot_kw={"projection": ccrs.PlateCarree()}
        )
        grid_object.plot(ax=ax, colorbar="right", cb_label="Power", show=False)
        fig.legend(loc="upper right")
        fig.savefig(file_name, dpi=400)
        plt.close(fig)
        return

    def plot_combined_spectra(self, spectra_series, file_name, title=None):
        """Plot multiple spectra on one axis and save to file.

        Args:
            spectra_series: Iterable of ``(label, clm_object)`` where each
                ``clm_object`` is a pyshtools SHCoeffs instance.
            file_name: Output image path.
            title: Optional figure title.
        """

        fig, ax = plt.subplots(1, 1, figsize=(8, 6))

        for label, clm_object in spectra_series:
            spectrum = clm_object.spectrum()
            degrees = np.arange(spectrum.shape[0])

            # Skip degree zero for log plotting; it dominates and hides detail.
            ax.semilogy(degrees, spectrum, label=label)

        ax.set_xlabel("Spherical harmonic degree (l)")
        ax.set_ylabel("Power")
        if title is not None:
            ax.set_title(title)
        ax.grid(True, which="both", linestyle=":", alpha=0.5)
        ax.legend(loc="best")

        fig.tight_layout()
        fig.savefig(file_name, dpi=400)
        plt.close(fig)
        return

    def reshape_field(self, lon_lat_field, depth, varname):
        """Extract and reshape a 2D field at a given depth for spherical harmonic expansion."""

        zmesh = lon_lat_field[varname].sel(depth=float(depth))
        zmesh = zmesh.sortby("longitude")
        # zmesh.assign_coords(latitude=zmesh.latitude[::-1])
        # zmesh = np.flipud(zmesh.values)
        return zmesh.values

    def merge(self):
        """Merge regional and global models across all depth slices.

        Iterates over every depth in the regional model, calling
        :meth:`process_slice` for each, then concatenates the results into a
        single xr.Dataset.

        Returns:
            Merged model wrapped in a :class:`Dataset` instance (GLOBAL type).
        """

        depths = self.regional_model.depths

        # On Windows, process-based multiprocessing uses "spawn" which re-imports modules.
        # Using a process Pool with bound methods/self can also hang due to pickling.
        # A thread pool avoids both issues and is reliable here.
        merged_arrays = []
        for depth in depths:
            depth_val = float(depth)
            merged_arrays.append(self.process_slice(depth_val))

        # concat returns a DataArray; convert to Dataset immediately
        merged_all_arrays = xr.concat(merged_arrays, dim="depth")
        self.merge_model = merged_all_arrays.to_dataset(name=self.global_variable)

        self.merge_model = self.merge_model.sortby("depth")
        self.merge_model = Dataset(
            "", Dataset.GLOBAL, xr_dataset=self.merge_model, depth_units="km"
        )

        return self.merge_model
