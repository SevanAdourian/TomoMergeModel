# Quickstart

This guide shows the shortest path to run a merge with one global and one regional model.

## Prerequisites

- Python environment from `environment.yaml` (recommended) or `requirements.txt`.
- One global NetCDF model.
- One regional NetCDF model.

## Minimal Merge Example

```python
from Dataset import Dataset
from MergeMethods import MergeMethods
from ConfigParams import ConfigParams


def main() -> None:
    depths = [150]

    global_mod = Dataset(
        file_path="Data/netcdf/semucb-2014-ucb-vs.nc",
        model_type=Dataset.GLOBAL,
        depth_units="km",
    )

    regional_mod = Dataset(
        file_path="Data/netcdf/alaska.nc",
        model_type=Dataset.REGIONAL,
        depths=depths,
        depth_units="km",
        global_model=global_mod,
    )

    conf = ConfigParams(
        reg_lmax=239,
        win_lmax=30,
        lon_bounds=(197.5, 230.5),
        lat_bounds=(52.65, 71.55),
        win_type="spherical",
        mask_mode="bounds",
    )

    merger = MergeMethods(
        model_one=regional_mod,
        model_two=global_mod,
        config_params=conf,
        regional_variable="Vs",
        global_variable="vsv",
    )

    merged = merger.merge()
    merged.save_model("merged_output.nc")


if __name__ == "__main__":
    main()
```

## What Happens Internally

1. The regional model is interpolated to match the global grid shape/resolution.
2. Both models are transformed into spherical harmonics.
3. A geographic window is generated from either bounds or continent geometry.
4. Spectral blending is applied using a logistic transition across harmonic degree.
5. The blended result is spatially composited with the global model.
6. Slices are concatenated into one merged output dataset.

## Common First-Run Errors

- Wrong variable names: ensure `regional_variable` and `global_variable` exist in each dataset.
- Missing depth unit conversion: set `depth_units` correctly as `"m"` or `"km"`.
- Incorrect mask bounds: in bounds mode, longitudes must be in [0, 360] and latitudes in [-90, 90].
