# Dataset API Reference

`Dataset` wraps an `xarray.Dataset` and standardizes coordinate naming, depth units, longitude convention, optional interpolation, plotting, and persistence.

## Constants

- `Dataset.REGIONAL = 0`
- `Dataset.GLOBAL = 1`
- `Dataset.RADIUS_IN_METERS = 6371000`

## Constructor

```python
Dataset(
    file_path: str,
    model_type: int,
    depth_units: str = None,
    xr_dataset: xr.Dataset = None,
    global_model: Dataset = None,
    depths: list[int] = None,
)
```

## Initialization Behavior

1. Validates `file_path` and `model_type`.
2. Loads from `xr_dataset` when `file_path==""`, otherwise opens NetCDF from disk.
3. Normalizes coordinate names (`assign_names`).
4. Optionally converts depth/velocity units (`convert_units`).
5. For regional models, computes interpolation targets using `global_model` and interpolates (`interpolate_model`).
6. Converts longitudes to [0, 360] and sorts (`convert_longitude`).

## Key Methods

### Metadata / Accessors

- `getFilePath() -> str`
- `setFilePath(path: str)`
- `getDataset() -> xr.Dataset`
- `getModelType() -> int`
- `setModelType(model_type: int)`

### Coordinate and Unit Normalization

- `assign_names()`
- `convert_units(depth_unit: str = "m")`
- `convert_longitude()`

### Interpolation

- `getInterpolationParameters(depths, global_model) -> tuple`
- `interpolate_slice(...)` (static)
- `interpolate_model(ds, target_lats, target_lons, target_depths=None)` (static)

### Output and Visualization

- `save_model(file_path: str)`
- `plot_variable(varname, depth=None, ...) -> (fig, ax)`
- `plot_all_variables(depth=None, save_dir=None, ...) -> list[(varname, path)]`

## Data Assumptions

- Coordinate aliases are accepted and renamed to `latitude`, `longitude`, `depth`.
- Longitude is normalized to [0, 360].
- For regional interpolation, values outside spatial coverage remain `NaN`.

## Common Failure Cases

- Missing required coordinates after alias scan.
- Missing `global_model` for regional datasets when `depths is None`.
- Invalid `depth_units` values (must be `"m"` or `"km"`).
- Invalid file path in `setFilePath`.

## Example

```python
global_mod = Dataset(
    file_path="Data/netcdf/semucb-2014-ucb-vs.nc",
    model_type=Dataset.GLOBAL,
    depth_units="km",
)

regional_mod = Dataset(
    file_path="Data/netcdf/alaska.nc",
    model_type=Dataset.REGIONAL,
    depths=[150],
    depth_units="km",
    global_model=global_mod,
)
```
