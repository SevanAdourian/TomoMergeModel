# Data Requirements

This page describes the coordinate, unit, and variable assumptions expected by the package.

## Required Coordinates

Input datasets must provide equivalent latitude, longitude, and depth coordinates.

Accepted aliases are normalized automatically by `Dataset.assign_names`:

- Latitude: `lat`, `latitude`, `y`
- Longitude: `lon`, `longitude`, `x`
- Depth: `depth`, `z`, `level`

If none of the aliases are found for a required coordinate, initialization fails.

## Longitude Convention

`Dataset.convert_longitude` remaps longitudes from [-180, 180] into [0, 360] and sorts them.

Implications:

- Bounds in `ConfigParams.lon_bounds` must use [0, 360].
- Regional and global datasets should be interpreted with the same convention after loading.

## Depth Units

Set the `depth_units` parameter on `Dataset`:

- `"m"`: depth is converted to km and velocity fields `VSV`/`VSH` are scaled by 1000.
- `"km"`: no conversion is performed.

Any other value raises a validation error.

## Regional Dataset Interpolation Behavior

For regional models, the class interpolates to the global target grid using:

- `RectSphereBivariateSpline` for spatial interpolation.
- Optional linear depth interpolation if explicit target depths are provided.

Outside the regional footprint, values remain `NaN` before downstream masking/blending.

## Variable Naming

`MergeMethods` validates both variables passed to the constructor:

- `regional_variable` must exist in the regional dataset.
- `global_variable` must exist in the global dataset.

Use dataset inspection before merging if variable names differ by model source.
