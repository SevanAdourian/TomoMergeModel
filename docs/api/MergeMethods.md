# MergeMethods API Reference

`MergeMethods` coordinates the end-to-end merge of one regional model into one global model, depth by depth.

## Constructor

```python
MergeMethods(
    model_one: Dataset,
    model_two: Dataset,
    config_params: ConfigParams,
    regional_variable: str,
    global_variable: str,
)
```

## Constructor Validation

The constructor enforces:

- Exactly one `Dataset.REGIONAL` and one `Dataset.GLOBAL` dataset.
- Non-null configuration object.
- `regional_variable` exists in the regional dataset.
- `global_variable` exists in the global dataset.

If validation fails, it raises `ValueError`.

## Core Merge Methods

### `merge()`

Runs the merge across all depths in `regional_model.depths`:

1. Calls `process_slice(depth)` for each depth.
2. Concatenates returned arrays along `depth`.
3. Converts to `xr.Dataset` named with `global_variable`.
4. Wraps and returns a new `Dataset` object (global type).

Notes:

- Current implementation processes slices sequentially.
- Output depth coordinate is sorted ascending.

### `process_slice(depth)`

Per-depth workflow:

1. Extract global field for depth.
2. Convert global field to SH grid/coefficients.
3. If depth is within regional coverage:
   - Extract regional field.
   - Convert regional field to SH grid/coefficients.
   - Build blending window with `create_window`.
   - Blend and composite with `apply_window`.
   - Plot spectra comparison.
4. Else, pass through global slice.
5. Serialize to `xr.DataArray` via `write_model`.

### `apply_window(global_grid, global_clm, reg_grid, reg_win_grid, lcut=60, delta=5)`

Two-stage blending implementation:

- Spectral blending in SH space using logistic degree weights.
- Spatial compositing with window grid.

Returns a tuple:

1. `merged_grid`
2. `reg_clm`
3. `blended_clm`
4. `merged_clm`

## Window and Mask Construction

### `create_window(reg_field, global_clm)`

Builds a normalized window grid from a binary geographic mask.

Window options:

- `spherical`: multi-taper SH window (smooth edges).
- `rectangular`: binary mask SH reconstruction (sharper edges).

### `build_binary_mask(reg_field, mask_mode="bounds")`

Mask modes:

- `bounds`: simple axis-aligned lon/lat bounds from `ConfigParams`.
- `continent`: polygon-based mask from Natural Earth via Cartopy/Shapely.

Continent mode supports:

- `mask_continents`
- `mask_target` (`land` or `ocean`)
- `mask_resolution`

## Utility and Plotting Methods

- `convert_to_spherical_harmonics(zmesh, reg_lmax)`
- `reshape_field(lon_lat_field, depth, varname)`
- `write_model(depth, grid)`
- `plot_map_and_spectra(grid_object, clm_object, file_name)`
- `plot_combined_spectra(spectra_series, file_name, title=None)`

## Typical Usage

```python
merger = MergeMethods(
    model_one=regional_mod,
    model_two=global_mod,
    config_params=conf,
    regional_variable="Vs",
    global_variable="vsv",
)
merged = merger.merge()
```

## Operational Tips

- Keep `regional_variable` and `global_variable` explicit; model naming conventions often differ.
- Start with a moderate `reg_lmax`/`win_lmax` to limit runtime during tuning.
- Validate masking settings first (bounds vs continent) before large depth batches.
