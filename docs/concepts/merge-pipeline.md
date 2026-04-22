# Merge Pipeline

This package uses a two-stage blending strategy at each depth slice.

## Stage 1: Spectral Blending

In `MergeMethods.apply_window`, regional and global spherical harmonic coefficients are blended per degree:

- Global-dominated at low degree (long wavelengths).
- Regional-dominated at high degree (short wavelengths).

The transition is controlled by a logistic weight:

- `lcut`: central transition degree.
- `delta`: transition sharpness.

The resulting coefficient field is expanded back to a spatial grid (`blended_grid`).

## Stage 2: Spatial Compositing

A regional window grid (`reg_win_grid`) controls where the blended result replaces the global result:

- Outside window: mostly global model.
- Inside window: mostly blended model.

Composite equation implemented in `apply_window`:

- `merged = global * (1 - window) + blended * window`

## Window Construction

`MergeMethods.create_window` supports:

1. `spherical` window:
   - Uses Slepian/multi-taper window functions (`pyshtools.SHWindow.from_mask`).
   - Produces smooth edge transitions.
2. `rectangular` window:
   - Binary mask transformed through SH expansion/reconstruction.
   - Produces sharper boundaries.

Window mask sources are configured with `ConfigParams.mask_mode`:

- `bounds`: axis-aligned lat/lon box.
- `continent`: Natural Earth country polygons filtered by continent set.

## Depth Processing Logic

`MergeMethods.process_slice` behavior:

- If depth is within regional depth coverage: run full merge.
- If depth is outside regional depth coverage: pass through global slice.

`MergeMethods.merge` loops through target depths and concatenates slices into one output dataset.
