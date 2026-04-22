# ConfigParams API Reference

`ConfigParams` encapsulates validated parameters used by `MergeMethods` to construct windows and control harmonic resolution.

## Constructor

```python
ConfigParams(
    reg_lmax: int,
    win_lmax: int,
    lon_bounds: tuple[float, float],
    lat_bounds: tuple[float, float],
    win_type: str,
    mask_mode: str = "bounds",
    mask_continents: list[str] | tuple[str, ...] | None = None,
    mask_target: str = "land",
    mask_resolution: str = "110m",
)
```

## Parameter Details

- `reg_lmax`: Maximum spherical harmonic degree used for regional expansion/padding. Must be non-negative.
- `win_lmax`: Harmonic bandwidth used in spherical window construction. Must be non-negative.
- `lon_bounds`: `(min_lon, max_lon)` for bounds mask mode. Must be length 2, sorted ascending, in [0, 360].
- `lat_bounds`: `(min_lat, max_lat)` for bounds mask mode. Must be length 2, sorted ascending, in [-90, 90].
- `win_type`: Either `"spherical"` or `"rectangular"`.
- `mask_mode`: Either `"bounds"` or `"continent"`.
- `mask_continents`: Optional continent names used when `mask_mode="continent"`.
- `mask_target`: `"land"` or `"ocean"` selection for continent masking.
- `mask_resolution`: Natural Earth resolution string (for example `"110m"`).

## Attributes Set by Constructor

Common attributes:

- `reg_lmax`
- `win_lmax`
- `mask_mode`
- `mask_target`
- `mask_resolution`
- `mask_continents`
- `win_type`

Bounds-specific attributes (set when `mask_mode="bounds"`):

- `lon_min_mask`
- `lon_max_mask`
- `lat_min_mask`
- `lat_max_mask`

Continent mode sets the four bounds attributes to `None`.

## Validation Rules

The constructor raises `ValueError` for invalid inputs, including:

- Negative `reg_lmax` or `win_lmax`.
- Invalid `mask_mode`, `mask_target`, or `win_type`.
- Invalid bounds tuple lengths/order/range when in bounds mode.

## Example

```python
conf = ConfigParams(
    reg_lmax=239,
    win_lmax=30,
    lon_bounds=(197.5, 230.5),
    lat_bounds=(52.65, 71.55),
    win_type="spherical",
    mask_mode="bounds",
)
```
