# TomoMergeModel Documentation

This documentation is organized from high-level concepts to low-level API details.

## Documentation Map

### High-level

- [Getting Started](getting-started/quickstart.md)
- [Data Requirements](getting-started/data-requirements.md)
- [Merge Pipeline](concepts/merge-pipeline.md)

### API Reference (Low-level)

- [Dataset](api/Dataset.md)
- [ConfigParams](api/ConfigParams.md)
- [MergeMethods](api/MergeMethods.md)

## Recommended Reading Order

1. Start with [Getting Started](getting-started/quickstart.md).
2. Read [Merge Pipeline](concepts/merge-pipeline.md) to understand how blending works.
3. Use the pages in [API Reference](api/Dataset.md) while implementing scripts.

## Package Architecture

The package revolves around three classes:

1. `Dataset` handles model loading, normalization, interpolation, plotting, and saving.
2. `ConfigParams` defines validated merge and mask settings.
3. `MergeMethods` performs spherical-harmonic conversion, window construction, and depth-by-depth merge execution.
