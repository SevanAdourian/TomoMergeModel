from __future__ import annotations

import multiprocessing
import cmcrameri.cm as cm

from Dataset import Dataset
from MergeMethods import MergeMethods
from ConfigParams import ConfigParams


def main() -> None:
    depths = [150]
    global_mod = Dataset(
        "Data/netcdf/glad-m25-vs-0.0-n4.nc",
        Dataset.GLOBAL,
        depthUnits="km",
    )
    regional_mod = Dataset(
        filePath="Data/netcdf/alaska.nc",
        modelType=Dataset.REGIONAL,
        depths=depths,
        depthUnits="km",
        globalModel=global_mod,
    )
    configParams = ConfigParams(
        239, 80, 80, (197.5, 230.5), (52.65, 71.55), "spherical"
    )

    merger = MergeMethods(regional_mod, global_mod, configParams, "Vs", "vsv")
    print("Created merge methods instance")

    print("Merging global and regional")
    merged_mod = merger.merge()

    print("Displaying Merged Maps")
    cmap_seismic = cm.vik
    for depth in merged_mod.getDataset().depth.values:
        merged_mod.plot_all_variables(
            depth=depth,
            cmap=cmap_seismic
        )
    print("End of Merged Map Displays")


if __name__ == "__main__":
    main()
