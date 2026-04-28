from __future__ import annotations

import cmcrameri.cm as cm

from geostitch.Dataset import Dataset
from geostitch.MergeMethods import MergeMethods
from geostitch.ConfigParams import ConfigParams


def main() -> None:
    """Run a merge with Alaska and display the results.

    Workflow:
    1. Load the global GLAD-M25 model.
    2. Merge the Alaska regional model into the global model at the target depths.
    4. Plot all variables of the final merged model at each depth.
    """
    depths = [150]
    global_mod = Dataset(
        "../data/glad-m25-vs-0.0-n4.nc",
        Dataset.GLOBAL,
        depth_units="km",
    )
    reg_alaska = Dataset(
        file_path="../data/alaska.nc",
        model_type=Dataset.REGIONAL,
        depths=depths,
        depth_units="km",
        global_model=global_mod,
    )
    configParams = ConfigParams(
        200, 60, (197.5, 230.5), (52.65, 71.55),
        win_type="spherical",
        # blend_mode="adaptive",
        preserve_global_low_lmax=30,
        reg_noise_floor=1e-11,
        glo_noise_floor=1e-11,
        blend_lcut=80,
        blend_delta=10.0
    )

    merger = MergeMethods(reg_alaska, global_mod, configParams, "Vs", "vsv")
    print("Merging Alaska regional model into global...")
    merged_mod = merger.merge()
   
    
    print("Displaying merged maps...")
    cmap_seismic = cm.vik_r
    for depth in merged_mod.getDataset().depth.values:
        merged_mod.plot_all_variables(
            depth=depth,
            cmap=cmap_seismic
        )
    print("Done.")


if __name__ == "__main__":
    main()
