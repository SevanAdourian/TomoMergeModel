from __future__ import annotations

import multiprocessing

from Dataset import Dataset
from MergeMethods import MergeMethods


def main() -> None:
	global_mod = Dataset(
		"Data/netcdf/glad-m25-vs-0.0-n4.nc",
		Dataset.GLOBAL,
		depthUnits="km",
	)
	regional_mod = Dataset(
		"Data/netcdf/alaska.nc",
		Dataset.REGIONAL,
		"Data/spline.par",
		"km",
		globalModel=global_mod,
	)

	merger = MergeMethods(regional_mod, global_mod, "conf.yaml")
	print("Created merge methods instance")

	print("Merging global and regional")
	merged_mod = merger.merge()

	print("Displaying Merged Maps")
	for depth in merged_mod.getDataset().depth.values:
		merged_mod.plot_all_variables(depth=depth, save_dir="alaska_plots")
	print("End of Merged Map Displays")


if __name__ == "__main__":
	main()
