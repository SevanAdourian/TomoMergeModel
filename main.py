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
		"Data/netcdf/CANVAS_15-60s_400km.nc",
		Dataset.REGIONAL,
		"Data/spline.par",
		"m",
		globalModel=global_mod,
	)

	merger = MergeMethods(regional_mod, global_mod, "conf.yaml")
	print("Created merge methods instance")

	print("Merging global and regional")
	merged_mod = merger.merge()

	print("Displaying Merged Maps")
	merged_mod.plot_all_variables()
	print("End of Merged Map Displays")


if __name__ == "__main__":
	# Required on Windows when using multiprocessing (directly or indirectly)
	multiprocessing.freeze_support()
	main()
