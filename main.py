from Dataset import Dataset
from MergeMethods import MergeMethods
regional_mod = Dataset("Data/netcdf/CANVAS_15-60s_400km.nc", Dataset.REGIONAL, "Data/spline.par_26", 'm')
global_mod = Dataset("Data/netcdf/semucb-2014-ucb-vs.nc", Dataset.GLOBAL, depthUnits='m')

# print("Displaying Regional Maps")
# regional_mod.plot_all_variables()
# print("End of Regional Map Displays")
# print("Displaying Global Maps")
# global_mod.plot_all_variables()
# print("End of Global Map Displays")

merger = MergeMethods(regional_mod, global_mod, "conf.yaml")
print("Created merge methods instance")

print("Merging global and regional")
merged_mod = merger.merge()

print("Displaying Merged Maps")
merged_mod.plot_all_variables()
print("End of Merged Map Displays")