import pdb

from Dataset import Dataset
from MergeMethods import MergeMethods
regional_mod = Dataset("Data/netcdf/CANVAS_15-60s_400km.nc", Dataset.REGIONAL, "Data/spline.par_26", 'm')
# print(regional_mod.getDataset())
# pdb.set_trace()
global_mod = Dataset("Data/netcdf/glad-m25-vs-0.0-n4.nc", Dataset.GLOBAL, depthUnits='km')

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
