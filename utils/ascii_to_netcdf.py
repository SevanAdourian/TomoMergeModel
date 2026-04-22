import pandas as pd

"""
Program to load in ascii file, let the user interactively rename the columns if needed, and convert to xarray on certain dimensions
"""
# load in the ascii file
df = pd.read_csv("alaska_files/joint_800.finmod.ascii", sep=r"\s+")

print(df.head())
print(df.columns)

# rename the columns
cols = df.columns
pairs = {}

for col in cols:
    print(col)
    rename = input("Rename col: ")
    if rename.strip() != "":
        pairs[col] = rename

if len(pairs) != 0:
    df = df.rename(columns=pairs)

print(df.columns)

# convert to xarray
ds = df.set_index(["depth", "latitude", "longitude"]).to_xarray()
print(ds)

# Add metadata


# write netcdf
encoding = {var: {"zlib": True, "complevel": 4} for var in ds.data_vars}

ds.to_netcdf("alaska.nc", encoding=encoding)
print("saved netcdf")
