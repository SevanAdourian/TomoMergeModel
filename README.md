################################
Author: Sevan Adourian
This package has the purpose to smoothly merge a global model within a global model.

# First you need to create the conda environment using the environment file
# Discard this part if you already have such an environment
conda env create -f environment.yaml

# Add the ucbpy libraries to what is read by python
cd /path/to/ucbpy/python3/libraries/
conda develop .

# If you have a PYTHONPATH defined by the bash_rc by default, discard it
PYTHONPATH=''

# You are now ready to use the package

# Download glad-m25 (to run the example)
wget https://ds.iris.edu/files/products/emc/emc-files/glad-m25-vp-0.0-n4.nc
mv glad-m25-vp-0.0-n4.nc ./Data/netcdf

##### First step, write asciis of your model at depths of interest #####
To run this code you need to edit conf.yaml and simply run
python 0_netcdf_to_ascii.py conf.yaml

If you have both tomography models in netcdf, great. However, this package is not
yet capable of detecting the exact tree and variable names of the netcdf packages, 
so you might have to tweak the code 0_netcdf_to_ascii.py (this feature will be coming soon)
This code outputs ascii files for the global and regional model as follows in 3 columns:
Longitude Latitude dV where dV is the velocity perturbation with respect to a 1D model.
If your models are not in netcdf, this is the kind of file you will need to reproduce.
The spacing in latitude and longitude does not need to be the same for global and regional.

For those of you using SEMUCB, I directly provide the ascii files (since the netcdf file is not
defined at enough depth knots to interpolate).
See in ./Data/ascii/semucb

##### Second step, perform the merging #####
If you followed step 1, you can go ahead and directly run 
python 1_TomoModelMerge.py conf.yaml

A good value for the nreg_win is 100, but this takes a bit more time to compute
Adjust lreg_win depending on the resolution of your regional model
The output is the merged ascii file at each depth. From there you can use the package ascii2A3d
to write your A3d model.
