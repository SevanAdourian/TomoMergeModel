#To merge a regional and global tomographic model
# Sevan Adourian & Dan Frost
# UC Berkeley
# Full Description:
# 1. Read in xyz file of regional model (expecting lon, lat, dv)
#    Converts to Spherical Harmoncs
#    Plots dV and spectrum
# 2. Read in xyz file of global model (expecting lon, lat, dv)
#    Converts to Spherical Harmoncs
#    Plots dV and spectrum
# 3. Construct mask using specified lon/lat region and number of spherical harmonc windows
#    Normalise mask from 0 to 1
# 4. Apply windows to regional model. Also apply to global model although this is technically irrelevant
#    Sum regional and global models
# 5. Write out summed model

import pyshtools
import multiprocessing
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import time
import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import yaml
import pdb
import cartopy.crs as ccrs
import logging

np.set_printoptions(threshold=sys.maxsize)

class ForkedPdb(pdb.Pdb):
    """A Pdb subclass that may be used
    from a forked multiprocessing child

    """
    def interaction(self, *args, **kwargs):
        _stdin = sys.stdin
        try:
            sys.stdin = open('/dev/stdin')
            pdb.Pdb.interaction(self, *args, **kwargs)
        finally:
            sys.stdin = _stdin

def setup(conf):
    # Parameters
    depth_file = conf['depth_file'] 
    reg_lwin = conf['reg_lwin'] #Lmax of regional model (in spherical harmonics)
    reg_nwin = conf['reg_nwin'] #Maximum number of windowing functions. higher = sharper
    
    #  - Mask parameters
    # Part of the regional model that will be preserved in the merging.
    # Read in the depth file
    conf['depth_knots_radius'] =  np.flipud(np.loadtxt(conf['depth_file'], skiprows=1))
    conf['depth_knots'] = conf['radius_in_meters']/1000.0 - conf['depth_knots_radius']
    
    conf['base_global_ascii_files'] = conf['path_to_ascii_files']+'/global_'
    conf['base_regional_ascii_files']    = conf['path_to_ascii_files']+'/regional_'
    conf['base_merged_ascii_files'] = conf['path_to_ascii_files']+'/merged_'

    return conf

def convert_longitude_0_360(lon_lat_field):
    ''' 
    Converts longitude from -180/180 to 0/360
    Deals with redundancy in -180 to put it in 360 (needed by pyshtools)
    Need to add checks whether to avoid random crashing.
    '''

    # Removing redundancy in 180 to put it in 0 360
    ind = np.where(lon_lat_field[:,0] == -180)
    lon_lat_field = np.delete(lon_lat_field, ind, 0)
    ind = np.where(lon_lat_field[:,0] == 0)
    redundant_values = np.squeeze(lon_lat_field[ind,:])
    redundant_values[:,0] = 360.
    ind = np.array(ind).flatten()
    lon_lat_field = np.insert(lon_lat_field, ind, redundant_values,0)
    
    # Convert lon -180/180 to 0/360
    ind =  np.where(lon_lat_field[:,0] < 0)
    lon_lat_field[ind,0] = lon_lat_field[ind,0] + 360.

    return lon_lat_field

def sorting_columns(lon_lat_field):
    ''' 
    Sort first on lon (col1 0to360) then lat (col2 from 90to-90) from NorthWest to SouthEast 
    '''
    
    ind = np.lexsort((lon_lat_field[:,0],-lon_lat_field[:,1]))
    lon_lat_field = lon_lat_field[ind]
    
    return lon_lat_field

def process_slice(conf, depth):
    # Reading in global tomography model
    path_to_global_ascii = conf['base_global_ascii_files'] + str(format(int(depth),'04d'))+".xyz"
    global_tomo = process_ascii_files(conf, path_to_global_ascii)
    zmesh_global = reshape_field(global_tomo)

    global_grid, global_clm = convert_to_spherical_harmonics(zmesh_global, conf["reg_lwin"])
    
    if depth < conf['max_depth_regional_model']:
        # Above where the regional model is defined in depth, actual merging
        # Reading in regional tomography model
        path_to_regional_ascii = conf['base_regional_ascii_files']+str(format(int(depth),'04d'))+".xyz"
        regional_tomo = process_ascii_files(conf, path_to_regional_ascii)
        zmesh_regional = reshape_field(regional_tomo)
        reg_grid, reg_clm = convert_to_spherical_harmonics(zmesh_regional, conf["reg_lwin"])
        
        # Doing mask windowing
        print("Creating spherical windowing", flush=True)
        reg_win_energy_grid = create_window(conf, regional_tomo, global_clm)
        print("Applying window to regional and global field", flush=True)
        summed_grid = apply_window(conf, global_grid, global_clm, reg_grid, reg_win_energy_grid, depth)
        # ForkedPdb().set_trace()
        print("Finished merging models, for depth "+str(depth), flush=True)
        write_model(conf, depth, summed_grid)
        
    else:
        # Below where the regional model is defined, we just write the global model
        write_model(conf, depth, global_grid)
        
    return

def write_model(conf, depth, grid):
    '''
    Write out masked, merged model as 3 columns: lon, lat, dv 
    '''
    #
    m_dv = grid.data
    c_lats = grid.lats()
    c_lons = grid.lons()
    m_lonlat = np.meshgrid(c_lons, c_lats)
    c_lonlat = np.transpose(np.reshape(m_lonlat,(2,len(c_lats)*len(c_lons))))
    c_dv = np.reshape(m_dv,(len(c_lats)*len(c_lons),1))
    c_lonlatdv = np.column_stack((c_lonlat,c_dv))
    
    model_out = conf['base_merged_ascii_files']+str(format(int(depth),'04d'))+".xyz"
    dataf = pd.DataFrame(data=c_lonlatdv, columns=('LON','LAT','DV') ) 
    dataf.to_csv(model_out, sep=' ', header=True, float_format='%.4f', index=False)

    return
    
def process_ascii_files(conf, path_to_ascii_file):
    ''' 
    Read the tomographic model in ascii
    The tomographic model in ascii should be stored in three columns:
    longitude latitude field (needs to be same format for regional and global) 
    '''
    #
    lon_lat_field = np.loadtxt(path_to_ascii_file) # tomographic model

    # Make sure it is defined between 0 to 360 in longitude
    if min(lon_lat_field[:,0]) < 0:
        print("Converting longitudes from -180/180 to 0/360 and sorting")
        lon_lat_field = convert_longitude_0_360(lon_lat_field)
        lon_lat_field = sorting_columns(lon_lat_field)
    else:
        print("Longitudes already in the 0/360 range, just sorting")
        sorting_columns(lon_lat_field)
        
    return lon_lat_field

def reshape_field(lon_lat_field):
    #  - Create mesh of z
    xvar = np.unique(lon_lat_field[:,0])
    xlen = (len(xvar))
    yvar = np.unique(lon_lat_field[:,1])
    ylen = (len(yvar))
    zmesh = (np.reshape(lon_lat_field[:,2],(ylen,xlen)))

    return zmesh


def convert_to_spherical_harmonics(zmesh, reg_lwin):
    #   - Convert to spherical harmonics with pySHtools
    grid = pyshtools.SHGrid.from_array(zmesh, grid= 'DH')
    clm = pyshtools.SHGrid.expand(grid)
    clm = clm.pad(reg_lwin)  #Pad to match global clm

    # Trim coefficients in SPH to match number used in regional model
    grid = pyshtools.SHCoeffs.expand(clm,lmax=reg_lwin)

    return grid, clm

def plot_map_and_spectra(conf, grid_object, clm_object, file_name):
    fig, (col1, col2) = plt.subplots(2, 1)
    grid_object.plot(ax=col1,colorbar='right', cb_label='Power', show=False)
    clm_object.plot_spectrum(ax=col2)
    fig.legend(loc = 'upper right')
    fig.savefig(file_name, dpi=400) #Original global data
    plt.close(fig)

    return

def create_window(conf, reg_field, global_clm):
    '''
    Create a mask of spherical harmonic windows using the range given above
    Sharpness of the edges of the mask are controlled by reg_nwin
    Normalise mask between 0 and 1 and then apply to the regional data
    '''
    
    #-----------
    xvar = np.unique(reg_field[:,0])
    xlen = (len(xvar))
    yvar = np.unique(reg_field[:,1])
    ylen = (len(yvar))

    lon_left=conf['lon_min_mask']
    lon_right=conf['lon_max_mask']
    lat_bottom = conf['lat_min_mask'] # possible range: -90, 90 deg
    lat_top    =  conf['lat_max_mask'] # possible range: -90, 90 deg

    #  - Regional mask - set 1s inside region and 0s outside
    reg_zmask= np.where( ( (reg_field[:,0] >= lon_left) & (reg_field[:,0] <= lon_right) )&\
                         ( (reg_field[:,1] >= lat_bottom) & (reg_field[:,1] <= lat_top) ),1,0)
    # 
    reg_zmesh_mask=(np.reshape(reg_zmask,(ylen,xlen)))
    reg_mask_xyz=np.column_stack((reg_field[:,0],reg_field[:,1],reg_zmask))
    
    #-----------
    if conf['win_type'] == 'spherical': # spherical or rectangular'
        #   - Construct spherical harmonic window function from mask
        reg_win=pyshtools.SHWindow.from_mask(reg_zmesh_mask,lwin=conf['reg_lwin'])
        reg_win_clm=pyshtools.SHWindow.to_shcoeffs(reg_win,0)
        reg_win_clm_pad=reg_win_clm.pad(global_clm.lmax)  #Pad to match global clm
        
        reg_win_energy = (reg_win.to_shgrid(0).to_array())**2
        for i in range(1,conf['reg_nwin']):
            reg_win_energy += (reg_win.to_shgrid(i).to_array())**2
            reg_win_energy_grid = pyshtools.SHGrid.from_array(reg_win_energy)
            
        #  - Normalise window (0 to 1) so that we can mask outside and inside
        valmax=(np.amax(reg_win_energy_grid.data))
        reg_win_energy_grid=reg_win_energy_grid/valmax
    
        reg_win_energy_clm = pyshtools.SHGrid.expand(reg_win_energy_grid)
        reg_win_energy_clm=reg_win_energy_clm.pad(conf['reg_lwin'])  #Pad to match global clm
        reg_win_energy_grid = pyshtools.SHCoeffs.expand(reg_win_energy_clm)  #Grid of mask

    elif conf['win_type'] == 'rectangular': # spherical or rectangular
        # Construct rectangular window from mask (no smoothing at edges, sharp window)
        reg_win_grid=pyshtools.SHGrid.from_array(reg_zmesh_mask)
        reg_win_clm = pyshtools.SHGrid.expand(reg_win_grid)
        reg_win_clm_pad=reg_win_clm.pad(global_clm.lmax)  #Pad to match reg clm
        reg_win_energy_grid=pyshtools.SHCoeffs.expand(reg_win_clm_pad)

        
    return reg_win_energy_grid


def apply_window(conf, global_grid, global_clm, reg_grid, reg_win_energy_grid, depth):

    # - Multiply grid by mask
    reg_grid_masked=reg_grid*reg_win_energy_grid  
    reg_clm_masked=pyshtools.SHGrid.expand(reg_grid_masked)
    
    # Plot reg map and spec (masked)
    plot_map_and_spectra(conf, reg_grid_masked, reg_clm_masked, \
                         conf['path_to_figures']+'/Regional_MaskedSpec_'+str(int(depth))+'km.png')    
    
    #   - Apply windows to global grid
    global_grid_masked=global_grid*reg_win_energy_grid
    global_clm_masked=pyshtools.SHGrid.expand(global_grid_masked)
    
    # Plot global map and spec (masked)
    plot_map_and_spectra(conf, global_grid_masked, global_clm_masked, \
                         conf['path_to_figures']+'/Global_MaskedSpec_'+str(int(depth))+'km.png')
    
    # Sum regional and global spectra inside box
    summed_clm_masked=reg_clm_masked+global_clm_masked
    summed_grid_masked=pyshtools.SHCoeffs.expand(summed_clm_masked)

    plot_map_and_spectra(conf, summed_grid_masked, summed_clm_masked, \
                         conf['path_to_figures']+'/Summed_DoubleMaskedSpec_'+str(int(depth))+'km.png')
    
    # Sum regional spectrum inside box (masked) with global spectra outside box (unmasked)
    summed_clm=reg_clm_masked+global_clm
    summed_grid=pyshtools.SHCoeffs.expand(summed_clm)
    
    # Plot summed map and spec (masked)
    plot_map_and_spectra(conf, summed_grid, summed_clm, \
                         conf['path_to_figures']+'/Summed_MaskedSpec_'+str(int(depth))+'km.png')
    
    return summed_grid

    
# -------- MAIN ----------#
if __name__ == '__main__':
    # Initiate computing slices on multiple processes
    multiprocessing.set_start_method('spawn')

    # Load yaml file
    for conf_file in sys.argv[1:]:
        print('using configuration loaded from: %s' % (conf_file))
        with open(conf_file) as f:
            conf = yaml.safe_load(f)
            conf = setup(conf)
        print('Loaded the input file and configured the run')

    # Sequential way:
    # Loop over depths:
    # for depth in conf['depth_knots'][0:1]:
    
    #     process_slice(conf, depth)
    
    # Loop over depths with multiprocessing:
    depths = conf['depth_knots']  # list of slices to process
    n_depths = len(depths)  # number of times to run process_slice
    
    t0 = time.time()
    with multiprocessing.Pool() as pool:
        for i in range(n_depths):
            print("Processing depth "+str(depths[i]))
            pool.apply_async(process_slice, args=(conf, depths[i]))
        
        # Wait for all processes to complete
        pool.close()
        pool.join()

    t1 = time.time()
    print("Generated the merged model in "+str(t1-t0)+" seconds")
    print("All done and executed normally")
