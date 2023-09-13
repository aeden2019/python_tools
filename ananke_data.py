"""
Ananke function to apply Vera Rubin parameters to FIRE data. 

Requirements:
-------------
gizmo_analysis



author: Andrew Eden
github: aeden2019

"""
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
import sys
import pynbody

sys.path.append('/home/jovyan/home/python_tools')
import gizmo_analysis as ga
import halo_analysis as halo
import nba

import ananke as an



def get_ananke(sim, snap, pos, vel, mask):
    """
    
        
    Parameters:
    ----------
    self - 
    snap
    pos - transformed position array input
    vel - transformed velocity array input
    

    Returns:
    --------
    pos_out - position array outputed from ananke
    vel_out - velocity array outputed from ananke
    
    """
    part = ga.io.Read.read_snapshots(['star', 'dark'], 'index', snap, sim, 
                              assign_hosts=True, sort_dark_by_id=True)
    
    # TESTING - LENGTH
    if len(part['star']['mass']) == len(mask):
        print("SAME LENGTH")
    else:
        print("DIFFERENT LENGTH")
    
    # Create p dictionary to store particle data
    p = {}
    # Already masked
    p['pos3'] = pos       # position in kpc
    p['vel3'] = vel       # velocity in km/s
    # Not yet masked
    # CHECK IF MASK IS SAME LENGTH, AND BOOLEAN
    p['mass'] = part['star']['mass'][mask]                   # mass in solar masses
    p['age'] = part['star'].prop('age')[mask]                  # log age in Gyr
    p['feh'] = part['star'].prop('metallicity.fe')[mask]       # [Fe/H]
    p['helium'] = part['star'].prop('metallicity.he')[mask]    # [He/H]
    p['carbon'] = part['star'].prop('metallicity.c')[mask]     # [C/H]
    p['nitrogen'] = part['star'].prop('metallicity.n')[mask]   # [N/H]
    p['neon'] = part['star'].prop('metallicity.ne')[mask]      # [Ne/H]
    p['magnesium'] = part['star'].prop('metallicity.mg')[mask] # [Mg/H]
    p['silicon'] = part['star'].prop('metallicity.si')[mask]   # [Si/H]
    p['sulphur'] = part['star'].prop('metallicity.s')[mask]    # [S/H]
    p['calcium'] = part['star'].prop('metallicity.ca')[mask]   # [Ca/H]
    p['oxygen'] = part['star'].prop('metallicity.o')[mask]     # [O/H]
    p['alpha'] = part['star'].prop('metallicity.mg - metallicity.fe')[mask]          # [Mg/Fe]
    p['parentid'] = part['star']['id'][mask]                   # indices of parent particles in snapshot
    p['dform'] = np.zeros(part['star']['position'].shape[0], dtype='float32')[mask]  # dummy variable for now
    
    # TESTING 
    print(f"Length of pos: {len(pos)}")
    print(f"Length of vel: {len(vel)}")
    print(f"Length of mass: {len(part['star']['mass'])}")
    print(f"Length of mask: {len(mask)}")
    
    print("ANANKE - Subsample")
    # Declare subsampled factor and subsampled dictionary
    subsample_factor = 1000
    p_subsampled = {}
    # Iterate through part dictionary and subsample
    for key, array in p.items():
        # TESTING
        print(f"Key: {key}, Length: {len(array)}")
        # Subsample the array
        array_subsampled = array[::subsample_factor]
        # Store subsampled array in new dictionary
        p_subsampled[key] = array_subsampled
        print(f"Key: {key}, Length: {len(array_subsampled)} (Subsampled)")
        
    
    # Run the ananke process
    name='sim'
    print("ANANKE - Call")
    ananke = an.Ananke(p_subsampled, name, photo_sys='padova/LSST', cmd_magnames='rmag,gmag-rmag'
                                                , app_mag_lim_lo=17, app_mag_lim_hi=27.5, abs_mag_lim_lo=-7.0, abs_mag_lim_hi=10.0
                                                , color_lim_lo=-1000, color_lim_hi=1000, r_max=1000)
    print("ANANKE - Run")
    ananke.run()
    
    # Make a survey using LSST (outputted as vaex data structure)
    print("ANANKE - Make Vaex Survey")
    survey = ananke._output
    df = survey._vaex
    
    # Combine the position and velocities
    print("ANANKE - Combine position and velocities")
    pos_out = np.column_stack((df['px'].values, df['py'].values, df['pz'].values))
    vel_out = np.column_stack((df['vx'].values, df['vy'].values, df['vz'].values))
    
    print("ANANKE - End")

    # TESTING - RUN WITHOUT ANANKE
#     pos_out = p_subsampled['pos3']
#     vel_out = p_subsampled['vel3']
    
    return pos_out, vel_out
    
    
    
    
    
    
    
    