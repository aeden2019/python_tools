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



def get_ananke(sim, snap, pos, vel):
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
    part = ga.io.Read.read_snapshots(['dark', 'star'], 'index', snap, sim, 
                              assign_hosts=True, particle_subsample_factor=0.1, sort_dark_by_id=True)
    
    # Create p dictionary to store particle data
    p = {}
    p['pos3'] = pos       # position in kpc
    p['vel3'] = vel       # velocity in km/s
    p['mass'] = part['star']['mass']                     # mass in solar masses
    p['age'] = part['star'].prop('age')                  # log age in Gyr
    p['feh'] = part['star'].prop('metallicity.fe')       # [Fe/H]
    p['helium'] = part['star'].prop('metallicity.he')    # [He/H]
    p['carbon'] = part['star'].prop('metallicity.c')     # [C/H]
    p['nitrogen'] = part['star'].prop('metallicity.n')   # [N/H]
    p['neon'] = part['star'].prop('metallicity.ne')      # [Ne/H]
    p['magnesium'] = part['star'].prop('metallicity.mg') # [Mg/H]
    p['silicon'] = part['star'].prop('metallicity.si')   # [Si/H]
    p['sulphur'] = part['star'].prop('metallicity.s')    # [S/H]
    p['calcium'] = part['star'].prop('metallicity.ca')   # [Ca/H]
    p['oxygen'] = part['star'].prop('metallicity.o')     # [O/H]
    p['alpha'] = part['star'].prop('metallicity.mg - metallicity.fe')          # [Mg/Fe]
    p['parentid'] = part['star']['id']                   # indices of parent particles in snapshot
    p['dform'] = np.zeros(part['star']['position'].shape[0], dtype='float32')  # dummy variable for now
    
    # Run the ananke process
    name='sim'
    ananke = an.Ananke(p, name, fsample=0.01, photo_sys='padova/LSST', cmd_magnames='rmag,gmag-rmag'
                                                , app_mag_lim_lo=17, app_mag_lim_hi=27.5, abs_mag_lim_lo=-7.0, abs_mag_lim_hi=10.0
                                                , color_lim_lo=-1000, color_lim_hi=1000, r_max=1000)
    ananke.run()
    
    # Make a survey using LSST (outputted as vaex data structure)
    survey = ananke._output
    df = survey._vaex
    
    # Combine the position and velocities
    pos_out = np.column_stack((df['px'].values, df['py'].values, df['pz'].values))
    vel_out = np.column_stack((df['vx'].values, df['vy'].values, df['vz'].values))
    
    return pos_out, vel_out
    
    
    
    
    
    
    
    