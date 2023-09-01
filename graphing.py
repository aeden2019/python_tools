"""
Graphing functions to rotate FIRE simulations and graph profiles. Does not include subhalo object removal. 

Requirements:
-------------
numpy
pynbody_routines
pynbody
astropy
matplotlib

author: Andrew Eden
github: aeden2019

"""

import numpy as np
import pynbody_routines  as pr 
import pynbody
from astropy import units as u
import matplotlib
import matplotlib.pyplot as plt
import nba
import plotting as pl
import json
from matplotlib.cm import get_cmap


def mollweide(pos, vel, rmin, rmax, bmin, bmax, sim, snap, matter):
    """
    Graph a mollweide plot
        
    Parameters:
    ----------
    pos: pynbody SimArray
    vel: pynbody SimArray
    rmin: float
    rmin: float
    bmin: int 
    bmax: int
    sim: string
    snap: string
    matter: string (star or dark)

    Returns:
    --------

    
    """
    
    # Cut stellar distribution based on radius limits
    dist = np.sqrt(np.sum(pos**2, axis=1))
    dist_cut1 = np.where((dist > rmin) & (dist< rmax)) 
        
    # Create kinematic profile
    kinematics1 = nba.kinematics.Kinematics(pos[dist_cut1],  vel[dist_cut1])
    
    # Convert pos and vel into galactic coordinate from kinematics 
    pos_galactic = kinematics1.pos_cartesian_to_galactic()
    vel_galactic = kinematics1.vel_cartesian_to_galactic()
    
    # Declare figure title
    fig_title = "{} satellite vr {}; {}-{} kpc; {} snapshot".format(sim, matter, rmin, rmax, snap)
    
    # Call mollweide projection function 
    pl.mollweide_projection(pos_galactic[0]*180/np.pi, pos_galactic[1]*180/np.pi, 0, 0, 
                            title=fig_title, bmin=bmin, bmax=bmax, nside=40, smooth=5, overdensity=True)
    
    return



def rotate_fire(part, matter):
    """
    Reposition and rotate the simulation to see the halo face on
        
    Parameters:
    ----------
    part: gizmo particle dictionary 
    matter: string

    Returns:
    --------
    hfaceon: pynbody SimSnap
    pos: pynbody SimArray
    vel: pynbody SimArray
    
    """
    
    # Create masks for stellar and dark matter particles
    npart = len(part['dark'].prop('mass'))
    mask_sub = np.ones(npart, dtype=bool)
    nparts = len(part['star'].prop('mass'))
    mask_subs = np.ones(nparts, dtype=bool)
    
    # Create pynbody object
    hfaceon = pr.pynbody_halo(part, mask=mask_sub, masks=mask_subs)
    
    # Apply rotation 
    pynbody.analysis.angmom.faceon(hfaceon, cen=(0,0,0))
    
    # Declare velocity unit conversion 
    f = 1* (u.km/u.s).to(u.kpc/u.Gyr)
    
    # Get positions and velocities from pynbody object
    if matter == 'star':
        pos = hfaceon.star['pos']
        vel = hfaceon.star['vel']*f
    elif matter == 'dark':
        pos = hfaceon.dark['pos']
        vel = hfaceon.dark['vel']*f
    else:
        print('Selected matter type is not defined, please use star or dark instead')
    
    return hfaceon, pos, vel



def graph_profiles(hfaceon, rmin, rmax, sim, snap):
    """
    Graph tangential velocity dispersion, radial velocity dispersion, and density profiles
        
    Parameters:
    ----------
    hfaceon: pynbody SimSnap
    rmin: float
    rmax: float
    sim: string
    snap: string

    Returns:
    --------
    
    """
    
    # Create star and dark matter profiles
    pStar = pynbody.analysis.profile.Profile(hfaceon.s,min=rmin,max=rmax,ndim=3)
    pDark = pynbody.analysis.profile.Profile(hfaceon.d,min=rmin,max=rmax,ndim=3)
    
    # Declare figure with three subplots
    fig, axs = plt.subplots(2,3,figsize=(21,12))
    
    # Add title 
    plt.title('Velocity and density profiles for {} snapshot {}'.format(sim, snap))
    
    # Stellar plots
    axs[0,0].plot(pStar['rbins'].in_units('kpc'),pStar['vr_disp'].in_units('km s^-1'))
    axs[0,0].set_title('Stellar Radial Velocity Dispersion')
    axs[0,0].set_xlabel('R [kpc]')
    axs[0,0].set_ylabel('$\sigma_{r}$')
    
    axs[0,1].plot(pStar['rbins'].in_units('kpc'),pStar['vt_disp'].in_units('km s^-1'))
    axs[0,1].set_title('Stellar Tangential Velocity Dispersion')
    axs[0,1].set_xlabel('R [kpc]')
    axs[0,1].set_ylabel('$\sigma_{r}$')
    
    axs[0,2].plot(pStar['rbins'],pStar['density'])
    axs[0,2].semilogy()
    axs[0,2].set_title('Stellar Density')
    axs[0,2].set_xlabel('R [kpc]')
    axs[0,2].set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
    
    # Dark matter plots   
    axs[1,0].plot(pDark['rbins'].in_units('kpc'),pDark['vr_disp'].in_units('km s^-1'))
    axs[1,0].set_title('Dark Matter Radial Velocity Dispersion')
    axs[1,0].set_xlabel('R [kpc]')
    axs[1,0].set_ylabel('$\sigma_{r}$')
    
    axs[1,1].plot(pDark['rbins'].in_units('kpc'),pDark['vt_disp'].in_units('km s^-1'))
    axs[1,1].set_title('Dark Matter Tangential Velocity Dispersion')
    axs[1,1].set_xlabel('R [kpc]')
    axs[1,1].set_ylabel('$\sigma_{r}$')
    
    axs[1,2].plot(pDark['rbins'],pDark['density'])
    axs[1,2].semilogy()
    axs[1,2].set_title('Dark Matter Density')
    axs[1,2].set_xlabel('R [kpc]')
    axs[1,2].set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-3}$]')
    
    #plt.close()
    
    return



def store_profile_data(hfaceon, rmin, rmax, sim, snap, file_path):
    
    """
    Extracts data from pynbody profiles, and stores them for graphing later
        
    Parameters:
    ----------
    hfaceon: pynbody SimSnap
    rmin: float
    rmax: float
    sim: string
    snap: string
    file_path

    Returns:
    --------
    
    """
    
    # Create star and dark matter profiles
    pStar = pynbody.analysis.profile.Profile(hfaceon.s,min=rmin,max=rmax,ndim=3)
    pDark = pynbody.analysis.profile.Profile(hfaceon.d,min=rmin,max=rmax,ndim=3)
    
    # Create a profile instance to store all the necessary data 
    instance_data = {
        'pStar': {
            'radial_bins_kpc': pStar['rbins'].in_units('kpc').tolist(),
            'radial_bins': pStar['rbins'].tolist(),
            'vr_disp': pStar['vr_disp'].in_units('km s^-1').tolist(),
            'vt_disp': pStar['vt_disp'].in_units('km s^-1').tolist(),
            'density': pStar['density'].tolist()
        },
        'pDark': {
            'radial_bins_kpc': pDark['rbins'].in_units('kpc').tolist(),
            'radial_bins': pDark['rbins'].tolist(),
            'vr_disp': pDark['vr_disp'].in_units('km s^-1').tolist(),
            'vt_disp': pDark['vt_disp'].in_units('km s^-1').tolist(),
            'density': pDark['density'].tolist()
        }
    }
        
    # Read the existing data from the file
    try:
        with open(file_path, 'r') as file:
            profiles = json.load(file)
    except FileNotFoundError:
        profiles = {}
    
    # Add the new instance data to the existing profiles
    profiles[f'snap_{snap}'] = instance_data
    
    # Write the updated data back to the file
    with open(file_path, 'w') as file:
        json.dump(profiles, file, indent=4)
        

