"""
Requirements:
-------------
numpy
pynbody_routines
pynbody
astropy
matplotlib

"""

import numpy as np
import pynbody_routines  as pr 
import pynbody
from astropy import units as u
import matplotlib
import matplotlib.pyplot as plt



def rotate_fire(part):
    """
    Reposition and rotate the simulation to see the halo face on
        
    Parameters:
    ----------
    part: gizmo particle dictionary 

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
    pos = hfaceon.star['pos']
    vel = hfaceon.star['vel']*f
    
    return hfaceon, pos, vel



def graph_profiles(hfaceon, rmin, rmax):
    """
    Graph tangential velocity dispersion, radial velocity dispersion, and density profiles
        
    Parameters:
    ----------
    hfaceon: pynbody SimSnap
    rmin: float
    rmax: float

    Returns:
    --------
    
    """
    
    # Create star and dark matter profiles
    pStar = pynbody.analysis.profile.Profile(hfaceon.s,min=rmin,max=rmax,ndim=3)
    pDark = pynbody.analysis.profile.Profile(hfaceon.d,min=rmin,max=rmax,ndim=3)
    
    # Declare figure with three subplots
    fig, axs = plt.subplots(2,3,figsize=(21,12))
    
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
    axs[0,2].set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-2}$]')
    
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
    axs[1,2].set_ylabel(r'$\Sigma_{\star}$ [M$_{\odot}$ kpc$^{-2}$]')
    
    return


