"""
This code compute quantities related to the anisotropy parameter
beta defined as:

beta = 1 - (sigma_t^2 / 2*sigma_r^2)

where sigma_t is the tangential velocity dispersion and sigma_r the
radial velocity dispersion.

the code compute beta as a function of radius beta(r), beta as a
function of r and z. and the velocity dispersions as a function of
radius.

Requirements:
-------------
Numpy.

"""

import numpy as np
#import mock_observations
from mass_estimators import octant_cut

def beta_r(pos, vel, n_bins, rmax):
    dr = np.linspace(0, rmax, n_bins)
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    beta_dr = np.zeros(len(dr)-1)
    for i in range(len(dr)-1):
        index = np.where((r<dr[i+1]) & (r>dr[i]))
        beta_dr[i] = beta(pos[index], vel[index])
    return dr, beta_dr


# def beta_octants(pos, vel, lmin, lmax, bmin, bmax, n_bins, rmax):
# 	l_lmc, b_lmc = mock_observations.galactic_coodinates(pos, vel)
# 	l_oct, b_oct, pos_oct, vel_oct = octant_cut(lmin, lmax, bmin, bmax, l_lmc, b_lmc, pos, vel)
# 	r_oct, beta_r_oct = beta_r(pos_oct, vel_oct, n_bins, rmax)
# 	return r_oct, beta_r_oct

def beta_r_phi_z(pos, vel, n_bins, phi_bins, rmax):
    minz = -200
    maxz = 200
    binz = 20
    dr = np.linspace(0, rmax, n_bins)
    dphi = np.linspace(-np.pi, np.pi, phi_bins)
    dz = np.linspace(minz, maxz, binz)
    z = pos[:,2]
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    phi = np.arctan2(pos[:,1], pos[:,0])
    Beta_dr = np.zeros((len(dz)-1, len(dr)-1, len(dphi)-1))
    for k in range(len(dz)-1):
        for i in range(len(dr)-1):
            for j in range(len(dphi)-1):
                index = np.where((r<dr[i+1]) & (r>dr[i]) & \
                                 (phi<dphi[j+1]) & (phi>dphi[j])\
                                & (z<dz[k+1]) & (z>dz[k]))
                Beta_dr[k][i][j] = beta(pos[index], vel[index])
    return Beta_dr


def beta(xyz, vxyz):
    """
    Computes the anisotropy $\beta$.

    Parameters:
    -----------
    xyz : 3d array with the postions of the particles.
    vxyz : 3d array with the velcities of the particles.

    Returns:
    --------

    Beta : double
        The value of the anisotropy parameter.
    """

    sigma_r, sigma_theta, sigma_phi = velocity_dispersion(xyz, vxyz)
    sigma_t = ((sigma_theta**2 + sigma_phi**2))**0.5
    Beta = 1 - sigma_t**2.0/(2.0*sigma_r**2.0)
    return Beta


def velocity_dispersion(pos, vel):
    """
    Computes the velocity dispersions in spherical coordinates.

    """

    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    theta = np.arccos(pos[:,2]/r)
    phi = np.arctan2(pos[:,1], pos[:,0])

    vr = np.sin(theta)*np.cos(phi)*vel[:,0] +\
         np.sin(theta)*np.sin(phi)*vel[:,1] +\
         np.cos(theta)*vel[:,2]

    v_theta = np.cos(theta)*np.cos(phi)*vel[:,0] +\
              np.cos(theta)*np.sin(phi)*vel[:,1] - \
              np.sin(theta)*vel[:,2]

    v_phi = -np.sin(phi)*vel[:,0] + np.cos(phi)*vel[:,1]

    vr_disp = np.std(vr)
    vtheta_disp = np.std(v_theta)
    vphi_disp = np.std(v_phi)
    v_tan_disp = np.sqrt(vtheta_disp**2 + vphi_disp**2)

    return vr_disp, vtheta_disp, vphi_disp




def velocities_means(pos, vel):
    """
    Computes the velocity dispersions in spherical coordinates.

    """

    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    theta = np.arccos(pos[:,2]/r)
    phi = np.arctan2(pos[:,1], pos[:,0])

    vr = np.sin(theta)*np.cos(phi)*vel[:,0] +\
         np.sin(theta)*np.sin(phi)*vel[:,1] +\
         np.cos(theta)*vel[:,2]

    v_theta = np.cos(theta)*np.cos(phi)*vel[:,0] +\
              np.cos(theta)*np.sin(phi)*vel[:,1] - \
              np.sin(theta)*vel[:,2]

    v_phi = -np.sin(phi)*vel[:,0] + np.cos(phi)*vel[:,1]
    v_tan = np.sqrt(v_phi**2 + v_theta**2)
    return np.mean(vr), np.mean(v_phi), np.mean(v_theta)

def velocity_dispersions_r(pos, vel, n_bins, rmax):
    dr = np.linspace(0, rmax, n_bins)
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**0.5
    vr_disp_r = np.zeros(len(dr)-1)
    vtheta_disp_r = np.zeros(len(dr)-1)
    vphi_disp_r = np.zeros(len(dr)-1)

    for i in range(len(dr)-1):
        index = np.where((r<dr[i+1]) & (r>dr[i]))
        vr_disp_r[i], vtheta_disp_r[i], vphi_disp_r[i] = velocity_dispersion(pos[index], vel[index])

    return vr_disp_r, vtheta_disp_r, vphi_disp_r
