import numpy as np
from astropy import constants
from astropy import units
# import mock_observations

G = constants.G

def radial_velocity(pos, vel):
    # transforming data to spherical coordinates.
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**(0.5)
    theta = np.arccos(pos[:,2]/r)
    phi = np.arctan2(pos[:,1], pos[:,0])

    # Radial velocity:
    vr = np.sin(theta)*np.cos(phi)*vel[:,0] + \
         np.sin(theta)*np.sin(phi)*vel[:,1] + \
         np.cos(theta)*vel[:,2]
    
    return vr


def watkins11(alpha, beta, gamma, pos, vel, r_out):
    """

    MW mass estimator. Eqns 16 & 24 From Watkins et al 2011
    M = C1/G <v^2 r^{alpha}>
    M = C2/G <vr^2 r^{alpha}>
    
    C1 = (alpha + gamma - 2 beta)/ (3 - 2 beta) * r_out^{1-alpha}
    C2 = (alpha + gamma - 2 beta) * r_out^{1-alpha}
    
    Where:
    gamma = - dlogrho / dr
    alpha = 
    
    Parameters:
    -----------

    alpha : float
        Potential powerlaw expoenent [-1, 1]
        
    beta :
        anisotropy profile
    gamma : float
        Density powerlaw exponent.
    pos : 3d numpy array
        positions of the halo particles in cartessian coordinate system.
        in [kpc]
    vel : 3d numpy array
        velocities of the halo particles in cartessian coordinate system.
        in [kpc]
    r_out : float
        maximum radius of the halo in [kpc].

    return:
    -------

    Mv : enclosed mass using the full velocity vector in units of Msun/1E10.
    Mvr : enclosed mass using the LOS velocity vector in units of Msun/1E10.

    """
    r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**(0.5)
    v = (vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)**(0.5)

    
    vr = radial_velocity(pos, vel)
    
    # Quantities units.
    r = r*units.kpc
    v = v*units.km / units.s
    vr = vr*units.km / units.s
    r_out = r_out * units.kpc

    # Computing constants of eqns 16 & 24 and Masses.
    C1 = (alpha + gamma - 2*beta) / (3 - 2*beta) * r_out**(1-alpha)
    Mv = 1/G * np.mean(C1*v**2*r**(alpha))

    C2 = (alpha + gamma - 2*beta) * r_out**(1-alpha)
    Mvr = 1/G * np.mean(C2*vr**2*r**(alpha))

    return Mv.to(units.Msun).value/1E10, Mvr.to(units.Msun).value/1E10




def mass_estimator_profile(alpha, beta, gamma, pos, vel, r_in, r_out, nbins=30):
	"""
	Estimates the mass profile using Watkins+11 mass estimators.

	Parameters:
	-----------
	alpha : float
		power law exponent of the underlying DM potential
	beta : float or numpy.array
		beta parameter if input is beta profile it has to match
	gamma : float 
		Tracers density powe law
	pos : numpy.ndarray
		3d array of the distances of the DM particles in cartessian galactocentric coordinates 
	vel : numpy.ndarray
		3d array of the velocities of the DM particles in cartessian galactocentric coordinates
	r_in : float
		minimum distance to compute the enclosed mass
	r_out :
		maximum distance to compute the enclosed mass
	nbins :int
		number of means to compute the density profile

	Returns:
	--------
	r : numpy.ndarray
		galactocentric distance bins
	mass_v : numpy.ndarray
		enclosed mass profile computed with 3d velocities
	mass_vr : numpy.ndarray
		enclosed mass profile computed with radial velocities
		
	"""
	
	r_cuts = np.linspace(r_in, r_out, nbins)
	dr = (r_cuts[1]-r_cuts[0])/2.
	r = (pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)**(0.5)
	v = (vel[:,0]**2 + vel[:,1]**2 + vel[:,2]**2)**(0.5)
	Mv = np.zeros(len(r_cuts))
	Mvr = np.zeros(len(r_cuts))

	if type(beta)==np.ndarray:
		assert len(beta)==len(r_cuts), "Error: beta profile should match the desired distance bins"
	
	for i in range(1,len(r_cuts)-1):
		index_cut = np.where((r<r_cuts[i]) & (r>r_cuts[i-1]))[0]
		index_cumulative = np.where(r<r_cuts[i])[0]
		r_cut = r[index_cut]
		v_cut = v[index_cut]

		if ((type(beta)==float) | (type(beta)==int)):
			Mv[i], Mvr[i] = watkins11(alpha, beta, gamma, \
			                          pos[index_cut],\
			                          vel[index_cut], r_out)
	
		elif type(beta)==np.ndarray:
			Mv[i], Mvr[i] = watkins11(alpha, beta[i], gamma,\
			                          pos[index_cut], vel[index_cut],\
			                          r_out)

	return r_cuts[:-1], Mv, Mvr



def octant_cut(lmin, lmax, bmin, bmax, l, b, pos, vel):
    assert((lmin<lmax) & (bmin < bmax))
    oct_cut = np.where((l<lmax) & (l>lmin) & (b>bmin) & (b<bmax))
    return l[oct_cut], b[oct_cut], pos[oct_cut], vel[oct_cut]

# def mass_estimator_oct_profile(alpha, beta, gamma, pos, vel, r_in, r_out, lmin, lmax, bmin, bmax, nbins=30):
# 	"""
	
# 	Estimates the mass profile using Watkins+11 mass estimators in regions in the sky selected 
# 	in galactocentric coordinates.

# 	Parameters:
# 	-----------
# 	alpha : float
# 		power law exponent of the underlying DM potential
# 	beta : float or numpy.array
# 		beta parameter if input is beta profile it has to match
# 	gamma : float 
# 		Tracers density powe law
# 	pos : numpy.ndarray
# 		3d array of the distances of the DM particles in cartessian galactocentric coordinates 
# 	vel : numpy.ndarray
# 		3d array of the velocities of the DM particles in cartessian galactocentric coordinates
# 	r_in : float
# 		minimum distance to compute the enclosed mass
# 	r_out :
# 		maximum distance to compute the enclosed mass

# 	lmin : float
# 		minimum longitude in degrees
# 	lmax : float
# 		maximum longitude in degrees
# 	bmin : float
# 		minimum latitude in degrees
# 	bmax : float
# 		maximum latitude in degrees
# 	nbins :int
# 		number of means to compute the density profile

# 	Returns:
# 	--------
# 	r : numpy.ndarray
# 		galactocentric distance bins
# 	mass_v : numpy.ndarray
# 		enclosed mass profile computed with 3d velocities
# 	mass_vr : numpy.ndarray
# 		enclosed mass profile computed with radial velocities
		
		

# 	"""
# 	l_lmc, b_lmc = mock_observations.galactic_coodinates(pos, vel)
# 	l_oct, b_oct, pos_oct, vel_oct = octant_cut(lmin, lmax, bmin, bmax, l_lmc, b_lmc, pos, vel)
# 	r_cuts, Mv_oct, Mr_oct = mass_estimator_profile(alpha, beta, gamma, pos_oct, vel_oct, r_in, r_out)
# 	return r_cuts, Mv_oct, Mr_oct	
