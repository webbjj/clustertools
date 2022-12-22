"""Functions to calculate StarCluster properties

Designed to accept StarCluster instance as in put to
calculate key parameters

"""

__author__ = "Jeremy J Webb"
__all__ = [
    "find_centre",
    "find_centre_of_density",
    "find_centre_of_mass",
    "relaxation_time",
    "half_mass_relaxation_time",
    "core_relaxation_time",
    "energies",
    "closest_star",
    "rlagrange",
    "virial_radius",
    "virial_radius_inverse_distance",
    "virial_radius_critical_density",
    "mass_function",
    "tapered_mass_function",
    "eta_function",
    "meq_function",
    "ckin",
    "rcore",
    "rtidal",
    "rlimiting",
]

import numpy as np
import numba
try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion
from galpy import potential
from galpy.potential import rtide
from scipy.optimize import curve_fit
from scipy.spatial import cKDTree

from ..util.recipes import *
from ..util.constants import _get_grav
from ..util.plots import _plot,_lplot,_scatter
from ..util.units import _convert_length,_convert_time,_convert_velocity

try:
    import amuse.units.units as u
except:
    pass

import matplotlib.pyplot as plt


def find_centre(
    cluster,
    xstart=0.0,
    ystart=0.0,
    zstart=0.0,
    vxstart=0.0,
    vystart=0.0,
    vzstart=0.0,
    indx=None,
    nsigma=1.0,
    nsphere=100,
    density=True,
    rmin=0.1,
    rmax=None,
    nmax=100,
    method='harfst',
    nneighbour=6,
):
    """Find the cluster's centre

    - The default assumes the cluster's centre is the centre of density, calculated via the find_centre_of_density function.
    - For density=False, the routine first works to identify a sphere of nsphere stars around the centre in which to perform a centre of mass calculation (similar to NBODY6). Stars beyond nsigma standard deviations are removed from the calculation until only nsphere stars remain. This step prevents long tidal tails from affecting the calculation

    Parameters
    ----------
    cluster : class
        StarCluster
    xstart,ystart,zstart : float
        starting position for centre
    vxstart,vystart,vzstart :
        starting velocity for centre
    indx : bool
        subset of stars to use when finding center
    nsigma : int
        number of standard deviations to within which to keep stars
    nsphere : int
        number of stars in centre sphere (default:100)
    density : bool
        use Yohai Meiron's centre of density calculator instead (default: True)
    rmin : float
        minimum radius to start looking for stars
    rmax : float
        maxmimum radius of sphere around which to estimate density centre (default: None cluster.units, uses maximum r)
    nmax : int
        maximum number of iterations to find centre
    method : str
        method for finding the centre of density ('harfst' (default), 'casertano')

    if method=='casertano'
        nneighbour : int
            number of neighbours for calculation local densities

    Returns
    -------
    xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

    History
    -------
    2019 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if indx is None:
        indx = np.ones(cluster.ntot, bool)
    elif np.sum(indx) == 0.0:
        print("NO SUBSET OF STARS GIVEN")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    rhos=np.zeros(cluster.ntot)

    if density:
        xc,yc,zc,vxc,vyc,vzc=find_centre_of_density(
            cluster=cluster,
            xstart=xstart,
            ystart=ystart,
            zstart=zstart,
            vxstart=vxstart,
            vystart=vystart,
            vzstart=vzstart,
            indx=indx,
            nsphere=nsphere,
            rmin=rmin,
            rmax=rmax,
            nmax=nmax,
            method=method,
            nneighbour=nneighbour,
        )
    else:

        x = cluster.x[indx] - xstart
        y = cluster.y[indx] - ystart
        z = cluster.z[indx] - zstart
        r = np.sqrt(x ** 2.0 + y ** 2.0 + z ** 2.0)
        i_d = cluster.id[indx]

        while len(r) > nsphere:
            sigma = nsigma * np.std(r)
            indx = r < sigma

            if len(r[indx]) > nsphere:
                i_d = i_d[indx]
                x = x[indx] - np.mean(x[indx])
                y = y[indx] - np.mean(y[indx])
                z = z[indx] - np.mean(z[indx])
                r = np.sqrt(x * x + y * y + z * z)
            else:
                break

        # Find centre of mass and velocity of inner stars:
        indx = np.in1d(cluster.id, i_d)

        xc = np.sum(cluster.m[indx] * cluster.x[indx]) / np.sum(cluster.m[indx])
        yc = np.sum(cluster.m[indx] * cluster.y[indx]) / np.sum(cluster.m[indx])
        zc = np.sum(cluster.m[indx] * cluster.z[indx]) / np.sum(cluster.m[indx])

        vxc = np.sum(cluster.m[indx] * cluster.vx[indx]) / np.sum(cluster.m[indx])
        vyc = np.sum(cluster.m[indx] * cluster.vy[indx]) / np.sum(cluster.m[indx])
        vzc = np.sum(cluster.m[indx] * cluster.vz[indx]) / np.sum(cluster.m[indx])

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if cluster.units =='amuse':
        return xc | u.pc, yc | u.pc, zc | u.pc, vxc | u.kms, vyc | u.kms, vzc | u.kms
    else:
        return xc, yc, zc, vxc, vyc, vzc

def find_centre_of_density(
    cluster,
    xstart=0.0,
    ystart=0.0,
    zstart=0.0,
    vxstart=0.0,
    vystart=0.0,
    vzstart=0.0,
    indx=None,
    nsphere=100,
    rmin=0.1,
    rmax=None,
    nmax=100,
    method='harfst',
    nneighbour=6,
):
    """Find cluster's centre of density

    - Find cluster's centre of density:
        if method=='harfst' use (Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357) courtesy of Yohai Meiron
        if method=='casertano' use (Casertano, S., Hut, P. 1985, ApJ, 298, 80)

    Parameters
    ----------
    cluster : class
        StarCluster
    xstart,ystart,zstart : float
        starting position for centre (default: 0,0,0)
    vxstart,vystart,vzstart : float
        starting velocity for centre (default: 0,0,0)
    indx: bool
        subset of stars to perform centre of density calculation on (default: None)
    nsphere : int
        number of stars in centre sphere (default:100)
    rmin : float
        minimum radius of sphere around which to estimate density centre (default: 0.1 cluster.units)
    rmax : float
        maxmimum radius of sphere around which to estimate density centre (default: None cluster.units, uses maximum r)
    nmax : float
        maximum number of iterations (default:100)

    method : str
        method for finding the centre of density ('harfst' (default), 'casertano')

    if method=='casertano'
        nneighbour : int
            number of neighbours for calculation local densities

    Returns
    -------
    xc,yc,zc,vxc,vyc,vzc : float
        coordinates of centre of mass

    HISTORY
    -------
    2019 - Written - Webb (UofT) with Yohai Meiron (UofT)
    2022 - Written - Webb (UofT) - add method=='casertano'
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0
    if cluster.units=='amuse':
        cluster.to_pckms()

    if method=='harfst':
        xdc, ydc, zdc,vxdc, vydc, vzdc=find_centre_of_density_harfst(cluster,xstart=xstart,
            ystart=ystart,zstart=zstart,vxstart=vxstart,vystart=vystart,vzstart=vzstart,indx=indx,
            nsphere=nsphere,rmin=rmin,rmax=rmax,nmax=nmax)
    elif method=='casertano':
        xdc, ydc, zdc,vxdc, vydc, vzdc=find_centre_of_density_casertano(cluster,nneighbour=nneighbour)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if cluster.units =='amuse':
        return xdc | u.pc, ydc | u.pc, zdc | u.pc, vxdc | u.kms, vydc | u.kms, vzdc | u.kms
    else:
        return xdc, ydc, zdc,vxdc, vydc, vzdc

def find_centre_of_density_casertano(
    cluster,
    nneighbour=6,
):

    """Find cluster's centre of density

    - The motivation behind this piece of code comes from Casertano, S., Hut, P. 1985, ApJ, 298, 80


    Parameters
    ----------
    cluster : class
        StarCluster
    nneighbour : int
        number of neighbours for calculation local densities

    Returns
    -------
    xc,yc,zc,vxc,vyc,vzc : float
        coordinates of centre of mass

    HISTORY
    -------
    2022 - Written - Webb (UofT)
    """

    pos=np.column_stack([cluster.x,cluster.y,cluster.z])
    vel=np.column_stack([cluster.vx,cluster.vy,cluster.vz])

    tree=cKDTree(pos)
    dist, arg = tree.query(pos, k=nneighbour)

    #Sum the masses but ignore last neighbour
    mass=np.sum(cluster.m[arg],axis=1)-cluster.m[arg[:,-1]]
    #Set volume equal using cube of distance to outermost neighbour
    vol=(dist[:,-1]**3.)
    rhos=mass/vol
    rhototal=np.sum(rhos)

    xdc=np.sum(rhos*cluster.x/rhototal)
    ydc=np.sum(rhos*cluster.y/rhototal)
    zdc=np.sum(rhos*cluster.z/rhototal)
    vxdc=np.sum(rhos*cluster.vx/rhototal)
    vydc=np.sum(rhos*cluster.vy/rhototal)
    vzdc=np.sum(rhos*cluster.vz/rhototal)

    return xdc, ydc, zdc,vxdc, vydc, vzdc


def find_centre_of_density_harfst(
    cluster,
    xstart=0.0,
    ystart=0.0,
    zstart=0.0,
    vxstart=0.0,
    vystart=0.0,
    vzstart=0.0,
    indx=None,
    nsphere=100,
    rmin=0.1,
    rmax=None,
    nmax=100,
):
    """Find cluster's centre of density

    - The motivation behind this piece of code comes from phigrape (Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357) courtesy of Yohai Meiron
    - The routine first finds the centre of density of the whole system, and then works to identify a sphere stars around the centre in which to perform the final centre of density calculation. Stars with radii outside 80% of the maximum radius are removed from the calculation until the final subset of stars are enclosed within a radius rmin. The maximum size of the final subset is nmax. This step prevents long tidal tails from affecting the calculation

    Parameters
    ----------
    cluster : class
        StarCluster
    xstart,ystart,zstart : float
        starting position for centre (default: 0,0,0)
    vxstart,vystart,vzstart : float
        starting velocity for centre (default: 0,0,0)
    indx: bool
        subset of stars to perform centre of density calculation on (default: None)
    nsphere : int
        number of stars in centre sphere (default:100)
    rmin : float
        minimum radius of sphere around which to estimate density centre (default: 0.1 cluster.units)
    rmax : float
        maxmimum radius of sphere around which to estimate density centre (default: None cluster.units, uses maximum r)
    nmax : float
        maximum number of iterations (default:100)

    Returns
    -------
    xc,yc,zc,vxc,vyc,vzc : float
        coordinates of centre of mass

    HISTORY
    -------
    2019 - Written - Webb (UofT) with Yohai Meiron (UofT)
    """
    if indx is None:
        indx = np.ones(cluster.ntot, bool)

    m = cluster.m[indx]
    x = cluster.x[indx] - xstart
    y = cluster.y[indx] - ystart
    z = cluster.z[indx] - zstart
    vx = cluster.vx[indx] - vxstart
    vy = cluster.vy[indx] - vystart
    vz = cluster.vz[indx] - vzstart

    r = np.sqrt(x ** 2.0 + y ** 2.0 + z ** 2.0)

    if rmax is None:
        rlim = np.amax(r)
    else:
        rlim=rmax

    xdc, ydc, zdc = xstart, ystart, zstart
    vxdc, vydc, vzdc = vxstart, vystart, vzstart

    n = 0

    while (rlim > rmin) and (n < nmax):
        r2 = x ** 2.0 + y ** 2.0 + z ** 2.0
        indx = r2 < rlim ** 2
        nc = np.sum(indx)
        mc = np.sum(m[indx])

        if mc == 0:
            xc, yc, zc = 0.0, 0.0, 0.0
            vxc, vyc, vzc = 0.0, 0.0, 0.0
        else:

            xc = np.sum(m[indx] * x[indx]) / mc
            yc = np.sum(m[indx] * y[indx]) / mc
            zc = np.sum(m[indx] * z[indx]) / mc

            vxc = np.sum(m[indx] * vx[indx]) / mc
            vyc = np.sum(m[indx] * vy[indx]) / mc
            vzc = np.sum(m[indx] * vz[indx]) / mc

        if (mc > 0) and (nc > nsphere):
            x -= xc
            y -= yc
            z -= zc
            xdc += xc
            ydc += yc
            zdc += zc

            vx -= vxc
            vy -= vyc
            vz -= vzc
            vxdc += vxc
            vydc += vyc
            vzdc += vzc

        else:
            break
        rlim *= 0.8
        n += 1

    return xdc, ydc, zdc,vxdc, vydc, vzdc


def find_centre_of_mass(cluster):
    """ Find the centre of mass of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    xc,yc,zc,vxc,vyc,vzc : float
        coordinates of centre of mass

    HISTORY
    -------
    2018 - Written - Webb (UofT)
    """
    xc = np.sum(cluster.m * cluster.x) / np.sum(cluster.m)
    yc = np.sum(cluster.m * cluster.y) / np.sum(cluster.m)
    zc = np.sum(cluster.m * cluster.z) / np.sum(cluster.m)

    vxc = np.sum(cluster.m * cluster.vx) / np.sum(cluster.m)
    vyc = np.sum(cluster.m * cluster.vy) / np.sum(cluster.m)
    vzc = np.sum(cluster.m * cluster.vz) / np.sum(cluster.m)

    return xc, yc, zc,vxc, vyc, vzc


def relaxation_time(cluster, rad=None, coulomb=0.4, projected=False,method='spitzer'):
    """Calculate the relaxation time (Spitzer & Hart 1971) within a given radius of the cluster

    - Spitzer, L. Jr, Hart, M.H. 1971, ApJ, 164, 399 (Equation 5)
    - Need to adjust amplitude for different input units
    Parameters
    ----------
    cluster : class
      StarCluster
    rad : float
      radius within which to calculate the relaxation time (defult: cluster.rm)
    coulomb : float
      Coulomb parameter (default: 0.4)
    projected : bool
      use projected values (default: False)
    method : str
      choose between Spitzer & Hart 1971 and other methods (in development)

    Returns
    -------
       trelax : float
          relaxation time within radius rad

    History
    -------
    2020 - Written - Webb (UofT)

    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre()
    else:
        cluster.sortstars()

    cluster.to_pckms()
    grav=_get_grav(cluster)


    if rad is None and projected:
        rad=cluster.rmpro
        vel=cluster.vpro
    elif rad is None:
        rad=cluster.rm
        vel=cluster.v

    if projected:
        rindx=cluster.rpro <= rad
    else:
        rindx=cluster.r <= rad
        
    ntot=np.sum(rindx)
    mbar=np.mean(cluster.m[rindx])
    vol=4.0*np.pi*(rad**3.)/3.0
    rho=ntot/vol
    
    v2=np.mean(vel**2.)
    
    #v2=0.4*grav*np.sum(cluster.m)/rad
    
    lnlambda=np.log(coulomb*cluster.ntot)
    
    trelax=v2**(3./2.)/(15.4*grav**2.*mbar**2.*rho*lnlambda)

    # Units of Myr
    trelax*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    trelax=_convert_time(trelax,'pckms',cluster)

    return trelax

def half_mass_relaxation_time(cluster, coulomb=0.4, projected=False):
    """ Calculate the half-mass relaxation time (Spitzer 1987) of the cluster
    - Spitzer, L. 1987, Dynamical evolution of globular clusters
    - Need to adjust amplitude for different input units

    Parameters
    ----------
    cluster : class
      StarCluster

    coulomb : float
      Coulomb parameter (default: 0.4)

    projected : bool
      use projected values (default: False)

    Returns
    -------
       trh : float
          half-mass relaxation time within radius rm (in Myr)

    History
    -------
       2019 - Written - Webb (UofT)

    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre()
    else:
        cluster.sortstars()

    cluster.to_pckms()
    grav=_get_grav(cluster)

    mass=np.sum(cluster.m)
    ntot=float(cluster.ntot)
    mbar=mass/ntot
    lnlambda = np.log(coulomb*ntot)

    if projected:
        rm=cluster.rmpro
    else:
        rm=cluster.rm


    trh=0.138*(mass**0.5)*(rm**1.5)/(mbar*np.sqrt(grav)*lnlambda)
    # Units of Myr
    trh*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    trh=_convert_time(trh,'pckms',cluster)


    return trh


def core_relaxation_time(cluster, coulomb=0.4, projected=False):
    """ Calculate the core relaxation time (Stone & Ostriker 2015) of the cluster
    
    - Stone, N.C. & Ostriker, J.P. 2015, ApJ, 806, 28

    Parameters

    cluster : class
      StarCluster

    coulomb : float
      Coulomb parameter (default: 0.4)

    projected : bool
      use projected values (default: False)

    method : str
      choose between Stone & Ostriker 2015 and other methods (in development)

    Returns

     trc (in Myr)

    History
    -------
    2019 - Written - Webb (UofT)

    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre()
    else:
        cluster.sortstars()

    cluster.to_pckms()
    grav=_get_grav(cluster)

    lnlambda=np.log(coulomb*cluster.ntot)
    mtot=np.sum(cluster.m)
    mbar=np.mean(cluster.m)
    if projected:
        rc=cluster.rcore(projected=True)
        rh=cluster.rmpro
    else:
        rc=cluster.rcore()
        rh=cluster.rm

    trc=(0.39/lnlambda)*np.sqrt(rc**3./(grav*mtot))*(mtot/mbar)*np.sqrt(rc*rh)/(rc+rh)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    trc=_convert_time(trc,'pckms',cluster)


    return trc


def energies(cluster, specific=True, ids=None, projected=False, softening=0.0, full=True, parallel=False, **kwargs):
    """Calculate kinetic and potential energy of every star
    Parameters
    ----------
    cluster : class
      StarCluster instance
    specific : bool
      find specific energies (default: True)
    ids: bool or int
      if given, find the energues of a subset of stars defined either by an array of
      star ids, or a boolean array that can be used to slice the cluster. (default: None)
    projected : bool
      use projected values (default: False)
    softening : float
      Plummer softening length in cluster.units (default: 0.0)
    full : bool
      calculate distance of full array of stars at once with numbra (default: True)
    parallel : bool
      calculate distances in parallel (default: False)
    Returns
    -------
    kin,pot : float
      kinetic and potential energy of every star if the i_d argument is not used. If i_d
      argument is used, return an arrays with potential and kinetic energy in the same shape
      of i_d
    History
    -------
       2019 - Written - Webb (UofT)
       2022 - Updated with support for multiple ids or an idexing array - Erik Gillis (UofT)
    """

    if ids is None:
        ids=kwargs.get('i_d',None)

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_cluster(sortstars=False)
    if cluster.units=='amuse':
        cluster.to_pckms()

    grav=_get_grav(cluster)

    if projected:
        if specific:
            kin = 0.5 * (cluster.vpro ** 2.0)
        else:
            kin = 0.5 * cluster.m * (cluster.vpro ** 2.0)
    else:

        if specific:
            kin = 0.5 * (cluster.v ** 2.0)
        else:
            kin = 0.5 * cluster.m * (cluster.v ** 2.0)

    if ids is not None:
        
        # Convert ids to boolean array if given as an array of ids
        if isinstance(ids,int) or isinstance(ids,float) or isinstance(ids,np.int64) or isinstance(ids,np.float64):
            ids = cluster.id == ids
        elif isinstance(ids[0],int) or isinstance(ids[0],float) or isinstance(ids[0],np.int64) or isinstance(ids[0],np.float64):
            ids = np.in1d(cluster.id, ids)
    
        # Get gravitational constant
        grav = _get_grav(cluster)

        kin = 0.5 * cluster.m[ids] * cluster.v[ids]**2

        cluster_full = np.array([cluster.x, cluster.y, cluster.z, cluster.m]).T
        cluster_sub  = np.array([cluster.x[ids], cluster.y[ids], 
                                 cluster.z[ids], cluster.m[ids]]).T

        if parallel:
            pot = grav * np.array(_potential_energy_subset_parallel(cluster_sub, cluster_full,softening))
        else:
            pot = grav * np.array(_potential_energy_subset(cluster_sub, cluster_full,softening))

        if specific:
            pot /= cluster.m[ids]
            kin /= cluster.m[ids]
    
    elif full:
        if projected:
            x = np.array([cluster.x, cluster.y, np.zeros(len(cluster.x)), cluster.m]).T
        else:
            x = np.array([cluster.x, cluster.y, cluster.z, cluster.m]).T
        if parallel:
            pot = grav * np.array(_potential_energy_parallel(x,softening))
        else:
            pot = grav * np.array(_potential_energy(x,softening))

        if specific:
            pot /= cluster.m
    else:
        pot = []

        for i in range(0, cluster.ntot):
            dx = cluster.x[i] - cluster.x
            dy = cluster.y[i] - cluster.y
            dz = cluster.z[i] - cluster.z
            if specific:
                m = cluster.m
            else:
                m = cluter.m[i] * cluster.m

            if projected:
                dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + softening**2.) 
            else:
                dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + dz ** 2.0 + softening**2.)

            indx = dr != 0.0
            gmr = -grav * m[indx] / dr[indx]

            pot.append(np.sum(gmr))

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if cluster.units=='amuse' and specific:
        kin=kin | (u.kms*u.kms)
        pot=pot | (u.kms*u.kms)
    elif cluster.units=='amuse' and not specific:
        kin=kin | (u.MSun*u.kms*u.kms)
        pot=pot | (u.MSun*u.kms*u.kms)
    return kin, pot


@numba.njit
def _potential_energy(cluster,softening=0.0):
    """Find potential energy for each star in a cluster
    - uses numba

    Parameters
    ----------
    cluster : float
        positions and masses of stars within the StarCluster
    softening : float
      Plummer softening length in cluster.units (default: 0.0)
    Returns
    -------
        pot : float
            potential energy of every star

    History
    -------
       2019 - Written - Webb (UofT)
    """
    pot = [0.0] * len(cluster)
    for i in range(len(cluster) - 1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i], cluster[j])

            if r==0: r=np.nan

            m2 = cluster[i, 3] * cluster[j, 3]
            pot[i] += -m2 / np.sqrt(r**2.+softening**2.)
            pot[j] += -m2 / np.sqrt(r**2.+softening**2.)

    return pot


@numba.njit(parallel=True)
def _potential_energy_parallel(cluster,softening=0.0):
    """Find potential energy for each star in a cluster in parallel
    - uses numba

    Parameters
    ----------
    cluster : class
        positions and masses of stars within the StarCluster
    softening : float
      Plummer softening length in cluster.units (default: 0.0)
    Returns
    -------
        pot : float
            potential energy of every star

    History
    -------
       2019 - Written - Webb (UofT)
    """
    pot = [0.0] * len(cluster)
    for i in numba.prange(len(cluster) - 1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i], cluster[j])

            if r==0: r=np.nan

            m2 = cluster[i, 3] * cluster[j, 3]
            pot[i] += -m2 / np.sqrt(r**2.+softening**2.)
            pot[j] += -m2 / np.sqrt(r**2.+softening**2.)

    return pot


@numba.njit()
def _potential_energy_subset(cluster_sub, cluster_full,softening=0.0):
    """Find the potential energy for a subset of stars in a bigger cluster
    
    Parameters
    ----------
    cluster_sub : float
        2d numpy array with x, y, z and mass at each index of the sub cluster
    cluster_full : float
        2d numpy array with x, y, z and mass at each index of the cluster
    softening : float
      Plummer softening length in cluster.units (default: 0.0)
    Returns:
    --------
        potential : numpy array
            array of the potential energy of each star in the subcluster
            
    History
    -------
        2022 - Written - Erik Gillis (UofT)
    """
    pot = [0.0] * len(cluster_sub)
    
    for i in range(len(cluster_sub)):
        for j in range(len(cluster_full)):
            r = distance(cluster_sub[i], cluster_full[j])

            if r!=0:

                m = cluster_sub[i,3] * cluster_full[j,3]

                pot[i] += -m / np.sqrt(r**2.+softening**2.)

        
    return pot
    
@numba.njit()
def _potential_energy_subset_parallel(cluster_sub, cluster_full,softening=0.0):
    """Find the potential energy for a subset of stars in a bigger cluster
    
    Parameters
    ----------
    cluster_sub : float
        2d numpy array with x, y, z and mass at each index of the sub cluster
    cluster_full : float
        2d numpy array with x, y, z and mass at each index of the cluster
    softening : float
      Plummer softening length in cluster.units (default: 0.0)    
    Returns:
    --------
        potential : numpy array
            array of the potential energy of each star in the subcluster
            
    History
    -------
        2022 - Written - Erik Gillis (UofT)
    """
    pot = [0.0] * len(cluster_sub)
    
    for i in numba.prange(len(cluster_sub)):
        for j in range(len(cluster_full)):
            r = distance(cluster_sub[i], cluster_full[j])

            if r!=0:

                m = cluster_sub[i,3] * cluster_full[j,3]

                pot[i] += -m / np.sqrt(r**2.+softening**2.)
        
    return pot


def closest_star(cluster, projected=False, argument=False):
    """Find distance to closest star for each star
    - uses numba

    Parameters
    ----------
    cluster : class
        positions of stars within the StarCluster
    projected : bool
      use projected values (default: False)
    argument : bool
      return argument of closest star as well (default: False)

    Returns
    -------
        minimum_distance : float
            distance to closest star for each star

        if argument:
            arg : int
                argument of closest star for each star            

    History
    -------
       2019 - Written - Webb (UofT)
    """

    if projected:
        x=np.column_stack([cluster.x,cluster.y])
    else:
        x=np.column_stack([cluster.x,cluster.y,cluster.z])

    tree=cKDTree(x)
    dist, arg = tree.query(x, k=2)

    if argument:
        return dist[:,1],arg[:,1]
    else:
        return dist[:,1]

def rlagrange(cluster, nlagrange=10, mfrac=None, projected=False):
    """Calculate lagrange radii of the cluster by mass

    Parameters
    ----------
    cluster : class
        StarCluster
    nlagrange : int 
        number of lagrange radii bins (default: 10)
    mfrac : float
        Exact masss fraction to calculate radius. Will supersede nlagrange if not None (default : None)
    projected : bool
        calculate projected lagrange radii (default: False)

    Returns
    -------
    rn : float
        lagrange radii

    History
    -------
       2019 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre()
    else:
        cluster.sortstars()


    #Array for Lagrange radii
    rn = []

    if projected:
        rorder = cluster.rproorder
        r=cluster.rpro
    else:
        rorder = cluster.rorder
        r=cluster.r

    msum = np.cumsum(cluster.m[rorder])

    if mfrac is None:

        for i in range(1, nlagrange):
            indx = msum >= np.sum(cluster.m) * float(i) / float(nlagrange)
            rn.append(r[rorder[indx][0]])

        while len(rn) != nlagrange:
            rn.append(np.max(r))

    else:
        indx=(msum/cluster.mtot)>=mfrac
        rn=r[rorder][indx][0]

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if cluster.units=='amuse':
        rn=_convert_length(rn,'pckms',cluster)


    return rn

def virial_radius(cluster, method='inverse_distance',
    full=True,
    H=70.0,
    Om=0.3,
    overdens=200.0,
    projected=False,
    plot=False,
    **kwargs):
    """Calculate virial radius of the cluster
    - Virial radius is calculated using either:
    -- the average inverse distance between particles, weighted by their masses (default)
    -- the radius at which the density is equal to the critical density of the Universe at the redshift of the system, multiplied by an overdensity constant
    Parameters
    ----------
    cluster : class
        StarCluster
    method : str
        method for calculating virial radius (default: 'inverse_distance', alternative: 'critical_density')

    Returns
    -------
    rv : float
        virial radius

    Other Parameters
    ----------------
    full : bool
        Use Numba to calculate average inverse distance between stars (default:True)
    H : float
        Hubble constant
    Om : float
        density of matter
    overdens : float
        overdensity constant
    projected : bool
        calculate projected virial radius (default: False)
    plot : bool
        plot cluster density profile and illustrate virial radius calculation
    kwargs : str
        key word arguments for plotting function

    History
    -------
       2019 - Written - Webb (UofT)
    """

    if method=='inverse_distance':
        rv=virial_radius_inverse_distance(cluster,projected=projected,full=full)
    elif method=='critical_density':
        rv=virial_radius_critical_density(cluster,H,Om,overdens,projected,plot,**kwargs)

    return rv

def virial_radius_inverse_distance(cluster, projected=False, full=True):
    """ Calculate virial radius of the cluster 
    - Virial radius is defined as the inverse of the average inverse distance between particles, weighted by their masses
    -- Definition taken from AMUSE (www.amusecode.org)
    -- Portegies Zwart S., McMillan S., 2018, Astrophysical Recipes; The art ofAMUSE, doi:10.1088/978-0-7503-1320-9

    Parameters
    ----------
    cluster : class
        StarCluster
    projected : bool
        calculate projected virial radius (default: False)
    full : bool
        Use Numba to calculate average inverse distance between stars (default:True)

    Returns
    -------
    r_v : float
        virial radius

    History
    -------
    2019 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse': cluster.to_pckms()

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre(sortstars=False)


    if full:
        if projected:
            x = np.array([cluster.x, cluster.y, np.zeros(cluster.ntot), cluster.m]).T
        else:
            x = np.array([cluster.x, cluster.y, cluster.z, cluster.m]).T

        ms = cluster.m
        partial_sum = _weighted_inverse_distance_sum(x)

    else:
        partial_sum = 0.0

        ms = cluster.m
        xs = cluster.x
        ys = cluster.y
        if projected:
            zs = np.zeros(cluster.ntot)
        else:
            zs = cluster.z

        for i in range(cluster.ntot - 1):
            x = xs[i]
            y = ys[i]
            z = zs[i]
            dx = x - xs[i + 1 :]
            dy = y - ys[i + 1 :]
            dz = z - zs[i + 1 :]
            dr2 = (dx * dx) + (dy * dy) + (dz * dz)
            dr = np.sqrt(dr2)
            m_m = ms[i] * ms[i + 1 :]
            partial_sum += np.sum(m_m / dr)

    r_v=(np.sum(ms) ** 2) / (2 * partial_sum)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if cluster.units=='amuse':
        r_v=_convert_length(r_v,'pckms',cluster)

    return r_v

@numba.njit
def _weighted_inverse_distance_sum(cluster):
    """Find the sum of the mass weighted inverse distance for each star

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    weighted_sum : float
        sum of the mass weighted inverse distance for each star

    History
    -------
    2019 - Written - Webb (UofT)
    """
    weighted_sum = 0.0
    for i in range(len(cluster) - 1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i], cluster[j])
            m2 = cluster[i, 3] * cluster[j, 3]
            weighted_sum += m2 / r

    return weighted_sum


def virial_radius_critical_density(
    cluster,
    H=70.0,
    Om=0.3,
    overdens=200.0,
    projected=False,
    plot=False,
    **kwargs
):
    """Calculate virial radius of the cluster
    
    - Virial radius is defined as the radius at which the density is equal to the critical density of the Universe at the redshift of the system, multiplied by an overdensity constant
    - Note that this a quick method that is a bit of an approximation as it interpolates the cluster's density profile. A more accurate (but expensive)
    approach would be to subtract the product of the critical density and the overdensity constant from the density profile and find the root (in development)

    Parameters
    ----------
    cluster : class
        StarCluster
    H : float
        Hubble constant
    Om : float
        density of matter
    overdens : float
        overdensity constant
    projected : bool
        calculate projected virial radius (default: False)
    plot : bool
        plot cluster density profile and illustrate virial radius calculation

    Returns
    -------
    r_v : float
        virial radius

    Other Parameters
    ----------------
    kwargs : str
        key word arguments for plotting function

    History
    -------
    2019 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre()
    else:
        cluster.sortstars()

    cluster.to_pckms()

    rhocrit=conversion.dens_in_msolpc3(vo=cluster._vo,ro=cluster._ro)/conversion.dens_in_criticaldens(vo=cluster._vo,ro=cluster._ro,H=H)

    if projected:
        if not cluster.projected: cluster.analyze(sortstars=True,projected=True)
        indx = cluster.rproorder
    else:
        indx = cluster.rorder

    msum = np.cumsum(cluster.m[indx])

    if projected:
        vsum = (4.0 / 3.0) * np.pi * (cluster.rpro[indx] ** 3.0)
        pprof = msum / vsum
        rprof = cluster.rpro[indx]
    else:
        vsum = (4.0 / 3.0) * np.pi * (cluster.r[indx] ** 3.0)
        pprof = msum / vsum
        rprof = cluster.r[indx]

    rho_local = rhocrit * overdens

    rindx=np.argmin(np.fabs(pprof-rho_local))

    if rindx==len(pprof)-1 or rindx==0:
        r_v=rprof[rindx]
    elif pprof[rindx]==rho_local:
        r_v=rprof[rindx]
    else:
        if pprof[rindx] < rho_local:
            r1=rprof[rindx-1]
            r2=rprof[rindx]
            p1=pprof[rindx-1]
            p2=pprof[rindx]
        elif pprof[rindx] > rho_local:
            r1=rprof[rindx]
            r2=rprof[rindx+1]
            p1=pprof[rindx]
            p2=pprof[rindx+1]

        m=(p2-p1)/(r2-r1)
        b=p2-m*r2

        r_v=(rho_local-b)/m

    if plot:

        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        xunits = " (pc)"
        if projected:
            yunits = " Msun/pc^2"
        else:
            yunits = " Msun/pc^3"

        x, y = rprof, pprof
        _lplot(
            x,
            y,
            xlabel=r"$R" + xunits + "$",
            ylabel=r"$rho" + yunits + "$",
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )
        _lplot(x, np.ones(len(x)) * rho_local, "--", overplot=True)
        _lplot(np.ones(len(y)) * r_v, y, "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    r_v=_convert_length(r_v,'pckms',cluster)

    return r_v


def mass_function(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    npop=None,
    indx=None,
    projected=False,
    mcorr=None,
    plot=False,
    **kwargs
):
    """Find mass function over a given mass range

    - mass bins are set up so that there are an equal number of stars in each bin

    Parameters
    ----------
    cluster : class
        StarCluster instance
    mmin/mmax : float
        specific mass range
    nmass : 
        number of mass bins used to calculate alpha
    rmin/rmax : 
        specific radial range
    vmin/vmax : float
        specific velocity range
    emin/emax : float
        specific energy range
    kwmin/kwmax : int
        specific stellar evolution type range
    npop : int
        population number
    indx : bool 
        specific subset of stars
    projected : bool 
        use projected values (default: False)
    mcorr : float
        completeness correction for masses
    plot : bool 
        plot the mass function

    Returns
    -------
    m_mean : float
        mean mass in each bin
    m_hist : float
        number of stars in each bin
    dm : float
        dN/dm of each bin
    alpha : float
        power-law slope of the mass function (dN/dm ~ m^alpha)
    ealpha : float
        error in alpha
    yalpha : float
        y-intercept of fit to log(dN/dm) vs log(m)
    eyalpha : float
        error in yalpha

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    """
    if rmin == None:
        rmin = np.min(r)
    if rmax == None:
        rmax = np.max(r)
    if vmin == None:
        vmin = np.min(v)
    if vmax == None:
        vmax = np.max(v)
    if mmin == None:
        mmin = np.min(cluster.m)
    if mmax == None:
        mmax = np.max(cluster.m)

    if indx is None:
        indx = cluster.id > -1

    # Build subcluster containing only stars in the full radial and mass range:
    indx *= (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m < mmax)
        * (v >= vmin)
        * (v <= vmax)
    )

    if len(cluster.kw) > 0:
        indx *= (cluster.kw >= kwmin) * (cluster.kw <= kwmax)

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax
    """

    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)


    if mcorr is None: 
        mcorr=np.ones(cluster.ntot)
        return_error=False
    else:
        return_error=True


    if np.sum(indx) >= nmass:

        if kwargs.get('bintype','num')=='fix' or kwargs.get('mbintype','num')=='fix':
            m_lower, m_mean, m_upper, m_hist = binmaker(cluster.m[indx], nmass)
        else:
            m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)

        m_corr_hist = np.zeros(len(m_hist))
        for i in range(0, len(m_hist)):
            mindx = (cluster.m >= m_lower[i]) * (cluster.m < m_upper[i]) * indx
            m_hist[i]=np.sum(mindx)
            m_corr_hist[i] = np.sum(1.0 / mcorr[mindx])

        mbinerror = m_hist / m_corr_hist

        lm_mean = np.log10(m_mean)
        dm = m_corr_hist / (m_upper - m_lower)
        ldm = np.log10(dm)

        (alpha, yalpha), V = np.polyfit(lm_mean, ldm, 1, cov=True)
        ealpha = np.sqrt(V[0][0])
        eyalpha = np.sqrt(V[1][1])

        if plot:
            filename = kwargs.get("filename", None)
            _plot(m_mean, np.log10(dm), xlabel="M", ylabel="LOG(dN/dM)",**kwargs)
            mfit = np.linspace(np.min(m_mean), np.max(m_mean), nmass)
            dmfit = 10.0 ** (alpha * np.log10(mfit) + yalpha)
            _lplot(
                mfit, np.log10(dmfit), overplot=True, label=(r"$\alpha$ = %f" % alpha),
            )

            plt.legend()

            if filename != None:
                plt.savefig(filename)
        if return_error:
            return m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbinerror
        else:
            return m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha
    else:
        print("NOT ENOUGH STARS TO ESTIMATE MASS FUNCTION")

        if return_error:
            return (
                np.zeros(nmass),
                np.zeros(nmass),
                np.zeros(nmass),
                -1000.0,
                -1000.0,
                -1000.0,
                -1000.0,
                np.zeros(nmass),
            )
        else:

            return (
                np.zeros(nmass),
                np.zeros(nmass),
                np.zeros(nmass),
                -1000.0,
                -1000.0,
                -1000.0,
                -1000.0,
            )

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


def tapered_mass_function(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    npop=None,
    indx=None,
    projected=False,
    mcorr=None,
    plot=False,
    **kwargs
):
    """Find a tapered mass function over a given mass range

    - mass bins are set up so that there are an equal number of stars in each bin
    - functional form of the tapered mass function is taken from De Marchi, Paresce & Portegies Zwart 2010
    Parameters
    ----------
    cluster : class
        StarCluster instance
    mmin/mmax : float
        specific mass range
    nmass : 
        number of mass bins used to calculate alpha
    rmin/rmax : 
        specific radial range
    vmin/vmax : float
        specific velocity range
    emin/emax : float
        specific energy range
    kwmin/kwmax : int
        specific stellar evolution type range
    npop : int
        population number
    indx : bool 
        specific subset of stars
    projected : bool 
        use projected values (default: False)
    mcorr : float
        completeness correction for masses
    plot : bool 
        plot the mass function

    Returns
    -------
    m_mean : float
        mean mass in each bin
    m_hist : float
        number of stars in each bin
    dm : float
        dN/dm of each bin
    alpha : float
        power-law slope of the mass function (dN/dm ~ m^alpha)
    ealpha : float
        error in alpha
    yalpha : float
        y-intercept of fit to log(dN/dm) vs log(m)
    eyalpha : float
        error in yalpha

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    """
    if rmin == None:
        rmin = np.min(r)
    if rmax == None:
        rmax = np.max(r)
    if vmin == None:
        vmin = np.min(v)
    if vmax == None:
        vmax = np.max(v)
    if mmin == None:
        mmin = np.min(cluster.m)
    if mmax == None:
        mmax = np.max(cluster.m)

    if indx is None:
        indx = cluster.id > -1

    # Build subcluster containing only stars in the full radial and mass range:
    indx *= (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m < mmax)
        * (v >= vmin)
        * (v <= vmax)
    )

    if len(cluster.kw) > 0:
        indx *= (cluster.kw >= kwmin) * (cluster.kw <= kwmax)

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax
    """

    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)


    if mcorr is None: mcorr=np.ones(cluster.ntot)

    if np.sum(indx) >= nmass:

        m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)

        m_corr_hist = np.zeros(len(m_hist))
        for i in range(0, len(m_hist)):
            mindx = (cluster.m >= m_lower[i]) * (cluster.m < m_upper[i]) * indx
            m_hist[i]=np.sum(mindx)
            m_corr_hist[i] = np.sum(1.0 / mcorr[mindx])

        mbinerror = m_hist / m_corr_hist

        lm_mean = np.log10(m_mean)
        dm = m_corr_hist / (m_upper - m_lower)
        ldm = np.log10(dm)

        lower_bounds=kwargs.get('lower_bounds',[0.,-1.*np.inf,np.amin(cluster.m),-1.*np.inf])
        upper_bounds=kwargs.get('upper_bounds',[np.inf,np.inf,np.amax(cluster.m),np.inf])

        (A, alpha, mc, beta), V=curve_fit(tpl_func,10.0**np.array(lm_mean),10.0**np.array(ldm) ,bounds=(lower_bounds,upper_bounds))

        eA = np.sqrt(V[0][0])
        ealpha = np.sqrt(V[1][1])
        emc = np.sqrt(V[2][2])
        ebeta = np.sqrt(V[3][3])

        if plot:
            filename = kwargs.get("filename", None)
            _plot(m_mean, np.log10(dm), xlabel="M", ylabel="LOG(dN/dM)",**kwargs)
            mfit = np.linspace(np.min(m_mean), np.max(m_mean), nmass)

            dmfit = tpl_func(mfit,A,alpha,mc,beta)

            _lplot(
                mfit, np.log10(dmfit), overplot=True, label=(r"$\alpha = %f, mc = %f, \beta = %f$" % (alpha,mc,beta)),
            )

            plt.legend()

            if filename != None:
                plt.savefig(filename)

        return m_mean, m_hist, dm, A, eA, alpha, ealpha, mc, emc, beta, ebeta
    else:
        print("NOT ENOUGH STARS TO ESTIMATE MASS FUNCTION")
        return (
            np.zeros(nmass),
            np.zeros(nmass),
            np.zeros(nmass),
            -1000.0,
            -1000.0,
            -1000.0,
            -1000.0,
        )

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

def tpl_func(m,A,alpha,mc,beta):

    dm=A*(m**alpha)*(1.0-np.exp(-1.*(m/mc)**beta))

    return dm

def eta_function(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    npop=None,
    indx=None,
    projected=False,
    plot=False,
    meq=False,
    **kwargs
):
    """
    NAME: Find power slope of velocity dispersion versus mass
    
    - mass bins are set up so that there are an equal number of stars in each bin

    Parameters
    ----------
    cluster : class
        StarCluster instance
    mmin/mmax : float
        specific mass range
    nmass : 
        number of mass bins used to calculate alpha
    rmin/rmax : 
        specific radial range
    vmin/vmax : float
        specific velocity range
    emin/emax : float
        specific energy range
    kwmin/kwmax : int
        specific stellar evolution type range
    npop : int
        population number
    indx : bool 
        specific subset of stars
    projected : bool 
        use projected values (default: False)
    plot : bool 
        plot the mass function

    Returns
    -------
    m_mean : float
        mean mass in each bin
    sigvm : float
        velocity dispersion of stars in each bin
    eta : float
        power-law slope of (sigvm ~ m^eta)
    eeta : float
        error in eta
    yeta : float
        y-intercept of fit to log(sigvm) vs log(m)
    eeta : float
        error in yeta

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    """
    if rmin == None:
        rmin = np.min(r)
    if rmax == None:
        rmax = np.max(r)
    if vmin == None:
        vmin = np.min(v)
    if vmax == None:
        vmax = np.max(v)
    if mmin == None:
        mmin = np.min(cluster.m)
    if mmax == None:
        mmax = np.max(cluster.m)

    if indx is None:
        indx = cluster.id > -1

    # Build subcluster containing only stars in the full radial and mass range:
    indx *= (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m <= mmax)
        * (v >= vmin)
        * (v <= vmax)
    )

    if len(cluster.kw) > 0:
        indx *= (cluster.kw >= kwmin) * (cluster.kw <= kwmax)

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax
    """

    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)

    if np.sum(indx) >= 2 * nmass:

        m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)
        lm_mean = np.log10(m_mean)

        sigvm = []
        lsigvm = []
        for i in range(0, nmass):

            mindx = indx * (cluster.m >= m_lower[i]) * (cluster.m < m_upper[i])
            sigvm.append(np.std(v[mindx]))
            lsigvm.append(np.log10(sigvm[-1]))

        if meq:
            (eta, yeta), V=curve_fit(meq_func,10.0**np.array(lm_mean),10.0**np.array(lsigvm),bounds=([np.amin(cluster.m),0.],[np.amax(cluster.m),np.inf]))
        else:
            (eta, yeta), V = np.polyfit(lm_mean, lsigvm, 1, cov=True)

        eeta = np.sqrt(V[0][0])
        eyeta = np.sqrt(V[1][1])

        if plot:
            filename = kwargs.get("filename", None)
            _plot(m_mean, np.log10(sigvm), xlabel="M", ylabel=r"$\log_{10} \\ \sigma_v$", **kwargs)
            mfit = np.linspace(np.min(m_mean), np.max(m_mean), nmass)

            if meq:
                sigfit=meq_func(mfit,eta,yeta)
                _lplot(mfit, np.log10(sigfit), overplot=True, label=(r"$m_{eq}$ = %f" % eta))
            else:
                sigfit = 10.0 ** (eta * np.log10(mfit) + yeta)
                _lplot(mfit, np.log10(sigfit), overplot=True, label=(r"$\eta$ = %f" % eta))

            plt.legend()

            if filename != None:
                plt.savefig(filename)

        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        return m_mean, sigvm, eta, eeta, yeta, eyeta
    else:
        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        print("NOT ENOUGH STARS TO ESTIMATE SIGMA-MASS RELATION")
        return (
            np.zeros(nmass),
            np.zeros(nmass),
            -1000.0,
            -1000.0,
            -1000.0,
            -1000.0,
        )

def meq_func(m,meq,sigma0):
    sigma_eq=sigma0*np.exp(-0.5)
    if isinstance(m,float):
        if m<=meq:
            sigma=sigma0*np.exp(-0.5*m/meq)
        else:
            sigma=sigma_eq*((m/meq)**(-0.5))
    else:
        indx=m<=meq
        sigma=np.zeros(len(m))
        sigma[indx]=sigma0*np.exp(-0.5*m[indx]/meq)
        indx=m>meq

        sigma[indx]=sigma_eq*((m[indx]/meq)**(-0.5))

    return sigma

def meq_function(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    npop=1,
    indx=None,
    projected=False,
    plot=False,
    **kwargs
):
    """
    NAME: Find meq from velocity dispersion versus mass
    
    - mass bins are set up so that there are an equal number of stars in each bin
    - As per Bianchini, P. et al. 2016, MNRAS, 458, 3644, velocity dispersion 
      versus mass is fit with the following:
      sigma(m)= sigma e^(-1/2 m/meq) if m<= meq
              = sigma0 e^(-1/2) (m/meq)^-1/2 if m > meq

    Parameters
    ----------
    cluster : class
        StarCluster instance
    mmin/mmax : float
        specific mass range
    nmass : 
        number of mass bins used to calculate alpha
    rmin/rmax : 
        specific radial range
    vmin/vmax : float
        specific velocity range
    emin/emax : float
        specific energy range
    kwmin/kwmax : int
        specific stellar evolution type range
    npop : int
        population number
    indx : bool 
        specific subset of stars
    projected : bool 
        use projected values (default: False)
    plot : bool 
        plot the mass function

    Returns
    -------
    m_mean : float
        mean mass in each bin
    sigvm : float
        velocity dispersion of stars in each bin
    meq : float
        Bianchini fit to sigvm vs m
    emeq : float
        error in Bianchini fit to sigvm vs m
    sigma0 : float
        Bianchini fit to sigvm vs m
    esigma0 : float
        error in Bianchini fit to sigvm vs m

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2020
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    m_mean, sigvm, meq, emq, sigma0, esigma0 = eta_function(cluster,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rmin,
            rmax=rmax,
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            npop=npop,
            indx=indx,
            projected=projected,
            plot=plot,
            meq=True,
            **kwargs
        )

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    return m_mean, sigvm, meq, emq, sigma0, esigma0

def ckin(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    npop=None,
    indx=None,
    projected=False,
    **kwargs,
):
    """
    NAME: Find the kinematic concentration parameter ck
    
    - see Bianchini et al. 2018, MNRAS, 475, 96

    Parameters
    ----------
    cluster : class
        StarCluster instance
    mmin/mmax : float
        specific mass range
    nmass : 
        number of mass bins used to calculate alpha
    rmin/rmax : 
        specific radial range
    vmin/vmax : float
        specific velocity range
    emin/emax : float
        specific energy range
    kwmin/kwmax : int
        specific stellar evolution type range
    npop : int
        population number
    indx : bool 
        specific subset of stars
    projected : bool 
        use projected values (default: False)

    Returns
    -------
    ck : float
        kinematic concentration

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2020
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    rn=rlagrange(cluster, nlagrange=10, projected=projected)

    m_mean50, sigvm50, meq50, emq50, sigma050, esigma050 = eta_function(cluster,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rn[3],
            rmax=rn[5],
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            npop=npop,
            indx=indx,
            projected=projected,
            plot=False,
            meq=True,
            **kwargs,
        )

    m_mean, sigvm, meq, emq, sigma0, esigma0 = eta_function(cluster,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=0.,
            rmax=rn[4],
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            npop=npop,
            indx=indx,
            projected=projected,
            plot=False,
            meq=True,
            **kwargs,
        )

    ck=meq/meq50

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return ck

def rcore(
    cluster,
    method='casertano',
    nneighbour=6,
    mfrac=0.1,
    projected=False,
    plot=False,
    **kwargs
):
    """Calculate core radius of the cluster
    -- The default method (method='casertano') follows Casertano, S., Hut, P. 1985, ApJ, 298, 80 to find the core
    -- An alternative metrhod (method=='isothermal') assumes the cluster is an isothermal sphere the core radius is where density drops to 1/3 central value
    --- For projected core radius, the core radius is where the surface density profile drops to 1/2 the central value
    --- Note that the inner mass fraction of stars used to calculate central density is set by mfrac (default 0.1 = 10%)

    Parameters
    ----------
    cluster : class
        StarCluster instance
    method : string
        method of calculating the core radius of a star cluster (default 'casertano')
    if method =='casertano':
        nneighbour : int
            number of neighbours for calculation local densities
    if method=='isothermal':
        mfrac : float
            inner mass fraction to be used to establish the central density
    projected : bool
        use projected values (default: False)
    plot : bool
        plot the density profile and mark the core radius of the cluster (default: False)
    Returns
    -------
    rc : float
        core radius

    Other Parameters
    ----------------
    None

    History
    -------
    2021 - Written - Webb (UofT)
    2022 - Written - Webb (UofT) - add method='casertano'
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if cluster.origin0 != 'centre':
        if cluster.xc==0. and cluster.yc==0. and cluster.zc==0.0:
            if method=='casertano':
                cluster.find_centre(method='casertano')
            else:
                cluster.find_centre()

        cluster.to_centre(sortstars=True)

    if method=='casertano':
        
        if projected:
            r=cluster.rpro
            pos=np.column_stack([cluster.x,cluster.y])
        else:
            r=cluster.r
            pos=np.column_stack([cluster.x,cluster.y,cluster.z])

        tree=cKDTree(pos)
        dist, arg = tree.query(pos, k=nneighbour)

        mass=np.sum(cluster.m[arg],axis=1)-cluster.m[arg[:,-1]]
        if projected:
            area=dist[:,-1]**2.
            rhos=mass/area
        else:
            vol=dist[:,-1]**3.
            rhos=mass/vol

        rc2=np.sum((rhos**2.)*(r**2.))/np.sum((rhos**2.))
        rc=np.sqrt(rc2)

        if plot:
            nrad=int(np.ceil(1./mfrac)) 
            rprof,pprof,nprof=_rho_prof(cluster,nrad=nrad,projected=projected,plot=False,**kwargs)
            rcindx=(cluster.r<rc)
            mc=np.sum(cluster.m[rcindx])
            volc=(4./3.)*np.pi*(rc**3.)
            rho_c=mc/volc

        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)
        if cluster.units=='amuse': rc=_convert_length(rc,'pckms',cluster)

    else:

        mo=conversion.mass_in_msol(ro=cluster._ro,vo=cluster._vo)
        dens_in_msolpc2=(mo/cluster._ro**2.)/(1000.0**2.)

        cluster.to_pckms()

        if projected:
            r=cluster.rpro
            rorder=cluster.rproorder
            v=cluster.vpro
        else:
            r=cluster.r
            rorder=cluster.rorder
            v=cluster.v

        rcentral=rlagrange(cluster,mfrac=mfrac,projected=projected)


        nrad=int(np.ceil(1./mfrac))

        rprof,pprof,nprof=_rho_prof(cluster,nrad=nrad,projected=projected,plot=False,**kwargs)


        #interpolate

        if projected:
            rho_c=0.5*pprof[0]
        else:
            rho_c=pprof[0]/3.

        rindx=pprof < rho_c

        r1=rprof[np.invert(rindx)][-1]
        r2=rprof[rindx][0]

        p1=pprof[np.invert(rindx)][-1]
        p2=pprof[rindx][0]
       
        m=(p2-p1)/(r2-r1)
        b=p2-m*r2

        rc=(rho_c-b)/m


        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        rc=_convert_length(rc,'pckms',cluster)

    if plot:

        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            rprof /= cluster.rbar

            if projected:
                pprof *= (
                    (cluster.rbar ** 2.0)
                    / cluster.zmbar
                )

                rho_c *= (
                    (cluster.rbar ** 2.0)
                    / cluster.zmbar
                )

            else:
                pprof *= (
                    (cluster.rbar ** 3.0)
                    / cluster.zmbar
                )

                rho_c *= (
                    (cluster.rbar ** 3.0)
                    / cluster.zmbar
                )

            xunits = " (NBODY)"
            yunits = " (NBODY)"

        elif cluster.units == "kpckms" or cluster.units=='kpcgyr':
            rprof /=1000.0
            if projected:
                pprof *= (1000.0 ** 2.0)
                rho_c *= (1000.0 ** 2.0)
            else:
                pprof *= (1000.0 ** 3.0)
                rho_c *= (1000.0 ** 3.0)

            xunits = " (kpc)"
            if projected:
                yunits = " Msun/kpc^2"
            else:
                yunits = " Msun/kpc^3"

        elif cluster.units=='galpy':
            rprof /=(1000.0*ro)

            if projected:
                pprof /= dens_in_msolpc2
                rho_c /= dens_in_msolpc2
            else:
                pprof /= conversion.dens_in_msolpc3(ro=ro, vo=vo) 
                rho_c /= conversion.dens_in_msolpc3(ro=ro, vo=vo)               

            xunits = "(GALPY)"
            yunits = "(GALPY)"

        elif cluster.units == "pckms" or cluster.units == "pcmyr":
            xunits = " (pc)"
            if projected:
                yunits = " Msun/pc^2"
            else:
                yunits = " Msun/pc^3"

        else:
            xunits = ""
            yunits = ""

        x, y, n = rprof, pprof, nprof
        _lplot(
            x,
            y,
            xlabel=r"$R %s$" % (xunits),
            ylabel=r"$\rho %s$" % (yunits),
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )
        _lplot(x, np.ones(len(x)) * rho_c, "--", overplot=True)
        _lplot(np.ones(len(y)) * rc, y, "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    return rc

def rtidal(
    cluster,
    pot=None,
    rtiterate=0,
    rtconverge=0.9,
    indx=None,
    rgc=None,
    zgc=None,
    from_centre=False,
    plot=False,
    verbose=False,
    **kwargs,
):
    """Calculate tidal radius of the cluster
    - The calculation uses Galpy (Bovy 2015_, which takes the formalism of Bertin & Varri 2008 to calculate the tidal radius
    -- Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005
    -- Bovy J., 2015, ApJS, 216, 29
    - riterate = 0 corresponds to a single calculation of the tidal radius based on the cluster's mass (np.sum(cluster.m))
    -- Additional iterations take the mass within the previous iteration's calculation of the tidal radius and calculates the tidal
       radius again using the new mass until the change is less than 90%
    - for cases where the cluster's orbital parameters are not set, it is possible to manually set rgc which is assumed to be in kpc.

    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate tidal radius (default: None)
    rtiterate : int
        how many times to iterate on the calculation of r_t (default: 0)
    rtconverge : float
        criteria for tidal radius convergence within iterations (default 0.9)
    indx : bool
        subset of stars to use when calculate the tidal radius (default: None)
    rgc : float
        Manually set galactocentric distance in kpc at which the tidal radius is to be evaluated (default: None)
    zgc : float
        For non-spherically symmetric potentials, manually set distance in kpc above disk at which the tidal radius is to be evaluated. When set, rgc becomes radius in cylindrical coordinates (default: None)
    ro : float
        GALPY radius scaling parameter
    vo : float
        GALPY velocity scaling parameter
    from_centre : bool
        calculate tidal radius based on location of cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    plot : bool
        plot the x and y coordinates of stars and mark the tidal radius of the cluster (default: False)
    verbose : bool
        Print information about iterative calculation of rt

    Returns
    -------
    rt : float
        tidal radius

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2019 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre(sortstars=False)

    cluster.to_galpy()

    if rgc != None:
        R = rgc / ro
    else:
        if from_centre:
            R = np.sqrt((cluster.xgc+cluster.xc) ** 2.0 + (cluster.ygc+cluster.yc) ** 2.0)
        else:
            R = np.sqrt(cluster.xgc ** 2.0 + cluster.ygc ** 2.0)

    if zgc !=None:
        z = zgc/ ro
    else:
        if from_centre:
            z = cluster.zgc+cluster.zc
        else:
            z = cluster.zgc

    if indx is None:
        indx=np.ones(len(cluster.m),dtype=bool)

    # Calculate rtide
    rt = rtide(pot, R, z, M=np.sum(cluster.m[indx]),use_physical=False)
    nit = 0
    for i in range(0, rtiterate):
        msum = 0.0

        indx *= cluster.r < rt
        msum = np.sum(cluster.m[indx])

        rtnew = rtide(pot, R, z, M=msum,use_physical=False)

        if rtnew == 0.:
            print('RT DID NOT CONVERGE')
            rt=0.
            break

        if verbose:
            print(rt, rtnew, rt/rtnew, msum / np.sum(cluster.m))

        if rtnew / rt >= rtconverge:
            break
        rt = rtnew
        nit += 1

    if verbose:
        print(
            "FINAL RT: ",
            rt * ro * 1000.0,
            "pc after",
            nit,
            " of ",
            rtiterate,
            " iterations",
        )

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    rt=_convert_length(rt,'galpy',cluster)


    if plot:

        if cluster.units == "nbody":
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms" or cluster.units=="pcmyr":
            xunits = " (pc)"
            yunits = " (pc)"

        elif cluster.units == "kpckms" or cluster.units=="kpcgyr":
            xunits = " (kpc)"
            yunits = " (kpc)"
        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"
        else:
            xunits = ""
            yunits = ""


        _plot(
            cluster.x,
            cluster.y,
            xlabel=r"$x %s$" % xunits,
            ylabel=r"$y %s$" % yunits,
            title="Time = %f" % cluster.tphys,
            log=False,
            filename=None,
            **kwargs,
        )

        x=np.linspace(-rt,rt,100)
        y=np.sqrt(rt**2.-x**2.)

        x=np.append(x,x)
        y=np.append(y,-y)

        if cluster.origin0=='galaxy':
            if from_centre:
                x+=(cluster.xgc+cluster.xc)
                y+=(cluster.ygc+cluster.yc)
            else:
                x+=cluster.xgc
                y+=cluster.ygc
        elif cluster.origin0=='cluster' and from_centre:
            x+=cluster.xc
            y+=cluster.yc

        _lplot(x,y,overplot=True,linestyle='--',color='k')


    return rt


def rlimiting(
    cluster,
    pot=None,
    rgc=None,
    zgc=None,
    nrad=20,
    projected=False,
    plot=False,
    from_centre=False,
    **kwargs
):
    """Calculate limiting radius of the cluster
       
    - The limiting radius is defined to be where the cluster's density reaches the local background density of the host galaxy
    - for cases where the cluster's orbital parameters are not set, it is possible to manually set rgc which is assumed to be in kpc.

    Parameters
    ----------

    cluster : class
        StarCluster
    pot : class 
        GALPY potential used to calculate actions
    rgc : 
        Manually set galactocentric distance in kpc at which the tidal radius is to be evaluated (default: None)
    zgc : float
        For non-spherically symmetric potentials, manually set distance in kpc above disk at which the tidal radius is to be evaluated. When set, rgc becomes radius in cylindrical coordinates (default: None)
    nrad : int
        number of radial bins used to calculate density profile (Default: 20)
    projected : bool
        use projected values (default: False)
    plot : bool
        plot the density profile and mark the limiting radius of the cluster (default: False)
    from_centre : bool
        calculate tidal radius based on location of cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    verbose : bool
    Returns
    -------
        rl : float
            limiting radius

    Other Parameters
    ----------------
    kwargs : str
        key words for plotting

    History
    -------
    2019 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    mo=conversion.mass_in_msol(ro=ro,vo=vo)
    dens_in_msolpc2=(mo/ro**2.)/(1000.0**2.)

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre(sortstars=False)

    cluster.to_galpy()

    if rgc != None:
        R = rgc / ro
    else:
        if from_centre:
            R = np.sqrt((cluster.xgc+cluster.xc) ** 2.0 + (cluster.ygc+cluster.yc) ** 2.0)
        else:
            R = np.sqrt(cluster.xgc ** 2.0 + cluster.ygc ** 2.0)

    if zgc !=None:
        z = zgc/ ro
    else:
        if from_centre:
            z = cluster.zgc+cluster.zc
        else:
            z = cluster.zgc

    # Calculate local density:
    rho_local = potential.evaluateDensities(
        pot, R, z, ro=ro, vo=vo, use_physical=False
     )


    rprof, pprof, nprof = _rho_prof(cluster, nrad=nrad, projected=projected)

    #Approximate projected local density across entire area of cluster
    if projected: rho_local*=(4.*rprof[-1]/3)

    if pprof[-1] > rho_local:
        rl = rprof[-1]
    elif pprof[0] < rho_local:
        rl = 0.0
    else:
        indx = np.argwhere(pprof < rho_local)[0][0]
        r1 = (rprof[indx - 1], pprof[indx - 1])
        r2 = (rprof[indx], pprof[indx])

        rl = interpolate(r1, r2, y=rho_local)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    rl=_convert_length(rl,'galpy',cluster)

    if plot:

        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            rprof *= ro * 1000.0 / cluster.rbar

            if projected:
                pprof *= (
                    dens_in_msolpc2
                    * (cluster.rbar ** 2.0)
                    / cluster.zmbar
                )

                rho_local *= (
                    dens_in_msolpc2
                    * (cluster.rbar ** 2.0)
                    / cluster.zmbar
                )

            else:
                pprof *= (
                    conversion.dens_in_msolpc3(ro=ro, vo=vo)
                    * (cluster.rbar ** 3.0)
                    / cluster.zmbar
                )

                rho_local *= (
                    conversion.dens_in_msolpc3(ro=ro, vo=vo)
                    * (cluster.rbar ** 3.0)
                    / cluster.zmbar
                )

            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms" or cluster.units == "pcmyr":
            rprof *= ro * 1000.0

            if projected:
                pprof *= dens_in_msolpc2
                rho_local *= dens_in_msolpc2
            else:
                pprof *= conversion.dens_in_msolpc3(ro=ro, vo=vo)
                rho_local *= conversion.dens_in_msolpc3(ro=ro, vo=vo)

            xunits = " (pc)"
            if projected:
                yunits = " Msun/pc^2"
            else:
                yunits = " Msun/pc^3"
        elif cluster.units == "kpckms" or cluster.units == "kpcgyr":
            rprof *= ro
            if projected:
                pprof *= dens_in_msolpc2 * (1000.0 ** 2.0)
                rho_local *= dens_in_msolpc2 * (1000.0 ** 2.0)
            else:
                pprof *= conversion.dens_in_msolpc3(ro=ro, vo=vo) * (1000.0 ** 3.0)
                rho_local *= conversion.dens_in_msolpc3(ro=ro, vo=vo) * (1000.0 ** 3.0)

            xunits = " (kpc)"
            if projected:
                yunits = " Msun/kpc^2"
            else:
                yunits = " Msun/kpc^3"
        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"

        else:
            xunits = ""
            yunits = ""

        x, y, n = rprof, pprof, nprof
        _lplot(
            x,
            y,
            xlabel=r"$R %s$" % (xunits),
            ylabel=r"$\rho %s$" % (yunits),
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )
        _lplot(x, np.ones(len(x)) * rho_local, "--", overplot=True)
        _lplot(np.ones(len(y)) * rl, y, "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    return rl

def _rho_prof(
    cluster,
    mmin=None,
    mmax=None,
    rmin=None,
    rmax=None,
    nrad=20,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=15,
    npop=None,
    indx=None,
    bins=None,
    projected=False,
    normalize=False,
    plot=False,
    **kwargs
):
    """Measure the density profile of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    mmin/mmax : float
        minimum and maximum stellar mass
    rmin/rmax : float
        minimum and maximum stellar radii
    nrad : int
        number of radial bins
    vmin/vmax : float 
        minimum and maximum stellar velocity
    emin/emax : float
        minimum and maximum stellar energy
    kwmin/kwmax : float
        minimum and maximum stellar type (kw)
    npop : int
        population number
    indx : float
        user defined boolean array from which to extract the subset
    bins : float
        User defined bins in the form of (rlower,rmean,rupper) (default: None)
    projected : bool
        use projected values and constraints (default:False)
    normalize : bool
        normalize radial bins by cluster's half-mass radius (default: False)
    plot : bool 
        plot the density profile (default: False)

    Returns
    -------
    rprof : float
        radius bins
    pprof : float
        mass density in each bin
    nprof : float
        number of stars in each bin

    Other Parameters
    ----------------
    kwrags : str
        key word arguments for plotting

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre(sortstars=normalize)
    elif normalize:
        cluster.sortstars()

    rprof = np.array([])
    pprof = np.array([])
    nprof = np.array([])

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    """
    if rmin == None:
        rmin = np.min(r)
    if rmax == None:
        rmax = np.max(r)
    if vmin == None:
        vmin = np.min(v)
    if vmax == None:
        vmax = np.max(v)
    if mmin == None:
        mmin = np.min(cluster.m)
    if mmax == None:
        mmax = np.max(cluster.m)

    if indx is None:
        indx = cluster.id > -1

    # Build subcluster containing only stars in the full radial and mass range:
    indx *= (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m <= mmax)
        * (v >= vmin)
        * (v <= vmax)
    )

    if len(cluster.kw)>0:
        indx*=(cluster.kw >= kwmin) * (cluster.kw <= kwmax)

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax
    """
    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))
    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])
        rprof = np.append(rprof, r_mean[i])
        if projected:
            vol = np.pi * (r_upper[i] ** 2 - r_lower[i] ** 2.0)
        else:
            vol = (4.0 / 3.0) * np.pi * (r_upper[i] ** 3 - r_lower[i] ** 3.0)

        pprof = np.append(pprof, np.sum(cluster.m[rindx] / vol))
        nprof = np.append(nprof, np.sum(rindx))

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms" or cluster.units == "pcmyr":
            xunits = " (pc)"
            if projected:
                yunits = " Msun/pc^2"
            else:
                yunits = " Msun/pc^3"
        elif cluster.units == "kpckms" or cluster.units == "kpcgyr":
            xunits = " (kpc)"
            if projected:
                yunits = " Msun/kpc^2"
            else:
                yunits = " Msun/kpc^3"
        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"

        else:
            xunits = ""
            yunits = ""

        if projected:
            xlabel=r"$R \ %s$" % xunits
            ylabel=r"$\Sigma \ %s$" % yunits
        else:
            xlabel=r"$r \ %s$" % xunits
            ylabel=r"$\rho \ %s$" % yunits


        x, y, n = rprof, pprof, nprof

        if normalize:
            x/=cluster.rm

        _lplot(
            x,
            y,
            xlabel=xlabel,
            ylabel=ylabel,
            title="Time = %f" % cluster.tphys,
            log=kwargs.pop('log',True),
            overplot=overplot,
            filename=filename,
            **kwargs,
        )

        if filename != None:
            plt.savefig(filename)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if normalize:
        rprof/=cluster.rm

    return rprof, pprof, nprof
