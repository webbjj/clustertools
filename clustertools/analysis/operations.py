""" Perform an operation on a cluster and return a new cluster

"""

__author__ = "Jeremy J Webb"

__all__ = [
    "to_pckms",
    "to_kpckms",
    "to_nbody",
    "to_radec",
    "to_galpy",
    "to_units",
    "to_centre",
    "to_cluster",
    "to_galaxy",
    "to_sky",
    "to_origin",
    "save_cluster",
    "return_cluster",
    "reset_nbody_scale",
    "add_rotation",
    "virialize",
]

import numpy as np
from galpy.util import bovy_coords,bovy_conversion,_rotate_to_arbitrary_vector
from copy import copy

def to_pckms(cluster, do_key_params=False):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and km/s

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units == "galpy" or cluster.units == "radec":
        cluster.to_kpckms()

    if cluster.units == "nbody":
        cluster.m *= cluster.zmbar
        cluster.x *= cluster.rbar
        cluster.y *= cluster.rbar
        cluster.z *= cluster.rbar
        cluster.vx *= cluster.vstar
        cluster.vy *= cluster.vstar
        cluster.vz *= cluster.vstar

        cluster.xc *= cluster.rbar
        cluster.yc *= cluster.rbar
        cluster.zc *= cluster.rbar
        cluster.vxc *= cluster.vstar
        cluster.vyc *= cluster.vstar
        cluster.vzc *= cluster.vstar

        cluster.xgc *= cluster.rbar
        cluster.ygc *= cluster.rbar
        cluster.zgc *= cluster.rbar
        cluster.vxgc *= cluster.vstar
        cluster.vygc *= cluster.vstar
        cluster.vzgc *= cluster.vstar

        cluster.units = "pckms"

    elif cluster.units == "kpckms":
        cluster.x *= 1000.0
        cluster.y *= 1000.0
        cluster.z *= 1000.0

        cluster.xgc *= 1000.0
        cluster.ygc *= 1000.0
        cluster.zgc *= 1000.0

        cluster.xc *= 1000.0
        cluster.yc *= 1000.0
        cluster.zc *= 1000.0

        cluster.units = "pckms"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_kpckms(cluster, do_key_params=False, ro=8.0, vo=220.0):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and km/s

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    ro : float
        galpy radius scaling parameter
    vo : float
        galpy velocity scaling parameter

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units == "radec":
        from_radec(cluster)

    if cluster.units == "galpy":
        cluster.m *= bovy_conversion.mass_in_msol(ro=ro, vo=vo)
        cluster.x *= ro
        cluster.y *= ro
        cluster.z *= ro
        cluster.vx *= vo
        cluster.vy *= vo
        cluster.vz *= vo

        cluster.xc *= ro
        cluster.yc *= ro
        cluster.zc *= ro
        cluster.vxc *= vo
        cluster.vyc *= vo
        cluster.vzc *= vo

        cluster.xgc *= ro
        cluster.ygc *= ro
        cluster.zgc *= ro
        cluster.vxgc *= vo
        cluster.vygc *= vo
        cluster.vzgc *= vo

        cluster.units = "kpckms"

    elif cluster.units == "nbody":
        cluster.m *= cluster.zmbar
        cluster.x *= cluster.rbar / 1000.0
        cluster.y *= cluster.rbar / 1000.0
        cluster.z *= cluster.rbar / 1000.0
        cluster.vx *= cluster.vstar
        cluster.vy *= cluster.vstar
        cluster.vz *= cluster.vstar

        cluster.xc *= cluster.rbar / 1000.0
        cluster.yc *= cluster.rbar / 1000.0
        cluster.zc *= cluster.rbar / 1000.0
        cluster.vxc *= cluster.vstar
        cluster.vyc *= cluster.vstar
        cluster.vzc *= cluster.vstar

        cluster.xgc *= cluster.rbar / 1000.0
        cluster.ygc *= cluster.rbar / 1000.0
        cluster.zgc *= cluster.rbar / 1000.0
        cluster.vxgc *= cluster.vstar
        cluster.vygc *= cluster.vstar
        cluster.vzgc *= cluster.vstar

        cluster.units = "kpckms"

    elif cluster.units == "pckms":
        cluster.x /= 1000.0
        cluster.y /= 1000.0
        cluster.z /= 1000.0

        cluster.xgc /= 1000.0
        cluster.ygc /= 1000.0
        cluster.zgc /= 1000.0

        cluster.xc /= 1000.0
        cluster.yc /= 1000.0
        cluster.zc /= 1000.0

        cluster.units = "kpckms"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_nbody(cluster, do_key_params=False, ro=8.0, vo=220.0):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Nbody units
       
    - requires that cluster.zmbar, cluster.rbar, cluster.vstar are set (defaults are 1)

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    ro : float
        galpy radius scaling parameter
    vo : float
        galpy velocity scaling parameter

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units != "pckms":
        cluster.to_pckms(do_key_params=False)

    if cluster.units == "pckms":
        cluster.m /= cluster.zmbar
        cluster.x /= cluster.rbar
        cluster.y /= cluster.rbar
        cluster.z /= cluster.rbar
        cluster.vx /= cluster.vstar
        cluster.vy /= cluster.vstar
        cluster.vz /= cluster.vstar

        cluster.xc /= cluster.rbar
        cluster.yc /= cluster.rbar
        cluster.zc /= cluster.rbar
        cluster.vxc /= cluster.vstar
        cluster.vyc /= cluster.vstar
        cluster.vzc /= cluster.vstar

        cluster.xgc /= cluster.rbar
        cluster.ygc /= cluster.rbar
        cluster.zgc /= cluster.rbar
        cluster.vxgc /= cluster.vstar
        cluster.vygc /= cluster.vstar
        cluster.vzgc /= cluster.vstar

        cluster.units = "nbody"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_radec(cluster, do_key_params=False, ro=8.0, vo=220.0):
    """Convert to on-sky position, proper motion, and radial velocity of cluster
    
    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    ro : float
        galpy radius scaling parameter
    vo : float
        galpy velocity scaling parameter

    Returns
    -------
    None

    History:
    -------
   2018 - Written - Webb (UofT)
    """
    if len(cluster.ra)==len(cluster.x):
        cluster.x = copy(cluster.ra)
        cluster.y = copy(cluster.dec)
        cluster.z = copy(cluster.dist)
        cluster.vx = copy(cluster.pmra)
        cluster.vy = copy(cluster.pmdec)
        cluster.vz = copy(cluster.vlos)

        cluster.units = "radec"
        cluster.origin = "sky"

    else:

        units0, origin0 = cluster.units, cluster.origin

        cluster.to_galaxy()
        cluster.to_kpckms()

        x0, y0, z0 = bovy_coords.galcenrect_to_XYZ(
            cluster.x, cluster.y, cluster.z, Xsun=8.0, Zsun=0.025
        ).T

        cluster.dist = np.sqrt(x0 ** 2.0 + y0 ** 2.0 + z0 ** 2.0)

        vx0, vy0, vz0 = bovy_coords.galcenrect_to_vxvyvz(
            cluster.vx,
            cluster.vy,
            cluster.vz,
            Xsun=8.0,
            Zsun=0.025,
            vsun=[-11.1, 244.0, 7.25],
        ).T

        cluster.vlos = (vx0 * x0 + vy0 * y0 + vz0 * z0) / np.sqrt(
            x0 ** 2.0 + y0 ** 2.0 + z0 ** 2.0
        )

        l0, b0, cluster.dist = bovy_coords.XYZ_to_lbd(x0, y0, z0, degree=True).T
        cluster.ra, cluster.dec = bovy_coords.lb_to_radec(l0, b0, degree=True).T

        vr0, pmll0, pmbb0 = bovy_coords.vxvyvz_to_vrpmllpmbb(
            vx0, vy0, vz0, l0, b0, cluster.dist, degree=True
        ).T
        cluster.pmra, cluster.pmdec = bovy_coords.pmllpmbb_to_pmrapmdec(
            pmll0, pmbb0, l0, b0, degree=True
        ).T

        x0, y0, z0 = bovy_coords.galcenrect_to_XYZ(
            cluster.xgc, cluster.ygc, cluster.zgc, Xsun=8.0, Zsun=0.025
        )
        vx0, vy0, vz0 = bovy_coords.galcenrect_to_vxvyvz(
            cluster.vxgc,
            cluster.vygc,
            cluster.vzgc,
            Xsun=8.0,
            Zsun=0.025,
            vsun=[-11.1, 244.0, 7.25],
        )

        cluster.vlos_gc = (vx0 * x0 + vy0 * y0 + vz0 * z0) / np.sqrt(
            x0 ** 2.0 + y0 ** 2.0 + z0 ** 2.0
        )

        l0, b0, cluster.dist_gc = bovy_coords.XYZ_to_lbd(x0, y0, z0, degree=True)
        cluster.ra_gc, cluster.dec_gc = bovy_coords.lb_to_radec(l0, b0, degree=True)

        vr0, pmll0, pmbb0 = bovy_coords.vxvyvz_to_vrpmllpmbb(
            vx0, vy0, vz0, l0, b0, cluster.dist_gc, degree=True
        )
        cluster.pmra_gc, cluster.pmdec_gc = bovy_coords.pmllpmbb_to_pmrapmdec(
            pmll0, pmbb0, l0, b0, degree=True
        )

        cluster.x = copy(cluster.ra)
        cluster.y = copy(cluster.dec)
        cluster.z = copy(cluster.dist)
        cluster.vx = copy(cluster.pmra)
        cluster.vy = copy(cluster.pmdec)
        cluster.vz = copy(cluster.vlos)

        cluster.xgc = copy(cluster.ra_gc)
        cluster.ygc = copy(cluster.dec_gc)
        cluster.zgc = copy(cluster.dist_gc)
        cluster.vxgc = copy(cluster.pmra_gc)
        cluster.vygc = copy(cluster.pmdec_gc)
        cluster.vzgc = copy(cluster.vlos_gc)

        cluster.units = "radec"
        cluster.origin = "sky"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def from_radec(cluster, do_key_params=False):
    """Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units == "radec" and cluster.origin == "sky":

        origin0 = cluster.origin

        l, b = bovy_coords.radec_to_lb(cluster.ra, cluster.dec, degree=True).T
        x0, y0, z0 = bovy_coords.lbd_to_XYZ(l, b, cluster.dist, degree=True).T
        cluster.x, cluster.y, cluster.z = bovy_coords.XYZ_to_galcenrect(
            x0, y0, z0, Xsun=8.0, Zsun=0.025
        ).T

        pml, pmb = bovy_coords.pmrapmdec_to_pmllpmbb(
            cluster.pmra, cluster.pmdec, cluster.ra, cluster.dec, degree=True
        ).T
        vx0, vy0, vz0 = bovy_coords.vrpmllpmbb_to_vxvyvz(
            cluster.vlos, pml, pmb, l, b, cluster.dist, degree=True
        ).T
        cluster.vx, cluster.vy, cluster.vz = bovy_coords.vxvyvz_to_galcenrect(
            vx0,
            vy0,
            vz0,
            vsun=[0.0, 220.0, 0.0],
            Xsun=8.0,
            Zsun=0.025,
            _extra_rot=True,
        ).T

        l_gc, b_gc = bovy_coords.radec_to_lb(cluster.ra_gc, cluster.dec_gc, degree=True)
        x0_gc, y0_gc, z0_gc = bovy_coords.lbd_to_XYZ(
            l_gc, b_gc, cluster.dist_gc, degree=True
        )
        cluster.xgc, cluster.ygc, cluster.zgc = bovy_coords.XYZ_to_galcenrect(
            x0_gc, y0_gc, z0_gc, Xsun=8.0, Zsun=0.025
        )

        pml_gc, pmb_gc = bovy_coords.pmrapmdec_to_pmllpmbb(
            cluster.pmra_gc, cluster.pmdec_gc, cluster.ra_gc, cluster.dec_gc, degree=True
        )
        vx0_gc, vy0_gc, vz0_gc = bovy_coords.vrpmllpmbb_to_vxvyvz(
            cluster.vlos_gc, pml_gc, pmb_gc, l_gc, b_gc, cluster.dist_gc, degree=True
        )
        cluster.vx_gc, cluster.vy_gc, cluster.vz_gc = bovy_coords.vxvyvz_to_galcenrect(
            vx0_gc,
            vy0_gc,
            vz0_gc,
            vsun=[0.0, 220.0, 0.0],
            Xsun=8.0,
            Zsun=0.025,
            _extra_rot=True,
        )

        cluster.origin = "galaxy"
        cluster.units = "kpckms"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_galpy(cluster, do_key_params=False, ro=8.0, vo=220.0):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to galpy units
    
    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    ro : float
        galpy radius scaling parameter
    vo : float
        galpy velocity scaling parameter

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units != "kpckms" and cluster.units != "galpy":
        cluster.to_kpckms(do_key_params=False)

    if cluster.units == "kpckms":
        cluster.m = cluster.m / bovy_conversion.mass_in_msol(ro=ro, vo=vo)
        cluster.x /= ro
        cluster.y /= ro
        cluster.z /= ro
        cluster.vx /= vo
        cluster.vy /= vo
        cluster.vz /= vo

        cluster.xc /= ro
        cluster.yc /= ro
        cluster.zc /= ro
        cluster.vxc /= vo
        cluster.vyc /= vo
        cluster.vzc /= vo

        cluster.xgc /= ro
        cluster.ygc /= ro
        cluster.zgc /= ro
        cluster.vxgc /= vo
        cluster.vygc /= vo
        cluster.vzgc /= vo

        cluster.units = "galpy"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_units(cluster, units, do_key_params=False, ro=8.0, vo=220.0):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to user defined units

    Parameters
    ----------
    cluster : class
        StarCluster
    units : str
        units to be converted to (nbody,galpy,pckms,kpckms,radec)
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    ro : float
        galpy radius scaling parameter
    vo : float
        galpy velocity scaling parameter

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if units == "nbody":
        cluster.to_nbody(do_key_params=do_key_params)
    elif units == "galpy":
        cluster.to_galpy(do_key_params=do_key_params, ro=ro, vo=vo)
    elif units == "pckms":
        cluster.to_pckms(do_key_params=do_key_params)
    elif units == "kpckms":
        cluster.to_kpckms(do_key_params=do_key_params)
    elif units == "radec":
        origin0 = cluster.origin
        cluster.to_radec(do_key_params=do_key_params)
        cluster.to_origin(origin0, do_key_params=do_key_params)

def to_centre(cluster, do_key_params=False, centre_method=None):
    """Shift coordinates such that origin is the centre of mass

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.origin != "centre":

        if cluster.origin != "cluster":
            cluster.to_cluster(do_key_params=False, centre_method=centre_method)

        cluster.x -= cluster.xc
        cluster.y -= cluster.yc
        cluster.z -= cluster.zc
        cluster.vx -= cluster.vxc
        cluster.vy -= cluster.vyc
        cluster.vz -= cluster.vzc

        cluster.origin = "centre"

        cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_cluster(cluster, do_key_params=False, centre_method=None):
    """Shift coordinates to clustercentric reference frame

    - When units='radec' and origin='sky', the cluster will be shifted to clustercentric coordinates using either:
    --centre_method=None: angular distance between each star's RA/DEC and the RA/DEC of the cluster's centre with proper motions directly subtracted off
    --centre_method='orthographic' , positions and velocities changed to orthnormal coordinates (Helmi et al. 2018)
    -- centre_method='VandeVen' , positions and velocities changed to clustercentric coordinates using method outlined by Van de Ven et al. 2005

    - Note the the line of site positions and velocities will just have the galactocentric coordinates of the cluster
    subtracted off. Be sure to set projected=True when making any calculations to use only x and y coordinates

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates. (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if centre_method is not None:
        cluster.centre_method = centre_method

    if cluster.origin != "cluster":
        if cluster.units == "radec" and cluster.origin == "sky":
            ra = np.radians(cluster.x)
            dec = np.radians(cluster.y)
            pmra = np.radians(cluster.vx / (1000.0 * 3600.0))
            pmdec = np.radians(cluster.vy / (1000.0 * 3600.0))
            ra_gc = np.radians(cluster.xgc)
            dec_gc = np.radians(cluster.ygc)

            if cluster.centre_method == "orthographic":
                cluster.x = np.cos(dec) * np.sin(ra - ra_gc)
                cluster.y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
                    dec_gc
                ) * np.cos(ra - ra_gc)
                cluster.z = np.zeros(len(cluster.x))

                cluster.vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
                    ra - ra_gc
                )
                cluster.vy = pmra * np.sin(dec_gc) * np.sin(
                    ra - ra_gc
                ) + cluster.pmdec * (
                    np.cos(dec) * np.cos(dec_gc)
                    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
                )

                cluster.z -= cluster.zgc

            else:
                if cluster.centre_method == "VandeVen":
                    cluster.x = (
                        (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
                    )
                    cluster.y = (
                        (10800.0 / np.pi)
                        * (
                            np.sin(dec) * np.cos(dec_gc)
                            - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
                        )
                        / 60.0
                    )
                else:
                    cluster.x = (cluster.x - cluster.xgc) * np.cos(np.radians(cluster.ygc))
                    cluster.y = cluster.y - cluster.ygc

                cluster.z -= cluster.zgc

                cluster.vx -= -cluster.vxgc
                cluster.vy -= cluster.vygc
                cluster.vz -= cluster.vzgc

        elif cluster.origin == "centre":
            cluster.x += cluster.xc
            cluster.y += cluster.yc
            cluster.z += cluster.zc
            cluster.vx += cluster.vxc
            cluster.vy += cluster.vyc
            cluster.vz += cluster.vzc
        elif cluster.origin == "galaxy":
            cluster.x -= cluster.xgc
            cluster.y -= cluster.ygc
            cluster.z -= cluster.zgc
            cluster.vx -= cluster.vxgc
            cluster.vy -= cluster.vygc
            cluster.vz -= cluster.vzgc

        cluster.rv3d()

        cluster.origin = "cluster"
        if do_key_params:
            cluster.key_params()

def to_galaxy(cluster, do_key_params=False):
    """Shift coordinates to galactocentric reference frame

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units == "radec" and cluster.origin == "sky":
        cluster.from_radec(do_key_params=False)

    elif cluster.origin != "galaxy":
        if cluster.origin == "centre":
            cluster.to_cluster(do_key_params=False)

        cluster.x += cluster.xgc
        cluster.y += cluster.ygc
        cluster.z += cluster.zgc
        cluster.vx += cluster.vxgc
        cluster.vy += cluster.vygc
        cluster.vz += cluster.vzgc

        cluster.origin = "galaxy"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_sky(cluster, do_key_params=False):
    """Calculate on-sky position, proper motion, and radial velocity of cluster
        
    - Also changes units to radec

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)
    """
    cluster.to_radec(do_key_params=do_key_params)

def from_sky(cluster, do_key_params=False):
    """Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster

    - Also changes units to kpckms

    Parameters
    ----------
    cluster : class
        StarCluster
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    cluster.from_radec(do_key_params=do_key_params)

def to_origin(cluster, origin, do_key_params=False):
    """Shift cluster to origin as defined by keyword

    Parameters
    ----------
    cluster : class
        StarCluster
    origin : str
        origin of coordinate system to be shifted to (centre, cluster, galaxy, sky)
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if origin == "centre":
        cluster.to_centre(do_key_params=do_key_params)
    elif origin == "cluster":
        cluster.to_cluster(do_key_params=do_key_params)
    elif origin == "galaxy":
        cluster.to_galaxy(do_key_params=do_key_params)
    elif origin == "sky":
        cluster.to_sky(do_key_params=do_key_params)


def save_cluster(cluster):
    """Save cluster's units and origin

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    cluster.units, cluster.origin : str
        units and origin of StarCluster

    History:
    -------
    2018 - Written - Webb (UofT)
    """
    return cluster.units, cluster.origin


def return_cluster(cluster, units0, origin0, do_key_params=False):
    """ return cluster to a specific combination of units and origin

    Parameters
    ----------

    cluster : class
        StarCluster
    units0 : str
        units that StarCluster will be changed to
   origin0 : str
        origin that StarCluster will be changed to
    do_key_params : bool
        call key_params to calculate key parameters after unit change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)
    """
    if cluster.units != units0:
        cluster.to_units(units0)
    if cluster.origin != origin0:
        cluster.to_origin(origin0, do_key_params=do_key_params)

def reset_nbody_scale(cluster, mass=True, radii=True, rvirial=True, projected=False,**kwargs):
    """ Assign new conversions for real mass, size, and velocity to Nbody units
    
    - kwargs are passed to the virial_radius function. See the virial_radius documenation in functions.py

    Parameters
    ----------

    cluster : class
        StarCluster instance
    mass : bool
        find new mass conversion (default: True)
    radii : bool
        find new radius conversion (default: True)
    rvirial : bool (default: True)
        use virial radius to set conversion rate for radii as opposed to the approximation that rbar=4/3 rm
    projected : bool
        use projected values to calculate radius conversion (default: False)

    Returns
    -------

    zmbar : float
        mass conversion
    rbar : float
        radius conversion
    vstar : float
        velocity conversion
    tstar : float
        time conversion

    History:

    2018 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)

    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre(do_key_params=True)
    else:
        cluster._order_check()

    cluster.to_pckms(do_key_params=True)

    if mass:
        zmbar = cluster.mtot
    else:
        zmbar = cluster.zmbar

    if radii:
        if rvirial:

            method=kwargs.get("method", 'inverse_distance')
            full=kwargs.get("full", True)
            H = kwargs.get("H", 70.0)
            Om = kwargs.get("Om", 0.3)
            overdens = kwargs.get("overdens", 200.0)
            nrad = kwargs.get("nrad", 20.0)
            plot = kwargs.get("plot", False)

            cluster.virial_radius(method=method,full=full,
                H=H,Om=Om,overdens=overdens,nrad=nrad,projected=projected,
                plot=plot,**kwargs)

            rbar=cluster.rv
        else:
            rbar = 4.0 * cluster.rm / 3.0
    else:
        rbar = cluster.rbar

    vstar = 0.06557 * np.sqrt(zmbar / rbar)
    tstar = rbar / vstar

    return_cluster(cluster, units0, origin0)

    return zmbar,rbar,vstar,tstar

def virialize(cluster, specific=True, full=True, projected=False):
    """ Adjust stellar velocities so cluster is in virial equilibrium

    Parameters
    ----------
    cluster : class
        StarCluster
    specific : bool
        find specific energies (default: True)
    full: bool
        do full array of stars at once with numba (default: True)
    projected : bool
        use projected values when calculating energies (default: False)

    Returns
    -------
    qv : float
        scaling factor used to adjust velocities

    History
    -------
    2019 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)
    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre()

    try:
        qv = np.sqrt(np.abs(0.5 / cluster.qvir))
    except:
        print("NEED TO CALCULATE ENERGIES FIRST")
        cluster.energies(specific=specific, full=full, projected=projected)
        qv = np.sqrt(np.abs(0.5 / cluster.qvir))

    return_cluster(cluster, units0, origin0)

    return qv

def add_rotation(cluster, qrot):
    """Add a degree of rotation to an already generated StarCluster

    Parameters
    ----------
    cluster : class
        StarCluster
    qrot : float
        fraction of stars with v_phi < 0 that are switched to vphi > 0

    Returns
    -------
    x,y,z,vx,vy,vz : float
        stellar positions and velocities (now with rotation)

    History
    -------
    2019 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)

    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre()

    r, theta, z = bovy_coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
    vr, vtheta, vz = bovy_coords.rect_to_cyl_vec(
        cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
    )

    indx = vtheta < 0.0
    rindx = np.random.rand(cluster.ntot) < qrot

    vtheta[indx * rindx] = np.sqrt(vtheta[indx * rindx] * vtheta[indx * rindx])
    x,y,z = bovy_coords.cyl_to_rect(r, theta, z)
    vx,vy,vz = bovy_coords.cyl_to_rect_vec(
        vr, vtheta, vz, theta
    )

    return_cluster(cluster, units0, origin0)

    return x,y,z,vx,vy,vz
