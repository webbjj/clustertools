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
from galpy.util import _rotate_to_arbitrary_vector
try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion
    
from copy import copy
try:
    import numba
    full_default=True
except:
    full_default=False
    pass

def to_pckms(cluster):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and km/s

    Parameters
    ----------
    cluster : class
        StarCluster

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
        cluster.tphys *= cluster.tbar
        cluster.m *= cluster.zmbar
        cluster.x *= cluster.rbar
        cluster.y *= cluster.rbar
        cluster.z *= cluster.rbar
        cluster.vx *= cluster.vbar
        cluster.vy *= cluster.vbar
        cluster.vz *= cluster.vbar

        cluster.xc *= cluster.rbar
        cluster.yc *= cluster.rbar
        cluster.zc *= cluster.rbar
        cluster.vxc *= cluster.vbar
        cluster.vyc *= cluster.vbar
        cluster.vzc *= cluster.vbar

        cluster.xgc *= cluster.rbar
        cluster.ygc *= cluster.rbar
        cluster.zgc *= cluster.rbar
        cluster.vxgc *= cluster.vbar
        cluster.vygc *= cluster.vbar
        cluster.vzgc *= cluster.vbar

        cluster.units = "pckms"

    elif cluster.units == "kpckms":

        cluster.tphys *= 1000.0

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

    cluster.analyze()

def to_kpckms(cluster, ro=8.0, vo=220.0):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and km/s

    Parameters
    ----------
    cluster : class
        StarCluster
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
        cluster.tphys *= conversion.time_in_Gyr(ro=ro,vo=vo)
        cluster.m *= conversion.mass_in_msol(ro=ro, vo=vo)
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
        cluster.tphys *= cluster.tbar * 1000.0
        cluster.m *= cluster.zmbar
        cluster.x *= cluster.rbar / 1000.0
        cluster.y *= cluster.rbar / 1000.0
        cluster.z *= cluster.rbar / 1000.0
        cluster.vx *= cluster.vbar
        cluster.vy *= cluster.vbar
        cluster.vz *= cluster.vbar

        cluster.xc *= cluster.rbar / 1000.0
        cluster.yc *= cluster.rbar / 1000.0
        cluster.zc *= cluster.rbar / 1000.0
        cluster.vxc *= cluster.vbar
        cluster.vyc *= cluster.vbar
        cluster.vzc *= cluster.vbar

        cluster.xgc *= cluster.rbar / 1000.0
        cluster.ygc *= cluster.rbar / 1000.0
        cluster.zgc *= cluster.rbar / 1000.0
        cluster.vxgc *= cluster.vbar
        cluster.vygc *= cluster.vbar
        cluster.vzgc *= cluster.vbar

        cluster.units = "kpckms"

    elif cluster.units == "pckms":

        cluster.tphys /= 1000.0

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

    cluster.analyze()

def to_nbody(cluster, ro=8.0, vo=220.0):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Nbody units
       
    - requires that cluster.zmbar, cluster.rbar, cluster.vbar are set (defaults are 1)

    Parameters
    ----------
    cluster : class
        StarCluster
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
        cluster.to_pckms()

    if cluster.units == "pckms":
        cluster.tphys /= cluster.tbar
        cluster.m /= cluster.zmbar
        cluster.x /= cluster.rbar
        cluster.y /= cluster.rbar
        cluster.z /= cluster.rbar
        cluster.vx /= cluster.vbar
        cluster.vy /= cluster.vbar
        cluster.vz /= cluster.vbar

        cluster.xc /= cluster.rbar
        cluster.yc /= cluster.rbar
        cluster.zc /= cluster.rbar
        cluster.vxc /= cluster.vbar
        cluster.vyc /= cluster.vbar
        cluster.vzc /= cluster.vbar

        cluster.xgc /= cluster.rbar
        cluster.ygc /= cluster.rbar
        cluster.zgc /= cluster.rbar
        cluster.vxgc /= cluster.vbar
        cluster.vygc /= cluster.vbar
        cluster.vzgc /= cluster.vbar

        cluster.units = "nbody"

    cluster.analyze()

def to_radec(cluster, sortstars=True, from_centre=False, ro=8.0, vo=220.0):
    """Convert to on-sky position, proper motion, and radial velocity of cluster
    
    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
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

    if cluster.units != "radec" and cluster.origin != "sky":

        if cluster.orbit is None: cluster.initialize_orbit(from_centre,ro=ro,vo=vo)
        if cluster.orbits is None: cluster.initialize_orbits(ro=ro,vo=vo)

        cluster.ra=cluster.orbits.ra()
        cluster.dec=cluster.orbits.dec()
        cluster.dist=cluster.orbits.dist()
        cluster.pmra=cluster.orbits.pmra()
        cluster.pmdec=cluster.orbits.pmdec()
        cluster.vlos=cluster.orbits.vlos()

        cluster.ra_gc=cluster.orbit.ra()
        cluster.dec_gc=cluster.orbit.dec()
        cluster.dist_gc=cluster.orbit.dist()
        cluster.pmra_gc=cluster.orbit.pmra()
        cluster.pmdec_gc=cluster.orbit.pmdec()
        cluster.vlos_gc=cluster.orbit.vlos()

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

        cluster.analyze(sortstars=sortstars)

    elif cluster.units=='radec' and (cluster.origin=='cluster' or cluster.origin=='centre'):
        cluster.x+=cluster.xgc
        cluster.y+=cluster.ygc
        cluster.z+=cluster.zgc
        cluster.vx+=cluster.vxgc
        cluster.vy+=cluster.vygc
        cluster.vz+=cluster.vzgc
        cluster.origin = "sky"
        cluster.analyze(sortstars=sortstars)

def from_radec(cluster, sortstars=True, from_centre=False, ro=8.0, vo=220.0):
    """Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
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
    if cluster.units == "radec" and cluster.origin == "sky":

        if cluster.orbit is None: cluster.initialize_orbit(from_centre,ro=ro,vo=vo)
        if cluster.orbits is None: cluster.initialize_orbits(ro=ro,vo=vo)

        cluster.x = cluster.orbits.x()
        cluster.y = cluster.orbits.y()
        cluster.z = cluster.orbits.z()
        cluster.vx = cluster.orbits.vx()
        cluster.vy = cluster.orbits.vy()
        cluster.vz = cluster.orbits.vz()

        cluster.xgc = cluster.orbit.x()
        cluster.ygc = cluster.orbit.y()
        cluster.zgc = cluster.orbit.z()
        cluster.vxgc = cluster.orbit.vx()
        cluster.vygc = cluster.orbit.vy()
        cluster.vzgc = cluster.orbit.vz()

        cluster.origin = "galaxy"
        cluster.units = "kpckms"

def to_galpy(cluster, ro=8.0, vo=220.0):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to galpy units
    
    Parameters
    ----------
    cluster : class
        StarCluster
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
        cluster.to_kpckms()

    if cluster.units == "kpckms":
        cluster.tphys /= conversion.time_in_Gyr(ro=ro, vo=vo)
        cluster.m = cluster.m / conversion.mass_in_msol(ro=ro, vo=vo)
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

    cluster.analyze()

def to_units(cluster, units, ro=8.0, vo=220.0):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to user defined units

    Parameters
    ----------
    cluster : class
        StarCluster
    units : str
        units to be converted to (nbody,galpy,pckms,kpckms,radec)
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
        cluster.to_nbody()
    elif units == "galpy":
        cluster.to_galpy(ro=ro, vo=vo)
    elif units == "pckms":
        cluster.to_pckms()
    elif units == "kpckms":
        cluster.to_kpckms()
    elif units == "radec":
        origin0 = cluster.origin
        cluster.to_radec(sortstars=False)
        cluster.to_origin(origin0)

def convert_binary_units(cluster,param,from_units,to_units):
        """Convert units of binary star parameters

            - This function is a wip
        
        Parameters
        ----------

            cluster : class
                StarCluster instance
            param : str
                list of strings for parameters to be convered (pb,semi,mass)
            from_units : str
                list of strings for current parameter units ((days, years, nbody), (pc, su, au, nbody), (Msun,nbody))
            from_units : str (default: True)
                list of strings for converted parameter units ((days, years, nbody), (pc, su, au, nbody), (Msun,nbody))

        Returns
        -------

            None

        History:

           2020 - Written - Webb (UofT)
        """
        yrs = (self.rbar * 1296000.0 / (2.0 * np.pi)) ** 1.5 / np.sqrt(self.zmbar)
        days = 365.25 * yrs
        au = 1.49597870700e13
        pc = 1296000.0e0/(2.0*np.pi)*au
        rsun=6.960e10
        su=pc/rsun*self.rbar


        param=np.array(param)
        from_units=np.array(from_units)
        tp_units=np.array(to_units)

        for i in range(0,len(param)):
            p=param[i]

            if p=='pb':
                #Convert to nbody first
                if from_units[i]=='days':
                    self.pb/=days
                elif from_units[i]=='years':
                    self.pb/=yrs
                elif from_units[i]=='nbody':
                    pass
                else:
                    print('UNIT %s NOT FOUND' % from_units[i])

                if to_units[i]=='days':
                    self.pb*=days
                elif to_units[i]=='years':
                    self.pb*=yrs
                elif to_units[i]=='nbody':
                    pass
                else:
                    print('UNIT %s NOT FOUND' % from_units[i])

            elif p=='semi':
                #Convert to nbody first
                if from_units[i]=='pc':
                    self.semi/=self.rbar
                elif from_units[i]=='su':
                    self.semi/=su
                elif from_units[i]=='au':
                    self.semi/=(pc/au)*self.rbar
                elif from_units[i]=='nbody':
                    pass
                else:
                    print('UNIT %s NOT FOUND' % from_units[i])

                if to_units[i]=='pc':
                    self.semi*=self.rbar
                elif to_units[i]=='su':
                    self.semi*=su
                elif to_units[i]=='au':
                    self.semi*=(pc/au)*self.rbar
                elif to_units[i]=='nbody':
                    pass
                else:
                    print('UNIT %s NOT FOUND' % to_units[i])

            elif p=='mass':
                if from_units[i]=='Msun' or from_units[i]=='msun' :
                    self.m1/=self.zmbar
                    self.m2/=self.zmbar
                elif from_units[i]=='nbody':
                    pass

                if to_units=='Msun' or to_units[i]=='msun' :
                    self.m1*=self.zmbar
                    self.m2*=self.zmbar
                elif to_units[i]=='nbody':
                    pass

def to_centre(cluster, sortstars=True, centre_method=None):
    """Shift coordinates such that origin is the centre of mass

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
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
            cluster.to_cluster(centre_method=centre_method)

        cluster.x -= cluster.xc
        cluster.y -= cluster.yc
        cluster.z -= cluster.zc
        cluster.vx -= cluster.vxc
        cluster.vy -= cluster.vyc
        cluster.vz -= cluster.vzc

        cluster.origin = "centre"

        cluster.analyze(sortstars=sortstars)


def to_cluster(cluster, sortstars=True, centre_method=None):
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
    sortstars : bool
        sort star by radius after coordinate change (default: False)
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
                cluster.vz -= cluster.vzgc


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


        cluster.origin = "cluster"

        cluster.analyze(sortstars=sortstars)


def to_galaxy(cluster, sortstars=True):
    """Shift coordinates to galactocentric reference frame

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.units!='galaxy':
        if cluster.units == "radec" and cluster.origin == "sky":
            from_radec(cluster)

        elif cluster.origin != "galaxy":
            if cluster.origin == "centre":
                cluster.to_cluster(sortstars=False)

            cluster.x += cluster.xgc
            cluster.y += cluster.ygc
            cluster.z += cluster.zgc
            cluster.vx += cluster.vxgc
            cluster.vy += cluster.vygc
            cluster.vz += cluster.vzgc

            cluster.origin = "galaxy"

        cluster.analyze(sortstars=sortstars)


def to_sky(cluster, sortstars=True,):
    """Calculate on-sky position, proper motion, and radial velocity of cluster
        
    - Also changes units to radec

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)
    """
    cluster.to_radec()
    cluster.analyze(sortstars=sortstars)

def from_sky(cluster):
    """Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster

    - Also changes units to kpckms

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    cluster.from_radec()

def to_origin(cluster, origin, sortstars=True):
    """Shift cluster to origin as defined by keyword

    Parameters
    ----------
    cluster : class
        StarCluster
    origin : str
        origin of coordinate system to be shifted to (centre, cluster, galaxy, sky)
    sortstars : bool
        sort star by radius after coordinate change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if origin == "centre":
        cluster.to_centre(sortstars=sortstars)
    elif origin == "cluster":
        cluster.to_cluster(sortstars=sortstars)
    elif origin == "galaxy":
        cluster.to_galaxy(sortstars=sortstars)
    elif origin == "sky":
        cluster.to_sky(sortstars=sortstars)

    cluster.analyze(sortstars=sortstars)


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

    return cluster.units, cluster.origin, cluster.rorder, cluster.rorder_origin


def return_cluster(cluster, units0, origin0, rorder0, rorder_origin0):
    """ return cluster to a specific combination of units and origin

    Parameters
    ----------

    cluster : class
        StarCluster
    units0 : str
        units that StarCluster will be changed to
    origin0 : str
        origin that StarCluster will be changed to
    sortstars : bool
        sort star by radius after coordinate change (default: False)

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
        cluster.to_origin(origin0)

    cluster.rorder=rorder0
    cluster.rorder_origin=rorder_origin0
    cluster.analyze()

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
    vbar : float
        velocity conversion
    tbar : float
        time conversion

    History:

    2018 - Written - Webb (UofT)
    """
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)

    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre(sortstars=True)

    cluster.to_pckms()

    if mass:
        zmbar = np.sum(cluster.m)
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

    vbar = 0.06557 * np.sqrt(zmbar / rbar)
    tbar = rbar / vbar

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return zmbar,rbar,vbar,tbar

def virialize(cluster, specific=True, full=full_default, projected=False):
    """ Adjust stellar velocities so cluster is in virial equilibrium

    Parameters
    ----------
    cluster : class
        StarCluster
    specific : bool
        find specific energies (default: True)
    full: bool
        do full array of stars at once with numba (default: full_default)
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
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre()

    try:
        qv = np.sqrt(np.abs(0.5 / cluster.qvir))
    except:
        print("NEED TO CALCULATE ENERGIES FIRST")
        cluster.energies(specific=specific, full=full, projected=projected)
        qv = np.sqrt(np.abs(0.5 / cluster.qvir))

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

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
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)

    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre()

    r, theta, z = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
    vr, vtheta, vz = coords.rect_to_cyl_vec(
        cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
    )

    indx = vtheta < 0.0
    rindx = np.random.rand(cluster.ntot) < qrot

    vtheta[indx * rindx] = np.sqrt(vtheta[indx * rindx] * vtheta[indx * rindx])
    x,y,z = coords.cyl_to_rect(r, theta, z)
    vx,vy,vz = coords.cyl_to_rect_vec(
        vr, vtheta, vz, theta
    )

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return x,y,z,vx,vy,vz

def _get_grav(cluster):
    if cluster.units == "nbody":
        grav = 1.0
    elif cluster.units == "pckms":
        # G has units of pc (km/s)^2 / Msun
        grav = 4.302e-3
    elif cluster.units == "kpckms":
        # G has units of kpc (km/s)^2 / Msun
        grav = 4.302e-6
    else:
        grav = 1.0

    return grav
