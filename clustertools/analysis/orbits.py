""" Functions and Operations related to the cluster's orbit

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "initialize_orbit",
    "initialize_orbits",
    "interpolate_orbit",
    "interpolate_orbits",
    "orbit_interpolate",
    "orbits_interpolate",
    "orbital_path",
    "orbital_path_match",
]

from galpy.orbit import Orbit

try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion


from galpy import potential

import numpy as np
import matplotlib.pyplot as plt

from ..util.recipes import interpolate, binmaker
from ..util.plots import starplot,skyplot,_plot,_lplot,_scatter
from ..util.constants import *
from ..util.units import _convert_length,_convert_time,_convert_velocity

import astropy.coordinates as coord
import astropy.units as u

try:
    import amuse.units.units as u
except:
    pass

def initialize_orbit(cluster, from_centre=False, ro=None, vo=None, zo = None, solarmotion=None):
    """ Initialize a galpy orbit instance for the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    from_centre : bool
        intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)

    Returns
    -------
    orbit : class
        GALPY orbit

    History
    -------
    2018 - Written - Webb (UofT)
    """

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    if cluster.units == "radec":
        o = Orbit(
            [
                cluster.ra_gc,
                cluster.dec_gc,
                cluster.dist_gc,
                cluster.pmra_gc,
                cluster.pmdec_gc,
                cluster.vlos_gc,
            ],
            radec=True,
            ro=ro,
            vo=vo,
            zo=zo,
            solarmotion=solarmotion,
        )
    else:

        cluster.save_cluster()
        units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

        cluster.to_galpy()

        x, y, z = cluster.xgc, cluster.ygc, cluster.zgc
        vx, vy, vz = cluster.vxgc, cluster.vygc, cluster.vzgc

        if from_centre:
            x+=cluster.xc
            y+=cluster.yc
            z+=cluster.zc
            vx+=cluster.vxc
            vy+=cluster.vyc
            vz+=cluster.vzc


        R, phi, z = coords.rect_to_cyl(x, y, z)
        vR, vT, vz = coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)
        o = Orbit(
            [R, vR, vT, z, vz, phi], ro=ro, vo=vo, zo=zo, solarmotion=solarmotion
        )
    
        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return o


def initialize_orbits(cluster, ro=None, vo=None, zo = None, solarmotion=None):
    """Initialize a galpy orbit for every star in the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)

    Returns
    -------
    orbit : class
        GALPY orbit

    History
    -------
    2018 - Written - Webb (UofT)
    """

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    if cluster.units == "radec" and cluster.origin=='sky':

        vxvv=np.column_stack([cluster.ra,cluster.dec,cluster.dist,cluster.pmra,cluster.pmdec,cluster.vlos])

        os = Orbit(
            vxvv,
            radec=True,
            ro=ro,
            vo=vo,
            zo=zo,
            solarmotion=solarmotion,
        )

    elif cluster.units == "radec" and (cluster.origin=='cluster' or cluster.origin=='centre'):
        print('NO METHOD FOR INITALIZING STELLAR ORBITS WITH RESPET TO CENTRE OR CLUSTER WHEN UNITS==RADEC')
    else:

        cluster.save_cluster()
        units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

        cluster.to_galpy()
     
        x, y, z = cluster.x, cluster.y, cluster.z
        vx, vy, vz = cluster.vx, cluster.vy, cluster.vz

        R, phi, z = coords.rect_to_cyl(x, y, z)
        vR, vT, vz = coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)

        vxvv = np.column_stack([R, vR, vT, z, vz, phi])

        if cluster.origin=='cluster' or cluster.origin=='centre':
            os = Orbit(vxvv, ro=ro, vo=vo, zo=zo, solarmotion=[0,0,0])
        else:
            os = Orbit(vxvv, ro=ro, vo=vo, zo=zo, solarmotion=solarmotion)

        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return os


def _integrate_orbit(
    cluster, pot=None, tfinal=None, nt=1000, from_centre=False, ro=None, vo=None, zo=None, solarmotion=None, plot=False
):
    """Integrate a galpy orbit instance for the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential that orbit is to be integrate in (default: None)
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    from_centre : bool
        intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)
    plot : float
        show plot of cluster's orbit

    Returns
    -------
    ts : float
        timesteps
    o : class
        galpy orbit

    History
    -------
       2018 - Written - Webb (UofT)
    """

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    o = initialize_orbit(cluster,from_centre=from_centre,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)

    if tfinal is None:
        tfinal=12./conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='pckms' or cluster.units=='pcmyr':
        tfinal/=(1000.*conversion.time_in_Gyr(ro=ro, vo=vo))
    elif cluster.units=='amuse':
        tfinal=tfinal.value_in(u.Gyr)/conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='kpckms' or cluster.units=='radec' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
        tfinal/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        tfinal*=((cluster.tbar/1000.)/conversion.time_in_Gyr(ro=ro, vo=vo))

    ts = np.linspace(0, tfinal, nt)
    o.integrate(ts, pot)

    if plot:
        o.plot()

    return ts, o

def _integrate_orbits(
    cluster, pot=None, tfinal=None, nt=1000, ro=None, vo=None,zo = None, solarmotion=None, plot=False
):
    """Integrate a galpy orbit instance for each star

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential for host cluster that orbit is to be integrated in
        if None, assume a Plumme Potential
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)
    plot : float
        show plot of cluster's orbit

    Returns
    -------
    ts : float
        timesteps
    o : class
        galpy orbit

    History
    -------
       2018 - Written - Webb (UofT)
    """

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    if tfinal is None:
        tfinal=12./conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='pckms' or cluster.units=='pcmyr':
        tfinal/=(1000.*conversion.time_in_Gyr(ro=ro, vo=vo))
    elif cluster.units=='amuse':
        tfinal=tfinal.value_in(u.Gyr)/conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='kpckms' or cluster.units=='radec' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
        tfinal/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        tfinal*=((cluster.tbar/1000.)/conversion.time_in_Gyr(ro=ro, vo=vo))

    os = initialize_orbits(cluster, ro=ro,vo=vo, zo=zo, solarmotion=solarmotion)
    ts = np.linspace(0, tfinal, nt)

    #integrate orbits of stars in combined potential of GC and galaxy
    os.integrate(ts, pot)

    if plot:
        os.plot()

    return ts, os

def interpolate_orbit(
    cluster,
    pot=None,
    tfinal=None,
    nt=1000,
    from_centre=False,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None,
):
    """
    Interpolate past or future position of cluster and escaped stars

    Parameters
    ----------
    cluster : class
        StarCluster
    cluster_pot : class
        Galpy potential for host cluster that orbit is to be integrated in
        if None, assume a Plumme Potential
    pot : class
        galpy Potential that orbit is to be integrate in (default: None)
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    from_centre : bool
        intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)

    Returns
    -------
    x,y,z : float
        interpolated positions of each star
    vx,vy,vz : float
        interpolated velocities of each star    

    History
    -------
    2021 - Written - Webb (UofT)
    """

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()
        tfinal=tfinal.value_in(u.Myr)

    ts,o=_integrate_orbit(cluster, pot=pot, tfinal=tfinal, nt=nt, from_centre=from_centre, ro=ro, vo=vo,zo=zo,solarmotion=solarmotion, plot=False)

    if cluster.units=='radec':
        if cluster.origin!='sky':
            cluster.to_sky()
        xgc,ygc,zgc=o.ra(ts[-1]),o.dec(ts[-1]),o.dist(ts[-1])
        vxgc,vygc,vzgc=o.pmra(ts[-1]),o.pmdec(ts[-1]),o.vlos(ts[-1])  

        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    else:
        xgc,ygc,zgc=o.x(ts[-1]),o.y(ts[-1]),o.z(ts[-1])
        vxgc,vygc,vzgc=o.vx(ts[-1]),o.vy(ts[-1]),o.vz(ts[-1])

        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        xgc=_convert_length(xgc,'kpckms',cluster)
        ygc=_convert_length(ygc,'kpckms',cluster)
        zgc=_convert_length(zgc,'kpckms',cluster)
        vxgc=_convert_velocity(vxgc,'kpckms',cluster)
        vygc=_convert_velocity(vygc,'kpckms',cluster)
        vzgc=_convert_velocity(vzgc,'kpckms',cluster)

    return xgc,ygc,zgc,vxgc,vygc,vzgc

def interpolate_orbits(
    cluster,
    pot=None,
    tfinal=None,
    nt=1000,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None
):
    """
    Interpolate past or future position of stars within the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential for host cluster that orbit is to be integrated in
        if None, assume a Plumme Potential
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)

    Returns
    -------
    x,y,z : float
        interpolated positions of each star
    vx,vy,vz : float
        interpolated velocities of each star    

    History
    -------
    2021 - Written - Webb (UofT)
    """

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()
        tfinal=tfinal.value_in(u.Myr)

    if cluster.origin=='centre':
        ts,o=_integrate_orbits(cluster, pot=pot, tfinal=tfinal, nt=nt, ro=ro, vo=vo, zo=zo,solarmotion=[0,0,0], plot=False)
    elif cluster.origin=='cluster':
        ts,o=_integrate_orbits(cluster, pot=pot, tfinal=tfinal, nt=nt, ro=ro, vo=vo, zo=zo,solarmotion=[0,0,0], plot=False)
    else:
        ts,o=_integrate_orbits(cluster, pot=pot, tfinal=tfinal, nt=nt, ro=ro, vo=vo, zo=zo, solarmotion=solarmotion,plot=False)

    if cluster.units=='radec' and cluster.origin=='sky':
        x,y,z=o.ra(ts[-1]),o.dec(ts[-1]),o.dist(ts[-1])
        vx,vy,vz=o.pmra(ts[-1]),o.pmdec(ts[-1]),o.vlos(ts[-1])
        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    elif cluster.units=='radec' and (cluster.origin=='centre' or cluster.origin=='cluster'):
        print('CANT INTEGRATE ORBITS WITH FROM_CENTRE OR FROM_CLUSTER AND RETURN IN SKY COORDINATES')

    else:
        x,y,z=o.x(ts[-1]),o.y(ts[-1]),o.z(ts[-1])
        vx,vy,vz=o.vx(ts[-1]),o.vy(ts[-1]),o.vz(ts[-1])
        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        x=_convert_length(x,'kpckms',cluster)
        y=_convert_length(y,'kpckms',cluster)
        z=_convert_length(z,'kpckms',cluster)
        vx=_convert_velocity(vx,'kpckms',cluster)
        vy=_convert_velocity(vy,'kpckms',cluster)
        vz=_convert_velocity(vz,'kpckms',cluster)

    return x,y,z,vx,vy,vz

# Renamed to interpolate_orbit and interpolates_orbit, but will keep old names for legacy purposes
def orbit_interpolate(
    cluster,
    pot=None,
    tfinal=None,
    nt=1000,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None,
):
    """
    Interpolate past or future position of cluster and escaped stars

    - same as interpolate_orbit, but included for legacy purposes

    Parameters
    ----------
    cluster : class
        StarCluster
    cluster_pot : class
        Galpy potential for host cluster that orbit is to be integrated in
        if None, assume a Plumme Potential
    pot : class
        galpy Potential that orbit is to be integrate in (default: None)
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    from_centre : bool
        intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)

    Returns
    -------
    x,y,z : float
        interpolated positions of each star
    vx,vy,vz : float
        interpolated velocities of each star    

    History
    -------
    2021 - Written - Webb (UofT)
    """
    return interpolate_orbit(cluster, pot=pot, tfinal=tfinal, nt=nt, ro=ro, vo=vo,zo=zo,solarmotion=solarmotion,plot=False)

def orbits_interpolate(
    cluster,
    pot=None,
    tfinal=None,
    nt=1000,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None,
):
    """
    Interpolate past or future position of stars within the cluster

    - same as interpolate_orbits, but kept for legacy purposes

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential for host cluster that orbit is to be integrated in
        if None, assume a Plumme Potential
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)

    Returns
    -------
    x,y,z : float
        interpolated positions of each star
    vx,vy,vz : float
        interpolated velocities of each star    

    History
    -------
    2021 - Written - Webb (UofT)
    """
    return interpolate_orbits(cluster, cluster_pot=cluster_pot, pot=pot, tfinal=tfinal, nt=nt, ro=ro, vo=vo,zo=zo,solarmotion=solarmotion,plot=False)


def orbital_path(
    cluster,
    tfinal=0.1,
    nt=1000,
    pot=None,
    from_centre=False,
    skypath=False,
    initialize=False,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None,
    plot=False,
    **kwargs,
):
    """Calculate the cluster's orbital path

    Parameters
    ----------
    cluster : class
        StarCluster
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 0.1 Gyr)
    nt : int
        number of timesteps
    pot : class
        galpy Potential that orbit is to be integrate in (default: None)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    skypath : bool
        return sky coordinates instead of cartesian coordinates (default: False)
    initialize : bool
        Initialize and return Orbit (default: False)
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)
    plot : bool
        plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)

    Returns
    -------
    t : float
        times for which path is provided
    x,y,z : float
        orbit positions
    vx,vy,vz : float
        orbit velocity
    o : class
        galpy orbit (if initialize==True)
    History
    -------
    2018 - Written - Webb (UofT)
    """

    #Legacy - allow for dt to be given instead of tfinal
    tfinal=kwargs.get('dt',tfinal)

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    if tfinal is None:
        tfinal=0.1/conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='pckms' or cluster.units=='pcmyr':
        tfinal/=(1000.*conversion.time_in_Gyr(ro=ro, vo=vo))
    elif cluster.units=='amuse':
        tfinal=tfinal.value_in(u.Myr)/(1000.*conversion.time_in_Gyr(ro=ro, vo=vo))
    elif cluster.units=='kpckms' or cluster.units=='radec' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
        tfinal/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        tfinal*=((cluster.tbar/1000.)/conversion.time_in_Gyr(ro=ro, vo=vo))

    o = initialize_orbit(cluster, from_centre=from_centre, ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)

    ts = np.linspace(0, -1.0 * tfinal, nt)
    o.integrate(ts, pot)

    R, phi, z = coords.rect_to_cyl(o.x(ts[-1]), o.y(ts[-1]), o.z(ts[-1]))
    vR, vT, vz = coords.rect_to_cyl_vec(
        o.vx(ts[-1]), o.vy(ts[-1]), o.vz(ts[-1]), o.x(ts[-1]), o.y(ts[-1]), o.z(ts[-1])
    )
    o = Orbit(
        [R / ro, vR / vo, vT / vo, z / ro, vz / vo, phi],
        ro=ro,
        vo=vo,
        zo=zo,
        solarmotion=solarmotion,
    )
    ts = np.linspace(
        -1.0 * tfinal,
        tfinal,
        nt,
    )

    o.integrate(ts, pot)

    if skypath:
        ra = np.array(o.ra(ts))
        dec = np.array(o.dec(ts))
        dist = np.array(o.dist(ts))
        pmra = np.array(o.pmra(ts))
        pmdec = np.array(o.pmdec(ts))
        vlos = np.array(o.vlos(ts))

        if cluster.units == "pckms" or cluster.units=='pcmyr':
            t = ts * conversion.time_in_Gyr(ro=ro, vo=vo) * 1000.0
        elif cluster.units=='amuse':
            t = ts * conversion.time_in_Gyr(ro=ro, vo=vo) * 1000.0 | u.Myr
        elif cluster.units == "nbody":
            t = ts * conversion.time_in_Gyr(ro=ro, vo=vo) * 1000.0 / cluster.tbar
        elif cluster.units == "galpy":
            t = ts
        elif cluster.units=='kpckms' or cluster.units=='radec' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
            t = ts * conversion.time_in_Gyr(ro=ro, vo=vo)
        else:
            print('TIME RETURNED IN GALPY UNITS')
            t=ts

        if plot:
            filename = kwargs.pop("filename", None)
            overplot = kwargs.pop("overplot", False)
            skyplot(cluster)
            plt.plot(ra,dec)
            if filename != None:
                plt.savefig(filename)

        if initialize:
            return t, ra, dec, dist, pmra, pmdec, vlos, o
        else:
            return t, ra, dec, dist, pmra, pmdec, vlos
    else:
        x = np.array(o.x(ts))
        y = np.array(o.y(ts))
        z = np.array(o.z(ts))
        vx = np.array(o.vx(ts))
        vy = np.array(o.vy(ts))
        vz = np.array(o.vz(ts))

        x=_convert_length(x,'kpckms',cluster)
        y=_convert_length(y,'kpckms',cluster)
        z=_convert_length(z,'kpckms',cluster)
        vx=_convert_velocity(vx,'kpckms',cluster)
        vy=_convert_velocity(vy,'kpckms',cluster)
        vz=_convert_velocity(vz,'kpckms',cluster)
        t=_convert_time(ts,'galpy',cluster)

        if plot:
            filename = kwargs.pop("filename", None)
            overplot = kwargs.pop("overplot", False)
            starplot(cluster,coords='xy',overplot=overplot)
            _lplot(x,y,overplot=True)

            if filename != None:
                plt.savefig(filename)

        if initialize:
            return t, x, y, z, vx, vy, vz, o
        else:
            return t, x, y, z, vx, vy, vz

def orbital_path_match(
    cluster,
    tfinal=0.1,
    nt=1000,
    pot=None,
    path=None,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None,
    plot=False,
    projected=False,
    **kwargs,
):
    """Match stars to a position along the orbital path of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 0.1 Gyr)
    nt : int
        number of timesteps
    pot : class
        galpy Potential that orbit is to be integrate in (default: None)
    path : array
        array of (t,x,y,x,vx,vy,vz) corresponding to the tail path. If none path is calculated (default: None)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    skypath : bool
        return sky coordinates instead of cartesian coordinates (default: False)
        if True, projected is set to True
    to_path : bool
        measure distance to the path itself instead of distance to central point along the path (default: False)
    do_full : bool
        calculate dpath all at once in a single numpy array (can be memory intensive) (default:False)
    ro : float
        distance to the Galactic centre (Default: None)
    vo : float
        circular velocity at ro (Default: None)
    zo : float
        Sun's distance above the Galactic plane (default: None)
    solarmotion : float
        array representing U,V,W of Sun (default: None)
    plot : bool
        plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)
    projected : bool
        match to projected orbital path, which means matching just x and y coordinates or Ra and Dec coordinates (not z, or dist) (default:False)

    Returns
    -------
    tstar : float
        orbital time associated with star
    dprog : float
        distance along the orbit to the progenitor
    dpath : 
        distance to centre of the orbital path bin (Default) or the orbit path (to_path = True)

    History
    -------
    2018 - Written - Webb (UofT)
    """

    #Legacy - allow for dt to be given instead of tfinal
    tfinal=kwargs.get('dt',tfinal)

    if ro is None: ro=cluster._ro
    if vo is None: vo=cluster._vo
    if zo is None: zo=cluster._zo
    if solarmotion is None: solarmotion=cluster._solarmotion

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if skypath:
        cluster.to_radec()
        projected=True
    else:
        cluster.to_galaxy(sortstars=False)
        cluster.to_kpckms()

    if path is None:
        ts, x, y, z, vx, vy, vz = orbital_path(
            cluster,
            tfinal=tfinal,
            nt=nt,
            pot=pot,
            from_centre=from_centre,
            skypath=skypath,
            initialize=False,
            ro=ro,
            vo=vo,
            zo=zo,
            solarmotion=solarmotion,
        )
    else:
        ts, x, y, z, vx, vy, vz = path

    pindx = np.argmin(np.fabs(ts))

    dx = np.tile(np.array(x), cluster.ntot).reshape(
        cluster.ntot, len(ts)
    ) - np.repeat(cluster.x, len(ts)).reshape(cluster.ntot, len(ts))
    dy = np.tile(np.array(y), cluster.ntot).reshape(
        cluster.ntot, len(ts)
    ) - np.repeat(cluster.y, len(ts)).reshape(cluster.ntot, len(ts))
    dz = np.tile(np.array(z), cluster.ntot).reshape(
        cluster.ntot, len(ts)
    ) - np.repeat(cluster.z, len(ts)).reshape(cluster.ntot, len(ts))

    if projected:
        dr = np.sqrt(dx ** 2.0 + dy ** 2.0)
    else:
        dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + dz ** 2.0)

    indx = np.argmin(dr, axis=1)
    tstar = ts[indx]
    dpath = np.amin(dr, axis=1)

    dxo = x[1:] - x[0:-1]
    dyo = y[1:] - y[0:-1]
    dzo = z[1:] - z[0:-1]

    dprogx = np.cumsum(np.fabs(dxo))
    dprogy = np.cumsum(np.fabs(dyo))
    dprogz = np.cumsum(np.fabs(dzo))

    dprogx = np.insert(dprogx, 0, 0.0)
    dprogy = np.insert(dprogy, 0, 0.0)
    dprogz = np.insert(dprogz, 0, 0.0)

    if projected:
        dprogr = np.sqrt(dprogx ** 2.0 + dprogy ** 2.0)
    else:
        dprogr = np.sqrt(dprogx ** 2.0 + dprogy ** 2.0 + dprogz ** 2.0)

    dprog = dprogr[indx] - dprogr[pindx]

    # Find distance to path instead of to central point
    if to_path:
        dxo = np.append(dxo, dxo[-1])
        dyo = np.append(dyo, dyo[-1])
        dzo = np.append(dzo, dzo[-1])

        if do_full:
            # Typically it is too expensive to calculate dpath all at once, but will allow option via do_full

            if projected:
                ovec = np.column_stack([dxo, dyo])
                mag_ovec = np.sqrt(dxo ** 2.0 + dyo ** 2.0)
                svec = np.column_stack([dx[:, indx], dy[:, indx]])
                mag_svec = dr[:, indx]
                theta = np.arccos(np.dot(ovec[indx], svec) / (mag_ovec[indx] * mag_svec))
                dpath = mag_svec * np.sin(theta)
            else:
                ovec = np.column_stack([dxo, dyo, dzo])
                mag_ovec = np.sqrt(dxo ** 2.0 + dyo ** 2.0 + dzo ** 2.0)
                svec = np.column_stack([dx[:, indx], dy[:, indx], dz[:, indx]])
                mag_svec = dr[:, indx]
                theta = np.arccos(np.dot(ovec[indx], svec) / (mag_ovec[indx] * mag_svec))
                dpath = mag_svec * np.sin(theta)
        else:
            # Need to optimize this via numba
            dpath = np.array([])
            if projected:
                for i in range(0, cluster.ntot):
                        ovec = [dxo[indx[i]], dyo[indx[i]]]
                        mag_ovec = np.sqrt(
                            dxo[indx[i]] ** 2.0 + dyo[indx[i]] ** 2.0)

                        svec = [dx[i, indx[i]], dy[i, indx[i]]]
                        mag_svec = dr[i, indx[i]]

                        theta = np.arccos(
                            (ovec[0] * svec[0] + ovec[1] * svec[1])
                            / (mag_ovec * mag_svec)
                        )
                        dpath = np.append(dpath, mag_svec * np.sin(theta))

            else:
                for i in range(0, cluster.ntot):
                    ovec = [dxo[indx[i]], dyo[indx[i]], dzo[indx[i]]]
                    mag_ovec = np.sqrt(
                        dxo[indx[i]] ** 2.0 + dyo[indx[i]] ** 2.0 + dzo[indx[i]] ** 2.0
                    )

                    svec = [dx[i, indx[i]], dy[i, indx[i]], dz[i, indx[i]]]
                    mag_svec = dr[i, indx[i]]

                    theta = np.arccos(
                        (ovec[0] * svec[0] + ovec[1] * svec[1] + ovec[2] * svec[2])
                        / (mag_ovec * mag_svec)
                    )
                    dpath = np.append(dpath, mag_svec * np.sin(theta))

    # Assign negative to stars with position vectors in opposite direction as local angular momentum vector
    if not projected:
        rgc = np.column_stack([x[indx], y[indx], z[indx]])
        vgc = np.column_stack([vx[indx], vy[indx], vz[indx]])
        lz = np.cross(rgc, vgc)

        rstar = np.column_stack(
            [
                cluster.x - x[indx],
                cluster.y - y[indx],
                cluster.z - z[indx],
            ]
        )

        ldot = np.sum(rstar * lz, axis=1)

        dpath[ldot < 0] *= -1.0
    else:
        yindx = np.argmin(np.fabs(dy),axis=1)
        dpath*=np.sign(dy[np.arange(0,len(dy)),yindx])


    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if skypath:
            xlabel=r"$\rm D_{prog} (Degree)$"
            ylabel=r"$ \rm D_{path} (Degree)$"
        else:
            xlabel=r"$\rm D_{prog} (kpc)$"
            ylabel=r"$ \rm D_{path} (kpc)$"

        _scatter(dprog,dpath,xlabel=xlabel,ylabel=ylabel,overplot=overplot)

        if filename != None:
            plt.savefig(filename)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return np.array(tstar), np.array(dprog), np.array(dpath)