""" Functions and Operations that heavily use galpy and focus on the cluster's orbit

"""
__author__ = "Jeremy J Webb"

from galpy.orbit import Orbit
from galpy.util import bovy_coords, bovy_conversion
from galpy import potential
from galpy.potential import LogarithmicHaloPotential, MWPotential2014, rtide
from galpy.actionAngle import actionAngleStaeckel
from galpy.actionAngle.actionAngleIsochroneApprox import actionAngleIsochroneApprox

import numpy as np

from ..util.recipes import rotate, interpolate, binmaker
from .operations import save_cluster, return_cluster
from .profiles import rho_prof
from ..util.plots import *

import astropy.coordinates as coord
import astropy.units as u

def rtidal(
    cluster,
    pot=MWPotential2014,
    rtiterate=0,
    rtconverge=0.9,
    rgc=None,
    ro=8.0,
    vo=220.0,
    verbose=False,
):
    """Calculate tidal radius of the cluster
    - The calculation uses Galpy (Bovy 2015_, which takes the formalism of Bertin & Varri 2008 to calculate the tidal radius
    -- Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005
    -- Bovy J., 2015, ApJS, 216, 29
    - riterate = 0 corresponds to a single calculation of the tidal radius based on the cluster's mass (cluster.mtot)
    -- Additional iterations take the mass within the previous iteration's calculation of the tidal radius and calculates the tidal
       radius again using the new mass until the change is less than 90%
    - for cases where the cluster's orbital parameters are not set, it is possible to manually set rgc which is assumed to be in kpc.

    Parameters

    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate tidal radius (default: MWPotential2014)
    rtiterate : int
        how many times to iterate on the calculation of r_t (default: 0)
    rtconverge : float
        criteria for tidal radius convergence within iterations (default 0.9)
    rgc : float
        Manually set galactocentric distance in kpc at which the tidal radius is to be evaluated (default: None)
    ro : float
        GALPY radius scaling parameter
    vo : float
        GALPY velocity scaling parameter
    verbose : bool
        Print information about iterative calculation of rt

    Returns
    -------
    rt : float
        tidal radius

    History
    _______
    2019 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)

    cluster.to_centre()
    cluster.to_galpy()

    if rgc != None:
        R = rgc / ro
        z = 0.0
    else:
        R = np.sqrt(cluster.xgc ** 2.0 + cluster.ygc ** 2.0)
        z = cluster.zgc

    # Calculate rtide
    rt = rtide(pot, R, z, M=cluster.mtot,use_physical=False)
    nit = 0
    for i in range(0, rtiterate):
        msum = 0.0

        indx = cluster.r < rt
        msum = np.sum(cluster.m[indx])

        rtnew = rtide(pot, R, z, M=msum,use_physical=False)

        if verbose:
            print(rt, rtnew, rtnew / rt, msum / cluster.mtot)

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

    if units0 == "pckms":
        rt *= 1000.0 * ro
    elif units0 == "kpckms":
        rt *= ro
    elif units0 == "nbody":
        rt *= 1000.0 * ro / cluster.rbar

    return_cluster(cluster, units0, origin0)

    return rt


def rlimiting(
    cluster,
    pot=MWPotential2014,
    rgc=None,
    ro=8.0,
    vo=220.0,
    nrad=20,
    projected=False,
    plot=False,
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
    ro : float
        GALPY radius scaling parameter
    vo : float
        GALPY velocity scaling parameter
    nrad : int
        number of radial bins used to calculate density profile (Default: 20)
    projected : bool
        use projected values (default: False)
    plot : bool
        plot the density profile and mark the limiting radius of the cluster (default: False)

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
    units0, origin0 = save_cluster(cluster)

    cluster.to_centre()
    cluster.to_galpy()


    if rgc != None:
        R = rgc / ro
        z = 0.0
    else:
        R = np.sqrt(cluster.xgc ** 2.0 + cluster.ygc ** 2.0)
        z = cluster.zgc

    # Calculate local density:
    rho_local = potential.evaluateDensities(
        pot, R, z, ro=ro, vo=vo, use_physical=False
    ) / bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo)

    rprof, pprof, nprof = rho_prof(cluster, nrad=nrad, projected=projected)

    if pprof[-1] > rho_local:
        rl = rprof[-1]
    elif pprof[0] < rho_local:
        rl = 0.0
    else:
        indx = np.argwhere(pprof < rho_local)[0][0]
        r1 = (rprof[indx - 1], pprof[indx - 1])
        r2 = (rprof[indx], pprof[indx])

        rl = interpolate(r1, r2, y=rho_local)

    if verbose:
        print("FINAL RL: ", rl * ro * 1000.0, "pc")

    if units0 == "pckms":
        rl *= 1000.0 * ro
    elif units0 == "kpckms":
        rl *= ro
    elif units0 == "nbody":
        rl *= 1000.0 * ro / cluster.rbar

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    if plot:
        if verbose:
            print("LOCAL DENSITY = ", rho_local)

        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            rprof *= ro * 1000.0 / cluster.rbar
            pprof *= (
                bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo)
                * (cluster.rbar ** 3.0)
                / cluster.zmbar
            )
            rho_local *= (
                bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo)
                * (cluster.rbar ** 3.0)
                / cluster.zmbar
            )
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms":
            rprof *= ro * 1000.0
            pprof *= bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo)
            rho_local *= bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo)
            xunits = " (pc)"
            if projected:
                yunits = " Msun/pc^2"
            else:
                yunits = " Msun/pc^3"
        elif cluster.units == "kpckms":
            rprof *= ro
            pprof *= bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo) * (1000.0 ** 3.0)
            rho_local *= bovy_conversion.dens_in_msolpc3(ro=ro, vo=vo) * (1000.0 ** 3.0)

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
        nlplot(
            x,
            y,
            xlabel=r"$R %s$" % (xunits),
            ylabel=r"$\rho %s$" % (yunits),
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )
        nlplot(x, np.ones(len(x)) * rho_local, "--", overplot=True)
        nlplot(np.ones(len(y)) * rl, y, "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    return rl


def initialize_orbit(cluster, from_centre=False, ro=8.0, vo=220.0):
    """ Initialize a galpy orbit instance for the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    ro : float
        galpy distance scale (default: 8.)
    vo : float
        galpy velocity scale (default: 220.)

    Returns
    -------
    orbit : class
        GALPY orbit

    History
    -------
    2018 - Written - Webb (UofT)
    """
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
            solarmotion=[-11.1, 24.0, 7.25],
        )
    else:
        units0, origin0 = save_cluster(cluster)
        cluster.to_galpy()

        if from_centre:
            x, y, z = (
                cluster.xgc + cluster.xc,
                cluster.ygc + cluster.yc,
                cluster.zgc + cluster.zc,
            )
            vx, vy, vz = (
                cluster.vxgc + cluster.vxc,
                cluster.vygc + cluster.vyc,
                cluster.vzgc + cluster.vzc,
            )
        else:
            x, y, z = cluster.xgc, cluster.ygc, cluster.zgc
            vx, vy, vz = cluster.vxgc, cluster.vygc, cluster.vzgc

        R, phi, z = bovy_coords.rect_to_cyl(x, y, z)
        vR, vT, vz = bovy_coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)
        o = Orbit(
            [R, vR, vT, z, vz, phi], ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25]
        )
        
        return_cluster(cluster, units0, origin0)

    return o


def initialize_orbits(cluster, ro=8.0, vo=220.0):
    """Initialize a galpy orbit for every star in the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    ro : float
        galpy distance scale (default: 8.)
    vo : float
        galpy velocity scale (default: 220.)

    Returns
    -------
    orbit : class
        GALPY orbit

    History
    -------
    2018 - Written - Webb (UofT)
    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_galpy()

    x, y, z = cluster.x, cluster.y, cluster.z
    vx, vy, vz = cluster.vx, cluster.vy, cluster.vz

    R, phi, z = bovy_coords.rect_to_cyl(x, y, z)
    vR, vT, vz = bovy_coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)

    vxvv = np.column_stack([R, vR, vT, z, vz, phi])
    os = Orbit(vxvv, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

    return_cluster(cluster, units0, origin0)

    return os


def integrate_orbit(
    cluster, pot=MWPotential2014, tfinal=12.0, nt=1000, ro=8.0, vo=220.0, plot=False
):
    """Integrate a galpy orbit instance for the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential that orbit is to be integrate in (default: MWPotential2014)
    tfinal : float
        final time (in Gyr) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
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
    o = initialize_orbit(cluster)
    ts = np.linspace(0, tfinal / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), nt)
    o.integrate(ts, pot)

    if plot:
        o.plot()

    return ts, o

def integrate_orbits(
    cluster, pot=MWPotential2014, tfinal=12.0, nt=1000, ro=8.0, vo=220.0, plot=False
):
    """Integrate a galpy orbit instance for each star

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential that orbit is to be integrate in (default: MWPotential2014)
    tfinal : float
        final time (in Gyr) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
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
    os = initialize_orbits(cluster)
    ts = np.linspace(0, tfinal / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), nt)
    os.integrate(ts, pot)

    if plot:
        os.plot()

    return ts, os


def orbit_interpolate(
    cluster,
    dt,
    pot=MWPotential2014,
    from_centre=False,
    do_tails=False,
    rmin=None,
    rmax=None,
    emin=None,
    emax=None,
    indx=None,
    ro=8.0,
    vo=220.0,
):
    """
    NAME: Interpolate past or future position of cluster and escaped stars

    - When moving the cluster centre and stars backwards or forwards along their orbits, stars within the cluster are shifted with the cluster centre.
    Tail stars, identified using either rmin/rmax, emin/emax, or indx can be integrated separately in the potential 
    - Note that this function operates on the cluster and changes the positions and velocities of all stars

    Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    pot : class
        galpy Potential that orbit is to be integrate in (default: MWPotential2014)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    do_tails : bool
        interpolate the orbits of tail stars separately (default: False)
    rmin/rmax : float
        radial range corresponding to cluster (needed to identify tail stars)
    emin/emax : float
        energy range corresponding to cluster (needed to identify tail stars)
    indx : bool 
        specific subset of stars
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)

    Returns
    -------
    x,y,z : float
        interpolated positions of each star
    x,y,z : float
        interpolated velocities of each star    

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.tphys += dt
    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()

    if do_tails:

        cluster.to_cluster()
        if from_centre:
            cluster.to_centre()

        if rmin == None:
            rmin = np.min(cluster.r)
        if rmax == None:
            rmax = np.max(cluster.r)
        rindx = (cluster.r >= rmin) * (cluster.r <= rmax)

        if len(cluster.etot) == cluster.ntot:
            if emin == None:
                emin = np.min(cluster.etot)
            if emax == None:
                emax = np.max(cluster.etot)
            eindx = (cluster.etot >= emin) * (cluster.etot <= emax)
        else:
            eindx = cluster.id > -1

        if indx is None:
            indx = rindx * eindx
        else:
            indx*=(rindx*eindx)

        tindx = np.invert(indx)

        cluster.to_galaxy()

    else:
        indx = cluster.id > -1

    if cluster.orbit is None:
        cluster.orbit = initialize_orbit(cluster, from_centre)

    ts = np.linspace(0, dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), 10)

    cluster.orbit.integrate(ts, pot)

    cluster.to_kpckms()

    x,y,z=cluster.x,cluster.y,cluster.z
    vx,vy,vz=cluster.vx,cluster.vy,cluster.vz

    if from_centre:
        dx = cluster.orbit.x(ts[-1]) - cluster.xc - cluster.xgc
        dy = cluster.orbit.y(ts[-1]) - cluster.yc - cluster.ygc
        dz = cluster.orbit.z(ts[-1]) - cluster.zc - cluster.zgc
        dvx = cluster.orbit.vx(ts[-1]) - cluster.vxc - cluster.vxgc
        dvy = cluster.orbit.vy(ts[-1]) - cluster.vyc - cluster.vygc
        dvz = cluster.orbit.vz(ts[-1]) - cluster.vzc - cluster.vzgc
    else:
        dx = cluster.orbit.x(ts[-1]) - cluster.xgc
        dy = cluster.orbit.y(ts[-1]) - cluster.ygc
        dz = cluster.orbit.z(ts[-1]) - cluster.zgc
        dvx = cluster.orbit.vx(ts[-1]) - cluster.vxgc
        dvy = cluster.orbit.vy(ts[-1]) - cluster.vygc
        dvz = cluster.orbit.vz(ts[-1]) - cluster.vzgc

    x[indx] += dx
    y[indx] += dy
    z[indx] += dz
    vx[indx] += dvx
    vy[indx] += dvy
    vz[indx] += dvz

    if from_centre:
        xc, yc, zc = 0.0, 0.0, 0.0
        vxc, vyc, vzc = 0.0, 0.0, 0.0
    else:
        xc += dx
        yc += dy
        zc += dz
        vxc += dvx
        vyc += dvy
        vzc += dvz

    xgc, ygc, zgc = (
        cluster.orbit.x(ts[-1]),
        cluster.orbit.y(ts[-1]),
        cluster.orbit.z(ts[-1]),
    )
    vxgc, vygc, vzgc = (
        cluster.orbit.vx(ts[-1]),
        cluster.orbit.vy(ts[-1]),
        cluster.orbit.vz(ts[-1]),
    )

    if do_tails:
        cluster.to_galaxy()
        cluster.to_galpy()

        xt, yt, zt = cluster.x[tindx], cluster.y[tindx], cluster.z[tindx]
        vxt, vyt, vzt = cluster.vx[tindx], cluster.vy[tindx], cluster.vz[tindx]

        R, phi, z = bovy_coords.rect_to_cyl(x, y, z)
        vR, vT, vz = bovy_coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)

        vxvv = np.column_stack([R, vR, vT, z, vz, phi])
        otail = Orbit(vxvv, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

        cluster.to_kpckms()

        ts = np.linspace(0, dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), 10)

        otail.integrate(ts, pot)

        x[tindx] = np.array(otail.x(ts[-1]))
        y[tindx] = np.array(otail.y(ts[-1]))
        z[tindx] = np.array(otail.z(ts[-1]))
        vx[tindx] = np.array(otail.vx(ts[-1]))
        vy[tindx] = np.array(otail.vy(ts[-1]))
        vz[tindx] = np.array(otail.vz(ts[-1]))

    return_cluster(cluster, units0, origin0)

    return x,y,z,vx,vy,vz


def orbital_path(
    cluster,
    dt=0.1,
    nt=100,
    pot=MWPotential2014,
    from_centre=False,
    skypath=False,
    initialize=False,
    ro=8.0,
    vo=220.0,
    plot=False
):
    """Calculate the cluster's orbital path

    Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    nt : int
        number of timesteps
    pot : class
        galpy Potential that orbit is to be integrate in (default: MWPotential2014)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    sky_path : bool
        return sky coordinates instead of cartesian coordinates (default: False)
    initialize : bool
        Initialize and return Orbit (default: False)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
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
    o = initialize_orbit(cluster, from_centre=from_centre)

    ts = np.linspace(0, -1.0 * dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), nt)
    o.integrate(ts, pot)

    R, phi, z = bovy_coords.rect_to_cyl(o.x(ts[-1]), o.y(ts[-1]), o.z(ts[-1]))
    vR, vT, vz = bovy_coords.rect_to_cyl_vec(
        o.vx(ts[-1]), o.vy(ts[-1]), o.vz(ts[-1]), o.x(ts[-1]), o.y(ts[-1]), o.z(ts[-1])
    )
    o = Orbit(
        [R / ro, vR / vo, vT / vo, z / ro, vz / vo, phi],
        ro=ro,
        vo=vo,
        solarmotion=[-11.1, 24.0, 7.25],
    )
    ts = np.linspace(
        -1.0 * dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo),
        dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo),
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

        if cluster.units == "pckms":
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo)
        elif cluster.units == "nbody":
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo) / cluster.tstar
        elif cluster.units == "galpy":
            t = ts
        else:
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo)

        if plot:
            filename = kwargs.pop("filename", None)
            overplot = kwargs.pop("overplot", False)
            skyplot(cluster)
            plt.plot(ra,dec)
            if filename != None:
                plt.savefig(filename)


        if initialize:
            cluster.orbit = o
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

        if cluster.units == "pckms":
            x *= 1000.0
            y *= 1000.0
            z *= 1000.0
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo)
        elif cluster.units == "nbody":
            x *= 1000.0 / cluster.rbar
            y *= 1000.0 / cluster.rbar
            z *= 1000.0 / luster.rbar
            vx /= cluster.vstar
            vy /= cluster.vstar
            vz /= cluster.vstar
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo) / cluster.tstar

        elif cluster.units == "galpy":
            x /= ro
            y /= ro
            z /= ro
            vx /= vo
            vy /= vo
            vz /= vo
            t = ts
        else:
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo)

        if plot:
            filename = kwargs.pop("filename", None)
            overplot = kwargs.pop("overplot", False)
            starplot(cluster,coord='xy',overplot=overplot)
            nlplot(x,y,overplot=True)

            if filename != None:
                plt.savefig(filename)

        if initialize:
            return t, x, y, z, vx, vy, vz, o
        else:
            return t, x, y, z, vx, vy, vz

def orbital_path_match(
    cluster,
    dt=0.1,
    nt=100,
    pot=MWPotential2014,
    from_centre=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    plot=False,
):
    """Match stars to a position along the orbital path of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    nt : int
        number of timesteps
    pot : class
        galpy Potential that orbit is to be integrate in (default: MWPotential2014)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    to_path : bool
        measure distance to the path itself instead of distance to central point along the path (default: False)
    do_full : bool
        calculate dpath all at once in a single numpy array (can be memory intensive) (default:False)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    plot : bool
        plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)

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

    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_kpckms()

    t, x, y, z, vx, vy, vz, o = orbital_path(
        cluster,
        dt=dt,
        nt=nt,
        pot=pot,
        from_centre=from_centre,
        initialize=True,
        ro=ro,
        vo=vo,
    )

    ts = t / bovy_conversion.time_in_Gyr(ro=ro, vo=vo)

    x = o.x(ts)
    y = o.y(ts)
    z = o.z(ts)
    vx = o.vx(ts)
    vy = o.vy(ts)
    vz = o.vz(ts)

    pindx = np.argmin(np.fabs(ts))

    dx = np.tile(np.array(o.x(ts)), cluster.ntot).reshape(
        cluster.ntot, len(ts)
    ) - np.repeat(cluster.x, len(ts)).reshape(cluster.ntot, len(ts))
    dy = np.tile(np.array(o.y(ts)), cluster.ntot).reshape(
        cluster.ntot, len(ts)
    ) - np.repeat(cluster.y, len(ts)).reshape(cluster.ntot, len(ts))
    dz = np.tile(np.array(o.z(ts)), cluster.ntot).reshape(
        cluster.ntot, len(ts)
    ) - np.repeat(cluster.z, len(ts)).reshape(cluster.ntot, len(ts))
    dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + dz ** 2.0)

    indx = np.argmin(dr, axis=1)
    tstar = ts[indx] * bovy_conversion.time_in_Gyr(ro=ro, vo=vo)
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

    dprogr = np.sqrt(dprogx ** 2.0 + dprogy ** 2.0 + dprogz ** 2.0)
    dprog = dprogr[indx] - dprogr[pindx]

    # Find distance to path instead of to central point
    if to_path:
        dxo = np.append(dxo, dxo[-1])
        dyo = np.append(dyo, dyo[-1])
        dzo = np.append(dzo, dzo[-1])

        if do_full:
            # Typically it is too expensive to calculate dpath all at once, but will allow option via do_full

            ovec = np.column_stack([dxo, dyo, dzo])
            mag_ovec = np.sqrt(dxo ** 2.0 + dyo ** 2.0 + dzo ** 2.0)
            svec = np.column_stack([dx[:, indx], dy[:, indx], dz[:, indx]])
            mag_svec = dr[:, indx]
            theta = np.arccos(np.dot(ovec[indx], svec) / (mag_ovec[indx] * mag_svec))
            dpath = mag_svec * np.sin(theta)
        else:
            # Need to optimize this via numba
            dpath = np.array([])
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
    rgc = np.column_stack([o.x(ts[indx]), o.y(ts[indx]), o.z(ts[indx])])
    vgc = np.column_stack([o.vx(ts[indx]), o.vy(ts[indx]), o.vz(ts[indx])])
    lz = np.cross(rgc, vgc)

    rstar = np.column_stack(
        [
            cluster.x - o.x(ts[indx]),
            cluster.y - o.y(ts[indx]),
            cluster.z - o.z(ts[indx]),
        ]
    )

    ldot = np.sum(rstar * lz, axis=1)

    dpath[ldot < 0] *= -1.0

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)
        nscatter(dprog,dpath,xlabel="Dprog",ylabel="Dpath",overplot=overplot)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0)

    return np.array(tstar), np.array(dprog), np.array(dpath)


def tail_path(
    cluster, dt=0.1, nt=100, pot=MWPotential2014, from_centre=False, ro=8.0, vo=220.0,
    plot=False
):
    """Calculate tail path +/- dt Gyr around the cluster

        Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    nt : int
        number of timesteps
    pot : class
        galpy Potential that orbit is to be integrate in (default: MWPotential2014)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    plot : bool
        plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)

    Returns
    -------
    t : float
        times for which path is provided
    x,y,z : float
        tail path positions
    vx,vy,vz : float
        tail path velocities
    History
    -------
    2018 - Written - Webb (UofT)
    2019 - Implemented numpy array preallocation to minimize runtime - Nathaniel Starkman (UofT)
    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_kpckms()

    to, xo, yo, zo, vxo, vyo, vzo, o = orbital_path(
        cluster,
        dt=dt,
        nt=nt,
        pot=pot,
        from_centre=from_centre,
        initialize=True,
        ro=ro,
        vo=vo,
    )
    tstar, dprog, dpath = orbital_path_match(
        cluster=cluster, dt=dt, nt=nt, pot=pot, from_centre=from_centre, ro=ro, vo=vo
    )

    t_lower, t_mid, t_upper, t_hist = binmaker(to, nbin=nt)
    ttail = []
    xtail = []
    ytail = []
    ztail = []
    vxtail = []
    vytail = []
    vztail = []

    for i in range(0, len(t_mid)):
        indx = (tstar >= t_lower[i]) * (tstar <= t_upper[i])
        if np.sum(indx) > 0:
            ttail = np.append(ttail, t_mid[i])
            xtail = np.append(xtail, np.mean(cluster.x[indx]))
            ytail = np.append(ytail, np.mean(cluster.y[indx]))
            ztail = np.append(ztail, np.mean(cluster.z[indx]))
            vxtail = np.append(vxtail, np.mean(cluster.vx[indx]))
            vytail = np.append(vytail, np.mean(cluster.vy[indx]))
            vztail = np.append(vztail, np.mean(cluster.vz[indx]))

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)
        starplot(cluster,coord='xy',overplot=overplot)
        nlplot(xtail,ytail,overplot=True)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0)

    return ttail, xtail, ytail, ztail, vxtail, vytail, vztail


def tail_path_match(
    cluster,
    dt=0.1,
    nt=100,
    pot=MWPotential2014,
    from_centre=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    plot=False,
):
    """Match stars to a position along the tail path of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    nt : int
        number of timesteps
    pot : class
        galpy Potential that orbit is to be integrate in (default: MWPotential2014)
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    to_path : bool
        measure distance to the path itself instead of distance to central point along the path (default: False)
    do_full : bool
        calculate dpath all at once in a single numpy array (can be memory intensive) (default:False)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    plot : bool
        plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)

    Returns
    -------
    tstar : float
        orbital time associated with star
    dprog : float
        distance along the path to the progenitor
    dpath : 
        distance to centre of the tail path bin (default) or the tail path (to_path = True)

    History
    -------
    2018 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_kpckms()

    ts, x, y, z, vx, vy, vz = tail_path(
        cluster, dt=dt, nt=nt, pot=pot, from_centre=from_centre, ro=ro, vo=vo
    )
    pindx = np.argmin(np.fabs(ts))

    dx = np.tile(x, cluster.ntot).reshape(cluster.ntot, len(ts)) - np.repeat(
        cluster.x, len(ts)
    ).reshape(cluster.ntot, len(ts))
    dy = np.tile(y, cluster.ntot).reshape(cluster.ntot, len(ts)) - np.repeat(
        cluster.y, len(ts)
    ).reshape(cluster.ntot, len(ts))
    dz = np.tile(z, cluster.ntot).reshape(cluster.ntot, len(ts)) - np.repeat(
        cluster.z, len(ts)
    ).reshape(cluster.ntot, len(ts))
    dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + dz ** 2.0)

    indx = np.argmin(dr, axis=1)
    dpath = np.amin(dr, axis=1)
    tstar = ts[indx]  # *bovy_conversion.time_in_Gyr(ro=ro,vo=vo)

    dxo = x[1:] - x[0:-1]
    dyo = y[1:] - y[0:-1]
    dzo = z[1:] - z[0:-1]

    dprogx = np.cumsum(np.fabs(dxo))
    dprogy = np.cumsum(np.fabs(dyo))
    dprogz = np.cumsum(np.fabs(dzo))

    dprogx = np.insert(dprogx, 0, 0.0)
    dprogy = np.insert(dprogy, 0, 0.0)
    dprogz = np.insert(dprogz, 0, 0.0)

    dprogr = np.sqrt(dprogx ** 2.0 + dprogy ** 2.0 + dprogz ** 2.0)
    dprog = dprogr[indx] - dprogr[pindx]

    # Find distance to path instead of to central point
    if to_path:
        dxo = np.append(dxo, dxo[-1])
        dyo = np.append(dyo, dyo[-1])
        dzo = np.append(dzo, dzo[-1])

        if do_full:
            # Typically it is too expensive to calculate dpath all at once, but will allow option via do_full

            ovec = np.column_stack([dxo, dyo, dzo])
            mag_ovec = np.sqrt(dxo ** 2.0 + dyo ** 2.0 + dzo ** 2.0)
            svec = np.column_stack([dx[:, indx], dy[:, indx], dz[:, indx]])
            mag_svec = dr[:, indx]
            theta = np.arccos(np.dot(ovec[indx], svec) / (mag_ovec[indx] * mag_svec))
            dpath = mag_svec * np.sin(theta)
        else:
            # Need to optimize this via numba
            dpath = np.array([])
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
    rgc = np.column_stack([x[indx], y[indx], z[indx]])
    vgc = np.column_stack([vx[indx], vy[indx], vz[indx]])
    lz = np.cross(rgc, vgc)

    rstar = np.column_stack(
        [cluster.x - x[indx], cluster.y - y[indx], cluster.z - z[indx]]
    )

    ldot = np.sum(rstar * lz, axis=1)
    dpath[ldot < 0] *= -1

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)
        nscatter(dprog,dpath,xlabel="Dprog",ylabel="Dpath",overplot=overplot)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0)

    return np.array(tstar), np.array(dprog), np.array(dpath)


def get_cluster_orbit(gcname="mwglobularclusters",ro=8.0, vo=220.0):
    """Get the measured orbital parameters of a Galactic globular cluster
    - This is a simply wrapper for Orbit.from_name in galpy (Bovy 2015), which uses orbits measured by Vasiliev 2019 using Gaia DR2 via Galpy
    -- Bovy J., 2015, ApJS, 216, 29

    Parameters
    ----------
    gcname : str
        name of GC whose orbits is to be retrieved
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    Returns
    -------
    orbit : class
        galpy orbit

    History
    -------
    2019 - Written - Webb (UofT)

    """
    return Orbit.from_name(gcname,ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

def calc_actions(cluster, pot=MWPotential2014, ro=8.0, vo=220.0, **kwargs):
    """Calculate action angle values for each star 
    - This is a simple wrapper for calculating actions from an Orbit in galpy (Bovy 2015)
    -- Bovy J., 2015, ApJS, 216, 29

    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate actions
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)

    Returns
    -------
    JR,Jphi,Jz : float
        orbit actions
    OR,Ophi,Oz : float
        orbital frequencies
    TR,Tphi,Tz : float
        orbital periods

    Other Parameters
    ----------------
    type : str
        method for calculating actions (default: staeckel)
    delta : float
        focus for staeckel method (default: 0.45 - optimal for MWPotential2014)
    c : bool
        if True, always use C for calculations (default: True)
    kwargs : str
        key word arguments can be included for other action calculation methods in galpy
    History
    -------
    2019 - Written - Webb (UofT)
    """

    os = initialize_orbits(cluster, ro, vo)
    atype = kwargs.pop("type", "staeckel")
    delta = kwargs.pop("delta", 0.45)
    c = kwargs.pop("c", True)

    JR = os.jr(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    Jphi = os.jp(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    Jz = os.jz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    OR = os.Or(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    Ophi = os.Op(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    Oz = os.Oz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    TR = os.Tr(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    Tphi = os.Tp(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
    Tz = os.Tz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)

    return JR, Jphi, Jz, OR, Ophi, Oz, TR, Tphi, Tz

def ttensor(cluster, pot=MWPotential2014, ro=8.0, vo=220.0, eigenval=False):
    """Calculate the tidal tensor Tij=-d(Psi)(dxidxj)
    - This is a simple wrapper for calculating the tidal tensor in a potential in galpy (Bovy 2015)
    -- Bovy J., 2015, ApJS, 216, 29
    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate tidal tensor
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    eigenval : bool
        return eigenvalues if true (default; False)

    Returns
    -------
        Tidal Tensor
    History
    -------
    2018-03-21 - Written - Webb (UofT)
    """

    o = initialize_orbit(cluster, ro, vo)
    R=o.R()
    z=o.z()

    tij=potential.ttensor(R/ro,z/ro,eigenval=eigenval)

    return tij
