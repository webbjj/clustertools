""" Functions and Operations that heavily use galpy and focus on the cluster's orbit

"""
__author__ = "Jeremy J Webb"


__all__ = [
    "rtidal",
    "rlimiting",
    "initialize_orbit",
    "initialize_orbits",
    "integrate_orbit",
    "orbit_interpolate",
    "orbital_path",
    "orbital_path_match",
    "tail_path",
    "tail_path_match",
    "get_cluster_orbit",
    "calc_actions",
    "ttensor"
]

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
    """Calculate rtidal.
    NAME:

       rtidal

    PURPOSE:

       Calculate tidal radius of the cluster
       --> The calculation uses Galpy, which takes the formalism of Bertin & Varri 2008 to calculate the tidal radius
       --> riterate = 0 corresponds to a single calculation of the tidal radius based on the cluster's mass (cluster.mtot)
       --> More iterations take the mass within the previous iterations calculation of the tidal radius and calculates the tidal
           radius again until the change is less than 90%

    Parameters

       cluster - StarCluster instance

       pot - GALPY potential used to calculate actions

       rtiterate - how many times to iterate on the calculation of r_t

       rtconverge - criteria for tidal radius convergence within iterations

       rgc - Set galactocentric distance at which the tidal radius is to be evaluated

       ro,vo - GALPY scaling parameters

    Returns

        rt

    HISTORY:

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

    cluster.rt = rt

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
    verbose=False,
    **kwargs
):
    """
    NAME:

       rlimiting

    PURPOSE:

       Calculate limiting radius of the cluster
       --> The limiting radius is defined to be where the cluster's density reaches the local background density of the host galaxy

    Parameters

       cluster - StarCluster instance

       pot - GALPY potential used to calculate actions

       rgc - Set galactocentric distance at which the tidal radius is to be evaluated

       ro,vo - GALPY scaling parameters

       nrad - number of radial bins used to calculate density profile (Default: 20)

       projected - use projected values (Default: False)

       plot - plot the density profile and mark the limiting radius of the cluster (Default: False)

    KWARGS:

       Same as ..util.plots.nplot

    Returns

        rl

    HISTORY:

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

    cluster.rl = rl

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
    """
    NAME:

       initialize_orbit

    PURPOSE:

       Initialize a galpy orbit instance for the cluster

    Parameters

       cluster - StarCluster

       from_centre - genrate orbit from cluster's assigned galactocentric coordinates (default) or from its centre

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

    Returns

       orbit instance

    HISTORY:

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


    cluster.orbit = o


    return o


def initialize_orbits(cluster, ro=8.0, vo=220.0):
    """
    NAME:

       initialize_orbits

    PURPOSE:

       Initialize a galpy orbit for every stars in the cluster
       --> Note currently depends on version of galpy with Orbits, which is soon to be replace by a generic Orbit function

    Parameters

       cluster - StarCluster

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

    Returns

       orbit instance for each stars

    HISTORY:

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
    """
    NAME:

       integrate_orbit

    PURPOSE:

       Integrate a galpy orbit instance for the cluster

    Parameters

       cluster - StarCluster

       pot - Potential orbit is to be integrate in (Default: MWPotential2014)

       tfinal - final time (in Gyr) to integrate orbit to (Default: 12 Gyr)

       nt - number of timesteps

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

       plot - show plot of cluster's orbit

    Returns

       timesteps, orbit instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """
    cluster.orbit = initialize_orbit(cluster)
    ts = np.linspace(0, tfinal / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), nt)
    cluster.orbit.integrate(ts, pot)

    if plot:
        cluster.orbit.plot()

    return ts, cluster.orbit


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
    ro=8.0,
    vo=220.0,
):
    """
    NAME:

       orbit_interpolate

    PURPOSE:

       Move the cluster centre and stars backwards or forwards along its orbit assuming stars only feel a force from the galaxy

    Parameters

       cluster - StarCluster

       dt - timestep that StarCluster is to be moved to

       pot - Potential orbit is to be integrate in (Default: MWPotential2014)

       from_centre - interpolate from cluster's define position (default) or the measured centre of cluster 

       do_tails - interpolate the orbits of tail stars separately (Default: False)

       rmin/rmax - radial range corresponding to cluster (needed to identify tail stars)

       emin/emax - energy range corresponding to cluster (needed to identify tail stars)

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

    Returns

       None

    HISTORY:

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

        indx = rindx * eindx
        tindx = np.invert(indx)

        cluster.to_galaxy()

    else:
        indx = cluster.id > -1

    cluster.orbit = initialize_orbit(cluster, from_centre)
    ts = np.linspace(0, dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), 10)

    cluster.orbit.integrate(ts, pot)

    cluster.to_kpckms()

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

    cluster.x[indx] += dx
    cluster.y[indx] += dy
    cluster.z[indx] += dz
    cluster.vx[indx] += dvx
    cluster.vy[indx] += dvy
    cluster.vz[indx] += dvz

    if from_centre:
        cluster.xc, cluster.yc, cluster.zc = 0.0, 0.0, 0.0
        cluster.vxc, cluster.vyc, cluster.vzc = 0.0, 0.0, 0.0
    else:
        cluster.xc += dx
        cluster.yc += dy
        cluster.zc += dz
        cluster.vxc += dvx
        cluster.vyc += dvy
        cluster.vzc += dvz

    cluster.xgc, cluster.ygc, cluster.zgc = (
        cluster.orbit.x(ts[-1]),
        cluster.orbit.y(ts[-1]),
        cluster.orbit.z(ts[-1]),
    )
    cluster.vxgc, cluster.vygc, cluster.vzgc = (
        cluster.orbit.vx(ts[-1]),
        cluster.orbit.vy(ts[-1]),
        cluster.orbit.vz(ts[-1]),
    )

    if do_tails:
        cluster.to_galaxy()
        cluster.to_galpy()

        x, y, z = cluster.x[tindx], cluster.y[tindx], cluster.z[tindx]
        vx, vy, vz = cluster.vx[tindx], cluster.vy[tindx], cluster.vz[tindx]

        R, phi, z = bovy_coords.rect_to_cyl(x, y, z)
        vR, vT, vz = bovy_coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)

        vxvv = np.column_stack([R, vR, vT, z, vz, phi])
        otail = Orbit(vxvv, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

        cluster.to_kpckms()

        ts = np.linspace(0, dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), 10)

        otail.integrate(ts, pot)

        cluster.x[tindx] = np.array(otail.x(ts[-1]))
        cluster.y[tindx] = np.array(otail.y(ts[-1]))
        cluster.z[tindx] = np.array(otail.z(ts[-1]))

        cluster.vx[tindx] = np.array(otail.vx(ts[-1]))
        cluster.vy[tindx] = np.array(otail.vy(ts[-1]))
        cluster.vz[tindx] = np.array(otail.vz(ts[-1]))

    return_cluster(cluster, units0, origin0)


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
    """
    NAME:

       orbital_path

    PURPOSE:

       Calculate orbital path +/- dt Gyr around the cluster

    Parameters

       cluster - StarCluster

       dt - timestep that StarCluster is to be moved to

       nt - number of timesteps

       pot - Potential orbit is to be integrate in (Default: MWPotential2014)

       from_centre - interpolate from cluster's define position (default) or the measured centre of cluster 

       sky_path - return sky coordinates instead of cartesian coordinates (Default: False)

       initialize - Initialize and return Orbit (Default: False)

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

       plot - plot a snapshot of the cluster in galactocentric coordinates with the orbital path


    Returns

       t,x,y,z,vx,vy,vz

    HISTORY:

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
            cluster.orbit = o
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
    """
    NAME:

       orbital_path_match

    PURPOSE:

       Match stars to a position along the orbital path of the cluster

    Parameters

       cluster - StarCluster

       dt - timestep that StarCluster is to be moved to

       nt - number of timesteps

       pot - Potential orbit is to be integrate in (Default: MWPotential2014)

       from_centre - interpolate from cluster's define position (default) or the measured centre of cluster 

       to_path - measure distance to central point along the path (Default) or to the path itself 

       do_full - calculate dpath all at once in a single numpy array (can be memory intensive) (Default:False)

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

       plot - plot the distance of each star from the orbital path versus distance along the orbital path to the progenitor

    Returns

       tstar - orbital time associated with star

       dprog - distance along the orbit to the progenitor

       dpath - distance to centre of the orbital path bin (Default) or the orbit path (to_path = True)

    HISTORY:

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
    cluster, dt=0.1, nt=100, pot=MWPotential2014, from_centre=False, ro=8.0, vo=220.0
):
    """
    NAME:

       tail_path

    PURPOSE:

       Calculate tail path +/- dt Gyr around the cluster

    Parameters

       cluster - StarCluster

       dt - timestep that StarCluster is to be moved to

       nt - number of timesteps

       pot - Potential orbit is to be integrate in (Default: MWPotential2014)

       from_centre - interpolate from cluster's define position (default) or the measured centre of cluster 

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

    Returns

       t,x,y,z,vx,vy,vz
    HISTORY:

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
):
    """
    NAME:

       tail_path_match

    PURPOSE:

       Match stars to a position along the tail path of the cluster

    Parameters

       cluster - StarCluster

       dt - timestep that StarCluster is to be moved to

       nt - number of timesteps

       pot - Potential orbit is to be integrate in (Default: MWPotential2014)

       from_centre - interpolate from cluster's define position (default) or the measured centre of cluster 

       to_path - measure distance to central point along the path (Default) or to the path itself 

       ro - galpy distance scale (Default: 8.)

       vo - galpy velocity scale (Default: 220.)

    Returns

       tstar - tail time associated with star

       dprog - distance along the tail to the progenitor

       dpath - distance to centre of the tail path bin (Default) or the tail path (to_path = True)

    HISTORY:

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
    """
    NAME:

       get_cluster_orbit

    PURPOSE:

       Get the measured orbital parameters of a Galactic globular cluster has measured by Vasiliev 2019 using Gaia DR2 via Galpy

    Parameters

       gcname - name or list of GC names whose orbits are to be retrieved
        --> Note that list just returns the list to read, while 'all' gets the orbits for all the clusters

       pot - GALPY potential used to calculate actions

       ro,vo - GALPY scaling parameters

    Returns

        orbit

    HISTORY:

       2019 - Written - Webb (UofT)

    """
    return Orbit.from_name(gcname,ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

def calc_actions(cluster, pot=MWPotential2014, ro=8.0, vo=220.0, **kwargs):
    """
    NAME:

       calc_actions

    PURPOSE:

       Calculate action angle values for each star 

    Parameters

       cluster - StarCluster instance

       pot - GALPY potential used to calculate actions

       ro,vo - GALPY scaling parameters

    KWARGS:

       type - method for calculating actions (default: staeckel)

       delta - focus for staeckel method (default: 0.45 - optimal for MWPotential2014)

       c - if True, always use C for calculations

       Additional KWARGS can be included for other action angle calculation methods in galpy

    Returns

        JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz

    HISTORY:

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
    """

    NAME:

        ttensor

    PURPOSE:

        Calculate the tidal tensor Tij=-d(Psi)(dxidxj)

    Parameters

       cluster - StarCluster instance

       pot - GALPY potential used to calculate actions

       ro,vo - GALPY scaling parameters

       eigenval - return eigenvalues if true (optional; boolean)

    Returns

        Tidal Tensor

    HISTORY:

        2018-03-21 - Written - Webb (UofT)

    """

    o = initialize_orbit(cluster, ro, vo)
    R=o.R()
    z=o.z()

    tij=potential.ttensor(R/ro,z/ro,eigenval=eigenval)

    return tij
