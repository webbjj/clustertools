""" Functions and Operations related to the cluster's orbit

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "initialize_orbit",
    "initialize_orbits",
    "integrate_orbit",
    "integrate_orbits",
    "orbit_interpolate",
    "orbital_path",
    "orbital_path_match",
    "calc_actions",
    "ttensor",
]

from galpy.orbit import Orbit
from galpy.util import bovy_coords, bovy_conversion
from galpy import potential
from galpy.potential import MWPotential2014
from galpy.actionAngle import actionAngleStaeckel
from galpy.actionAngle.actionAngleIsochroneApprox import actionAngleIsochroneApprox

import numpy as np
import matplotlib.pyplot as plt

from ..util.recipes import interpolate, binmaker
from .operations import save_cluster, return_cluster
from .profiles import rho_prof
from ..util.plots import starplot,skyplot,_plot,_lplot,_scatter

import astropy.coordinates as coord
import astropy.units as u

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
        units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
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
        
        return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

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

    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_galaxy(sortstars=False)
    cluster.to_galpy()

    x, y, z = cluster.x, cluster.y, cluster.z
    vx, vy, vz = cluster.vx, cluster.vy, cluster.vz

    R, phi, z = bovy_coords.rect_to_cyl(x, y, z)
    vR, vT, vz = bovy_coords.rect_to_cyl_vec(vx, vy, vz, x, y, z)

    vxvv = np.column_stack([R, vR, vT, z, vz, phi])
    os = Orbit(vxvv, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

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
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_galaxy(sortstars=False)

    if do_tails:

        cluster.to_cluster(sortstars=False)
        if from_centre:
            cluster.to_centre(sortstars=False)

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

        cluster.to_galaxy(sortstars=False)

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
        xc = dx
        yc = dy
        zc = dz
        vxc = dvx
        vyc = dvy
        vzc = dvz

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
        cluster.to_galaxy(sortstars=False)
        cluster.to_galpy()

        xt, yt, zt = cluster.x[tindx], cluster.y[tindx], cluster.z[tindx]
        vxt, vyt, vzt = cluster.vx[tindx], cluster.vy[tindx], cluster.vz[tindx]

        Rt, phit, zt = bovy_coords.rect_to_cyl(xt, yt, zt)
        vRt, vTt, vzt = bovy_coords.rect_to_cyl_vec(vxt, vyt, vzt, xt, yt, zt)

        vxvvt = np.column_stack([Rt, vRt, vTt, zt, vzt, phit])
        otail = Orbit(vxvvt, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])

        cluster.to_kpckms()

        ts = np.linspace(0, dt / bovy_conversion.time_in_Gyr(ro=ro, vo=vo), 10)

        otail.integrate(ts, pot)

        x[tindx] = np.array(otail.x(ts[-1]))
        y[tindx] = np.array(otail.y(ts[-1]))
        z[tindx] = np.array(otail.z(ts[-1]))
        vx[tindx] = np.array(otail.vx(ts[-1]))
        vy[tindx] = np.array(otail.vy(ts[-1]))
        vz[tindx] = np.array(otail.vz(ts[-1]))

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

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
    plot=False,
    **kwargs,
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
    skypath : bool
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
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo) / cluster.tbar
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
            t = ts * bovy_conversion.time_in_Gyr(ro=ro, vo=vo) / cluster.tbar

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
    dt=0.1,
    nt=100,
    pot=MWPotential2014,
    from_centre=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    plot=False,
    **kwargs,
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

    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_galaxy(sortstars=False)
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
        _scatter(dprog,dpath,xlabel=r"$\rm D_{prog} (kpc)$",ylabel=r"$ \rm D_{path} (kpc)$",overplot=overplot)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return np.array(tstar), np.array(dprog), np.array(dpath)

def calc_actions(cluster, pot=MWPotential2014, ro=8.0, vo=220.0, full=False, **kwargs):
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
    full : bool
        return orbital frequencies and periods (default : False)
    Returns
    -------
    JR,Jphi,Jz : float
        orbit actions

    if full:
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
    if full:
        OR = os.Or(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
        Ophi = os.Op(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
        Oz = os.Oz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
        TR = os.Tr(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
        Tphi = os.Tp(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)
        Tz = os.Tz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo, **kwargs)

        return JR, Jphi, Jz, OR, Ophi, Oz, TR, Tphi, Tz
    else:
        return JR, Jphi, Jz


def ttensor(cluster, pot=MWPotential2014, ro=8.0, vo=220.0, eigenval=False, t=0.):
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
    time : float
        time to evaluate tidal tensor. Necessary if tidal field is time dependent (default: 0.)

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
    phi=o.phi()

    tij=potential.ttensor(pot,R/ro,z/ro,phi=phi,t=t/bovy_conversion.time_in_Gyr(ro=ro, vo=vo),eigenval=eigenval)

    return tij
