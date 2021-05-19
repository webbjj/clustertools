""" Functions and Operations for analysing a cluster's tidal tails

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "to_tail",
    "tail_path",
    "tail_path_match",
]

from galpy.util import bovy_conversion,_rotate_to_arbitrary_vector
from galpy import potential
from galpy.potential import MWPotential2014

import numpy as np
import matplotlib.pyplot as plt

from .orbit import orbital_path, orbital_path_match
from .operations import *
from ..util.recipes import binmaker
from ..util.coordinates import cart_to_sky

from ..util.plots import starplot,skyplot,_plot,_lplot,_scatter


def to_tail(cluster):
    """Calculate positions and velocities of stars when rotated such that clusters velocity vector
       points along x-axis

    - no change to coordinates in StarCluster

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    x_tail,y_tail,z_tail,vx_tail,vy_tail,vz_tail : float
        rotated coordinates with cluster's velocity vector point along x-axis

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre(sortstars=False)


    v_vec = np.array([cluster.vxgc, cluster.vygc, cluster.vzgc])
    new_v_vec = np.array([1.0, 0.0, 0.0])

    rot = _rotate_to_arbitrary_vector(
        np.atleast_2d(v_vec), new_v_vec, inv=False, _dontcutsmall=False
    )

    x_tail = (
        cluster.x * rot[:, 0, 0] + cluster.y * rot[:, 1, 0] + cluster.z * rot[:, 2, 0]
    )
    y_tail = (
        cluster.x * rot[:, 0, 1] + cluster.y * rot[:, 1, 1] + cluster.z * rot[:, 2, 1]
    )
    z_tail = (
        cluster.x * rot[:, 0, 2] + cluster.y * rot[:, 1, 2] + cluster.z * rot[:, 2, 2]
    )
    vx_tail = (
        cluster.vx * rot[:, 0, 0] + cluster.vy * rot[:, 1, 0] + cluster.vz * rot[:, 2, 0]
    )
    vy_tail = (
        cluster.vx * rot[:, 0, 1] + cluster.vy * rot[:, 1, 1] + cluster.vz * rot[:, 2, 1]
    )
    vz_tail = (
        cluster.vx * rot[:, 0, 2] + cluster.vy * rot[:, 1, 2] + cluster.vz * rot[:, 2, 2]
    )

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return x_tail,y_tail,z_tail,vx_tail,vy_tail,vz_tail

def tail_path(
    cluster, dt=0.1, nt=100, pot=MWPotential2014, from_centre=False, skypath=False, ro=8.0, vo=220.0,
    plot=False,**kwargs,
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
    skypath : bool
        return sky coordinates instead of cartesian coordinates (default: False)
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

    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_galaxy(sortstars=False)
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
        starplot(cluster,coords='xy',overplot=overplot)
        _lplot(xtail,ytail,overplot=True)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    if skypath:
        ratail,dectail,disttail,pmratail,pmdectail,vlostail=cart_to_sky(xtail, ytail, ztail, vxtail, vytail, vztail)
        return ttail,ratail,dectail,disttail,pmratail,pmdectail,vlostail
    else:
        return ttail, xtail, ytail, ztail, vxtail, vytail, vztail


def tail_path_match(
    cluster,
    dt=0.1,
    nt=100,
    pot=MWPotential2014,
    tail_path=None,
    from_centre=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    plot=False,
    **kwargs,
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
    tail_path : array
        array of (t,x,y,x,vx,vy,vz) corresponding to the tail path. If none path is calculated (default: None)
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
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_galaxy(sortstars=False)
    cluster.to_kpckms()

    if tail_path is None:
        ts, x, y, z, vx, vy, vz = tail_path(
            cluster, dt=dt, nt=nt, pot=pot, from_centre=from_centre, ro=ro, vo=vo
        )
    else:
        ts,x,y,z,vx,vy,vz=tail_path

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
        _scatter(dprog,dpath,xlabel=r"$\rm D_{prog} (kpc)$",ylabel=r"$ \rm D_{path} (kpc)$",overplot=overplot)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return np.array(tstar), np.array(dprog), np.array(dpath)