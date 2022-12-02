""" Functions and Operations for analysing a cluster's tidal tails

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "to_tail",
    "tail_path",
    "tail_path_match",
]

try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

from galpy.util import _rotate_to_arbitrary_vector
from galpy import potential

import numpy as np

from ..analysis.orbits import orbital_path, orbital_path_match
from ..cluster.operations import *
from ..util.recipes import binmaker,nbinmaker,roaming_binmaker,roaming_nbinmaker
from ..util.coordinates import cart_to_sky
from ..util.constants import *

from ..util.plots import starplot,skyplot,_plot,_lplot,_scatter

import matplotlib.pyplot as plt


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
    cluster, dt=0.1, no=1000, nt=100, ntail=100, pot=None, dmax=None, bintype = 'fix', from_centre=False, skypath=False, 
    to_path=False,
    do_full=False,
    ro=None,
    vo=None,
    zo=None,
    solarmotion=None,
    plot=False,projected=False,
    **kwargs,
):
    """Calculate tail path +/- dt Gyr around the cluster

        Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    no : int
        number of timesteps for orbit integration (default:1000)
    nt : int
        number of points along the tail to set the tail spacing (default: 100)
    ntail : int
        number of points along the tail with roaming average (default: 1000)
    pot : class
        galpy Potential that orbit is to be integrate in (default: None)
    dmax : float
        maximum distance (assumed to be same units as cluster) from orbital path to be included in generating tail path (default: None)
    bintype : str
        type of binning for tail stars (default : 'fix')
    from_centre : bool
        genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
    skypath : bool
        return sky coordinates instead of cartesian coordinates (default: False)
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

    #dmax is assumed to have same units as cluster
    if dmax is not None:
        if units0=='nbody':
            dmax*=cluster.rbar/1000.0
        elif units0=='pckms':
            dmax/=1000.
        elif units0=='galpy':
            dmax*=cluster._ro
        elif units0=='radec' and not skypath:
            dist=np.sqrt(cluster.xgc**2.+cluster.ygc**2.+cluster.zgc**2.)
            dmax=dist*np.tan(dmax)
        elif units0=='kpckms' and skypath:
            dist=np.sqrt(cluster.xgc**2.+cluster.ygc**2.+cluster.zgc**2.)
            dmax=np.arctan(dmax/dist)

    top, xop, yop, zop, vxop, vyop, vzop = orbital_path(
        cluster,
        dt=dt,
        nt=no,
        pot=pot,
        from_centre=from_centre,
        skypath=skypath,
        initialize=False,
        ro=ro,
        vo=vo,
        zo=zo,
        solarmotion=solarmotion,
    )

    path=(top, xop, yop, zop, vxop, vyop, vzop)


    if bintype=='fix':
        if ntail > nt:
            t_lower, t_mid, t_upper, t_hist = roaming_binmaker(top, nbin=nt,ntot=ntail)
        else:
            t_lower, t_mid, t_upper, t_hist = binmaker(top, nbin=nt)
    elif bintype=='num':
        if ntail>nt:
            t_lower, t_mid, t_upper, t_hist = roaming_nbinmaker(top, nbin=nt,ntot=ntail)
        else:
            t_lower, t_mid, t_upper, t_hist = nbinmaker(top, nbin=nt)


    tstar, dprog, dpath = orbital_path_match(
        cluster=cluster, dt=dt, nt=no, pot=pot, path=path, from_centre=from_centre, skypath=skypath, to_path=to_path,do_full=do_full, ro=ro, vo=vo, zo=zo, solarmotion=solarmotion,projected=projected
    )

    if dmax is None:
        dindx=np.ones(len(tstar),dtype=bool)
    else:
        dindx = (np.fabs(dpath) <= dmax)

    ttail = np.array([])
    xtail = np.array([])
    ytail = np.array([])
    ztail = np.array([])
    vxtail = np.array([])
    vytail = np.array([])
    vztail = np.array([])

    for i in range(0, len(t_mid)):
        indx = (tstar >= t_lower[i]) * (tstar <= t_upper[i]) * dindx
        if np.sum(indx) > 0:
            ttail = np.append(ttail, t_mid[i])
            xtail = np.append(xtail, np.mean(cluster.x[indx]))
            ytail = np.append(ytail, np.mean(cluster.y[indx]))
            ztail = np.append(ztail, np.mean(cluster.z[indx]))
            vxtail = np.append(vxtail, np.mean(cluster.vx[indx]))
            vytail = np.append(vytail, np.mean(cluster.vy[indx]))
            vztail = np.append(vztail, np.mean(cluster.vz[indx]))

    if skypath:
        ratail,dectail,disttail,pmratail,pmdectail,vlostail=cart_to_sky(xtail, ytail, ztail, vxtail, vytail, vztail)


    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)
        if skypath:
            skyplot(cluster,coords='radec',overplot=overplot)
            _lplot(ratail,dectail,overplot=True)
        else:
            starplot(cluster,coords='xy',overplot=overplot)
            _lplot(xtail,ytail,overplot=True)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    if skypath:
        return ttail,ratail,dectail,disttail,pmratail,pmdectail,vlostail
    else:
        return ttail, xtail, ytail, ztail, vxtail, vytail, vztail


def tail_path_match(
    cluster,
    dt=0.1,
    no=1000,
    nt=100,
    ntail=100,
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
    """Match stars to a position along the tail path of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    dt : float
        timestep that StarCluster is to be moved to
    no : int
        number of timesteps for orbit integration (default:1000)
    nt : int
        number of points along the tail to set the tail spacing (default: 100)
    ntail : int
        number of points along the tail with roaming average (default: 1000)
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
        distance along the path to the progenitor
    dpath : 
        distance to centre of the tail path bin (default) or the tail path (to_path = True)

    History
    -------
    2018 - Written - Webb (UofT)
    """
 
    if path is None:
        path = tail_path(
            cluster, dt=dt, no=no, nt=nt, ntail=ntail, pot=pot, from_centre=from_centre, skypath=skypath, ro=ro, vo=vo,zo=zo,solarmotion=solarmotion,
        )

    return orbital_path_match(cluster=cluster,dt=dt,nt=no,pot=pot,path=path,from_centre=from_centre,
        skypath=skypath,to_path=to_path,do_full=do_full,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion,plot=plot,projected=projected,**kwargs)