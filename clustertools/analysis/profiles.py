""" Determine radial profiles of key properties

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "rho_prof",
    "m_prof",
    "alpha_prof",
    "sigv_prof",
    "beta_prof",
    "v_prof",
    "v2_prof",
    "eta_prof",
    "vcirc_prof",
]

import numpy as np

try:
    from galpy.util import coords
except:
    import galpy.util.bovy_coords as coords

from ..util.constants import *
from ..util.recipes import *
from ..util.constants import _get_grav
from ..util.plots import _lplot,_plot
from ..util.coordinates import sphere_coords
from .functions import mass_function, eta_function
from ..util.units import _convert_length,_convert_time,_convert_velocity,_convert_density,_convert_mass,_convert_square_velocity

import matplotlib.pyplot as plt

def rho_prof(
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
        elif cluster.units == "kpckms" or cluster.units == "kpcgyr" or cluster.units== "WDunits":
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


    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    rprof=_convert_length(rprof,units,cluster)
    pprof=_convert_density(pprof,units,cluster,projected=projected)

    if normalize:
        rprof/=cluster.rm

    return rprof, pprof, nprof


def m_prof(
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
    cumulative=False,
    plot=False,
    **kwargs
):
    """ Measure the mass profile of the cluster

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
    cumalitive : bool
        determine the cumulative mass profile instead (default: False)
    plot : bool 
        plot the mass profile (default: False)

    Returns
    -------
    rprof : float
        radius bins
    mprof : float
        mass within radial bin (or enclosed mass if cumalitve==True)
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
    mprof = np.array([])
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
        if cumulative:
            rindx = indx * (r < r_upper[i])
        else:
            rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])
        rprof=np.append(rprof,r_mean[i])

        mprof=np.append(mprof,np.sum(cluster.m[rindx]))
        nprof=np.append(nprof,np.sum(rindx))

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms" or cluster.units == "pcmyr":
            xunits = " (pc)"
            yunits = " Msun"
        elif cluster.units == "kpckms" or cluster.units == "kpcgyr" or cluster.units== "WDunits":
            xunits = " (kpc)"
            yunits = " Msun"
        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"
        else:
            xunits = ""
            yunits = ""

        if normalize:
            x/=cluster.rm

        x, y, n = rprof, mprof, nprof
        _lplot(
            x,
            y,
            xlabel=r"$R %s $" % xunits,
            ylabel=r"$M %s $" % yunits,
            title="Time = %f" % cluster.tphys,
            log=kwargs.pop('log',True),
            overplot=overplot,
            filename=filename,
            **kwargs,
        )

        if filename != None:
            plt.savefig(filename)

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    rprof=_convert_length(rprof,units,cluster)
    mprof=_convert_mass(mprof,units,cluster)

    if normalize:
        rprof/=cluster.rm

    return rprof, mprof, nprof

def alpha_prof(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    nrad=20,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    npop=None,
    indx=None,
    bins=None,
    projected=False,
    normalize=False,
    r_lower=None,
    r_upper=None,
    aerror=False,
    mcorr=None,
    plot=False,
    **kwargs
):
    """Measure the radial variation in the mass function

    - measure the delta_alpha parameter from Webb & Vesperini 2016
    - Webb, J.J. & Vesperini, E. 2016, MNRAS, 463, 2383

    Parameters
    ----------
    cluster : class
        StarCluster
    mmin/mmax : float
        minimum and maximum stellar mass
    nmass : int
        number of mass bins (default: 10)
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
        normalize radial bins by cluster's half-mass radius (default: True)
    r_lower : float
        lower limits to radial bins
    r_upper : float
        upper limits to radial bins
    aerror : bool
        return error in alpha calculations (default:True)
    mcorr : float
        completeness correction for masses (default: None)
    plot : bool 
        plot the alpha profile (default: False)

    Returns
    -------
    rprofn : float
        natural log of each radius bin (normalized by half-mass radius)
    aprof : float
        slope of the mass function in each bin
    dalpha : float
        radial variation in alpha calculated as delta_alpha = d(alpha)/d(ln(r/rm) 
    edalpha : float
        error in dalpha
    ydalpha : float
        y-intercept in fit to alpha vs ln(r/rm)
    eydalpha : float
        error in ydalpha
    rbinerror : float
        if mcorr is not None, output lowest corrected mass fraction at each radius

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


    if mcorr is None:        
        mcorr = np.ones(cluster.ntot)
        return_error=False
    else:
        return_error=True


    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))
        r_mean=np.zeros(len(r_mean))

        for i in range(0,len(r_lower)):
            rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])
            r_mean[i]=np.mean(r[rindx])
            r_hist[i]=np.sum(rindx)

    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    rbinerror=np.zeros(len(r_mean))
    rprofn=np.zeros(len(r_mean))
    aprof=np.zeros(len(r_mean))
    eaprof=np.zeros(len(r_mean))

    if normalize:
        if projected:
            rprofn=r_mean / cluster.rmpro
        else:
            rprofn=r_mean / cluster.rm
    else:
        rprofn=r_mean

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])


        if return_error:
            m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbinerror = mass_function(cluster,nmass=nmass,indx=rindx,projected=projected,mcorr=mcorr,plot=False,**kwargs)
            rbinerror[i]=np.amin(mbinerror)
        else:
            m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(cluster,nmass=nmass,indx=rindx,projected=projected,mcorr=None,plot=False,**kwargs)
            rbinerror[i]=1.

        if alpha > -100:
            aprof[i]=alpha
            eaprof[i]=ealpha

    if len(rprofn) > 3:
        if projected:
            (dalpha, ydalpha), V = np.polyfit(np.log(r_mean/cluster.rmpro), aprof, 1, cov=True)
        else:
            (dalpha, ydalpha), V = np.polyfit(np.log(r_mean/cluster.rm), aprof, 1, cov=True)

        edalpha = np.sqrt(V[0][0])
        eydalpha = np.sqrt(V[1][1])
    else:
        dalpha = -100.0
        ydalpha = 0.0
        edalpha = 0.0
        eydalpha = 0.0

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if normalize:
            xlabel=r"$\ln(r/r_m)$"
        else:
            xlabel=r"$\ln(r)$"

        _plot(
            np.log(rprofn),
            aprof,
            xlabel=xlabel,
            ylabel=r"$\alpha$",
            overplot=overplot,
            **kwargs
        )
        if projected:
            rfit=r_mean/cluster.rmpro
        else:
            rfit=r_mean/cluster.rm

        afit = dalpha * np.log(rfit) + ydalpha
        _lplot(np.log(rprofn), afit, overplot=True, label=(r"d$\alpha$ = %f" % dalpha))
        plt.legend()

        if filename != None:
            plt.savefig(filename)

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)
    rprofn=_convert_length(rprofn,units,cluster)


    if aerror:
        return rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha, eaprof
    else:
        return rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha


def sigv_prof(
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
    coord=None,
    normalize=False,
    plot=False,
    **kwargs,
):
    """Measure the radial variation in the velocity dispersion

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
    coord : str
        choose what coordinate the velocity dispersion profile is to be returned in (default None returns (sigx**2.+sigy**2.+sigz**2.)^1/2).
        Alternatively can ask for 'r', 'phi', or 'theta' for spherical coordinate velocity dispersions.
    normalize : bool
        normalize radial bins by cluster's half-mass radius (default: False)
    plot : bool 
        plot the velocity disperions profile (default: False)

    Returns
    -------
    rprofn : float
        natural log of radius (normalized by half-mass radius)
    sigvprof : float
        velocity dispersion

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

    # Build subcluster containing only stars in the full radial and mass range:
    indx = (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m <= mmax)
        * (v >= vmin)
        * (v <= vmax)
    )

    if kwmin is not None and len(cluster.kw) > 0:
        indx*=(cluster.kw >= kwmin)
    if kwmax is not None and len(cluster.kw) > 0:
        indx*=(cluster.kw <= kwmax)

    if emin is not None:
        indx *= cluster.etot >= emin
    if emin is not None:
        indx *= cluster.etot <= emax
    """

    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)

    if coord is not None:

        if projected:
            r, phi, z = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
            vr, vp, vz = coords.rect_to_cyl_vec(
                cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
            )
        else:
            r, phi, theta, vr, vp, vt = sphere_coords(cluster)
            vz=cluster.vz

        if coord =='r' or coord=='vr':
            v=vr
            ylabel=r"$\sigma_{v_r}$"
        elif coord=='phi' or coord=='vp' or coord=='vphi':
            v=vp
            ylabel=r"$\sigma_{v_p}$"
        elif coord=='theta' or coord=='vt' or coord=='vtheta' :
            v=vt
            ylabel=r"$\sigma_{v_t}$"
        elif coord=='x' or coord=='vx':
            v=cluster.vx
            ylabel=r"$\sigma_{v_x}$"
        elif coord=='y' or coord=='vy':
            v=cluster.vy
            ylabel=r"$\sigma_{v_y}$"
        elif coord=='z' or coord=='vz':
            v=vz
            ylabel=r"$\sigma_{v_z}$"
        elif coord=='v':
            if projected:
                v=cluster.vpro
            else:
                v=cluster.v
    else:
        if projected:
            ylabel=r"$\sigma_{v_{pro}}$"
        else:
            ylabel=r"$\sigma_v$"

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))

    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    if normalize:
        if projected:
            rprofn=r_mean / cluster.rmpro
        else:
            rprofn=r_mean / cluster.rm
    else:
        rprofn=r_mean

    sigvprof = np.zeros(len(rprofn))


    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])
        if bins is not None: r_hist[i]=np.sum(rindx)

        if np.sum(rindx) > 3.0:

            sigv = np.std(v[rindx])

            sigvprof[i]=sigv

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if normalize:
            xlabel=r"$r/r_m$"
        else:
            xlabel=r"$r$"

        _plot(
            rprofn,
            sigvprof,
            xlabel=xlabel,
            ylabel=ylabel,
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    rprofn=_convert_length(rprofn,units,cluster)
    sigvprof=_convert_velocity(sigvprof,units,cluster)

    return rprofn, sigvprof

def beta_prof(
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
    **kwargs,
):
    """Measure the anisotropy profile of the cluster

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
    rprofn : float
        natural log of radius (normalized by half-mass radius)
    betaprof : float
        orbital anisotropy parameter beta

    Other Parameters
    ----------------
    kwrags : str
        key word arguments for plotting

    History
    -------
    2020 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre(sortstars=normalize)
    elif normalize:
        cluster.sortstars()

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

    # Build subcluster containing only stars in the full radial and mass range:
    indx = (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m <= mmax)
        * (v >= vmin)
        * (v <= vmax)
    )

    if kwmin is not None and len(cluster.kw) > 0:
        indx*=(cluster.kw >= kwmin)
    if kwmax is not None and len(cluster.kw) > 0:
        indx*=(cluster.kw <= kwmax)

    if emin is not None:
        indx *= cluster.etot >= emin
    if emin is not None:
        indx *= cluster.etot <= emax

    """
    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)


    if projected:
        r, phi, z = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
        vr, vp, vz = coords.rect_to_cyl_vec(
            cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
        )
        v=cluster.vpro
    else:
        r, phi, theta, vr, vp, vt = sphere_coords(cluster)
        v=cluster.v

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))
    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    if normalize:
        if projected:
            rprofn=r_mean / cluster.rmpro
        else:
            rprofn=r_mean / cluster.rm
    else:
        rprofn=r_mean

    betaprof=np.zeros(len(rprofn))

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])
        if bins is not None: r_hist[i]=np.sum(rindx)

        if np.sum(rindx) > 3.0:

            sigr = np.std(vr[rindx])
            sigp = np.std(vp[rindx])

            if projected:
                sigt = np.zeros(len(vr))
                beta = sigp / sigr - 1.0
            else:
                sigt = np.std(vt[rindx])
                beta = 1.0 - (sigt ** 2.0 + sigp ** 2.0) / (2.0 * (sigr ** 2.0))

            betaprof[i]=beta

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if projected:
            if normalize:
                xlabel=r"$R/R_m$"
            else:
                xlabel=r"$R$"
        else:
            if normalize:
                xlabel=r"$r/r_m$"
            else:
                xlabel=r"$r$"

        _lplot(
            rprofn,
            betaprof,
            xlabel=xlabel,
            ylabel=r"$\beta$",
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    if not normalize: rprofn=_convert_length(rprofn,units,cluster)

    return rprofn, betaprof


def v_prof(
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
    coord=None,
    normalize=False,
    plot=False,
    **kwargs,
):
    """Measure the radial variation in the mean velocity 

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
    coord : str
        choose what coordinate the mean velocity profile is to be returned in (default None returns (vx**2.+vy**2.+vz**2.)^1/2).
        Alternatively can ask for 'vr', 'vphi', or 'vtheta' for spherical coordinate velocity dispersions.
    normalize : bool
        normalize radial bins by cluster's half-mass radius (default: False)
    plot : bool 
        plot the velocity disperions profile (default: False)

    Returns
    -------
    rprofn : float
        natural log of radius (normalized by half-mass radius)
    vprof : float
        mean velocity

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

    if coord is not None:

        if projected:
            r, phi, z = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
            vr, vp, vz = coords.rect_to_cyl_vec(
                cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
            )
        else:
            r, phi, theta, vr, vp, vt = sphere_coords(cluster)

        if coord =='r':
            v=vr
            ylabel=r"$<v_r>$"
        elif coord=='phi':
            v=vp
            ylabel=r"$<v_p>$"
        elif coord=='theta':
            v=vt
            ylabel=r"$<v_t>$"
        elif coord=='x' or coord=='vx':
            v=cluster.vx
            ylabel=r"$\sigma_{v_x}$"
        elif coord=='y' or coord=='vy':
            v=cluster.vy
            ylabel=r"$\sigma_{v_y}$"
        elif coord=='z' or coord=='vz':
            v=vz
            ylabel=r"$\sigma_{v_z}$"
        elif coord=='v':
            if projected:
                v=cluster.vpro
            else:
                v=cluster.v
    else:
        if projected:
            ylabel=r"$<v_{pro}>$"
        else:
            ylabel=r"$<v>$"

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))
    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)


    if normalize:
        if projected:
            rprofn=r_mean / cluster.rmpro
        else:
            rprofn=r_mean / cluster.rm
    else:
        rprofn=r_mean

    vprof=np.zeros(len(rprofn))

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if np.sum(rindx) > 3.0:

            vmean = np.mean(v[rindx])
            vprof[i]=vmean

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if normalize:
            xlabel=r"$r/r_m$"
        else:
            xlabel=r"$r$"

        _plot(
            rprofn,
            vprof,
            xlabel=xlabel,
            ylabel=ylabel,
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    if not normalize: rprofn=_convert_length(rprofn,units,cluster)
    vprof=_convert_velocity(vprof,units,cluster)

    return rprofn, vprof


def v2_prof(
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
    coord=None,
    normalize=False,
    plot=False,
    **kwargs,
):
    """Measure the radial variation in the mean squared velocity 

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
    coord : str
        choose what coordinate the mean velocity profile is to be returned in (default None returns (vx**2.+vy**2.+vz**2.)^1/2).
        Alternatively can ask for 'vr', 'vphi', or 'vtheta' for spherical coordinate velocity dispersions.
    normalize : bool
        normalize radial bins by cluster's half-mass radius (default: False)
    plot : bool 
        plot the velocity disperions profile (default: False)

    Returns
    -------
    rprofn : float
        natural log of radius (normalized by half-mass radius)
    vprof : float
        mean velocity

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

    vprof = np.array([])

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

    if coord is not None:

        if projected:
            r, phi, z = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
            vr, vp, vz = coords.rect_to_cyl_vec(
                cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
            )
        else:
            r, phi, theta, vr, vp, vt = sphere_coords(cluster)

        if coord =='r':
            v=vr
            ylabel=r"$<v_r^2>$"
        elif coord=='phi':
            v=vp
            ylabel=r"$<v_p^2>$"
        elif coord=='theta':
            v=vt
            ylabel=r"$<v_t^2>$"
        elif coord=='z':
            v=vz
            ylabel=r"$<v_z^2>$"
    else:
        if projected:
            ylabel=r"$<v_{pro}^2>$"
        else:
            ylabel=r"$<v^2>$"

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))
    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    if normalize:
        if projected:
            rprofn=r_mean / cluster.rmpro
        else:
            rprofn=r_mean / cluster.rm
    else:
        rprofn=r_mean

    vprof=np.zeros(len(rprofn))

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if np.sum(rindx) > 3.0:

            vmean = np.mean(v[rindx]**2.)
            vprof[i]=vmean

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if normalize:
            xlabel=r"$r/r_m$"
        else:
            xlabel=r"$r$"

        _plot(
            rprofn,
            vprof,
            xlabel=xlabel,
            ylabel=ylabel,
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    if not normalize: rprofn=_convert_length(rprofn,units,cluster)
    vprof=_convert_square_velocity(vprof,units,cluster)

    return rprofn, vprof

def eta_prof(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
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
    meq=False,
    **kwargs,
):
    """Measure the radial variation in eta

    Parameters
    ----------
    cluster : class
        StarCluster
    mmin/mmax : float
        minimum and maximum stellar mass
    nmass : int
        number of mass bins (default: 10)
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
        normalize radial bins by cluster's half-mass radius (default: True)
    plot : bool 
        plot the alpha profile (default: False)

    Returns
    -------
    rprofn : float
        natural log of each radius bin (normalized by half-mass radius)
    eprof : float
        slope of the sigma_v-mass function
    deta : float
        radial variation in eta calculated as deta = d(eta)/d(ln(r/rm)
    edeta : float
        error in deta
    ydeta : float
        y-intercept in fit to eta vs ln(r/rm)
    eydeta : float
        error in ydeta

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

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        r_hist=np.zeros(len(r_mean))
    elif kwargs.pop('bintype','num')=='fix':
        r_lower, r_mean, r_upper, r_hist = binmaker(cluster.r[indx], nrad)
    else:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r[indx], nrad)

    if normalize:
        if projected:
            rprofn=r_mean / cluster.rmpro
        else:
            rprofn=r_mean / cluster.rm
    else:
        rprofn=r_mean

    eprof=np.zeros(len(rprofn))

    for i in range(0, len(r_mean)):

        m_mean, sigvm, eta, eeta, yeta, eyeta = eta_function(
            cluster,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=r_lower[i],
            rmax=r_upper[i],
            vmin=vmin,
            vmax=vmax,
            kwmin=kwmin,
            kwmax=kwmax,
            projected=projected,
            meq=meq,
            **kwargs,
        )

        if eta > -100:
            eprof[i]=eta

    if len(rprofn) > 3:
        if projected:
            (deta, ydeta), V = np.polyfit(np.log(r_mean/cluster.rmpro), eprof, 1, cov=True)
        else:
            (deta, ydeta), V = np.polyfit(np.log(r_mean/cluster.rm), eprof, 1, cov=True)
        edeta = np.sqrt(V[0][0])
        eydeta = np.sqrt(V[1][1])
    else:
        deta = -100.0
        ydeta = 0.0
        edeta = 0.0
        eydeta = 0.0

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if normalize:
            xlabel=r"$\ln(r/r_m)$"
        else:
            xlabel=r"$\ln(r)$"

        if meq:
            ylabel=r"$m_{eq}$"
        else:
            ylabel=r"$\eta$"


        _plot(
            np.log(rprofn),
            eprof,
            xlabel=xlabel,
            ylabel=ylabel,
            overplot=overplot,
            **kwargs
        )

        if projected:
            rfit=r_mean/cluster.rmpro
        else:
            rfit=r_mean/cluster.rm

        efit = deta * np.log(rfit) + ydeta
        if meq:
            _lplot(np.log(rprofn), efit, overplot=True, label=(r"d$m_{eq}$ = %f" % deta))
        else:
            _lplot(np.log(rprofn), efit, overplot=True, label=(r"d$\eta$ = %f" % deta))
        plt.legend()

        if filename != None:
            plt.savefig(filename)

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if not normalize: rprofn=_convert_length(rprofn,units,cluster)

    return rprofn, eprof, deta, edeta, ydeta, eydeta

def meq_prof(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
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
    meq=False,
    **kwargs,
):
    """Measure the radial variation in meq

    Parameters
    ----------
    cluster : class
        StarCluster
    mmin/mmax : float
        minimum and maximum stellar mass
    nmass : int
        number of mass bins (default: 10)
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
        normalize radial bins by cluster's half-mass radius (default: True)
    plot : bool 
        plot the alpha profile (default: False)

    Returns
    -------
    rprofn : float
        natural log of each radius bin (normalized by half-mass radius)
    eprof : float
        slope of the sigma_v-mass function
    deta : float
        radial variation in eta calculated as deta = d(eta)/d(ln(r/rm)
    edeta : float
        error in deta
    ydeta : float
        y-intercept in fit to eta vs ln(r/rm)
    eydeta : float
        error in ydeta

    Other Parameters
    ----------------
    kwrags : str
        key word arguments for plotting

    History
    -------
    2020 - Written - Webb (UofT)
    """

    rprofn, meqprof, dmeq, edmeq, ydmeq, eymeq=eta_prof(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=nmass,
        rmin=rmin,
        rmax=rmax,
        nrad=nrad,
        vmin=vmin,
        vmax=vmax,
        emin=emin,
        emax=emax,
        kwmin=kwmin,
        kwmax=kwmax,
        npop=npop,
        indx=indx,
        bins=bins,
        projected=projected,
        normalize=normalize,
        plot=plot,
        meq=True,
        **kwargs,
    )

    return rprofn, meqprof, dmeq, edmeq, ydmeq, eymeq
 

def vcirc_prof(
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
    """
    NAME:

       vcirc_prof

    PURPOSE:

       Measure the circulr velocity profile of the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    mmin/mmax : float
        minimum and maximum stellar mass
    nmass : int
        number of mass bins (default: 10)
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
        plot the alpha profile (default: False)

    Returns
    -------
    rprof : float
        radius bins
    vprof : float
        circular velocity
    rvmax : float
        radius of maximum circular velocity
    vmax : float
        maximum circular velocity
    
    Other Parameters
    ----------------
    kwrags : str
        key word arguments for plotting

    History
    -------
    2019 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if cluster.units=='amuse':
        cluster.to_pckms()

    if cluster.origin0 != 'cluster' and cluster.origin0 != 'centre':
        cluster.to_centre(sortstars=normalize)
    else:
        cluster.sortstars()

    grav=_get_grav(cluster)

    if projected:
        r = cluster.rpro[cluster.rproorder]
        v = cluster.vpro[cluster.rproorder]
        m = cluster.m[cluster.rproorder]
    else:
        r = cluster.r[cluster.rorder]
        v = cluster.v[cluster.rorder]
        m = cluster.m[cluster.rorder]


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

    r = r[indx]
    v = v[indx]
    m = m[indx]

    msum = np.cumsum(m)
    vcirc = np.sqrt(grav * msum / r)
    vmax = np.amax(vcirc)
    rvmax = r[np.argmax(vcirc)]

    if bins is not None:
        r_lower, r_mean, r_upper=bins[0],bins[1],bins[2]
        args=np.zeros(len(r_mean))
        for i in range(0,len(args)):
            args[i]=np.argmin(np.fabs(r-r_mean[i]))
    elif kwargs.get('bintype','num')=='fix' and nrad is not None:
        r_lower, r_mean, r_upper, r_hist = binmaker(cluster.r[indx], nrad)
        args=np.zeros(len(r_mean))
        for i in range(0,len(args)):
            args[i]=np.argmin(np.fabs(r-r_mean[i]))
    elif kwargs.get('bintype','num')=='num' and nrad is not None:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r[indx], nrad)
        args=np.zeros(len(r_mean))
        for i in range(0,len(args)):
            args[i]=np.argmin(np.fabs(r-r_mean[i]))
    elif nrad is None:
        args=np.arange(0,len(r),1)


    args=args.astype(int)

    rprof = r[args]
    vcprof = vcirc[args]

    if normalize:
        if projected:
            rprof/=cluster.rmpro
        else:
            rprof/=cluster.rm

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms" or cluster.units == "pcmyr":
            xunits = " (pc)"
            yunits = " km/s"
        elif cluster.units == "kpckms" or cluster.units == "kpcgyr" or cluster.units== "WDunits":
            xunits = " (kpc)"
            yunits = " km/s"

        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"

        else:
            xunits = ""
            yunits = ""

        if normalize:
            xlabel=r"$\ln(r/r_m)$"
        else:
            xlabel=r"$\ln(r)$"

        x, y = rprof, vcprof
        _lplot(
            x,
            y,
            xlabel=r"$R \ %s$" % xunits,
            ylabel=r"$v_c \ %s $" % yunits,
            title="Time = %f" % cluster.tphys,
            log=kwargs.pop('log',True),
            overplot=overplot,
            filename=filename,
            **kwargs,
        )
        _lplot([rvmax, rvmax], [np.amin(y), np.amax(y)], "--", overplot=True)
        _lplot([np.amin(x), np.amax(x)], [vmax, vmax], "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    units=cluster.units
    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if not normalize: rprof=_convert_length(rprof,units,cluster)
    vcprof=_convert_velocity(vcprof,units,cluster)

    rvmax=_convert_length(rvmax,units,cluster)
    vmax=_convert_velocity(vmax,units,cluster)

    return rprof, vcprof, rvmax, vmax
