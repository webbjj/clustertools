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
    "eta_prof",
    "vcirc_prof"
]

import numpy as np
from galpy.util import bovy_coords

from ..util.constants import *
from ..util.recipes import *
from .operations import *
from ..util.plots import *
from ..util.coordinates import sphere_coords
from .functions import mass_function


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
    indx=None,
    projected=False,
    plot=False,
    **kwargs
):
    """
    NAME:

       rho_prof

    PURPOSE:

       Measure the density profile of the cluster

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       plot - plot the density profile (Default: False)

    KWARGS:

        Same as for ..util.plot.nplot

    Returns

        rprof,pprof,nprof (radius, density, number of stars)

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    rprof = np.array([])
    pprof = np.array([])
    nprof = np.array([])

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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
        * (cluster.kw >= kwmin)
        * (cluster.kw <= kwmax)
    )

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax

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
        elif cluster.units == "pckms":
            xunits = " (pc)"
            if projected:
                yunits = " Msun/pc^2"
            else:
                yunits = " Msun/pc^3"
        elif cluster.units == "kpckms":
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
            xlabel=r"$R %s$" % xunits,
            ylabel=r"$\rho %s$" % yunits,
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

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
    indx=None,
    projected=False,
    cumulative=False,
    plot=False,
    **kwargs
):
    """
    NAME:

       m_prof

    PURPOSE:

       Measure the mass profile of the cluster

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       cumalitive - determine the cumulative mass profile instead (Default: False)

       plot - plot the mean mass profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    Returns

        rprof,mprof,nprof (radius, mass, number of stars)

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    rprof = []
    mprof = []
    nprof = []

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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
        * (cluster.kw >= kwmin)
        * (cluster.kw <= kwmax)
    )

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax

    r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    for i in range(0, len(r_mean)):
        if cumulative:
            rindx = indx * (r < r_upper[i])
        else:
            rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])
        rprof.append(r_mean[i])

        mprof.append(np.sum(cluster.m[rindx]))
        nprof.append(np.sum(rindx))

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms":
            xunits = " (pc)"
            yunits = " Msun"
        elif cluster.units == "kpckms":
            xunits = " (kpc)"
            yunits = " Msun"
        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"
        else:
            xunits = ""
            yunits = ""

        x, y, n = rprof, mprof, nprof
        nlplot(
            x,
            y,
            xlabel=r"$R %s $" % xunits,
            ylabel=r"$M %s $" % yunits,
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

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
    indx=None,
    projected=False,
    mcorr=None,
    omask=None,
    plot=False,
    **kwargs
):
    """
    NAME:

       alpha_prof

    PURPOSE:

       Measure the radial variation in the mass function

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       nmass - number of mass bins to calculate slope of mass function

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       mcorr - correction for masses
   
       omask - place a mask over the dataset to mimic observed data

       plot - plot the alpha profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    Returns

        lrprofn - natural log of radius (normalized by half-mass radius)

        aprof - slope of the mass function

        dalpha - delta_alpha = d(alpha)/d(ln(r/rm) 

        edalpha - error in dalpha

        ydalpha,eydalpha - y-intercept and error in fit to alpha vs ln(r/rm)

        rbinerror - if mcorr is not None, output lowest corrected mass fraction at each radius

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    if mcorr is None:
        if omask is not None:
            try:
                mcorr = omask.mcorr
                return_error=True
            except:
                mcorr = np.ones(cluster.ntot)
                return_error=False
        else:
            mcorr = np.ones(cluster.ntot)
            return_error=False
    else:
        return_error=True

    lrprofn = []
    aprof = []

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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
        * (cluster.kw >= kwmin)
        * (cluster.kw <= kwmax)
    )

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax

    if omask is None:
        r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)
    else:
        try:
            r_lower, r_mean, r_upper = omask.r_lower, omask.r_mean, omask.r_upper
        except:
            r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    rbinerror=np.zeros(len(r_mean))

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if return_error:
            m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbinerror = mass_function(cluster,nmass=nmass,indx=rindx,projected=projected,mcorr=mcorr,omask=omask,plot=False,**kwargs)
            rbinerror[i]=np.amin(mbinerror)
        else:
            m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(cluster,nmass=nmass,indx=rindx,projected=projected,mcorr=mcorr,omask=omask,plot=False,**kwargs)
            rbinerror[i]=1.

        if alpha > -100:
            if projected:
                lrprofn.append(np.log(r_mean[i] / cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i] / cluster.rm))

            aprof.append(alpha)

    if len(lrprofn) > 3:
        (dalpha, ydalpha), V = np.polyfit(lrprofn, aprof, 1, cov=True)
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

        nplot(
            lrprofn,
            aprof,
            xlabel=r"$\ln(r/r_m)$",
            ylabel=r"$\alpha$",
            overplot=overplot,
            **kwargs
        )
        rfit = np.linspace(np.min(lrprofn), np.max(lrprofn), nrad)
        afit = dalpha * rfit + ydalpha
        nlplot(rfit, afit, overplot=True, label=(r"d$\alpha$ = %f" % dalpha))
        plt.legend()

        if filename != None:
            plt.savefig(filename)

    cluster.dalpha = dalpha

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    if return_error:
        return lrprofn, aprof, dalpha, edalpha, ydalpha, eydalpha, rbinerror
    else:
        return lrprofn, aprof, dalpha, edalpha, ydalpha, eydalpha


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
    kwmin=None,
    kwmax=None,
    indx=None,
    projected=False,
    rcustom=None,
    coord=None,
    normalize=True,
    plot=False,
):
    """
    NAME:

       sigv_prof

    PURPOSE:

       Measure the radial variation in the velocity dispersion

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       rcustom - use custom radius to bin data (not well implemented)

       coord - choose what coordinate the velocity dispersion profile is to be returned in (default None returns (sigx**2.+sigy**2.+sigz**2.)^1/2)

       normalize - normalize radial bins by cluster's half-mass radius (default: True)

       plot - plot the velocity dispersion profile (Default: False)


    KWARGS:

       Same as for ..util.plot.nplot

    Returns

        lrprofn - natural log of radius (normalized by half-mass radius)

        sigvprof - velocity dispersion

        betaprof - anisotropy parameter 

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    lrprofn = []
    sigvprof = []

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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

    if kwmin is not None:
        indx*=(cluster.kw >= kwmin)
    if kwmax is not None:
        indx*=(cluster.kw <= kwmax)

    if emin is not None:
        indx *= cluster.etot >= emin
    if emin is not None:
        indx *= cluster.etot <= emax

    # Convert to cylindrical or spherical coordinates:
    if projected:
        r, theta, z = bovy_coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
        vr, vtheta, vz = bovy_coords.rect_to_cyl_vec(
            cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
        )
    else:
        r, phi, theta, vr, vp, vt = sphere_coords(cluster)

    if rcustom is not None:
        r=rcustom

    r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if np.sum(rindx) > 3.0:

            sigr = np.std(vr[rindx])
            sigt = np.std(vt[rindx])

            if projected:
                sigp = np.zeros(len(vr))
            else:
                sigp = np.std(vp[rindx])

            if coord is None:
                sigv = np.sqrt(sigr ** 2.0 + sigt ** 2.0 + sigp ** 2.0)
            elif coord=='r':
                sigv=sigr
            elif coord=='phi':
                sigv=sigp
            elif coord=='theta':
                sigv=sigt 

            if normalize:
                if projected:
                    lrprofn.append(np.log(r_mean[i] / cluster.rmpro))
                else:
                    lrprofn.append(np.log(r_mean[i] / cluster.rm))
            else:
                lrprofn.append(np.log(r_mean[i]))

            sigvprof.append(sigv)

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        nplot(
            lrprofn,
            sigvprof,
            xlabel=r"$\ln(r/r_m)$",
            ylabel=r"$\sigma_v$",
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    return lrprofn, sigvprof

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
    kwmin=None,
    kwmax=None,
    indx=None,
    projected=False,
    rcustom=None,
    coord=None,
    normalize=True,
    plot=False,
):
    """
    NAME:

       beta_prof

    PURPOSE:

       Measure the anisotropy profile of the cluster

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       rcustom - use custom radius to bin data (not well implemented)

       coord - choose what coordinate the velocity dispersion profile is to be returned in (default None returns (sigx**2.+sigy**2.+sigz**2.)^1/2)

       normalize - normalize radial bins by cluster's half-mass radius (default: True)

       plot - plot the orbital anisotropy profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    Returns

        lrprofn - natural log of radius (normalized by half-mass radius)

        betaprof - anisotropy parameter 

    HISTORY:

       2020 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    lrprofn = []
    betaprof = []

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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

    if kwmin is not None:
        indx*=(cluster.kw >= kwmin)
    if kwmax is not None:
        indx*=(cluster.kw <= kwmax)

    if emin is not None:
        indx *= cluster.etot >= emin
    if emin is not None:
        indx *= cluster.etot <= emax

    # Convert to cylindrical or spherical coordinates:
    if projected:
        r, theta, z = bovy_coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
        vr, vtheta, vz = bovy_coords.rect_to_cyl_vec(
            cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
        )
    else:
        r, phi, theta, vr, vp, vt = sphere_coords(cluster)

    if rcustom is not None:
        r=rcustom

    r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if np.sum(rindx) > 3.0:

            sigr = np.std(vr[rindx])
            sigt = np.std(vt[rindx])

            if projected:
                sigp = np.zeros(len(vr))
                beta = sigt / sigr - 1.0
            else:
                sigp = np.std(vp[rindx])
                beta = 1.0 - (sigt ** 2.0 + sigp ** 2.0) / (2.0 * (sigr ** 2.0))

            if normalize:
                if projected:
                    lrprofn.append(np.log(r_mean[i] / cluster.rmpro))
                else:
                    lrprofn.append(np.log(r_mean[i] / cluster.rm))
            else:
                lrprofn.append(np.log(r_mean[i]))

            betaprof.append(beta)

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        nplot(
            lrprofn,
            sigvprof,
            xlabel=r"$\ln(r/r_m)$",
            ylabel=r"$\beta$",
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    return lrprofn, betaprof


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
    indx=None,
    projected=False,
    plot=False,
):
    """
    NAME:

       v_prof

    PURPOSE:

       Measure the radial variation in the mean velocity 

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       plot - plot the mean velocity profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    Returns

        lrprofn - natural log of radius (normalized by half-mass radius)

        vprof - mean velocity

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    lrprofn = []
    sigvprof = []

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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
        * (cluster.kw >= kwmin)
        * (cluster.kw <= kwmax)
    )

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax

    # Convert to cylindrical or spherical coordinates:
    if projected:
        r, theta, z = bovy_coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
        vr, vtheta, vz = bovy_coords.rect_to_cyl_vec(
            cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
        )
    else:
        r, phi, theta, vr, vp, vt = sphere_coords(cluster)

    r_lower, r_mean, r_upper, r_hist = nbinmaker(r[indx], nrad)

    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if np.sum(rindx) > 3.0:

            vrmean = np.mean(vr[rindx])
            vtmean = np.mean(vt[rindx])

            if projected:
                vpmean = np.zeros(len(vr))
            else:
                vpmean = np.mean(vp[rindx])

            vmean = np.sqrt(vrmean ** 2.0 + vtmean ** 2.0 + vpmean ** 2.0)

            if projected:
                lrprofn.append(np.log(r_mean[i] / cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i] / cluster.rm))

            vprof.append(vmean)

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        nplot(
            lrprofn,
            vprof,
            xlabel=r"$\ln(r/r_m)$",
            ylabel=r"$<v>$",
            overplot=overplot,
            **kwargs
        )

        if filename != None:
            plt.savefig(filename)

    return lrprofn, vprof


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
    kwmax=1,
    indx=None,
    projected=False,
):
    """
    NAME:

       eta_prof

    PURPOSE:

       Measure the radial variation in eta

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       nmass - number of mass bins over which to measure eta

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       plot - plot the eta profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    Returns

        lrprofn - natural log of radius (normalized by half-mass radius)

        eprof - slope of the sigma_v-mass function

        deta - deta = d(eta)/d(ln(r/rm) 

        edeta - error in deta

        ydeta,eydeta - y-intercept and error in fit to eta vs ln(r/rm)

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    lrprofn = []
    eprof = []

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

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
        * (cluster.kw >= kwmin)
        * (cluster.kw <= kwmax)
    )

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax

    r_lower, r_mean, r_upper, r_hist = nbinmaker(cluster.r[indx], nrad)

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
        )

        if alpha > -100:
            if projected:
                lrprofn.append(np.log(r_mean[i] / cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i] / cluster.rm))

            eprof.append(eta)

    if len(lrprofn) > 3:
        (deta, ydeta), V = np.polyfit(lrprofn, eprof, 1, cov=True)
        edeta = np.sqrt(V[0][0])
        eydeta = np.sqrt(V[1][1])
    else:
        deta = -100.0
        ydeta = 0.0
        edeta = 0.0
        eydeta = 0.0

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    return lrprofn, eprof, deta, edeta, ydeta, eydeta


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
    indx=None,
    projected=False,
    plot=False,
    **kwargs
):
    """
    NAME:

       vcirc_prof

    PURPOSE:

       Measure the circulr velocity profile of the cluster

    Parameters

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       plot - plot the circular velocity profile (Default: False)

    KWARGS:

        Same as for ..util.plot.nplot

    Returns

        rprof,vprof,rvmax,vmax (radius, circular velocity, radius of maximum virc, maximum vcirc)

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

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

    rprof = np.array([])
    vcprof = np.array([])

    if projected:
        if cluster.rproorder is None:
            cluster.key_params(do_order=True)
        r = cluster.rpro[cluster.rproorder]
        v = cluster.vpro[cluster.rproorder]
        m = cluster.m[cluster.rproorder]
    else:
        if cluster.rorder is None:
            cluster.key_params(do_order=True)
        r = cluster.r[cluster.rorder]
        v = cluster.v[cluster.rorder]
        m = cluster.m[cluster.rorder]

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
        * (cluster.kw >= kwmin)
        * (cluster.kw <= kwmax)
    )

    if emin != None:
        indx *= cluster.etot >= emin
    if emin != None:
        indx *= cluster.etot <= emax

    r = r[indx]
    v = v[indx]
    m = m[indx]

    msum = np.cumsum(m)
    vcirc = np.sqrt(grav * msum / r)
    vmax = np.amax(vcirc)
    rvmax = r[np.argmax(vcirc)]

    rprof = r
    vcprof = vcirc

    if plot:
        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        if cluster.units == "nbody":
            xunits = " (NBODY)"
            yunits = " (NBODY)"
        elif cluster.units == "pckms":
            xunits = " (pc)"
            yunits = " km/s"
        elif cluster.units == "kpckms":
            xunits = " (kpc)"
            yunits = " km/s"

        elif cluster.units == "galpy":
            xunits = " (GALPY)"
            yunits = " (GALPY)"

        else:
            xunits = ""
            yunits = ""

        x, y = rprof, vcprof
        nlplot(
            x,
            y,
            xlabel=r"$R %s$" % xunits,
            ylabel=r"$vc %s $" % yunits,
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )
        nlplot([rvmax, rvmax], [np.amin(y), np.amax(y)], "--", overplot=True)
        nlplot([np.amin(x), np.amax(x)], [vmax, vmax], "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0, do_order=True, do_key_params=True)

    return rprof, vcprof, rvmax, vmax
