# Routines for analysing Nbody models as if they were Observations
import numpy as np
from ..main.operations import save_cluster, return_cluster
from ..util.recipes import *
from ..util.plots import *


def obs_mass_function(
    cluster,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    indx=None,
    mcorr=None,
    projected=True,
    omask=None,
    plot=False,
    **kwargs
):
    """
    NAME:

       obs_mass_function

    PURPOSE:

       Find mass function over a given mass range using nmass bins containing an equal number of stars.
       obs_mass_function differs from mass_function as it incorporates a mass completeness correction term and returns the completeness correction per bin

    INPUT:

       cluster - StarCluster instance
       mmin/mmax - specific mass range
       nmass - number of mass bins used to calculate alpha
       rmin/rmax - specific radial range
       vmin/vmax - specific velocity range
       emin/emax - specific energy range
       kwmin/kwmax - specific stellar evolution type range
       indx - specific subset of stars
       mcorr - correction for masses
       projected - use projected values
       omask - place a mask over the dataset and (WIP)
       plot - plot the mass function
       **kwargs - key words for plotting

    OUTPUT:

       m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    if mcorr is None:
        if omask is not None:
            try:
                mcorr = omask.mcorr
            except:
                mcorr = np.ones(cluster.ntot)
        else:
            mcorr = np.ones(cluster.ntot)

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
        mmin = np.amin(cluster.m)
    if mmax == None:
        mmax = np.amax(cluster.m)
    if kwmin == None:
        kwmin = np.amin(cluster.kw)
    if kwmax == None:
        kwmax = np.amax(cluster.kw)

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

    if np.sum(indx) >= nmass:

        if omask is None:
            m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)
        else:
            try:
                m_lower, m_mean, m_upper = omask.m_lower, omask.m_mean, omask.m_upper
                nmass = len(m_mean)
                m_hist = np.zeros(nmass)
            except:
                m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)

        m_corr_hist = np.zeros(len(m_hist))
        for i in range(0, len(m_hist)):
            mindx = (cluster.m >= m_lower[i]) * (cluster.m < m_upper[i]) * indx
            m_corr_hist[i] = np.sum(1.0 / mcorr[mindx])
            m_hist[i] = np.sum(mindx)

        mbincorr = m_hist / m_corr_hist

        lm_mean = np.log10(m_mean)
        dm = m_corr_hist / (m_upper - m_lower)
        ldm = np.log10(dm)

        (alpha, yalpha), V = np.polyfit(lm_mean, ldm, 1, cov=True)
        ealpha = np.sqrt(V[0][0])
        eyalpha = np.sqrt(V[1][1])

        if plot:
            filename = kwargs.get("filename", None)
            nplot(m_mean, np.log10(dm), xlabel="M", ylabel="LOG(dN/dM)", **kwargs)
            mfit = np.linspace(np.min(m_mean), np.max(m_mean), nmass)
            dmfit = 10.0 ** (alpha * np.log10(mfit) + yalpha)
            nlplot(
                mfit, np.log10(dmfit), overplot=True, label=(r"$\alpha$ = %f" % alpha)
            )

            plt.legend()

            if filename != None:
                plt.savefig(filename)

        return m_mean, m_corr_hist, dm, alpha, ealpha, yalpha, eyalpha, mbincorr
    else:
        print("NOT ENOUGH STARS TO ESTIMATE MASS FUNCTION")
        return (
            np.zeros(nmass),
            np.zeros(nmass),
            np.zeros(nmass),
            -1000.0,
            -1000.0,
            -1000.0,
            -1000.0,
            np.zeros(nmass),
        )


def obs_alpha_prof(
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
    mcorr=None,
    projected=True,
    omask=None,
    plot=False,
    **kwargs
):
    """
    NAME:

       alpha_prof

    PURPOSE:

       Measure the radial variation in the mass function

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       nmass - number of mass bins to calculate slope of mass function

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       mcorr - correction function for masses

       projected - use projected values and constraints (Default:False)

       omask - place a mask over the dataset and (WIP)

       plot - plot the density profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    OUTPUT:

        lrprofn - natural log of radius (normalized by half-mass radius)

        aprof - slope of the mass function

        dalpha - delta_alpha = d(alpha)/d(ln(r/rm) 

        edalpha - error in dalpha

        ydalpha,eydalpha - y-intercept and error in fit to alpha vs ln(r/rm)

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre(do_order=True, do_key_params=True)

    if mcorr is None:
        if omask is not None:
            try:
                mcorr = np.array(omask.mcorr)
            except:
                mcorr = np.ones(cluster.ntot)
        else:
            mcorr = np.ones(cluster.ntot)

    lrprofn = np.array([])
    aprof = np.array([])

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
    if kwmin == None:
        kwmin = np.min(cluster.kw)
    if kwmax == None:
        kwmax = np.max(cluster.kw)

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

    mbincorr=np.zeros(len(r_mean))


    for i in range(0, len(r_mean)):
        rindx = indx * (r >= r_lower[i]) * (r < r_upper[i])

        if omask is None:
            m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbincorr[i] = dx_corr_function(
                cluster.m[rindx], nmass, mcorr[rindx]
            )
        else:

            try:
                m_lower, m_mean, m_upper = omask.m_lower, omask.m_mean, omask.m_upper

                m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbincorr[i]= dx_corr_function(
                    cluster.m[rindx],
                    nmass,
                    mcorr[rindx],
                    x_lower=m_lower,
                    x_mean=m_mean,
                    x_upper=m_upper
                )
            except:
                m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbincorr[i]= dx_corr_function(
                    cluster.m[rindx],
                    nmass,
                    mcorr[rindx]
                )

        if projected:
            lrprofn=np.append(lrprofn,np.log(r_mean[i] / cluster.rmpro))
        else:
            lrprofn=np.append(lrprofn,np.log(r_mean[i] / cluster.rm))
        if alpha > -100:
            aprof=np.append(aprof,alpha)
        else:
            aprof=np.append(aprof,-100.)


    aindx=(aprof>-100.)

    if np.sum(aindx) > 3:
        (dalpha, ydalpha), V = np.polyfit(lrprofn[aindx], aprof[aindx], 1, cov=True)
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

    return lrprofn, aprof, dalpha, edalpha, ydalpha, eydalpha, mbincorr


def dx_corr_function(
    x, 
    nx=10, 
    xcorr=None, 
    bintype="num", 
    x_lower=None, 
    x_mean=None, 
    x_upper=None
):
    """
  NAME:

     dx_function

  PURPOSE:

     Find distribution function using nx bins containing an equal number of points
  INPUT:

     x - input array

     nx - number of bins

     xcorr - correction function for x (default:None)

     bintype - bin with equal number of stars per bin (bin) or evenly in x (fix) (default: num)

     omask - place a mask over the dataset and (WIP)

  OUTPUT:

     x_mean,x_hist,dx,alpha,ealpha,yalpha,eyalpha


  HISTORY:

     2018 - Written - Webb (UofT)

  """

    if x_lower is None:
        if bintype == "num":
            x_lower, x_mean, x_upper, x_hist = nbinmaker(x, nx)
        else:
            x_lower, x_mean, x_upper, x_hist = binmaker(x, nx)
    else:
        nx = len(x_lower)

    if xcorr is not None:
        x_hist_corr = np.zeros(nx)
        x_hist = np.zeros(nx)

        for i in range(0, len(x_hist_corr)):
            indx = (x >= x_lower[i]) * (x < x_upper[i])
            x_hist_corr[i] = np.sum(1.0 / xcorr[indx])
            x_hist[i]=np.sum(indx)
    else:
        x_hist_corr=x_hist

    xbincorr=np.amin(x_hist/x_hist_corr)

    lx_mean = np.log10(x_mean)
    dx = x_hist_corr / (x_upper - x_lower)
    ldx = np.log10(dx)

    (alpha, yalpha), V = np.polyfit(lx_mean, ldx, 1, cov=True)
    ealpha = np.sqrt(V[0][0])
    eyalpha = np.sqrt(V[1][1])

    return x_mean, x_hist, dx, alpha, ealpha, yalpha, eyalpha, xbincorr
