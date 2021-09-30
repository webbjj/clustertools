""" A matplotlib.pyplot wrapper for making key StarCluster plots

"""
__author__ = "Jeremy J Webb"
__all__ = [
  "starplot",
  "skyplot",
]

import matplotlib.pyplot as plt
import numpy as np
try:
    import galpy.util.plot as gplot
except:
    import galpy.util.bovy_plot as gplot
import seaborn as sns
import os
from scipy.ndimage import gaussian_filter
import matplotlib.colors as colors

try:
    gplot.start_print(axes_labelsize=18.0, xtick_labelsize=14.0, ytick_labelsize=14.0)
except:
    gplot.bovy_print(axes_labelsize=18.0, xtick_labelsize=14.0, ytick_labelsize=14.0)

current_palette = sns.color_palette()

def _scatter(
    x,
    y,
    z=None,
    xlabel="",
    ylabel="",
    legend=False,
    title="",
    xlim=None,
    ylim=None,
    scale=True,
    filename=None,
    overplot=False,
    **kwargs
):
    """Wrapper for matplotlib.pyplot.scatter 
      - allows for most pyplot commands to be assigned when calling the function

    Parameters
    ----------
    x,y : float
      points to be plotted
    z : float
      value for colour-coding
    xlabel : str
      X-axis label
    ylabel : str
      Y-axis label
    legend : bool
      Plot legend (default: False)
    title : str
      Title of plot
    xlim : float 
      X-axis limits (list)
    ylim : float list
      Y-axis limits (list)
    scale : bool
      if only xlim or ylim is given, set limits of axis to match points in given range (default: True)
    filename : str
      filename to save figure to
    overplot : bool
      overplot onto existing figure (default: False)

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.scatter

    History
    -------
    2018 - Written - Webb (UofT)

    """
    x = np.asarray(x)
    y = np.asarray(y)

    if not overplot:
        plt.figure()

    alpha = kwargs.pop("alpha", 1.0)

    if z is None:
        out = plt.scatter(x, y, alpha=alpha, **kwargs)
    else:
        z = np.asarray(z)
        out = plt.scatter(x, y, c=z, alpha=alpha, **kwargs)
        plt.colorbar(out)

    if overplot:
        return out
    else:
        if xlim != None:
            plt.xlim(xlim)
            if ylim == None and scale:
                xindx = (x >= xlim[0]) * (x <= xlim[1])
                ymin = np.min(y[xindx])
                ymax = np.max(y[xindx])
                plt.ylim((ymin, ymax))

        if ylim != None:
            plt.ylim(ylim)
            if xlim == None and scale:
                yindx = (y >= ylim[0]) * (y <= ylim[1])
                xmin = np.min(x[yindx])
                xmax = np.max(x[yindx])
                plt.xlim((xmin, xmax))

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        if legend:
            plt.legend()

        plt.tight_layout()

        if filename != None:
            plt.savefig(filename)

        return out


def _plot(
    x,
    y,
    ptype=".",
    xlabel="",
    ylabel="",
    legend=False,
    title="",
    xlim=None,
    ylim=None,
    scale=True,
    log=False,
    logx=False,
    logy=False,
    filename=None,
    overplot=False,
    **kwargs
):
    """Wrapper for plotting points with matplotlib.pyplot.plot
      - allows for most pyplot commands to be assigned when calling the function

    Parameters
    ----------
    x,y : float
      points to be plotted
    ptype : str
      point type
    xlabel : str
      X-axis label
    ylabel : str
      Y-axis label
    legend : bool
      Plot legend (default: False)
    title : str
      Title of plot
    xlim : float 
      X-axis limits (list)
    ylim : float list
      Y-axis limits (list)
    scale : bool
      if only xlim or ylim is given, set limits of axis to match points in given range (default: True)
    log : bool
      use log axis (default: False)
    logx : bool
      use log x-axis (default: False)
    logy : bool
      use log y-axis (default: False)
    filename : str
      filename to save figure to
    overplot : bool
      overplot onto existing figure (default: False)

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.plot

    History
    -------
    2018 - Written - Webb (UofT)
    """
    if not overplot:
        plt.figure()

    alpha = kwargs.pop("alpha", 1.0)

    if log or (logx and logy):
        out = plt.loglog(x, y, ptype, alpha=alpha, **kwargs)
    elif logx:
        out = plt.semilogx(x, y, ptype, alpha=alpha, **kwargs)
    elif logy:
        out = plt.semilogy(x, y, ptype, alpha=alpha, **kwargs)
    else:
        out = plt.plot(x, y, ptype, alpha=alpha, **kwargs)

    if overplot:
        return out
    else:
        if xlim != None:
            plt.xlim(xlim)
            if ylim == None and scale:
                xindx = (x >= xlim[0]) * (x <= xlim[1])
                ymin = np.amin(y[xindx])
                ymax = np.amax(y[xindx])
                plt.ylim((ymin, ymax))

        if ylim != None:
            plt.ylim(ylim)
            if xlim == None and scale:
                yindx = (y >= ylim[0]) * (y <= ylim[1])
                xmin = np.amin(x[yindx])
                xmax = np.amax(x[yindx])
                plt.xlim((xmin, xmax))

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        if legend:
            plt.legend()

        plt.tight_layout()

        if filename != None:
            plt.savefig(filename)

        return out


def _lplot(
    x,
    y,
    ltype="-",
    xlabel="",
    ylabel="",
    legend=False,
    title="",
    xlim=None,
    ylim=None,
    scale=True,
    log=False,
    logx=False,
    logy=False,
    filename=None,
    overplot=False,
    **kwargs
):
    """ Wrapper for plotting lines with matplotlib.pyplot.plot 
      - allows for most pyplot commands to be assigned when calling the function

    Parameters
    ----------
    x,y : float
      points to be plotted
    ltype : str
      line type
    xlabel : str
      X-axis label
    ylabel : str
      Y-axis label
    legend : bool
      Plot legend (default: False)
    title : str
      Title of plot
    xlim : float 
      X-axis limits (list)
    ylim : float list
      Y-axis limits (list)
    scale : bool
      if only xlim or ylim is given, set limits of axis to match points in given range (default: True)
    log : bool
      use log axis (default: False)
    logx : bool
      use log x-axis (default: False)
    logy : bool
      use log y-axis (default: False)
    filename : str
      filename to save figure to
    overplot : bool
      overplot onto existing figure (default: False)

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.plot

    History
    -------
    2018 
    """
    return _plot(
        x,
        y,
        ptype=ltype,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=legend,
        title=title,
        xlim=xlim,
        ylim=ylim,
        scale=scale,
        log=log,
        logx=logx,
        logy=logy,
        filename=filename,
        overplot=overplot,
        **kwargs
    )


def _hist(
    x,
    nbin=10,
    xlabel="",
    legend=False,
    title="",
    xlim=None,
    fill=False,
    filename=None,
    overplot=False,
    **kwargs
):
    """Wrapper for matplotlib.pyplot.hist 
      - allows for most pyplot commands to be assigned when calling the function

    Parameters
    ----------
    x : float
      points to make histogram with
    nbin : int
      number of bins
    xlabel : str
      X-axis label
    legend : bool
      Plot legend (default: False)
    title : str
      Title of plot
    xlim : float 
      X-axis limits (list)
    fill : bool
      Fill histogram bars (default: False)
    filename : str
      filename to save figure to
    overplot : bool
      overplot onto existing figure (default: False)

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.plot

    History
    -------
    2018 - Written - Webb (UofT)
    """
    nbin = kwargs.pop("bins", nbin)

    if not overplot:
        plt.figure()

    if fill:
        out = plt.hist(x, bins=nbin, **kwargs)
    else:
        histtype = kwargs.pop("histtype", "step")
        out = plt.hist(x, bins=nbin, histtype=histtype, **kwargs)

    if overplot:
        return out
    else:
        if xlim != None:
            plt.xlim(xlim)

        plt.xlabel(xlabel)
        plt.ylabel("N")

        plt.title(title)

        if legend:
            plt.legend()

        if filename != None:
            plt.savefig(filename)

        return out


def _hist2d(
    x,
    y,
    nbin=10,
    xlabel="",
    ylabel="",
    clabel="",
    title="",
    xlim=None,
    ylim=None,
    filename=None,
    **kwargs
):
    """Wrapper for matplotlib.pyplot._hist2d
      - allows for most pyplot commands to be assigned when calling the function

    Parameters
    ----------
    x,y : float
      points to make histogram with
    nbin : int
      number of bins
    xlabel : str
      X-axis label
    ylabel : str
      Y-axis label
    legend : bool
      Plot legend (default: False)
    title : str
      Title of plot
    xlim : float 
      X-axis limits (list)
    ylim : float 
      Y-axis limits (list)
    filename : str
      filename to save figure to
    clabel : str
       label for colorbar

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.hist2d

    History
    -------
    2018 - Written - Webb (UofT)
    """
    if xlim == None:
        xlim = [np.min(x), np.max(x)]
    if ylim == None:
        ylim = [np.min(y), np.max(y)]

    cmap = kwargs.pop("cmap", plt.cm.jet)

    out = plt.hist2d(x, y, bins=nbin, range=[xlim, ylim], cmap=cmap, **kwargs)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.title(title)
    plt.colorbar(label=clabel)

    if filename != None:
        plt.savefig(filename)

    return out


def _dens(
    x,
    y,
    z=None,
    normed=False,
    xstep=1.0,
    ystep=1.0,
    sigx=4.0,
    sigy=4.0,
    order=0,
    extent=None,
    interpolation="nearest",
    cmap="Spectral",
):
    """Wrapper for matplotlib.pyplot.histogram2d 
      - allows for most pyplot commands to be assigned when calling the function

    Parameters
    ----------
    x,y : float
      points to make histogram with
    z : float
      value for colour-coding
    normed : bool
      normalize (default: False)
    xstep : float
      step size between xbins
    ystep : float
      step size between ybins
    sigx/sigy : float
      dispersion for Gaussian filter
    order : float
      order of Gaussian filter (default: 0)
    extent : float
      axis range for imshow (default: None)
    interpolation : str
      interpolation type for imshow (default: nearest)
    cmap : str
      colour map for imshow (default: Spectral)

    Returns
    -------
    pyplot.figure()

    History
    -------
    2018 - Written - Webb (UofT)
    """
    xbins = np.arange(x.min(), x.max() + xstep, xstep)
    ybins = np.arange(y.min(), y.max() + ystep, ystep)

    if z is None:
        H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins])
    else:
        H, xedges, yedges = np.histogram2d(
            x, y, bins=[xbins, ybins], normed=normed, weights=z
        )

    fimg = gaussian_filter(H, sigma=(sigx, sigy), order=order)

    if extent is None:
        extent = [x.min(), x.max(), y.min(), y.max()]

    out = plt.imshow(
        fimg,
        norm=colors.LogNorm(),
        interpolation=interpolation,
        cmap=cmap,
        extent=extent,
    )
    plt.colorbar()

    return out


def starplot(
    cluster,
    coords="xyz",
    xlim=None,
    ylim=None,
    legend=False,
    filename=None,
    overplot=False,
    do_centre=False,
    **kwargs
):
    """Plot the xy/xz/yz coordinates of each star in the cluster

    Parameters
    ----------
    cluster : class
      StarCluster
    coords : str
      coordinates to be plotted ('xy','xz','yz','xyz') (default: 'xyz' returns a two panel plot)
    xlim : float 
      X-axis limits (list)
    ylim : float 
      Y-axis limits (list)
    legend : bool
      Plot legend (default: False)
    filename : str
      filename to save figure to
    overplot : bool
      overplot onto existing figure (default: False)
    do_centre : bool
      plot centre of cluster (default:False)

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.plot

    History
    -------
    2018 - Written - Webb (UofT)
    """

    alpha = kwargs.pop("alpha", 0.1)

    if coords != "xyz":

        if not overplot:
            plt.figure()

        if coords == "xz":
            x = cluster.x
            y = cluster.z
            if do_centre:
                xgc = cluster.xgc
                ygc = cluster.zgc
        elif coords == "yz":
            x = cluster.y
            y = cluster.z
            if do_centre:
                xgc = cluster.ygc
                ygc = cluster.zgc
        else:
            x = cluster.x
            y = cluster.y
            if do_centre:
                xgc = cluster.xgc
                ygc = cluster.ygc

        out = plt.plot(x, y, ".", alpha=alpha, **kwargs)

        if overplot:
            return out
        else:
            if do_centre:
                plt.plot(xgc, ygc, ".", alpha=1.0, label="COM", **kwargs)
            if cluster.units == "nbody":
                units = "(NBODY)"
            elif cluster.units == "pckms":
                units = "(pc)"
            elif cluster.units == "kpckms":
                units = "(kpc)"
            elif cluster.units == "galpy":
                units = "(GALPY)"
            elif cluster.units == 'radec':
                units = "(degree)"
            else:
                units = ""

            if coords == "xz":
                x = cluster.x
                y = cluster.z
                plt.xlabel("X " + units)
                plt.ylabel("Z " + units)
            elif coords == "yz":
                x = cluster.y
                y = cluster.z
                plt.xlabel("Y " + units)
                plt.ylabel("Z " + units)
            else:
                x = cluster.x
                y = cluster.y
                plt.xlabel("X " + units)
                plt.ylabel("Y " + units)

            plt.title("Time = %f" % cluster.tphys)

            if xlim != None:
                plt.xlim(xlim)
            if ylim != None:
                plt.ylim(ylim)

            plt.tight_layout()

            if filename != None:
                plt.savefig(filename)

            return out

    else:

        if cluster.units == "nbody":
            units = "(NBODY)"
        elif cluster.units == "pckms":
            units = "(pc)"
        elif cluster.units == "kpckms":
            units = "(kpc)"
        elif cluster.units == "galpy":
            units = "(GALPY)"
        elif cluster.units == 'radec':
            units = "(degree)"
        else:
            units = ""

        plt.subplot(1, 2, 1)

        plt.plot(cluster.x, cluster.z, ".", alpha=alpha, **kwargs)

        if cluster.origin == "galaxy" and do_centre:
            plt.plot(cluster.xgc, cluster.zgc, "o", label="Center")
            plt.plot(
                cluster.xgc + cluster.xc, cluster.zgc + cluster.zc, "o", label="COM"
            )
        elif cluster.origin != "centre" and do_centre:
            plt.plot(cluster.xc, cluster.zc, "o", label="COM")

        if xlim != None:
            plt.xlim(xlim)
        if ylim != None:
            plt.ylim(ylim)

        plt.xlabel("X " + units)
        plt.ylabel("Z " + units)

        plt.subplot(1, 2, 2)

        plt.plot(cluster.x, cluster.y, ".", alpha=alpha, **kwargs)

        if cluster.origin == "galaxy" and do_centre:
            plt.plot(cluster.xgc, cluster.ygc, "o", label="Center")
            plt.plot(
                cluster.xgc + cluster.xc, cluster.ygc + cluster.yc, "o", label="COM"
            )
        elif cluster.origin != "centre" and do_centre:
            plt.plot(cluster.xc, cluster.yc, "o", label="COM")

        if xlim != None:
            plt.xlim(xlim)
        if ylim != None:
            plt.ylim(ylim)

        plt.xlabel("X " + units)
        plt.ylabel("Y " + units)

        if legend:
            plt.legend()

        plt.tight_layout()

        if filename != None:
            plt.savefig(filename)
            return 0
        else:
            return 0

    return 0


def skyplot(
    cluster,
    coords="radec",
    xlim=None,
    ylim=None,
    legend=False,
    filename=None,
    overplot=False,
    do_centre=False,
    **kwargs
):
    """Plot the ra/dec or pmra/pmdec coordinates of each star in the cluster

    Parameters
    ----------
    cluster : class
      StarCluster
    coords : str
      coordinates to be plotted ('radec','pm') (default: 'radec')
    xlim : float 
      X-axis limits (list)
    ylim : float 
      Y-axis limits (list)
    legend : bool
      Plot legend (default: False)
    filename : str
      filename to save figure to
    overplot : bool
      overplot onto existing figure (default: False)
    do_centre : bool
      plot centre of cluster (default:False)

    Returns
    -------
    pyplot.figure()

    Other Parameters
    ----------------
    kwargs : str
      key word arguments for matplotlib.pyplot.plot

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    cluster.to_radec(sortstars=False)

    alpha = kwargs.pop("alpha", 0.1)

    if not overplot:
        plt.figure()

    if coords == "radec":
        x = cluster.ra
        y = cluster.dec
    elif coords == "pm":
        x = cluster.pmra
        y = cluster.pmdec

    out = plt.plot(x, y, ".", alpha=alpha, **kwargs)

    if overplot:
        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        return out
    else:
        if do_centre:
            xgc, ygc = cluster.ragc, cluster.decgc
            plt.plot(xgc, ygc, ".", alpha=1.0, label="COM", **kwargs)

        cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        if coords == "radec":
            plt.xlabel("Ra (degree)")
            plt.ylabel("Dec (degree)")
        elif coords == "pm":
            plt.xlabel(r"$\mu_{Ra}$ (mas $yr^{-1}$)")
            plt.ylabel(r"$\mu_{Dec}$ (mas $yr^{-1}$)")

        plt.title("Time = %f" % cluster.tphys)

        if xlim != None:
            plt.xlim(xlim)
        if ylim != None:
            plt.ylim(ylim)

        plt.tight_layout()

        if filename != None:
            plt.savefig(filename)

        return out
