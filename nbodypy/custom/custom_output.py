# Custum output files

import numpy as np

import matplotlib.pyplot as plt
from galpy.util import bovy_plot
import os

from ..main.cluster import sub_cluster
from ..util.coordinates import sky_coords
from ..main.functions import *
from ..main.profiles import *
from ..main.operations import *
from ..observations.observations import *


def p_prof_out(cluster, fileout, nrad=20, projected=False):
    # Write density profile (pprof.npy)
    fileout.write("%f " % (cluster.tphys))
    if cluster.rn == None or len(cluster.rn) != nrad:
        rn = rlagrange(cluster, nlagrange=nrad, projected=projected)
    mn = cluster.mtot / float(nrad)
    p_prof = []

    for r in rn:
        fileout.write("%f " % (r))

    for i in range(0, len(rn)):
        if i == 0:
            rmin = 0.0
            rmax = rn[i]
            vol = 4.0 * np.pi * (rmax ** 3.0) / 3.0
        else:
            rmin = rn[i - 1]
            rmax = rn[i]
            vol = 4.0 * np.pi * (rmax ** 3.0) / 3.0 - 4.0 * np.pi * (rmin ** 3.0) / 3.0

        p_prof.append(mn / vol)
        fileout.write("%f " % (p_prof[-1]))

    fileout.write("\n")


def alpha_prof_out(
    cluster,
    fileout,
    mmin=None,
    mmax=None,
    rmin=None,
    rmax=None,
    kwmin=None,
    kwmax=None,
    projected=False,
):
    # Write alpha_profile and dalpha for a given mass and radius range (alpha_prof.npy)
    m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=10,
        rmin=rmin,
        rmax=rmax,
        kwmin=kwmin,
        kwmax=kwmax,
        projected=projected,
    )
    lrprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = alpha_prof(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=10,
        kwmin=kwmin,
        kwmax=kwmax,
        projected=projected,
    )

    fileout.write("%f %f %f %f %f " % (cluster.tphys, alpha, ealpha, yalpha, eyalpha))
    for i in range(0, len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0, len(dm)):
        fileout.write("%f " % dm[i])
    for i in range(0, len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0, len(aprof)):
        fileout.write("%f " % aprof[i])

    fileout.write("%f %f %f %f\n" % (dalpha, edalpha, ydalpha, eydalpha))


def obs_alpha_prof_out(
    cluster,
    fileout,
    mmin=None,
    mmax=None,
    nmass=10,
    rmin=None,
    rmax=None,
    kwmin=None,
    kwmax=None,
    indx=None,
    projected=False,
    omask=None,
    **kwargs
):

    #Get custom values of provided:
    mtot=kwargs.pop('mtot',cluster.mtot)
    rm=kwargs.pop('rm',cluster.rmpro)
    trh=kwargs.pop('trh',relaxation_time(cluster,projected=True))

    # Write alpha_profile and dalpha for a given mass and radius range (alpha_prof.npy)
    m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbincorr = obs_mass_function(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=nmass,
        rmin=rmin,
        rmax=rmax,
        kwmin=kwmin,
        kwmax=kwmax,
        indx=indx,
        projected=projected,
        omask=omask,
        **kwargs,
    )

    if omask is None:
        rn=rlagrange(cluster,projected=projected)
        r50lower=rn[3]
        r50upper=rn[5]
    else:
        try:
            if projected:
                r50lower=omask.rmlower*(cluster.rmpro/omask.rm)
                r50upper=omask.rmupper*(cluster.rmpro/omask.rm) 
            else:
                r50lower=omask.rmlower*(cluster.rm/omask.rm)
                r50upper=omask.rmupper*(cluster.rm/omask.rm)  
        except:
            r50lower=None
            r50upper=None

        if r50lower==None:
            rn=rlagrange(cluster,projected=projected)
            r50lower=rn[3]
            r50upper=rn[5] 


    if projected:
        print('ALPHA50PRO: ',r50lower,cluster.rmpro,r50upper)
    else:
        print('ALPHA50', r50lower,cluster.rm,r50upper)

      
    lrprofn, aprof, dalpha, edalpha, ydalpha, eydalpha, mbincorr = obs_alpha_prof(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=nmass,
        rmin=rmin,
        rmax=rmax,
        kwmin=kwmin,
        kwmax=kwmax,
        indx=indx,
        projected=projected,
        omask=omask,
        **kwargs,
    )

    m_mean50, m_hist50, dm50, alpha50, ealpha50, yalpha50, eyalpha50, mbincorr50 = obs_mass_function(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=nmass,
        rmin=r50lower,
        rmax=r50upper,
        kwmin=kwmin,
        kwmax=kwmax,
        indx=indx,
        projected=projected,
        omask=omask,
        **kwargs,
    )

    print('a_g \t R/rm \t',lrprofn,'\t da/dlogr \t rm')
    print(alpha,aprof,dalpha,cluster.rm)
    print(mbincorr < 0.5)

    fileout.write("%f %f %f %f " % (cluster.tphys, mtot, rm, trh))
    fileout.write("%f %f %f %f " % (alpha50, ealpha50, yalpha50, eyalpha50))
    fileout.write("%f %f %f %f " % (alpha, ealpha, yalpha, eyalpha))

    for i in range(0, len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0, len(dm)):
        fileout.write("%f " % dm[i])
    for i in range(0, len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0, len(aprof)):
        fileout.write("%f " % aprof[i])

    fileout.write("%f %f %f %f\n" % (dalpha, edalpha, ydalpha, eydalpha))

def dalpha_out(
    cluster,
    fileout,
    mmin=[0.1, 0.3, 0.5],
    mmax=[0.5, 0.8, 0.8],
    rmin=None,
    rmax=None,
    kwmin=0,
    kwmax=1,
    projected=False,
):
    # Output alpha and dalpha for a range of values (dalpha_prof.npy)

    fileout.write("%f %f " % (cluster.tphys, cluster.mtot))

    for i in range(0, len(mmin)):
        m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(
            cluster,
            mmin=mmin[i],
            mmax=mmax[i],
            nmass=10,
            rmin=rmin,
            rmax=rmax,
            kwmin=kwmin,
            kwmax=kwmax,
            projected=projected,
        )
        fileout.write("%f " % alpha)

    for i in range(0, len(mmin)):
        lrprofn, aprof, da, eda, yda, eyda = alpha_prof(
            cluster, mmin=mmin[i], mmax=mmax[i], nmass=10, projected=projected
        )
        fileout.write("%f %f %f\n " % (da, yda, cluster.rm))


def sigv_out(cluster, fileout, projected=False):
    # Output velocity dispersion profile and anisotropy profile (dvprof.npy)

    fileout.write("%f %f " % (cluster.tphys, cluster.mtot))

    lrprofn, sigvprof, betaprof = sigv_prof(cluster, projected=projected)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for sig in sigvprof:
        fileout.write("%f " % sig)
    for beta in betaprof:
        fileout.write("%f " % beta)

    fileout.write("%f\n" % (cluster.rm))


def eta_prof_out(
    cluster,
    fileout,
    mmin=0.3,
    mmax=0.8,
    rmin=None,
    rmax=None,
    kwmin=0,
    kwmax=1,
    projected=False,
):
    # output eta profile (eta_prof.npy)

    m_mean, sigvm, eta, eeta, yeta, eyeta = eta_function(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=10,
        rmin=rmin,
        rmax=rmax,
        kwmin=kwmin,
        kwmax=kwmax,
        projected=projected,
    )
    lrprofn, eprof, deta, edeta, ydeta, eydeta = eta_prof(
        cluster,
        mmin=mmin,
        mmax=mmax,
        nmass=10,
        kwmin=kwmin,
        kwmax=kwmax,
        projected=projected,
    )

    fileout.write("%f %f %f %f %f " % (cluster.tphys, eta, eeta, yeta, eyeta))
    for i in range(0, len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0, len(sigvm)):
        fileout.write("%f " % sigvm[i])
    for i in range(0, len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0, len(eprof)):
        fileout.write("%f " % eprof[i])

    fileout.write("%f %f %f %f\n" % (deta, edeta, ydeta, eydeta))


def eta_out(
    cluster,
    fileout,
    mmin=[0.1, 0.3, 0.5],
    mmax=[0.5, 0.8, 0.8],
    rmin=None,
    rmax=None,
    kwmin=0,
    kwmax=1,
    projected=False,
):
    # Output eta and deta for a range of values (deta_prof.npy)

    fileout.write("%f %f " % (cluster.tphys, cluster.mtot))

    for i in range(0, len(mmin)):
        m_mean, sigvm, eta, eeta, yeta, eyeta = eta_function(
            cluster,
            mmin=mmin[i],
            mmax=mmax[i],
            nmass=10,
            rmin=rmin,
            rmax=rmax,
            kwmin=kwmin,
            kwmax=kwmax,
            projected=projected,
        )
        fileout.write("%f " % eta)

    for i in range(0, len(mmin)):
        lrprofn, eprof, deta, edeta, ydeta, eydeta = eta_prof(
            cluster,
            mmin=mmin[i],
            mmax=mmax[i],
            nmass=10,
            kwmin=kwmin,
            kwmax=kwmax,
            projected=projected,
        )
        fileout.write("%f %f %f\n " % deta, ydeta, cluster.rm)


def v_out(cluster, fileout, coord=None, projected=False):
    # Output mean velocity profile (vprof.npy)

    fileout.write("%f %f " % (cluster.tphys, cluster.mtot))

    lrprofn, vprof = v_prof(cluster, coord=coord, projected=projected)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for v in vprof:
        fileout.write("%f " % v)

    fileout.write("%f\n" % (cluster.rm))
