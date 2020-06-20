"""Functions to calculate StarCluster properties

Designed to accept StarCluster instance as in put to
calculate key parameters

"""

__author__ = "Jeremy J Webb"


__all__ = [
  "relaxation_time",
  "half_mass_relaxation_time",
  "core_relaxation_time",
  "energies",
  "closest_star",
  "virialize",
  "virial_radius",
  "virial_radius_inverse_distance",
  "virial_radius_critical_density",
  "mass_function",
  "eta_function",
  "surface_area",
]

import numpy as np
import numba
from galpy.util import bovy_coords

from ..util.constants import *
from ..util.recipes import *
from .operations import *
from ..util.plots import *

def relaxation_time(cluster, rad=None, multimass=True, projected=False,method='spitzer'):
    """
    NAME:

       relaxation_time

    PURPOSE:

       Calculate the relaxation time (Spitzer & Hart 1971) within a given radius of the cluster

    INPUT:

       cluster - StarCluster instance

       rad - radius within which to calculate the relaxation time

       multimass - use multimass (True) or single mass (False) value for ln lambda (default: True)

       projected - use projected values (default: False)

       method - choose between Spitzer & Hart 1971 and other methods to be added later

    OUTPUT:

       trelax

    HISTORY:

       2020 - Written - Webb (UofT)

    """

    if rad is None and projected:
        rad=cluster.rmpro
    elif rad is None:
        rad=cluster.rm

    grav=4.302e-3

    if projected:
        rindx=cluster.rpro < rad
    else:
        rindx=cluster.r < rad
        
    ntot=np.sum(rindx)
    mbar=np.mean(cluster.m[rindx])
    vol=4.0*np.pi*(rad**3.)/3.0
    rho=ntot/vol
    
    v2=np.mean(cluster.v**2.)
    
    #v2=0.4*grav*cluster.mtot/rad
    
    lnlambda=np.log(0.4*cluster.ntot)
    
    trelax=v2**(3./2.)/(15.4*grav**2.*mbar**2.*rho*lnlambda)

    # Units of Myr
    trelax*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

    return trelax

def half_mass_relaxation_time(cluster, multimass=True, projected=False):
    """
    NAME:

       relaxation_time

    PURPOSE:

       Calculate the half-mass relaxation time (Spitzer 1987) of the cluster

    INPUT:

       cluster - StarCluster instance

       multimass - use multimass (True) or single mass (False) value for ln lambda (default: True)

       projected - use projected values (default: False)

    OUTPUT:

       trelax

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    grav=4.302e-3
    mass=cluster.mtot
    ntot=float(cluster.ntot)
    mbar=mass/ntot
    lnlambda = np.log(0.4*ntot)

    if projected:
        rm=cluster.rmpro
    else:
        rm=cluster.rm


    trelax=0.138*(mass**0.5)*(rm**1.5)/(mbar*np.sqrt(grav)*lnlambda)
    # Units of Myr
    trelax*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

    return trelax


def core_relaxation_time(cluster, multimass=True, projected=False):
    """
    NAME:

       core_relaxation_time

    PURPOSE:

       Calculate the core relaxation time (Stone & Ostriker 2015) of the cluster

    INPUT:

       cluster - StarCluster instance

       multimass - use multimass (True) or single mass (False) value for ln lambda (default: True)

       projected - use projected values (default: False)

       method - choose between Spitzer 1987 and other methods to be added later

    OUTPUT:

       trelax

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    lnlambda=np.log(0.4*cluster.ntot)
    mtot=cluster.mtot
    mbar=np.mean(cluster.m)
    rc=cluster.r10
    rh=cluster.rm
    grav=4.302e-3

    trelax=(0.39/lnlambda)*np.sqrt(rc**3./(grav*mtot))*(mtot/mbar)*np.sqrt(rc*rh)/(rc+rh)

    return trelax


def energies(cluster, specific=True, i_d=None, full=True, parallel=False):
    """
    NAME:

       energies

    PURPOSE:

       Calculate kinetic and potential energy of every star

    INPUT:

       cluster - StarCluster instance

       specific - find specific energies (default: True)

       i_d - find energies for a specific star

       full - calculate distance of full array of stars at once with numbra (default: True)

       parallel - calculate distances in parallel if True (default: False)

    OUTPUT:

       ek,pot,etot

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)
    cluster.to_centre()

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

    if specific:
        ek = 0.5 * (cluster.v ** 2.0)
    else:
        ek = 0.5 * cluster.m * (cluster.v ** 2.0)

    if i_d != None:
        indx = cluster.id == i_d

        dx = cluster.x[indx] - cluster.x
        dy = cluster.y[indx] - cluster.y
        dz = cluster.z[indx] - cluster.z

        if specific:
            m = cluster.m
        else:
            m = cluter.m[indx] * cluster.m

        dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + dz ** 2.0)
        rindx = dr != 0.0
        gmr = -grav * m[rindx] / dr[rindx]

        pot = np.sum(gmr)
        ek = ek[indx]
        etot = ek + pot

    elif full:
        x = np.array([cluster.x, cluster.y, cluster.z, cluster.m]).T
        if parallel:
            pot = grav * np.array(potential_energy_parallel(x))
        else:
            pot = grav * np.array(potential_energy(x))

        if specific:
            pot /= cluster.m

        etot = ek + pot
        cluster.add_energies(ek, pot, etot)
    else:
        pot = []

        for i in range(0, cluster.ntot):
            dx = cluster.x[i] - cluster.x
            dy = cluster.y[i] - cluster.y
            dz = cluster.z[i] - cluster.z
            if specific:
                m = cluster.m
            else:
                m = cluter.m[i] * cluster.m

            dr = np.sqrt(dx ** 2.0 + dy ** 2.0 + dz ** 2.0)
            indx = dr != 0.0
            gmr = -grav * m[indx] / dr[indx]

            pot.append(np.sum(gmr))

        etot = ek + pot
        cluster.add_energies(ek, pot, etot)

    return_cluster(cluster, units0, origin0)

    return ek, pot, etot


@numba.njit
def potential_energy(cluster):
    """
    NAME:

       potential_energy

    PURPOSE:

       Find potential energy for each star in a cluster

    INPUT:

       cluster=[x,y,z,m].T

    OUTPUT:

        energy

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    energy = [0.0] * len(cluster)
    for i in range(len(cluster) - 1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i], cluster[j])
            m2 = cluster[i, 3] * cluster[j, 3]
            energy[i] += -m2 / r
            energy[j] += -m2 / r

    return energy


@numba.njit(parallel=True)
def potential_energy_parallel(cluster):
    """
    NAME:

       potential_energy

    PURPOSE:

       Find potential energy for each star in a cluster (done in parallel)

    INPUT:

       cluster=[x,y,z,m].T

    OUTPUT:

        energy

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    energy = [0.0] * len(cluster)
    for i in numba.prange(len(cluster) - 1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i], cluster[j])
            m2 = cluster[i, 3] * cluster[j, 3]
            energy[i] += -m2 / r
            energy[j] += -m2 / r

    return energy


def closest_star(cluster, projected=False):

    if projected:
        z = np.zeros(cluster.ntot)
        x = np.array([cluster.x, cluster.y, z]).T
    else:
        x = np.array([cluster.x, cluster.y, cluster.z]).T
    return minimum_distance(x)


def virialize(cluster, specific=True, full=True):
    """
    NAME:

       virialize

    PURPOSE:

       Adjust stellar velocities so cluster is in virial equilibrium

    INPUT:

       cluster - StarCluster instance
       specific - find specific energies (default: True)
       full - do full array of stars at once with numbra (default: True)

    OUTPUT:

       qv

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre()

    try:
        qv = np.sqrt(np.abs(0.5 / cluster.qvir))
    except:
        print("NEED TO CALCULATE ENERGIES FIRST")
        energies(cluster, specific=specific, full=full)
        qv = np.sqrt(np.abs(0.5 / cluster.qvir))

    cluster.vx *= qv
    cluster.vy *= qv
    cluster.vz *= qv
    cluster.key_params()

    energies(cluster, specific=specific, full=full)
    qv = np.sqrt(np.abs(0.5 / cluster.qvir))

    return_cluster(cluster, units0, origin0)

    return qv


def rlagrange(cluster, nlagrange=10, projected=False):
    """
    NAME:

       rlagrange

    PURPOSE:

       Calculate lagrange radii of the cluster by mass
       --> Note units of lagrange radii will be equal to cluster.units

    INPUT:

       cluster - StarCluster instance
       nlagrange - number of lagrange radii bins (default: 10)
       projected - calculate projected lagrange radii (default: False)

    OUTPUT:

       rn

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre()

    # Radially order the stars
    msum = 0.0
    nfrac = 1
    rn = []

    if projected:
        if cluster.rproorder is None:
            rorder = np.argsort(cluster.rpro)
        else:
            rorder = cluster.rproorder
    else:
        if cluster.rorder is None:
            rorder = np.argsort(cluster.r)
        else:
            rorder = cluster.rorder

    for i in range(0, cluster.ntot):
        mfrac = cluster.mtot * float(nfrac) / float(nlagrange)
        msum += cluster.m[rorder[i]]
        if msum >= mfrac:
            rn.append(cluster.r[rorder[i]])
            nfrac += 1

    while len(rn) != nlagrange:
        rn.append(np.max(cluster.r))

    return_cluster(cluster, units0, origin0)

    return rn

def virial_radius(cluster, method='inverse_distance',
    full=True,
    H=70.0,
    Om=0.3,
    overdens=200.0,
    nrad=20,
    projected=False,
    plot=False,
    **kwargs):

    if method=='inverse_distance':
        rv=virial_radius_inverse_distance(cluster,projected=projected,full=full)
    else:
        rv=virial_radius_critical_density(cluster,H,Om,overdens,nrad,projected,plot,**kwargs)

    return rv

def virial_radius_inverse_distance(cluster, projected=False, full=True):
    """
    NAME:

       virial_radius

    PURPOSE:

       Calculate virial radius of the cluster as the inverse of the 
                 average inverse distance between particles, weighted by their masses
       --> Definition taken from AMUSE (www.amusecode.org)

    INPUT:

       cluster - StarCluster instance

       projected - calculate projected virial radius (default: False)

       full - Use Numba (default:True)

    OUTPUT:

       rv

    HISTORY:

       2019 - Written - Webb (UofT)


    """

    if full:
        if projected:
            x = np.array([cluster.x, cluster.y, np.zeros(cluster.ntot), cluster.m]).T
        else:
            x = np.array([cluster.x, cluster.y, cluster.z, cluster.m]).T

        ms = cluster.m
        partial_sum = weighted_inverse_distance_sum(x)

    else:
        partial_sum = 0.0

        ms = cluster.m
        xs = cluster.x
        ys = cluster.y
        if projected:
            zs = np.zeros(cluster.ntot)
        else:
            zs = cluster.z

        for i in range(cluster.ntot - 1):
            x = xs[i]
            y = ys[i]
            z = zs[i]
            dx = x - xs[i + 1 :]
            dy = y - ys[i + 1 :]
            dz = z - zs[i + 1 :]
            dr2 = (dx * dx) + (dy * dy) + (dz * dz)
            dr = np.sqrt(dr2)
            m_m = ms[i] * ms[i + 1 :]
            partial_sum += np.sum(m_m / dr)

    return (np.sum(ms) ** 2) / (2 * partial_sum)


@numba.njit
def weighted_inverse_distance_sum(cluster):
    """
    NAME:

       weighted_inverse_distance_sum

    PURPOSE:

       Find the sum of the mass weighted inverse distance for each star

    INPUT:

       cluster=[x,y,z,m].T

    OUTPUT:

        weighted inverse distance sum

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    weighted_sum = 0.0
    for i in range(len(cluster) - 1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i], cluster[j])
            m2 = cluster[i, 3] * cluster[j, 3]
            weighted_sum += m2 / r

    return weighted_sum


def virial_radius_critical_density(
    cluster,
    H=70.0,
    Om=0.3,
    overdens=200.0,
    nrad=20,
    projected=False,
    plot=False,
    **kwargs
):
    """
    NAME:

       rvirial

    PURPOSE:

       Calculate virial radius of the cluster as the radius at which the density is equal to the critical 
           density of the Universe at the redshift of the system, multiplied by an overdensity constant

    INPUT:

       cluster - StarCluster instance

       H - Hubble constant

       Om - density of matter

       overdens - overdensity constant

       nrad - number of radial bins used to calculate cluster density profile

       projected - calculate projected virial radius (default: False)

    OUTPUT:

       rv

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_pckms()
    cluster.to_centre()

    H /= 1000000.0  # (km/s) / pc
    Grav = 4.302e-3  # pc (km/s)^2 / Msun

    rhocrit = 3.0 * (H ** 2.0) / (8.0 * np.pi * Grav)  # Msun/pc^3

    if projected:
        indx - cluster.rproorder
    else:
        indx = cluster.rorder

    msum = np.cumsum(cluster.m[indx])

    if projected:
        vsum = (4.0 / 3.0) * np.pi * (cluster.rpro[indx] ** 3.0)
        pprof = msum / vsum
        rprof = cluster.rpro[indx]
    else:
        vsum = (4.0 / 3.0) * np.pi * (cluster.r[indx] ** 3.0)
        pprof = msum / vsum
        rprof = cluster.r[indx]

    # Find radius where maxium density occurs
    rindx = np.argmax(pprof)
    rmax = rprof[rindx]

    indx1 = (rprof > rmax) * (pprof > rhocrit * overdens)
    indx2 = (rprof > rmax) * (pprof < rhocrit * overdens)

    if np.sum(indx2) == 0.0:
        print("SYSTEM IS NOT VIRIALIZED")
        r_v = -1.0
    else:
        r1 = rprof[indx1][-1]
        r2 = rprof[indx2][0]

        rho1 = pprof[indx1][-1]
        rho2 = pprof[indx2][0]

        print(r1, r2, rho1, rho2, rhocrit * overdens)
        r_v = interpolate([r1, rho1], [r2, rho2], y=rhocrit * overdens)

    if plot:
        rho_local = rhocrit * overdens

        filename = kwargs.pop("filename", None)
        overplot = kwargs.pop("overplot", False)

        xunits = " (pc)"
        if projected:
            yunits = " Msun/pc^2"
        else:
            yunits = " Msun/pc^3"

        x, y = rprof, pprof
        nlplot(
            x,
            y,
            xlabel=r"$R" + xunits + "$",
            ylabel=r"$rho" + yunits + "$",
            title="Time = %f" % cluster.tphys,
            log=True,
            overplot=overplot,
            filename=filename,
        )
        nlplot(x, np.ones(len(x)) * rho_local, "--", overplot=True)
        nlplot(np.ones(len(y)) * r_v, y, "--", overplot=True)

        if filename != None:
            plt.savefig(filename)

    return_cluster(cluster, units0, origin0)

    return r_v


def new_mass_function(
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
    projected=False,
    mcorr=None,
    omask=None,
    plot=False,
    **kwargs
):
    """
    NAME:

       mass_function

    PURPOSE:

       Find mass function over a given mass range using nmass bins containing an equal number of stars

    INPUT:

       cluster - StarCluster instance
       mmin/mmax - specific mass range
       nmass - number of mass bins used to calculate alpha
       rmin/rmax - specific radial range
       vmin/vmax - specific velocity range
       emin/emax - specific energy range
       kwmin/kwmax - specific stellar evolution type range
       indx - specific subset of stars
       projected - use projected values
       mcorr - correction for masses
       omask - place a mask over the dataset to mimic observed data
       plot - plot the mass function
       **kwargs - key words for plotting

    OUTPUT:

       if omask is None: m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha
       if omask is not None: m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha, mbinerror

    HISTORY:

       2018 - Written - Webb (UofT)
    """

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
        * (cluster.m < mmax)
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
            m_corr_hist=m_hist
            mbinerror=np.ones(nmass)
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
                m_hist[i]=np.sum(mindx)
                m_corr_hist[i] = np.sum(1.0 / mcorr[mindx])

            mbinerror = m_hist / m_corr_hist

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

        if return_error:
            return m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha, mbinerror
        else:
            return m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha
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
        )

def mass_function(
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
    projected=False,
    plot=False,
    **kwargs
):
    """
    NAME:

       mass_function

    PURPOSE:

       Find mass function over a given mass range using nmass bins containing an equal number of stars

    INPUT:

       cluster - StarCluster instance
       mmin/mmax - specific mass range
       nmass - number of mass bins used to calculate alpha
       rmin/rmax - specific radial range
       vmin/vmax - specific velocity range
       emin/emax - specific energy range
       kwmin/kwmax - specific stellar evolution type range
       indx - specific subset of stars
       projected - use projected values
       plot - plot the mass function
       **kwargs - key words for plotting

    OUTPUT:

       if omask is None: m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha
       if omask is not None: m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha, binerror

    HISTORY:

       2018 - Written - Webb (UofT)
    """

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

    if np.sum(indx) >= nmass:

        m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)

        lm_mean = np.log10(m_mean)
        dm = m_hist / (m_upper - m_lower)
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

        return m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha
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
        )


def eta_function(
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
    projected=False,
    plot=False,
    **kwargs
):
    """
    NAME:

       eta_function

    PURPOSE:

       Find eta over a given mass range using nmass bins containing an equal number of stars

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - specific mass range

       nmass - number of mass bins used to calculate eta

       rmin/rmax - specific radial range

       vmin/vmax - specific velocity range

       emin/emax - specific energy range

       kwmin/kwmax - specific stellar evolution type range

       indx - specific subset of stars

       projected - use projected values

       plot - plot the mass function
       
       **kwargs - key words for plotting

    OUTPUT:

       m_mean,sigvm,eta,eeta,yeta,eyeta

    HISTORY:

       2018 - Written - Webb (UofT)
    """
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

    if np.sum(indx) >= 2 * nmass:

        m_lower, m_mean, m_upper, m_hist = nbinmaker(cluster.m[indx], nmass)
        lm_mean = np.log10(m_mean)

        sigvm = []
        lsigvm = []
        for i in range(0, nmass):

            mindx = indx * (cluster.m >= m_lower[i]) * (cluster.m < m_upper[i])
            sigvm.append(np.std(v[mindx]))
            lsigvm.append(np.log10(sigvm[-1]))

        (eta, yeta), V = np.polyfit(lm_mean, lsigvm, 1, cov=True)
        eeta = np.sqrt(V[0][0])
        eyeta = np.sqrt(V[1][1])

        if plot:
            filename = kwargs.get("filename", None)
            nplot(m_mean, np.log10(sigvm), xlabel="M", ylabel=r"$\sigma_v$", **kwargs)
            mfit = np.linspace(np.min(m_mean), np.max(m_mean), nmass)
            sigfit = 10.0 ** (eta * np.log10(mfit) + yeta)
            nlplot(mfit, np.log10(sigfit), overplot=True, label=(r"$\eta$ = %f" % eta))
            plt.legend()

            if filename != None:
                plt.savefig(filename)

        return m_mean, sigvm, eta, eeta, yeta, eyeta
    else:
        print("NOT ENOUGH STARS TO ESTIMATE SIGMA-MASS RELATION")
        return (
            np.zeros(nmass),
            np.zeros(nmass),
            np.zeros(nmass),
            -1000.0,
            -1000.0,
            -1000.0,
            -1000.0,
        )


def surface_area(
    cluster,
    mmin=None,
    mmax=None,
    rmin=None,
    rmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=1,
    indx=None,
    projected=False,
    coords="xy",
    thresh=None,
    nrand=1000,
    method="sphere",
    full=False,
    plot=False,
):
    """
    NAME:

     surface_area

    PURPOSE:

     calculate surface area enclosed by cluster by finding what fraction of a random distribution overlaps with points

    INPUT:

     cluster - StarCluster instance
     
     mmin/mmax - specific mass range
     rmin/rmax - specific radial range
     vmin/vmax - specific velocity range
     emin/emax - specific energy range
     kwmin/kwmax - specific stellar evolution type range
     indx - specific subset of stars
     projected - use projected values
     
     coords - choose axis to project cluster (Default: xy)
     
     thresh - threshold for overlap between random distribution and points (Default: None - use maximum nearest neighbour)

     nrand - number of random points to be generated in uniform distribution (Default: 1000)

     method - generate spherical or rectangular distribution of random points (Default: sphere)
     
     plot - plot overlap (Default: False)

    OUTPUT:

     area

    HISTORY:

     2019 - Written - Webb (UofT)

    """

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

    if coords == "xy":
        x = cluster.x
        y = cluster.y
    elif coords == "xz":
        x = cluster.x
        y = cluster.z
    elif coords == "yz":
        x = cluster.y
        y = cluster.z

    return area_enclosed(
        x, y, thresh=thresh, nrand=nrand, method=method, full=full, plot=plot
    )
