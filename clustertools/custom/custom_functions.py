from galpy.orbit import Orbit

try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion


from galpy import potential
from galpy.potential import MWPotential2014
from galpy.actionAngle import actionAngleStaeckel
from galpy.actionAngle.actionAngleIsochroneApprox import actionAngleIsochroneApprox

import numpy as np
import matplotlib.pyplot as plt

from ..util.recipes import interpolate, binmaker
from ..util.plots import starplot,skyplot,_plot,_lplot,_scatter
from ..util.constants import *

import astropy.coordinates as coord
import astropy.units as u

def integrate_orbit(
    cluster, pot=MWPotential2014, tfinal=None, nt=1000, from_centre=False, ro=solar_ro, vo=solar_vo, solarmotion=solar_motion, plot=False
):
    """Integrate a galpy orbit instance for the cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential that orbit is to be integrate in (default: MWPotential2014)
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    from_centre : bool
        intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    solarmotion : float
        array representing U,V,W of Sun (default: solarmotion=solar_motion)
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
    o = initialize_orbit(cluster,from_centre=from_centre,solarmotion=solarmotion)

    if tfinal is None:
        tfinal=12./conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='pckms':
        tfinal/=1000.
    elif cluster.units=='kpckms':
        tfinal/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        tfinal*=(cluster.tbar/1000.)

    ts = np.linspace(0, tfinal, nt)
    o.integrate(ts, pot)

    if plot:
        o.plot()

    return ts, o

def integrate_orbits(
    cluster, pot=None, tfinal=None, nt=1000, ro=solar_ro, vo=solar_vo,solarmotion=solar_motion, plot=False
):
    """Integrate a galpy orbit instance for each star

    Parameters
    ----------
    cluster : class
        StarCluster
    pot : class
        Galpy potential for host cluster that orbit is to be integrated in
        if None, assume a Plumme Potential
    tfinal : float
        final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
    nt : int
        number of timesteps
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    solarmotion : float
        array representing U,V,W of Sun (default: solarmotion=solar_motion)
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

    if tfinal is None:
        tfinal=12./conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='pckms':
        tfinal/=(1000.*conversion.time_in_Gyr(ro=ro, vo=vo))
    elif cluster.units=='kpckms':
        tfinal/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        tfinal*=((cluster.tbar/1000.)/conversion.time_in_Gyr(ro=ro, vo=vo))

    if pot is None:

        if cluster.origin=='cluster' or cluster.origin=='centre':
            cluster.save_cluster()
            units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0
            cluster.to_galpy()
            pot=potential.PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
            cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

        else:
            pot=MWPotential2014

    os = initialize_orbits(cluster, ro=ro,vo=vo, solarmotion=solarmotion)
    ts = np.linspace(0, tfinal, nt)

    #integrate orbits of stars in combined potential of GC and galaxy
    os.integrate(ts, pot)

    if plot:
        os.plot()

    return ts, os

def calc_action(cluster, pot=None, ro=solar_ro, vo=solar_vo,solarmotion=solar_motion,full=False, **kwargs):
    """Calculate action angle values for cluster

    - This is a simple wrapper for calculating actions from an Orbit in galpy (Bovy 2015)
    -- Bovy J., 2015, ApJS, 216, 29
    -- 

    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate actions (default: MWPotential)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    solarmotion : float
        array representing U,V,W of Sun (default: solarmotion=solar_motion)
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

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if pot is None:
        pot=MWPotential2014

    os = initialize_orbit(cluster, ro=ro, vo=vo, solarmotion=solarmotion)

    if pot==MWPotential2014:
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

    else:
        ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
        os.integrate(ts,pot)

        JR=os.jr()
        Jphi=os.jp()
        Jz=os.jz()

        if full:
            OR = os.Or()
            Ophi = os.Op()
            Oz = os.Oz()
            TR = os.Tr()
            Tphi = os.Tp()
            Tz = os.Tz()

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if full:
        return JR, Jphi, Jz, OR, Ophi, Oz, TR, Tphi, Tz
    else:
        return JR, Jphi, Jz

def calc_actions(cluster, pot=None, ro=solar_ro, vo=solar_vo,solarmotion=solar_motion,full=False, **kwargs):
    """Calculate action angle values for each star

    - This is a simple wrapper for calculating actions from an Orbit in galpy (Bovy 2015)
    -- Bovy J., 2015, ApJS, 216, 29
    -- 

    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate actions (default: MWPotential for galaxy or Plummer for cluster)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    solarmotion : float
        array representing U,V,W of Sun (default: solarmotion=solar_motion)
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

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0


    if pot is None:
        if cluster.origin=='cluster' or cluster.origin=='centre':
            cluster.to_galpy()
            pot=potential.PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
            cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)
        else:
            pot=MWPotential2014

    os = initialize_orbits(cluster, ro, vo, solarmotion=solarmotion)

    if pot==MWPotential2014:
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

    else:
        ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
        os.integrate(ts,pot)

        JR=os.jr()
        Jphi=os.jp()
        Jz=os.jz()

        if full:
            OR = os.Or()
            Ophi = os.Op()
            Oz = os.Oz()
            TR = os.Tr()
            Tphi = os.Tp()
            Tz = os.Tz()

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    if full:
        return JR, Jphi, Jz, OR, Ophi, Oz, TR, Tphi, Tz
    else:
        return JR, Jphi, Jz


def ttensor(cluster, pot=None, ro=solar_ro, vo=solar_vo,solarmotion=solar_motion, eigenval=False, t=0.):
    """Calculate the tidal tensor Tij=-d(Psi)(dxidxj)
    
    - This is a simple wrapper for calculating the tidal tensor in a potential in galpy (Bovy 2015)
    -- Bovy J., 2015, ApJS, 216, 29
    -- Webb, J.J., Bovy, J., Carlberg, R.G., Gieles, M. 2019, MNRAS, 448, 4
    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate actions (default: MWPotential)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    solarmotion : float
        array representing U,V,W of Sun (default: solarmotion=solar_motion)
    eigenval : bool
        return eigenvalues if true (default; False)
    time : float
        time (in cluster.units) to evaluate tidal tensor. Necessary if tidal field is time dependent (default: 0.)

    Returns
    -------
        Tidal Tensor
    History
    -------
    2018-03-21 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if pot is None:
        pot=MWPotential2014

    o = initialize_orbit(cluster, ro=ro, vo=vo, solarmotion=solarmotion)
    R=o.R()
    z=o.z()
    phi=o.phi()

    if cluster.units=='pckms':
        t/=1000.
    elif cluster.units=='kpckms':
        t/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        t*=(cluster.tbar/1000.)

    tij=potential.ttensor(pot,R/ro,z/ro,phi=phi,t=t/conversion.time_in_Gyr(ro=ro, vo=vo),eigenval=eigenval)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return tij

def ttensors(cluster, pot=None, ro=solar_ro, vo=solar_vo,solarmotion=solar_motion, eigenval=False, t=0.):
    """Calculate the tidal tensor Tij=-d(Psi)(dxidxj) acting on all stars in the cluster

    - This is a simple wrapper for calculating the tidal tensor in a potential in galpy (Bovy 2015)
    -- Bovy J., 2015, ApJS, 216, 29
    -- Webb, J.J., Bovy, J., Carlberg, R.G., Gieles, M. 2019, MNRAS, 448, 4
    Parameters
    ----------
    cluster : class
        StarCluster instance
    pot : class 
        GALPY potential used to calculate actions (default: MWPotential for galaxy or Plummer for cluster)
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    solarmotion : float
        array representing U,V,W of Sun (default: solarmotion=solar_motion)
    eigenval : bool
        return eigenvalues if true (default; False)
    time : float
        time (in cluster.units) to evaluate tidal tensor. Necessary if tidal field is time dependent (default: 0.)

    Returns
    -------
        Tidal Tensor
    History
    -------
    2018-03-21 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if pot is None:
        if cluster.origin=='cluster' or cluster.origin=='centre':
            cluster.to_galpy()
            pot=potential.PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
            cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)
        else:
            pot=MWPotential2014

    o = initialize_orbits(cluster, ro=ro, vo=vo, solarmotion=solarmotion)
    R=o.R()
    z=o.z()
    phi=o.phi()

    if cluster.units=='pckms':
        t/=1000.
    elif cluster.units=='kpckms':
        t/=conversion.time_in_Gyr(ro=ro, vo=vo)
    elif cluster.units=='nbody':
        t*=(cluster.tbar/1000.)

    tij=potential.ttensor(pot,R/ro,z/ro,phi=phi,t=t/conversion.time_in_Gyr(ro=ro, vo=vo),eigenval=eigenval)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    return tij


