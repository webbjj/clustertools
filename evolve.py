import numpy as np
from galpy.potential import LogarithmicHaloPotential,MWPotential2014,rtide
from operations import *
from orbit import *

def evolve(cluster,dt=1.0,tfinal=12000.,pot=MWPotential2014,multimass=False,specific=True,do_batch=True,rtiterate=0,rgc=None,r0=8.,v0=220.):
    """
    NAME:

       evolve

    PURPOSE:

       Evolve a Star Cluster
       Notes:
        - Gieles et al. (2014) (EMACS)
        - Ebrahimi et al. 2019
    INPUT:

       cluster - StarCluster instance

    OUTPUT:

       time evolution of cluster mass and size

    HISTORY:

       2019 - Written - Webb (UofT)

    """    
    units0,origin0,center0=save_cluster(cluster)
    cluster.to_center()
    cluster.to_nbody()

    #Find constants
    #Gravitational Constant (pc km/s^2 / Msun)
    grav=1.0

    if multimass:
        lnlambda=np.log(0.02*cluster.ntot)
    else:
        lnlambda=np.log(0.11*cluster.ntot)

    mbar=cluster.mtot/float(cluster.ntot)

    trh=0.138*np.sqrt(float(cluster.ntot))*((cluster.rm)**(3./2.))/(np.sqrt(grav*mbar)*lnlambda)

    try:
        etot=cluster.ektot+cluster.ptot
    except:
        energies(cluster,specific=specific,do_batch=do_batch)
        etot=cluster.ektot+cluster.ptot

    rh=cluster.rm
    #Assume plummer:
    rv=1.3*rh
    rc=0.4*rh
    mgc=cluster.mtot

    initialize_orbit(cluster)
    ts=np.linspace(cluster.tphys/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),tfinal/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),(tfinal-cluster.tphys)/dt)
    cluster.orbit.integrate(ts,pot)
    rgc=cluster.orbit.r(ts)

    for t in ts:
    	Rhj=rh/rtidal(cluster,pot=pot,rtiterate=rtiterate,rgc=cluster.orbit.r(t),r0=r0,v0=v0)
    	Rvj=rv/rtidal(cluster,pot=pot,rtiterate=rtiterate,rgc=cluster.orbit.r(t),r0=r0,v0=v0)
    	Rcj=cluster.r10/rtidal(cluster,pot=pot,rtiterate=rtiterate,rgc=cluster.orbit.r(t),r0=r0,v0=v0)







