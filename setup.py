import numpy as np
from galpy.util import bovy_conversion
import os
from cluster import StarCluster
import limepy
from limepy import limepy
from limepy import sample

from profiles import m_prof

def setup_cluster(ctype,units0='realpc',origin0='cluster',orbit=None,**kwargs):
    #setup_cluster is the main function for setting up clusters, built around LIMEPY
   
    #ctype = limepy - use limepy to setup initial cluster conditions (limepy parameters read in via **kwargs)

    #Default units are set to realpc and default origin is set to cluster. These are really only used when reading in nbodypy snapshots, as other codes have their own default units and origin

    #orbit - allows for a galpy orbit to be used to set the cluster's orbit

    #kwargs:
    #Limepy parameters
   

    #Generate cluster using limepy
    if ctype=='limepy':
        g=kwargs.pop('g')
        cluster=get_limepy(g=g,**kwargs)
    elif ctype=='woolley':
        g=kwargs.pop('g',0)
        cluster=get_limepy(g=g,**kwargs)        
    elif ctype=='king':
        g=kwargs.pop('g',1)
        cluster=get_limepy(g=g,**kwargs)
    elif ctype=='wilson':
        g=kwargs.pop('g',2)
        cluster=get_limepy(g=g,**kwargs)


   #Add galpy orbit if given
    if orbit!=None:
        cluster.orbit=orbit
        t=(cluster.tphys/1000.)/bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
        cluster.add_orbit(orbit.x(t),orbit.y(t),orbit.z(t),orbit.vx(t),orbit.vy(t),orbit.vz(t),'realkpc')

    cluster.key_params()

    return cluster

def get_limepy(g=1,**kwargs):
    phi0=float(kwargs.get('phi0'))
    project=bool(kwargs.get('project',False))

    if 'M' in kwargs:
        units='realpc'
        M=float(kwargs.get('M'))
        if 'rt' in kwargs:
            rt=float(kwargs.get('rt'))  
            lmodel=limepy(phi0,g,M=M,rt=rt,project=project)
        elif 'rv' in kwargs:
            rv=float(kwargs.get('rv'))      
            lmodel=limepy(phi0,g,M=M,rv=rv,project=project)
        elif 'rh' in kwargs:
            rh=float(kwargs.get('rh'))      
            lmodel=limepy(phi0,g,M=M,rv=rh,project=project)
        elif 'r0' in kwargs:
            r0=float(kwargs.get('r0'))      
            lmodel=limepy(phi0,g,M=M,r0=r0,project=project)
        else:
            lmodel=limepy(phi0,g,M=M,r0=1.,project=project)
    else:
        units='nbody'
        lmodel=limepy(phi0,g,G=1, M=1, rv=1,project=project)

    N=int(kwargs.get('N',1000))

    ldata=sample(lmodel,N=N)

    cluster=StarCluster(N,units=units,origin='cluster')
    cluster.add_stars(np.linspace(1,N,N,dtype=int),ldata.m,ldata.x,ldata.y,ldata.z,ldata.vx,ldata.vy,ldata.vz)
    cluster.find_center()
    cluster.key_params()

    return cluster

def get_imf(mlimits = [0.1, 50.0], alphas = [-1.35], ntot=1, mtot=None, do_random=True):
    """
    NAME:

       get_imf

    PURPOSE:

       Generate stellar masses from a defined IMF
       Notes:
        -- This has been heaviliy borrowed (stolen) from AMUSE (amusecode.org), specifically
           the routine /src/amuse/ic/brokenimf.py

    INPUT:

       mlimits - mass limits (default: 0.1-50.0)

       alphas - slopes of the mass function between mass limits (default: -1.35 - Salpeter)

       ntot - number of stars to be generated (default: 1)

       mtot - total mass of system (overrides ntot) (default: None)

       do_random - use randomly generated masses or evenly distributed masses (default: True)

    OUTPUT:

       trelax

    HISTORY:

       2019 - Written - Webb (UofT)

    """ 
    
    mlimits = np.array(mlimits)
    alphas = np.array(alphas)

    #If mtot is given, initially assume all stars are equal to the lower mass limit
    if mtot!=None:
        ntot=int(mtot/mlimits[0])
        print(ntot)

    if do_random:
        random = np.random.random(ntot)
    else:
        random = np.linspace(0.0, 1.0, ntot)

    nbins = len(alphas)
    
    #Calculate fraction per bin
    nbin = []

    for i in range(0,nbins):
        alpha=alphas[i]
        if alpha == -1:
            factor = np.log(mlimits[i+1] / mlimits[i])
        else:
            factor = (mlimits[i+1]**(alpha+1) - mlimits[i]**(alpha+1)) / (alpha+1)

        for j in range(0,nbins - i - 1):
            factor *= mlimits[-j-2]**(alphas[-j-1]-alphas[-j-2])

        nbin.append(factor)
    total = sum(nbin, 0.0)
    fbin=np.array(nbin)/total 

    #Calculate cumultive fraction per bin, factors, and inverse alphas 
    cumulative_fractions = np.array([sum(fbin[:i]) for i in range(nbins+1)])
    cumulative_fractions[-1] = 1.0     # In case of round-off errors
    factors = pow(mlimits[1:] / mlimits[:-1], alphas + 1.0) - 1.0
    inv_alpha1s = np.array([np.inf if alpha==-1 else (1.0 / (alpha + 1.0)) for alpha in alphas])

    #Calculate masses
    indices = np.searchsorted(cumulative_fractions[:-1], random, 'right') - 1
    scaled = ((random - cumulative_fractions[indices]) / 
        (cumulative_fractions[indices+1] - cumulative_fractions[indices]))

    result = np.empty_like(random)
    zerodiv = alphas[indices]==-1
    normal = np.logical_not(zerodiv)
    result[zerodiv] = pow(mlimits[1:][indices[zerodiv]] / mlimits[:-1][indices[zerodiv]], scaled[zerodiv])
    result[normal] = pow(1.0 + factors[indices[normal]] * scaled[normal], inv_alpha1s[indices[normal]])

    masses=mlimits[:-1][indices] * result

    #If mtot is given, only return enough stars so that total mass = mtot
    if mtot!=None:
        ntot=int(mtot/np.mean(masses))
        masses=masses[0:ntot]

    return np.array(masses)

def get_segregated_imf(cluster,delta_alpha,alpha50, mlimits = [0.1, 50.0],do_random=True):
    #TO DO - maintain number profile or mass profile?

    if False:
        rprof,mprof,nprof=m_prof(cluster,cumulative=True,nrad=50)
        rprof=np.array(rprof)
        mprof=np.array(mprof)

        m=[]
        msum=0.
        ns=0
        while ns < cluster.ntot:
            indx=cluster.rorder[ns]
            if cluster.r[indx] <= rprof[0]:
                rindx=0
            elif cluster.r[indx] >= rprof[-1]:
                rindx=-1
            else:
                rindx=np.argwhere(cluster.r[indx]<=rprof)[-1]

            #Find local alpha
            ydalpha=delta_alpha*np.log(cluster.r[indx]/cluster.rm)+alpha50
            mtemp=get_imf(mlimits = mlimits, alphas = [ydalpha], ntot=1,do_random=do_random)[0]
            
            if ns==0:
                m.append(mtemp)
                msum+=m[-1]
                ns+=1
            else:
                if msum+mtemp <= mprof[rindx] or (mprof[rindx]-msum) <= mlimits[0]:
                    m.append(mtemp)
                    msum+=m[-1]
                    ns+=1
    else:
        m=[]
        for i in range(0,cluster.ntot):
            #Find local alpha
            ydalpha=delta_alpha*np.log(cluster.r[i]/cluster.rm)+alpha50
            m.append(get_imf(mlimits = mlimits, alphas = [ydalpha], ntot=1,do_random=do_random)[0])

    return np.array(m)



