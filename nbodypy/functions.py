"""
FUNCTIONS

Designed to accept StarCluster instance as in put to
calculate key parameters

"""
import numpy as np
import numba
from galpy.util import bovy_coords

from ..util.constants import *
from ..util.recipes import *
from .operations import *
from ..util.plots import *

def relaxation_time(cluster,local=False,multimass=True,projected=False):
    """
    NAME:

       relaxation_time

    PURPOSE:

       Calculate the relaxation time (Spitzer 1958) of the cluster (Default is half mass relaxation time)

    INPUT:

       cluster - StarCluster instance

       local - calcuate relaxation time at each lagrange bin (True) or the half mass radius (False) (default: False)

       multimass - use multimass (True) or single mass (False) value for ln lambda (default: True)

       projected - use projected values (default: False)

    OUTPUT:

       trelax

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 
    #Gravitational Constant (pc km/s^2 / Msun)
    grav=4.302E-3
    
    #Find Density within Half-Mass Radius
    if local:
        if projected:
            vol=(4.0/3.0)*np.pi*(np.max(cluster.rpro)**3.0-np.min(cluster.rpro)**3.0)
        else:
            vol=(4.0/3.0)*np.pi*(np.max(cluster.r)**3.0-np.min(cluster.r)**3.0)
        p50=cluster.mtot/vol
    else:
        if projected:
            p50=3.0*(0.5*cluster.mtot)/(4.0*np.pi*(cluster.rmpro**3.0))
        else:
            p50=3.0*(0.5*cluster.mtot)/(4.0*np.pi*(cluster.rm**3.0))

    #Find 1D Global Velocity Dispersion
    if projected:
        sigv_1d=np.std(cluster.vx)
    else:
        sigv_1d=np.sqrt((np.std(cluster.vx)**2.0+np.std(cluster.vy)**2.0+np.std(cluster.vz)**2.0)/3.0)

    #Mean stellar mass
    mbar=np.mean(cluster.m)
    
    if multimass:
        lnlambda=np.log(0.02*cluster.ntot)
    else:
        lnlambda=np.log(0.11*cluster.ntot)

    #Units of seconds * (pc/km)
    trelax=0.34*(sigv_1d**3.0)/((grav**2.0)*mbar*p50*lnlambda)

    #Units of Myr
    trelax=trelax*3.086e13/(3600.0*24.0*365.0*1000000.0)

    return trelax

def energies(cluster,specific=True,i_d=None,full=True,parallel=False):
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

       parllel - calculate distances in parallel if True (default: False)

    OUTPUT:

       ek,pot,etot

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    if cluster.units=='nbody':
        grav=1.0
    elif cluster.units=='realpc':
        #G has units of pc (km/s)^2 / Msun
        grav=4.302e-3
    elif cluster.units=='realkpc':
        #G has units of kpc (km/s)^2 / Msun
        grav=4.302e-6
    else:
        grav=1.0
    
    if specific:
        ek=0.5*(cluster.v**2.0)
    else:
        ek=0.5*cluster.m*(cluster.v**2.0)

    if i_d!=None:
        indx=(cluster.id==i_d)

        dx=cluster.x[indx]-cluster.x
        dy=cluster.y[indx]-cluster.y
        dz=cluster.z[indx]-cluster.z

        if specific:
            m=cluster.m
        else:
            m=cluter.m[indx]*cluster.m

        dr=np.sqrt(dx**2.+dy**2.+dz**2.)
        rindx=(dr!=0.0)
        gmr=-grav*m[rindx]/dr[rindx]

        pot=np.sum(gmr)
        ek=ek[indx]
        etot=ek+pot

    elif full:
        x=np.array([cluster.x,cluster.y,cluster.z,cluster.m]).T
        if parallel:
            pot=grav*np.array(potential_energy_parallel(x))
        else:
            pot=grav*np.array(potential_energy(x))


        if specific:
            pot/=cluster.m

        etot=ek+pot
        cluster.add_energies(ek,pot,etot)
    else:
        pot=[]

        for i in range(0,cluster.ntot):
            dx=cluster.x[i]-cluster.x
            dy=cluster.y[i]-cluster.y
            dz=cluster.z[i]-cluster.z
            if specific:
                m=cluster.m
            else:
                m=cluter.m[i]*cluster.m

            dr=np.sqrt(dx**2.+dy**2.+dz**2.)
            indx=(dr!=0.0)
            gmr=-grav*m[indx]/dr[indx]

            pot.append(np.sum(gmr))

        etot=ek+pot
        cluster.add_energies(ek,pot,etot)

    return_cluster(cluster,units0,origin0)

    return ek,pot,etot

@numba.njit
def distance(star1, star2):
    """
    NAME:

       distance

    PURPOSE:

       Find distance between two stars (made for use with numba)

    INPUT:

       star1=[x,y,z]
       star2=[x,y,z]

    OUTPUT:

       distance

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    dx = star2[0] - star1[0]
    dy = star2[1] - star1[1]
    dz = star2[2] - star1[2]

    r = (dx * dx + dy * dy + dz * dz) ** 0.5

    return r

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
    energy = [0.]*len(cluster)
    for i in range(len(cluster)-1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i],cluster[j])
            m2=cluster[i,3]*cluster[j,3]
            energy[i]+=-m2/r
            energy[j]+=-m2/r

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

    energy = [0.]*len(cluster)
    for i in numba.prange(len(cluster)-1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i],cluster[j])
            m2=cluster[i,3]*cluster[j,3]
            energy[i]+=-m2/r
            energy[j]+=-m2/r

    return energy

def closest_star(cluster,projected=False):

    if projected:
        z=np.zeros(cluster.ntot)
        x=np.array([cluster.x,cluster.y,z]).T
    else:
        x=np.array([cluster.x,cluster.y,cluster.z]).T
    return minimum_distance(x)

@numba.njit
def minimum_distance(cluster):
    """
    NAME:

       stellar_distances

    PURPOSE:

       Find distances between each star

    INPUT:

       cluster=[x,y,z].T

    OUTPUT:

       distances

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    min_distance=[-1.]*len(cluster)
    for i in range(len(cluster)-1):
        for j in range(i + 1, len(cluster)):
            r = distance(cluster[i],cluster[j])
            if min_distance[i]<0:
                min_distance[i]=r
            else:
                min_distance[i]=np.minimum(min_distance[i],r)

            if min_distance[j]<0:
                min_distance[j]=r
            else:
                min_distance[j]=np.minimum(min_distance[j],r)              

    return min_distance


def virialize(cluster,specific=True,full=True):
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

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    try:
        qv=np.sqrt(np.abs(0.5/cluster.qvir))
    except:
        print('NEED TO CALCULATE ENERGIES FIRST')
        energies(cluster,specific=specific,full=full)
        qv=np.sqrt(np.abs(0.5/cluster.qvir))

    print('QV = ',qv)

    cluster.vx*=qv
    cluster.vy*=qv
    cluster.vz*=qv
    cluster.key_params()

    energies(cluster,specific=specific,full=full)
    qv=np.sqrt(np.abs(0.5/cluster.qvir))

    return_cluster(cluster,units0,origin0)

    return qv

def rlagrange(cluster,nlagrange=10,projected=False):
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
    
    units0,origin0=save_cluster(cluster)
    cluster.to_cluster()

    #Radially order the stars
    msum=0.0
    nfrac=1
    rn=[]
    
    if projected:
        rorder=sorted(range(0,self.ntot),key=lambda k:self.rpro[k])
    else:
        rorder=cluster.rorder
    
    for i in range(0,cluster.ntot):
        mfrac=cluster.mtot*float(nfrac)/float(nlagrange)
        msum+=cluster.m[rorder[i]]
        if msum >= mfrac:
            rn.append(cluster.r[rorder[i]])
            nfrac+=1
        if nfrac>nlagrange:
            break
    if len(rn) != nlagrange:
        rn.append(np.max(cluster.r))

    return_cluster(cluster,units0,origin0)

    return rn

def rvirial(cluster,H=70.0,Om=0.3,overdens=200.,nrad=20,projected=False,plot=False,**kwargs):
    """
    NAME:

       rvirial

    PURPOSE:

       Calculate virial radius of the cluster
       --> Note: virial radius is the radius at which the density is equal to the critical 
           density of the Universe at the redshift of the system, multiplied by an overdensity constant

    INPUT:

       cluster - StarCluster instance
       H - Hubble constant
       Om - density of matter
       overdens - overdensity constant
       nrad - number of radial bins used to calculate cluster density profile
       projected - calculate projected lagrange radii (default: False)

    OUTPUT:

       rv

    HISTORY:

       2019 - Written - Webb (UofT)
    """    

    units0,origin0=save_cluster(cluster)
    cluster.to_realpc()
    cluster.to_centre()
 
    H/=(1000000.0) #(km/s) / pc
    Grav=4.302e-3 #pc (km/s)^2 / Msun
    
    rhocrit=3.0*(H**2.)/(8.0*np.pi*Grav) # Msun/pc^3
 

    indx=cluster.rorder
    msum=np.cumsum(cluster.m[indx])
    vsum=(4./3.)*np.pi*(cluster.r[indx]**3.)
    pprof=msum/vsum
    rprof=cluster.r[indx]


    #Find radius where maxium density occurs
    rindx=np.argmax(pprof)
    rmax=rprof[rindx]

    print('RMAX AT ',rindx,rmax)

    indx1=(rprof > rmax) * (pprof > rhocrit*overdens)
    indx2=(rprof > rmax) * (pprof < rhocrit*overdens)

    if np.sum(indx2)==0.:
        print('SYSTEM IS NOT VIRIALIZED')
        r_v=-1.
    else:
        r1=rprof[indx1][-1]
        r2=rprof[indx2][0]

        rho1=pprof[indx1][-1]
        rho2=pprof[indx2][0]    

        print(r1,r2,rho1,rho2,rhocrit*overdens)
        r_v=interpolate([r1,rho1],[r2,rho2],y=rhocrit*overdens)

    if plot:
        print('OVERDENSITY = ',rhocrit*overdens)    
        rho_local=rhocrit*overdens

        filename=kwargs.pop('filename',None)   
        overplot=kwargs.pop('overplot',False)        
     
        xunits=' (pc)'
        if projected:
            yunits=' Msun/pc^2'
        else:
            yunits=' Msun/pc^3'

        x,y=rprof,pprof
        nlplot(x,y,xlabel=r'$R'+xunits+'$',ylabel=r'$rho'+yunits+'$',title='Time = %f' % cluster.tphys,log=True,overplot=overplot,filename=filename)
        nlplot(x,np.ones(len(x))*rho_local,'--',overplot=True)
        nlplot(np.ones(len(y))*r_v,y,'--',overplot=True)

        if filename!=None:
            plt.savefig(filename)

    return_cluster(cluster,units0,origin0)
        
    return r_v


def mass_function(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
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
       obs_cut - trim data to match an observed dataset (WIP)
       plot - plot the mass function
       **kwargs - key words for plotting

    OUTPUT:

       m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    if np.sum(indx) >= nmass:

        m_lower,m_mean,m_upper,m_hist=nbinmaker(cluster.m[indx],nmass)
       
        lm_mean=np.log10(m_mean)
        dm=m_hist/(m_upper-m_lower)
        ldm=np.log10(dm)

        (alpha,yalpha),V=np.polyfit(lm_mean,ldm,1,cov=True)
        ealpha=np.sqrt(V[0][0])
        eyalpha=np.sqrt(V[1][1])

        if plot:
            filename=kwargs.get('filename',None)
            nplot(m_mean,np.log10(dm),xlabel='M',ylabel='LOG(dN/dM)',**kwargs)
            mfit=np.linspace(np.min(m_mean),np.max(m_mean),nmass)
            dmfit=10.0**(alpha*np.log10(mfit)+yalpha)
            nlplot(mfit,np.log10(dmfit),overplot=True,label=(r'$\alpha$ = %f' % alpha))

            plt.legend()

            if filename!=None:
                plt.savefig(filename)


        return m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha
    else:
        print('NOT ENOUGH STARS TO ESTIMATE MASS FUNCTION')
        return np.zeros(nmass),np.zeros(nmass),np.zeros(nmass),-1000.,-1000.,-1000.,-1000.

def eta_function(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
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

       obs_cut - trim data to match an observed dataset (WIP)

       plot - plot the mass function
       
       **kwargs - key words for plotting

    OUTPUT:

       m_mean,sigvm,eta,eeta,yeta,eyeta

    HISTORY:

       2018 - Written - Webb (UofT)
    """ 
    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    if np.sum(indx)>=2*nmass:

        m_lower,m_mean,m_upper,m_hist=nbinmaker(cluster.m[indx],nmass)
        lm_mean=np.log10(m_mean)

        sigvm=[]
        lsigvm=[]
        for i in range(0,nmass):

            mindx=indx * (cluster.m >=m_lower[i]) * (cluster.m<=m_upper[i])
            sigvm.append(np.std(v[mindx]))
            lsigvm.append(np.log10(sigvm[-1]))


        (eta,yeta),V=np.polyfit(lm_mean,lsigvm,1,cov=True)
        eeta=np.sqrt(V[0][0])
        eyeta=np.sqrt(V[1][1])

        if plot:
            filename=kwargs.get('filename',None)
            nplot(m_mean,np.log10(sigvm),xlabel='M',ylabel=r'$\sigma_v$',**kwargs)
            mfit=np.linspace(np.min(m_mean),np.max(m_mean),nmass)
            sigfit=10.0**(eta*np.log10(mfit)+yeta)
            nlplot(mfit,np.log10(sigfit),overplot=True,label=(r'$\eta$ = %f' % eta))
            plt.legend()

            if filename!=None:
                plt.savefig(filename)

        return m_mean,sigvm,eta,eeta,yeta,eyeta
    else:
        print('NOT ENOUGH STARS TO ESTIMATE SIGMA-MASS RELATION')
        return np.zeros(nmass),np.zeros(nmass),np.zeros(nmass),-1000.,-1000.,-1000.,-1000.
