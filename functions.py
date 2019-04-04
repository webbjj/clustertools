"""
FUNCTIONS

Designed to accept StarCluster instance as in put to
calculate key parameters

"""
import numpy as np
from galpy.util import bovy_coords

from constants import *
from recipes import *
from operations import *
from plots import *

def relaxation_time(cluster,local=False,multimass=True,projected=False):
    """
    NAME:

       relaxation_time

    PURPOSE:

       Calculate the relaxation time of the cluster (Default is half mass relaxation time)

    INPUT:

       cluster - StarCluster insane

       local - calcuate relaxation time at each lagrange bin (True) or the half mass radius (False) (default: False)

       multimass - use multimass (True) or single mass (False) value for ln lambda (default: True)

       projected - use projected values (default: False)

    OUTPUT:

       trelax

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 
    #Relaxation time - Spitzer (1987)
    #Default is to be calculated at the half-mass radius, but can be computed locally
    #multimass determines value of lnlambda

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

    cluster.trh=trelax

    return trelax

def energies(cluster,specific=True,do_batch=True,i_d=None):
    #Calculate energy on each star or just one star (i_d)
    units0,origin0,center0=save_cluster(cluster)
    cluster.to_center()

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

        pottot=np.sum(gmr)
        ektot=ek[indx]

    elif do_batch:

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

    else:

        dx=np.repeat(cluster.x,len(cluster.x))-np.tile(cluster.x,len(cluster.x))
        dy=np.repeat(cluster.y,len(cluster.x))-np.tile(cluster.y,len(cluster.x))
        dz=np.repeat(cluster.z,len(cluster.x))-np.tile(cluster.z,len(cluster.x))

        if specific:
            m=np.repeat(cluster.m,len(cluster.m))
        else:
            m=np.repeat(cluster.m,len(cluster.m))*np.tile(cluster.m,len(cluster.m))

        dr=np.sqrt(dx**2.+dy**2.+dz**2.)
        gmr=np.reshape(-grav*m/dr,(len(cluster.r),len(cluster.r)))

        indx=np.logical_or(np.isinf(gmr),np.isnan(gmr))
        gmr[indx]=0.0

        pot=np.sum(gmr,axis=1)


    etot=ek+pot

    cluster.add_energies(ek,pot,etot)

    return_cluster(cluster,units0,origin0,center0)

    return ek,pot,etot

def virialize(cluster,specific=True,do_batch=True):
    units0,origin0,center0=save_cluster(cluster)
    cluster.to_center()

    try:
        qv=np.sqrt(np.abs(0.5/cluster.qvir))
    except:
        print('NEED TO CALCULATE ENERGIES FIRST')
        energies(cluster,specific=specific,do_batch=do_batch)
        qv=np.sqrt(np.abs(0.5/cluster.qvir))

    print('QV = ',qv)

    cluster.vx*=qv
    cluster.vy*=qv
    cluster.vz*=qv
    cluster.key_params()

    energies(cluster,specific=specific,do_batch=do_batch)
    qv=np.sqrt(np.abs(0.5/cluster.qvir))

    return_cluster(cluster,units0,origin0,center0)

    return qv

def rlagrange(cluster,nlagrange=10,projected=False):
    #Calculate lagrange radii
    #Units will be whatever units the cluster is currently in
    
    units0,origin0,center0=save_cluster(cluster)
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

    cluster.rn=rn

    return_cluster(cluster,units0,origin0,center0)

    return rn

def rvirial(cluster,H=70.0,Om=0.3,overdens=200.,nrad=10,projected=False):
    #Calculate the virial radius of the cluster - 
    # radius at which the density is equal to the critical density of the Universe at the redshift of the system, multiplied by an overdensity constant

    units0,origin0,center0=save_cluster(cluster)
    cluster.to_realpc()
    cluster.to_center()
 
    H/=(1000000.0) #(km/s) / pc
    Grav=4.302e-3 #pc (km/s)^2 / Msun
    
    rhocrit=3.0*(H**2.)/(8.0*np.pi*Grav) # Msun/pc^3
 

    indx=cluster.rorder
    msum=np.cumsum(cluster.m[indx])
    vsum=(4./3.)*np.pi*(cluster.r[indx]**3.)
    pprof=msum/vsum


    #Find radius where maxium density occurs
    rindx=np.argwhere(pprof==np.max(pprof))[0][0]
    rmax=rprof[rindx]

    print('RMAX AT ',rindx,rmax)

    indx1=(rprof > rmax) * (pprof < rhocrit*overdens)
    indx2=(rprof > rmax) * (pprof > rhocrit*overdens)

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

    cluster.rv=r_v

    return_cluster(cluster,units0,origin0,center0)
        
    return r_v


def mass_function(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
    #Find mass function over a given mass range using nmass bins containing an equal number of stars
    #Radial range optional
    #Stellar evolution range (kw type) optional (default is MS stars)
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

    if indx==None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

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

    cluster.alpha=alpha

    return m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

#Find eta over a given mass range using nmass bins containing an equal number of stars
#Mass range optional (default is between 0.1 and 1.8 Msun)
#Radial range optional
#Stellar evolution range optional
def eta_function(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
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

    if indx==None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

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

    cluster.eta=eta

    return m_mean,sigvm,eta,eeta,yeta,eyeta
