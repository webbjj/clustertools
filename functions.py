#Key calculations (constantly adding to this file)
import math
import numpy as np
from cluster import *
from observations import *
from operations import *
from constants import *
from galpy.util import bovy_coords
from coordinates import *
from recipes import *

#Relaxation time - Spitzer (1987)
def relaxation_time(cluster,local=False,multimass=True):
    #Gravitational Constant (pc km/s^2 / Msun)
    grav=4.302E-3
    #Find Density within Half-Mass Radius
    if local:
        vol=(4.0/3.0)*math.pi*(np.max(cluster.r)**3.0-np.min(cluster.r)**3.0)
        p50=cluster.mtot/vol
    else:
        p50=3.0*(0.5*cluster.mtot)/(4.0*math.pi*(cluster.rm**3.0))

    #Find 1D Global Velocity Dispersion
    sigv_1d=math.sqrt((np.std(cluster.vx)**2.0+np.std(cluster.vy)**2.0+np.std(cluster.vz)**2.0)/3.0)
    #Mean stellar mass
    mbar=np.mean(cluster.m)
    
    if multimass:
        lnlambda=math.log(0.02*cluster.ntot)
    else:
        lnlambda=math.log(0.11*cluster.ntot)


    #Units of seconds * (pc/km)
    trelax=0.34*(sigv_1d**3.0)/((grav**2.0)*mbar*p50*lnlambda)

    #Units of Myr
    trelax=trelax*3.086e13/(3600.0*24.0*365.0*1000000.0)

    return trelax

def energies(cluster,specific=True,do_batch=True):
    origin0=cluster.origin
    center0=cluster.center

    
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

    if origin0!='cluster':
        cluster.to_cluster()
    if not center0:
        cluster.to_center()
    
    if specific:
        ek=0.5*(cluster.v**2.0)
    else:
        ek=0.5*cluster.m*(cluster.v**2.0)

    ektot=np.sum(ek)

    if do_batch:

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

    pottot=np.sum(pot)
    etot=ek+pot

    cluster.add_energies(ek,pot,etot)

    if not center0:
        cluster.from_center()
    if origin0!='cluster':
        cluster.to_galaxy()

    return ektot,pottot

#Calculate lagrange radii
#Units will be whatever units the cluster is currently in
def rlagrange(cluster,nlagrange=10):
    
    if cluster.origin!='cluster':
        print('CANNOT COMPUTE LAGRANGE RADII IF CLUSTER.ORIGIN=',cluster.origin)
        return 0

    #Radially order the stars
    msum=0.0
    nfrac=1
    rn=[]
    for i in range(0,cluster.ntot):
        mfrac=cluster.mtot*float(nfrac)/float(nlagrange)
        msum+=cluster.m[cluster.rorder[i]]
        if msum >= mfrac:
            rn.append(cluster.r[cluster.rorder[i]])
            nfrac+=1
        if nfrac>nlagrange:
            break
    if len(rn) != nlagrange:
        rn.append(np.max(cluster.r))

    cluster.rn=rn

    return rn

#Find mass function over a given mass range using nmass bins containing an equal number of stars
#Radial range optional
#Stellar evolution range (kw type) optional (default is MS stars)
def mass_function(cluster,mmin=0.1,mmax=1.0,nmass=10,rmin=None,rmax=None,se_min=0,se_max=1,projected=False):
    
    if projected:
        rad=cluster.rpro
    else:
        rad=cluster.r

    indx=(cluster.m>=mmin) * (cluster.m<=mmax) * (cluster.kw>=se_min) * (cluster.kw<=se_max)
    if rmin!=None:
        indx=indx * (rad >= rmin) * (rad <= rmax)

    msub=cluster.m[indx]

    if len(msub)>2.0*nmass:
        m_lower,m_mean,m_upper,m_hist=nbinmaker(msub,nmass)
    else:
        m_mean=np.zeros(nmass)
        m_hist=np.zeros(nmass)
        dm=np.zeros(nmass)
        alpha=-100.0
        ealpha=0.0
        yalpha=0.0
        eyalpha=0.0
        return m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

    indx=(m_upper!=m_lower)
    lm_mean=math.log(m_mean[indx],10.0)
    dm=m_hist[indx]/(m_upper[indx]-m_lower[indx])
    ldm=math.log(dm[indx],10.0)
    m_mean=m_mean[indx]
    m_hist=m_hist[indx]

    if len(lm_mean)>=3:
        (alpha,yalpha),V=np.polyfit(lm_mean,ldm,1,cov=True)
        ealpha=np.sqrt(V[0][0])
        eyalpha=np.sqrt(V[1][1])
    else:
        alpha=-100.0
        yalpha=0.0
        ealpha=0.0
        eyalpha=0.0

    if len(m_mean)!=nmass:
        for i in range(0,nmass-len(m_mean)):
            m_mean.append(0.0)
            m_hist.append(0.0)
            dm.append(0.0)

    cluster.alpha=alpha
    cluster.ealpha=ealpha

    return m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

#Measure the radial variation in the stellar mass function
#Mass range optional
#Radial range optional
#Stellar evolution range (kw type) optional (default is MS stars)
def alpha_prof(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,nrad=10,se_min=0,se_max=1,projected=False,obs_cut=None):

    stars=[]
    rprof=[]
    lrprofn=[]
    nskip=0

    if rmin==None: rmin=np.min(cluster.r)
    if rmax==None: rmax=np.max(cluster.r)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    #Build subcluster containing only stars in the full radial and mass range:
    subcluster=sub_cluster(cluster,rmin,rmax,mmin,mmax,se_min,se_max,projected=projected)

    #Make radial bins
    if obs_cut!=None and projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.rpro,cluster.rmpro,obs_cut)
        nrad=len(r_mean)
    elif obs_cut!=None and not projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.r,cluster.rm,obs_cut)
        nrad=len(r_mean)
    elif projected:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.rpro,nrad)
    else:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.r,nrad)


    rcluster=sub_cluster(subcluster,rmin=r_lower[i],rmax=r_upper[i],projected=projected)
    m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(rcluster,mmin,mmax,nmass)
     
    indx=(alpha!=-100.0)
    rprof=(r_mean[indx])
    if projected:
        lrprofn.append(math.log(rprof[indx]/cluster.rmpro))
    else:
        lrprofn.append(math.log(rprof[indx]/cluster.rm))

    if len(lrprofn)>3:
        (dalpha,ydalpha),V=np.polyfit(lrprofn,alpha,1,cov=True)
        edalpha=np.sqrt(V[0][0])
        eydalpha=np.sqrt(V[1][1])
    else:
        dalpha=-100.0
        ydalpha=0.0
        edalpha=0.0
        eydalpha=0.0

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
            aprof.append(0.0)
    
    cluster.dalpha=dalpha
    cluster.edalpha=dalpha

    return lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha

#Find eta over a given mass range using nmass bins containing an equal number of stars
#Mass range optional (default is between 0.1 and 1.8 Msun)
#Radial range optional
#Stellar evolution range optional
def eta_function(cluster,mmin=0.1,mmax=1.8,nmass=10,rmin=None,rmax=None,se_min=0,se_max=100,projected=False):

    if projected:
        rad=cluster.rpro
    else:
        rad=cluster.r

    indx=(cluster.m>=mmin) * (cluster.m<=mmax) * (cluster.kw>=se_min) * (cluster.kw<=se_max)
    if rmin!=None:
        indx=indx * (rad >= rmin) * (rad <= rmax)
    
    msub=cluster.m[indx]
    m_lower,m_mean,m_upper,m_hist=nbinmaker(msub,nmass)


    lm_mean=math.log(m_mean,10.0)

    sigvm=[]
    lsigvm=[]
    for i in range(0,nmass):
        indx= (cluster.m >=m_lower[i]) * (cluster.m<=m_upper) * (cluster.kw>=se_min) * (cluster.kw<=se_max)
        if rmin!=None:
            indx=indx * (cluster.r[i]>=rmin) * (cluster.r[i]<=rmax)

        sigx=np.std(cluster.vx[indx])
        sigy=np.std(cluster.vy[indx])
        sigz=np.std(cluster.vz[indx])

        sigvm.append(math.sqrt(sigvx**2.0+sigvy**2.0+sigvz**2.0))
        lsigvm.append(math.log(sigvm,10.0))


    (eta,yeta),V=np.polyfit(lm_mean,lsigvm,1,cov=True)
    eeta=np.sqrt(V[0][0])
    eyeta=np.sqrt(V[1][1])

    return m_mean,sigvm,eta,eeta,yeta,eyeta

#Measure the radial variation in the velocity dispersion
#Mass range optional
#Radial range optional
#Stellar evolution range (kw type) optional (default is all stars)
def sigv_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=10,se_min=0,se_max=15,projected=False,obs_cut=None):

    stars=[]
    rprof=[]
    lrprofn=[]
    sigvprof=[]
    nskip=0

    if rmin==None: rmin=np.min(cluster.r)
    if rmax==None: rmax=np.max(cluster.r)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    #Build subcluster containing only stars in the full radial and mass range:
    subcluster=sub_cluster(cluster,rmin,rmax,mmin,mmax,se_min,se_max,projected=projected)

    #Make radial bins
    if obs_cut!=None and projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.rpro,cluster.rmpro,obs_cut)
        nrad=len(r_mean)
    elif obs_cut!=None and not projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.r,cluster.rm,obs_cut)
        nrad=len(r_mean)
    elif projected:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.rpro,nrad)
    else:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.r,nrad)



    rcluster=sub_cluster(subcluster,rmin=r_lower,rmax=r_upper,projected=projected)
    sigx,sigy,sigz,sigv,lsigv=sigv_function(rcluster)

    if projected:
        lrprofn=math.log(rmean/cluster.rmpro)
    else:
        lrprofn=math.log(rmean/cluster.rm)

    sigvprof(sigv)

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
            sigv.append(0.0)
    
    return lrprofn,sigv


def beta_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=10,se_min=0,se_max=15,projected=False,obs_cut=None):

    stars=[]
    rprof=[]
    lrprofn=[]
    betaprof=[]
    nskip=0

    if rmin==None: rmin=np.min(cluster.r)
    if rmax==None: rmax=np.max(cluster.r)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    #Build subcluster containing only stars in the full radial and mass range:
    subcluster=sub_cluster(cluster,rmin,rmax,mmin,mmax,se_min,se_max,projected=projected)

    #Make radial bins
    if obs_cut!=None and projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.rpro,cluster.rmpro,obs_cut)
        nrad=len(r_mean)
    elif obs_cut!=None and not projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.r,cluster.rm,obs_cut)
        nrad=len(r_mean)
    elif projected:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.rpro,nrad)
    else:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.r,nrad)


    rcluster=sub_cluster(subcluster,rmin=r_lower,rmax=r_upper,projected=projected)

    sigr,sigt,sigp,sigv,lsigv=sigv_function(rcluster,spherical=True)

    if projected:
        lrprofn=math.log(rprof/cluster.rmpro)
    else:
        lrprofn=math.log(rprof/cluster.rm)

    betaprof=1.0-(sigt**2.0+sigp**2.0)/(2.*(sigr**2.0))

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
    
    return lrprofn,betaprof

def sigv_function(cluster,spherical=False):
    vx=[]
    vy=[]
    vz=[]
    
    if spherical:
        r,theta,phi=rect_to_sphere(cluster.x,cluster.y,cluster.z)
        vx,vy,vz=rect_to_sphere_vel(cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz)

    else:
        vx=cluster.vx
        vy=cluster.vy
        vz=cluster.vz

    sigx=np.std(vx)
    sigy=np.std(vy)
    sigz=np.std(vz)
    
    sigv=(math.sqrt(sigx**2.0+sigy**2.0+sigz**2.0))
    lsigv=(math.log(sigv,10.0))

    return sigx,sigy,sigz,sigv,lsigv

def v_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=10,se_min=0,se_max=15,coord=None,projected=False,obs_cut=None):

    stars=[]
    rprof=[]
    lrprofn=[]
    vprof=[]
    nskip=0

    if rmin==None: rmin=np.min(cluster.r)
    if rmax==None: rmax=np.max(cluster.r)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    #Build subcluster containing only stars in the full radial and mass range:
    subcluster=sub_cluster(cluster,rmin,rmax,mmin,mmax,se_min,se_max,projected=projected)

    #Make radial bins
    if obs_cut!=None and projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.rpro,cluster.rmpro,obs_cut)
        nrad=len(r_mean)
    elif obs_cut!=None and not projected:
        r_lower,r_mean,r_upper,r_hist=obsrbinmaker(subcluster.r,cluster.rm,obs_cut)
        nrad=len(r_mean)
    elif projected:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.rpro,nrad)
    else:
        r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.r,nrad)


    rcluster=sub_cluster(subcluster,rmin=r_lower,rmax=r_upper,projected=projected)

    if coord==None:
        vmean=np.mean(rcluster.v)
    elif coord=='x':
        vmean=np.mean(rcluster.vx)
    elif coord=='y':
        vmean=np.mean(rcluster.vy)
    elif coord=='z':
        vmean=np.mean(rcluster.vz)
    elif coord=='r' or coord=='theta' or coord=='phi':
        vr,vt,vp=rect_to_sphere_vel(cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz)

        if coord=='r':
            vmean=np.mean(vr)
        if coord=='theta':
            vmean=np.mean(vt)
        if coord=='phi':
            vmean=np.mean(vp)

    if projected:
        lrprofn=math.log(r_mean/cluster.rmpro)
    else:
        lrprofn=math.log(r_mean/cluster.rm)

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
            vmean.append(0.0)
    
    return lrprofn,vprof


