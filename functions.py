#Key calculations (constantly adding to this file)
import math
import numpy as np
from cluster import *
from observations import *
from constants import *
from galpy.util import bovy_coords
from coordinates import *
from recipes import *

#Relaxation time - Spitzer (1987)
def Relaxation_Time(cluster,local=False,multimass=True):
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

    #print('DEBUG: ',sigv_1d,mbar,p50,cluster.rm,np.min(cluster.r),np.max(cluster.r),trelax)


    return trelax

#Calculate Force Acting on Specific Star or All Stars
#Output units depends on input units
def get_forces(cluster,id=None):
    
    if id==None:
    #Calculate force acting on each star
        s2sd=np.zeros((cluster.ntot,cluster.ntot))
        for i in range(0,cluster.ntot):
            for j in range(i,cluster.ntot):
                d=math.sqrt((cluster.x[i]-cluster.x[j])**2.0+(cluster.y[i]-cluster.y[j])**2.0+(cluster.z[i]-cluster.z[j])**2.0)
                s2sd[i][j]=d
                s2sd[j][i]=d
    else:
        #Calculate force acting on specific star
        indx=cluster.id.index(id)
        fx=0.0
        fy=0.0
        fz=0.0
        pot=0.0

        for i in range(0,cluster.ntot):
            if cluster.id[i]!=id:       
                dx=cluster.x[indx]-cluster.x[i]
                dy=cluster.y[indx]-cluster.y[i]
                dz=cluster.z[indx]-cluster.z[i]

                if dx==0: dx=cluster.x[indx]-0.99*cluster.x[i]
                if dy==0: dy=cluster.y[indx]-0.99*cluster.y[i]
                if dz==0: dz=cluster.z[indx]-0.99*cluster.z[i]

                dr=math.sqrt(dx**2.0+dy**2.0+dz**2.0)

                xhat=dx/abs(dx)
                yhat=dy/abs(dy)
                zhat=dz/abs(dz)

                fx+=xhat*cluster.m[i]/(dx**2.0)
                fy+=yhat*cluster.m[i]/(dy**2.0)
                fz+=zhat*cluster.m[i]/(dz**2.0)
                pot-=cluster.m[i]/dr

    ek=0.5*math.sqrt((cluster.vx[indx]-np.mean(cluster.vx))**2.0+(cluster.vy[indx]-np.mean(cluster.vy))**2.0+(cluster.vz[indx]-np.mean(cluster.vz))**2.0)
    
    if cluster.units=='nbody':
        grav=1.0
        fscale=1.0
    elif cluster.units=='realpc':
        #G has units of pc (km/s)^2 / Msun
        grav=4.302e-3
        #Scale forces to units of (km/s) / Myr
        fscale=spermyr/kmperpc
    elif cluster.units=='realkpc':
        #G has units of kpc (km/s)^2 / Msun
        grav=4.302e-6
        #Scale forces to units of (km/s) / Myr
        fscale=spermyr/kmperkpc
    else:
        grav=1.0
        fscale=1.0

    fx=fx*grav*fscale
    fy=fy*grav*fscale
    fz=fz*grav*fscale
    pot=pot*grav

    return fx,fy,fz,pot,ek

def get_energies(cluster,specific=True):
    
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
    ektot=np.sum(ek)

    pot=np.zeros(cluster.ntot)
    etot=[]
    pottot=0.0

    for i in range(0,cluster.ntot):
        dx=cluster.x-cluster.x[i]
        dy=cluster.y-cluster.y[i]
        dz=cluster.z-cluster.z[i]
        dr=np.sqrt(dx**2.0+dy**2.0+dz**2.0)
        indx=dr>0

        if specific:
            pot[i]=np.sum(-grav*cluster.m[indx]/dr[indx])
        else:
            pot[i]+=np.sum(-grav*cluster.m[indx]/dr[indx])*cluster.m[i]

    pottot=np.sum(pot)
    etot=ek+pot
    cluster.add_energies(ek,pot,etot)

    return pottot,ektot


#Calculate lagrange radii
#Units will be whatever units the cluster is currently in
def rlagrange(cluster,nlagrange=10):
    
    #Shift to clustercentric origin if not so already:
    origin0=cluster.origin
    units0=cluster.units
    if origin0 !='cluster':
        xvshift(cluster,-cluster.xgc,-cluster.ygc,-cluster.zgc,-cluster.vxgc,-cluster.vygc,-cluster.vzgc,'cluster')

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

    #Shift back to original origin
    if origin0 == 'galaxy' and cluster.origin=='cluster':
        xvshift(cluster,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,'galaxy')

    return rn

#Find mass function over a given mass range using nmass bins containing an equal number of stars
#Radial range optional
#Stellar evolution range (kw type) optional (default is MS stars)
def mass_function(cluster,mmin=0.1,mmax=1.0,nmass=10,rmin=None,rmax=None,se_min=0,se_max=1,projected=False):

    msub=[]
    
    for i in range(0,cluster.ntot):
        if projected:
            rad=cluster.rpro[i]
        else:
            rad=cluster.r[i]
        if cluster.m[i]>=mmin and cluster.m[i]<=mmax and cluster.kw[i]>=se_min and cluster.kw[i]<=se_max:
            if rmin!=None:
                if rad>=rmin and rad<=rmax:
                    msub.append(cluster.m[i])
            else:
                msub.append(cluster.m[i])

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

    lm_mean=[]
    dm=[]
    ldm=[]
    m_mean_temp=[]
    m_hist_temp=[]
    #print('DEBUG ',m_lower,m_mean,m_upper,m_hist)
    for i in range(0,len(m_mean)):
        if m_upper[i]!=m_lower[i]:
            lm_mean.append(math.log(m_mean[i],10.0))
            dm.append(m_hist[i]/(m_upper[i]-m_lower[i]))
            ldm.append(math.log(dm[-1],10.0))
            m_mean_temp.append(m_mean[i])
            m_hist_temp.append(m_hist[i])

    m_mean=m_mean_temp
    m_hist=m_hist_temp

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



    return m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

#Measure the radial variation in the stellar mass function
#Mass range optional
#Radial range optional
#Stellar evolution range (kw type) optional (default is MS stars)
def alpha_prof(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,nrad=10,se_min=0,se_max=1,projected=False,obs_cut=None):

    stars=[]
    rprof=[]
    lrprofn=[]
    aprof=[]
    eaprof=[]
    nskip=0

    if rmin==None: rmin=np.min(cluster.r)
    if rmax==None: rmax=np.max(cluster.r)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    #Build subcluster containing only stars in the full radial and mass range:
    subcluster=sub_cluster(cluster,rmin,rmax,mmin,mmax,se_min,se_max,projected=projected)


    if not subcluster.keyparams:
        subcluster.keyparams=True
        subcluster.key_params()
    else:
        subcluster.key_params()

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

    for i in range(0,len(r_lower)):

        rcluster=sub_cluster(subcluster,rmin=r_lower[i],rmax=r_upper[i],projected=projected)

        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(rcluster,mmin,mmax,nmass)
     
        if alpha!=-100.0:
            rprof.append(r_mean[i])
            if projected:
                lrprofn.append(math.log(rprof[-1]/cluster.rmpro))
            else:
                lrprofn.append(math.log(rprof[-1]/cluster.rm))
            aprof.append(alpha)
            eaprof.append(ealpha)

    if len(lrprofn)>=3:
        (dalpha,ydalpha),V=np.polyfit(lrprofn,aprof,1,cov=True)
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
    
    return lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha

#Find eta over a given mass range using nmass bins containing an equal number of stars
#Mass range optional (default is between 0.1 and 1.8 Msun)
#Radial range optional
#Stellar evolution range optional
def eta_function(cluster,mmin=0.1,mmax=1.8,nmass=10,rmin=None,rmax=None,se_min=0,se_max=100,projected=False):

    msub=[]
    
    for i in range(0,cluster.ntot):
        if cluster.m[i]>=mmin and cluster.m[i]<=mmax and cluster.kw[i]>=se_min and cluster.kw[i]<=se_max:
            if rmin!=None:
                if cluster.r[i]>=rmin and cluster.r[i]<=rmax:
                    msub.append(cluster.m[i])
            else:
                msub.append(cluster.m[i])

    m_lower,m_mean,m_upper,m_hist=nbinmaker(msub,nmass)

    lm_mean=[]
    sigvm=[]
    lsigvm=[]
    for i in range(0,nmass):
        lm_mean.append(math.log(m_mean[i],10.0))
        vx=[]
        vy=[]
        vz=[]

        for j in range(0,cluster.ntot):
            if cluster.m[j]>=m_lower[i] and cluster.m[j]<=m_upper[i] and cluster.kw[i]>=se_min and cluster.kw[i]<=se_max:
                if rmin!=None:
                    if cluster.r[i]>=rmin and cluster.r[i]<=rmax:
                        vx.append(cluster.vx)
                        vy.append(cluster.vy)
                        vz.append(cluster.vz)
                else:
                    vx.append(cluster.vx)
                    vy.append(cluster.vy)
                    vz.append(cluster.vz)

        sigx=np.std(vx)
        sigy=np.std(vy)
        sigz=np.std(vz)

        sigvm.append(math.sqrt(sigvx**2.0+sigvy**2.0+sigvz**2.0))
        lsigvm.append(math.log(sigvm[-1],10.0))


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


    if not subcluster.keyparams:
        subcluster.keyparams=True
        subcluster.key_params()
    else:
        subcluster.key_params()

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

    for i in range(0,len(r_lower)):

        rcluster=sub_cluster(subcluster,rmin=r_lower[i],rmax=r_upper[i],projected=projected)

        sigx,sigy,sigz,sigv,lsigv=sigv_function(rcluster)
     
        rprof.append(r_mean[i])
        if projected:
            lrprofn.append(math.log(rprof[-1]/cluster.rmpro))
        else:
                lrprofn.append(math.log(rprof[-1]/cluster.rm))
        sigvprof.append(sigv)

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
            sigvprof.append(0.0)
    
    return lrprofn,sigvprof


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


    if not subcluster.keyparams:
        subcluster.keyparams=True
        subcluster.key_params()
    else:
        subcluster.key_params()

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

    for i in range(0,len(r_lower)):

        rcluster=sub_cluster(subcluster,rmin=r_lower[i],rmax=r_upper[i],projected=projected)

        sigr,sigt,sigp,sigv,lsigv=sigv_function(rcluster,spherical=True)
     
        rprof.append(r_mean[i])
        if projected:
            lrprofn.append(math.log(rprof[-1]/cluster.rmpro))
        else:
            lrprofn.append(math.log(rprof[-1]/cluster.rm))

        betaprof.append(1.0-(sigt**2.0+sigp**2.0)/(2.*(sigr**2.0)))

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
    
    return lrprofn,betaprof

def sigv_function(cluster,spherical=False):
    vx=[]
    vy=[]
    vz=[]
    
    for j in range(0,cluster.ntot):
        
        if spherical:
            r,theta,phi=rect_to_sphere(cluster.x[j],cluster.y[j],cluster.z[j])
            vr,vt,vp=rect_to_sphere_vel(cluster.x[j],cluster.y[j],cluster.z[j],cluster.vx[j],cluster.vy[j],cluster.vz[j])
            vx.append(vr)
            vy.append(vt)
            vz.append(vp)
        else:
            vx.append(cluster.vx[j])
            vy.append(cluster.vy[j])
            vz.append(cluster.vz[j])

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


    if not subcluster.keyparams:
        subcluster.keyparams=True
        subcluster.key_params()
    else:
        subcluster.key_params()

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

    for i in range(0,len(r_lower)):

        rcluster=sub_cluster(subcluster,rmin=r_lower[i],rmax=r_upper[i],projected=projected)

        if coord==None:
            vmean=np.mean(rcluster.v)
        elif coord=='x':
            vmean=np.mean(rcluster.vx)
        elif coord=='y':
            vmean=np.mean(rcluster.vy)
        elif coord=='z':
            vmean=np.mean(rcluster.vz)
        elif coord=='r' or coord=='theta' or coord=='phi':
            vr=[]
            vt=[]
            vp=[]
            for j in range(0,rcluster.ntot):
                vrad,vtheta,vphi=rect_to_sphere_vel(cluster.x[j],cluster.y[j],cluster.z[j],cluster.vx[j],cluster.vy[j],cluster.vz[j])
                vr.append(vrad)
                vt.append(vtheta)
                vp.append(vphi)
            if coord=='r':
                vmean=np.mean(vr)
            if coord=='theta':
                vmean=np.mean(vt)
            if coord=='phi':
                vmean=np.mean(vp)

        rprof.append(r_mean[i])
        if projected:
            lrprofn.append(math.log(rprof[-1]/cluster.rmpro))
        else:
                lrprofn.append(math.log(rprof[-1]/cluster.rm))
        
        vprof.append(vmean)

    #Add empty cells to array when rcluster.ntot<int(2.0*float(nmass))
    if len(lrprofn)!=nrad:
        for i in range(0,nrad-len(lrprofn)):
            lrprofn.append(0.0)
            vprof.append(0.0)
    
    return lrprofn,vprof


