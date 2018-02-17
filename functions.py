#Key calculations (constantly adding to this file)
import math
import numpy as np
from get_nbody import get_star
from get_nbody import get_starcluster
from cluster import *

#Half-Mass Relaxation Time (Units = Myr)
#Assumes virialized cluster
#TODO - ADD REFERENCE

def Half_Mass_Relaxation_Time(cluster):

    p50=3.0*(0.5*cluster.mtot)/(4.0*math.pi*(cluster.rm**3.0))
    #TODO - NEED TO EDIT DEFINITION OF SIG50
    sig50=np.std(cluster.v)
    trh=0.34*(sigv50**3.0)/(((4.302e-3)**2.0)*(1.0/3.0)*p50*math.log(0.4*cluster.mtot/(1.0/3.0)))
    trh=trh*3.086e13/(3600.0*24.0*365.0*1000000.0)

    return rh

#Relaxation time
#No assumptions of virialization
#TODO ADD REFERENCE
def Relaxation_Time(cluster):
    
    p50=3.0*(0.5*cluster.mtot)/(4.0*math.pi*(cluster.rm**3.0))

    trelax = float(cluster.ntot)/math.log(0.4*float(cluster.ntot),10)
    trelax=0.858*trelax*math.sqrt(math.pow(cluster.rm,3.0)/cluster.mtot)

    return trelax

#Calculate distances between stars in the cluster
def Inter_Stellar_Distances(cluster):
    s2sd=np.zeros((cluster.ntot,cluster.ntot))
    for i in range(0,cluster.ntot):
        for j in range(i,cluster.ntot):
            d=math.sqrt((cluster.x[i]-cluster.x[j])**2.0+(cluster.y[i]-cluster.y[j])**2.0+(cluster.z[i]-cluster.z[j])**2.0)
            s2sd[i][j]=d
            s2sd[j][i]=d

    return s2sd


#Calculate lagrange radii
#Units will be whatever units the cluster is currently in
def rlagrange(cluster,nlagrange=10):
    
    #Shift to clustercentric origin if not so already:
    origin0=cluster.origin
    units0=cluster.units
    if origin0 !='cluster':
        xvshift(cluster,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,'cluster')

    #Radially order the stars
    msum=0.0
    nfrac=1
    cluster.rn=[]
    for i in range(0,cluster.ntot):
        mfrac=cluster.mgc*float(nfrac)/float(nlagrange)
        msum+=cluster.m[cluster.rorder[i]]
        if msum >= mfrac:
            cluster.rn.append(cluster.r[cluster.rorder[i]])
            nfrac+=1
        if nfrac>nlagrange:
            break

    #Shift back to original origin
    if origin0 == 'galaxy' and cluster.origin=='cluster':
        xvshift(cluster,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,'galaxy')

#Split an array into nbin bin's of equal number elements
def nbinmaker(x,nbin):
    
    xorder=sorted(range(0,len(x)),key=lambda k:x[k])

    x_lower=[]
    x_upper=[]
    x_hist=[]
    x_sum=[]

    for i in range(0,nbin):
        indx=int(float(i)*float(len(x))/float(nbin))
        x_lower.append(x[xorder[indx]])
        indx=int(float(i+1)*float(len(x))/float(nbin))-1
        x_upper.append(x[xorder[indx]])
        x_hist.append(0.0)
        x_sum.append(0.0)

    for i in range(0,len(x)):
        for j in range(0,nbin):
            if x[i]>=x_lower[j] and x[i]<=x_upper[j]:
                x_hist[j]+=1.0
                x_sum[j]+=x[i]


    x_mean=[]
    for i in range(0,nbin):
        x_mean.append(x_sum[i]/x_hist[i])

    return x_lower,x_mean,x_upper,x_hist

#Find mass function over a given mass range using nmass bins containing an equal number of stars
def mass_function(cluster,mmin=0.1,mmax=1.0,nmass=10):

    msub=[]
    
    for i in range(0,cluster.ntot):
        if cluster.m[i]>=mmin and cluster.m[i]<=mmax:
            msub.append(cluster.m[i])

    m_lower,m_mean,m_upper,m_hist=nbinmaker(msub,nmass)

    lm_mean=[]
    dm=[]
    ldm=[]
    for i in range(0,nmass):
        lm_mean.append(math.log(m_mean[i],10.0))
        dm.append(m_hist[i]/(m_upper[i]-m_lower[i]))
        ldm.append(math.log(dm[i],10.0))


    (alpha,yalpha),V=np.polyfit(lm_mean,ldm,1,cov=True)
    ealpha=np.sqrt(V[0][0])
    eyalpha=np.sqrt(V[1][1])

    return m_mean,m_hist,alpha,ealpha,yalpha,eyalpha


def alpha_prof(cluster,mmin=0.1,mmax=100.0,nmass=10,rmin=0.0,rmax=100.0,nrad=10):

    stars=[]
    rprof=[]
    lrprofn=[]
    aprof=[]
    eaprof=[]

    
    #Build subcluster containing only stars in the full radial and mass range:
    for i in range(0,cluster.ntot):
        if cluster.r[i]>=rmin and cluster.r[i]<=rmax and cluster.m[i]>=mmin and cluster.m[i]<=mmax:
            stars.append(get_star(cluster,cluster.id[i]))

    subcluster=get_starcluster(stars)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(subcluster.r,nrad)

    for i in range(0,nrad):
        stars=[]
        for j in range(0,subcluster.ntot):
            if subcluster.r[j] >= r_lower[i] and subcluster.r[j] <= r_upper[i]:
                stars.append(get_star(subcluster,subcluster.id[j]))

        rcluster=get_starcluster(stars)

        m_mean,m_hist,alpha,ealpha,yalpha,eyalpha=mass_function(rcluster,mmin,mmax,nmass)

        rprof.append(r_mean[i])
        lrprofn.append(math.log(rprof[-1]/cluster.rm))
        aprof.append(alpha)
        eaprof.append(ealpha)

    (dalpha,ydalpha),V=np.polyfit(lrprofn,aprof,1,cov=True)
    edalpha=np.sqrt(V[0][0])
    eydalpha=np.sqrt(V[1][1])
    
    return lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha



        




