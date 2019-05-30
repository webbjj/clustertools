#Custum output files

import numpy as np

import matplotlib.pyplot as plt
from galpy.util import bovy_plot
import os

from ..nbodypy.cluster import sub_cluster
from ..util.coordinates import sky_coords
from ..nbodypy.functions import *
from ..nbodypy.profiles import *
from ..nbodypy.operations import *


from .custom_plots import dvplot, bplot, aplot

def trelax_prof_out(cluster,fileout,multimass=True,projected=False):
    #Write relaxation time profile (trhprof.npy)
    trelax=relaxation_time(cluster,local=False,multimass=multimass,projected=projected)
    fileout.write("%f %f " % (cluster.tphys,trelax))
    if cluster.rn==None:
        rn=rlagrange(cluster,nlagrange=10,projected=projected)
    trelax_prof=[]

    for r in rn:
        fileout.write("%f " % (r))

    for i in range(0,len(rn)):
        if i==0:
            rmin=0.0
            rmax=rn[i]
        else:
            rmin=rn[i-1]
            rmax=rn[i]

        rcluster=sub_cluster(cluster,rmin=rmin,rmax=rmax,projected=projected)
        trelax_prof.append(relaxation_time(rcluster,local=True,multimass=multimass,projected=projected))
        fileout.write("%f " % (trelax_prof[-1]))

    fileout.write("%f %f " % (cluster.rmpro,cluster.rhpro))

    fileout.write("\n")

def p_prof_out(cluster,fileout,nrad=20,projected=False):
    #Write density profile (pprof.npy)
    fileout.write("%f " % (cluster.tphys))
    if cluster.rn==None or len(cluster.rn)!=nrad:
        rn=rlagrange(cluster,nlagrange=nrad,projected=projected)
    mn=cluster.mtot/float(nrad)
    p_prof=[]

    for r in rn:
        fileout.write("%f " % (r))

    for i in range(0,len(rn)):
        if i==0:
            rmin=0.0
            rmax=rn[i]
            vol=4.0*np.pi*(rmax**3.0)/3.0
        else:
            rmin=rn[i-1]
            rmax=rn[i]
            vol=4.0*np.pi*(rmax**3.0)/3.0-4.0*np.pi*(rmin**3.0)/3.0

        p_prof.append(mn/vol)
        fileout.write("%f " % (p_prof[-1]))

    fileout.write("\n")

def alpha_prof_out(cluster,fileout,mmin=0.3,mmax=0.8,rmin=None,rmax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None):
    #Write alpha_profile and dalpha for a given mass and radius range (alpha_prof.npy)
    m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(cluster,mmin=mmin,mmax=mmax,nmass=10,rmin=rmin,rmax=rmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)
    lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha=alpha_prof(cluster,mmin=mmin,mmax=mmax,nmass=10,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)

    fileout.write("%f %f %f %f %f " % (cluster.tphys,alpha,ealpha,yalpha,eyalpha))
    for i in range(0,len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0,len(dm)):
        fileout.write("%f " % dm[i])
    for i in range(0,len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0,len(aprof)):
        fileout.write("%f " % aprof[i])

    fileout.write("%f %f %f %f\n" % (dalpha,edalpha,ydalpha,eydalpha))

def dalpha_out(cluster,fileout,mmin=[0.1,0.3,0.5],mmax=[0.5,0.8,0.8],rmin=None,rmax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None):
    #Output alpha and dalpha for a range of values (dalpha_prof.npy)

    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    for i in range(0,len(mmin)):
        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(cluster,mmin=mmin[i],mmax=mmax[i],nmass=10,rmin=rmin,rmax=rmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)
        fileout.write("%f " % alpha)

    for i in range(0,len(mmin)):
        lrprofn,aprof,da,eda,yda,eyda=alpha_prof(cluster,mmin=mmin[i],mmax=mmax[i],nmass=10,projected=projected)
        fileout.write("%f %f %f\n " % (da,yda,cluster.rm))

def sigv_out(cluster,fileout,projected=False):
    #Output velocity dispersion profile and anisotropy profile (dvprof.npy)

    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    lrprofn,sigvprof,betaprof=sigv_prof(cluster,projected=projected)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for sig in sigvprof:
        fileout.write("%f " % sig)
    for beta in betaprof:
        fileout.write("%f " % beta)

    fileout.write("%f\n" % (cluster.rm))

def eta_prof_out(cluster,fileout,mmin=0.3,mmax=0.8,rmin=None,rmax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None):
    #output eta profile (eta_prof.npy)


    m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=mmin,mmax=mmax,nmass=10,rmin=rmin,rmax=rmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)
    lrprofn,eprof,deta,edeta,ydeta,eydeta=eta_prof(cluster,mmin=mmin,mmax=mmax,nmass=10,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)

    fileout.write("%f %f %f %f %f " % (cluster.tphys,eta,eeta,yeta,eyeta))
    for i in range(0,len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0,len(sigvm)):
        fileout.write("%f " % sigvm[i])
    for i in range(0,len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0,len(eprof)):
        fileout.write("%f " % eprof[i])

    fileout.write("%f %f %f %f\n" % (deta,edeta,ydeta,eydeta))

def eta_out(cluster,fileout,mmin=[0.1,0.3,0.5],mmax=[0.5,0.8,0.8],rmin=None,rmax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None):
    #Output eta and deta for a range of values (deta_prof.npy)


    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    for i in range(0,len(mmin)):
        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=mmin[i],mmax=mmax[i],nmass=10,rmin=rmin,rmax=rmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)
        fileout.write("%f " % eta)

    for i in range(0,len(mmin)):
        lrprofn,eprof,deta,edeta,ydeta,eydeta=eta_prof(cluster,mmin=mmin[i],mmax=mmax[i],nmass=10,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)
        fileout.write("%f %f %f\n " % deta,ydeta,cluster.rm)


def v_out(cluster,fileout,coord=None,projected=False):
    #Output mean velocity profile (vprof.npy)

    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    lrprofn,vprof=v_prof(cluster,coord=coord,projected=projected)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for v in vprof:
        fileout.write("%f " % v)

    fileout.write("%f\n" % (cluster.rm))

def dvanimate(data,nsnap=0,tsnap=None,prefix='',nrad=20,save=True,**kwargs):
    """
    NAME:

       dvanimate

    PURPOSE:

       Plot velocity dispersion profiles using dvplot over a range of timesteps
       To Do:
        - use data from StarCluster as opposed to a written file

    INPUT:

       data - array from custom_output.sigv_out (t,m,rm,r[0:nrad],sig[0:nrad],beta[0:nrad])

       nsnap - starting snapshot (default:0)

       tsnap - starting time step (default: 0 Myr, overwrites nsnap)

       prefix - string noting the begining of directory to save images

       nrad - number of radial bins in profile

       save - option to save images or simply just show them (default: True)

    KWARGS:

       kwargs - for passing to plots.dvplot

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    t=data[:,0]
    m=data[:,1]
    rm=data[:,-1]

    
    if not os.path.exists('./%sdvmovie/' % prefix):
        os.makedirs('./%sdvmovie/' % prefix)

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    for i in range(nsnap,len(t)+1):

        if save:
            filename='./%sdvmovie/%sdvplot%s.png' % (prefix,prefix,str(i))
        else:
            filename=None

        npy.dvplot(data,nsnap,tsnap,prefix='',nrad=nrad,nsnap=i,filename=filename,**kwargs)

    return 0

def banimate(data,nsnap=0,tsnap=None,prefix='',nrad=20,save=True,**kwargs):
    """
    NAME:

       banimate

    PURPOSE:

       Plot anisotropy profiles using bplot over a range of timesteps
       To Do:
        - use data from StarCluster as opposed to a written file

    INPUT:

       data - array from custom_output.sigv_out (t,m,rm,r[0:nrad],sig[0:nrad],beta[0:nrad])

       nsnap - starting snapshot (default:0)

       tsnap - starting time step (default: 0 Myr, overwrites nsnap)

       prefix - string noting the begining of directory to save images

       nrad - number of radial bins in profile

       save - option to save images or simply just show them (default: True)

    KWARGS:

       kwargs - for passing to plots.bplot

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    t=data[:,0]
    m=data[:,1]
    rm=data[:,-1]

    
    if not os.path.exists('./%sbmovie/' % prefix):
        os.makedirs('./%sbmovie/' % prefix)

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    for i in range(nsnap,len(t)+1):

        if save:
            filename='./%sbmovie/%sbplot%s.png' % (prefix,prefix,str(i))
        else:
            filename=None

        npy.bplot(data,i,tsnap,prefix='',nrad=nrad,filename=filename,**kwargs)

    return 0

def aanimate(data,nsnap=0,tsnap=None,prefix='',nmass=10,nrad=20,save=True,**kwargs):
    """
    NAME:

       aanimate

    PURPOSE:

       Plot alpha profiles using aplot over a range of timesteps
       To Do:
        - use data from StarCluster as opposed to a written file

    INPUT:

       data - array from custom_output.alpha_prof_out (t,alpha,ealpha,yalpha,eyalpha,mmean,dm,lr,ar,dalpha,edalpha,ydalpha,eydalpha,)

       nsnap - starting snapshot (default:0)

       tsnap - starting time step (default: 0 Myr, overwrites nsnap)

       prefix - string noting the begining of directory to save images

       nrad - number of radial bins in profile

       save - option to save images or simply just show them (default: True)

    KWARGS:

       kwargs - for passing to plots.aplot

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    t=data[:,0]
    
    if not os.path.exists('./%samovie/' % prefix):
        os.makedirs('./%samovie/' % prefix)

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    for i in range(nsnap,len(t)+1):

        if save:
            filename='./%samovie/%saplot%s.png' % (prefix,prefix,str(i))
        else:
            filename=None

        aplot(data,i,tsnap,prefix='',nmass=nmass,nrad=nrad,filename=filename,**kwargs)
        plt.close()

    return 0
