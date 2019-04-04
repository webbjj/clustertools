#Output files containing key cluster properties 
#Only functions and profiles should be called here

from functions import *
from profiles import *
import numpy as np
from coordinates import sky_coords
from cluster import sub_cluster

def extrct_out(cluster,fileout,projected=False):
    #Write key properties (extrct.npy)
    
    units0=cluster.units
    origin0=cluster.origin

    if cluster.ntot==0:
        nb=0
        cluster.mtot=0.0
        trh=0.0
        rn=np.zeros(10)
    else:
        if origin0=='galaxy':
            cluster.to_cluster()
        if units0!='realpc':
            cluster.to_realpc()

        if cluster.nb>0:
            nb=len(cluster.m2)
        else:
            nb=0

        trh=relaxation_time(cluster,local=False,multimass=True,projected=projected)

        if cluster.ntot > 10:    
            if cluster.rn==None:
                rn=rlagrange(cluster,nlagrange=10,projected=projected)
        else:
            rn=np.zeros(10)
                
    fileout.write("%i %i %f %f %f " % (cluster.ntot,nb,cluster.tphys,trh,cluster.mtot))

    for i in range(0,len(rn)):
        fileout.write("%f " % rn[i])

    fileout.write("%f " % cluster.rmpro)

    if len(cluster.logl)>0:
        fileout.write("%f " % cluster.rhpro)

    fileout.write("\n")

    if units0=='realkpc':
        cluster.to_realkpc()
    elif units0=='nbody':
        cluster.to_nbody()
    elif units0=='galpy':
        cluster.to_galpy()
    if origin0=='galaxy':
        cluster.to_galaxy()

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

def p_prof_out(cluster,fileout,nrad=10,projected=False):
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
        fileout.write("%f %f %f\n " % da,yda,cluster.rm)

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

def snapout(cluster,filename,energies=False,observations=False):
    #Output a snapshot in nbodypy format

    if observations:

        ra,dec,d0,pmra,pmdec,vr0=sky_coords(cluster)

        if energies:
            np.savetxt(filename,np.olumn_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,cluster.kin,cluster.pot,cluster.etot,ra,dec,d0,pmra,pmdec,vr0]))
        else:
            np.savetxt(filename,np.column_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,ra,dec,d0,pmra,pmdec,vr0]))

    else:
        if energies:
            np.savetxt(filename,np.olumn_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,cluster.kin,cluster.pot,cluster.etot]))
        else:
            np.savetxt(filename,np.column_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw]))

    return 0
