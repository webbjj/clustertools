#Output files containing key cluster properties 
from functions import *
from coordinates import *
import numpy as np

def extrct_out(cluster,fileout):

    if cluster.ntot==0:
        origin0=cluster.origin
        units0=cluster.units
        nb=0
        cluster.mtot=0.0
        trh=0.0
        rn=np.zeros(10)
    
    else:
        if cluster.origin=='galaxy':
            origin0=cluster.origin
            xvshift(cluster,-cluster.xgc,-cluster.ygc,-cluster.zgc,-cluster.vxgc,-cluster.vygc,-cluster.vzgc,'cluster')
        else:
            origin0=cluster.origin

        if cluster.units=='realkpc':
            units0=cluster.units
            kpctopc(cluster)
        else:
            units0=cluster.units

        if cluster.nb>0:
            nb=len(cluster.m2)
        else:
            nb=0

        trh=relaxation_time(cluster,local=False,multimass=True)

        if cluster.ntot > 10:    
            rn=rlagrange(cluster,10)
        else:
            rn=np.zeros(10)
                
    fileout.write("%i %i %f %f %f " % (cluster.ntot,nb,cluster.tphys,trh,cluster.mtot))

    for i in range(0,len(rn)):
        fileout.write("%f " % rn[i])

    fileout.write("%f " % cluster.rmpro)

    if len(cluster.logl)>0:
        fileout.write("%f " % cluster.rhpro)

    fileout.write("\n")

    if units0!=cluster.units and units0=='realkpc':
        pctokpc(cluster)
    if origin0!=cluster.origin and origin0=='galaxy':
        xvshift(cluster,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,'galaxy')

def trelax_prof_out(cluster,fileout,multimass=True):
    trelax=relaxation_time(cluster,local=False,multimass=multimass)
    fileout.write("%f %f " % (cluster.tphys,trelax))
    rn=rlagrange(cluster,nlagrange=10)
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

        rcluster=sub_cluster(cluster,rmin=rmin,rmax=rmax)
        trelax_prof.append(relaxation_time(rcluster,local=True,multimass=multimass))
        fileout.write("%f " % (trelax_prof[-1]))

    fileout.write("%f %f " % (cluster.rmpro,cluster.rhpro))

    fileout.write("\n")

def p_prof_out(cluster,fileout,nrad=10):
    fileout.write("%f " % (cluster.tphys))
    rn=rlagrange(cluster,nlagrange=nrad)
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

#Output dvprof.dat (WIP)
def eta_out(cluster,fileout,m_lower=[0.1,0.1,0.1],m_upper=[50.0,0.5,1.8]):
    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    yeta_all=[]

    rn=rlagrange(cluster,nlagrange=10)
    
    for i in range(0,len(m_lower)):
        for j in range(0,len(rn)):
            if j==0:
                m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10,rmin=0.0,rmax=rn[i])
            else:
                m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10,rmin=rn[i-1],rmax=rn[i])
            fileout.write("%f " % eta)
            yeta_all.append(yeta)
        
        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10,rmin=0.0,rmax=cluster.rm)
        fileout.write("%f " % eta)
        yeta_all.append(yeta)

        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10)
        fileout.write("%f " % eta)
        yeta_all.append(yeta)

    for yeta in yeta_all:
        fileout.write("%f " % yeta)

    fileout.write("\n")

#Output dalpha_prof.dat (TODO - NEED TO ADD PROJECTED VALUES)
def dalpha_out(cluster,fileout,m_lower=[0.1,0.3,0.5,0.1,0.3,0.5],m_upper=[0.5,0.8,0.8,0.5,0.8,0.8]):


    r_lower=[0.0,0.0,0.0,None,None,None]
    r_upper=[cluster.rm,cluster.rm,cluster.rm,None,None,None]


    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    alphag=[]
    dalpha=[]
    ydalpha=[]
   
    cluster.key_params()

    for i in range(0,len(m_lower)):

        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(cluster,m_lower[i],m_upper[i],10,rmin=r_lower[i],rmax=r_upper[i])
        alphag.append(alpha)

        lrprofn,aprof,da,eda,yda,eyda=alpha_prof(cluster,m_lower[i],m_upper[i],10,rmin=r_lower[i],rmax=r_upper[i],nrad=10)
        dalpha.append(da)
        ydalpha.append(yda)

    for ag in alphag:
        fileout.write("%f " % ag)
    for da in dalpha:
        fileout.write("%f " % da)
    for yda in ydalpha:
        fileout.write("%f " % yda)
    fileout.write("%f\n" % cluster.rm)

def alpha_prof_out(cluster,fileout,mmin=0.3,mmax=0.8,rmin=None,rmax=None,se_min=0,se_max=1,projected=False,obs_cut=None):

    m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(cluster,mmin=mmin,mmax=mmax,nmass=10,rmin=rmin,rmax=rmax,se_min=se_min,se_max=se_max,projected=projected)
    lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha=alpha_prof(cluster,mmin=mmin,mmax=mmax,nmass=10,se_min=se_min,se_max=se_max,projected=projected,obs_cut=obs_cut)


    print(cluster.tphys,cluster.ntot,cluster.rm,cluster.rmpro,cluster.mtot,np.min(cluster.kw),np.max(cluster.kw),alpha)
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

#Output a snapshot
def snapout(cluster,filename,energies=False):

    if energies:
        np.savetxt(filename,np.olumn_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,cluster.kin,cluster.pot,cluster.etot]))
    else:
        np.savetxt(filename,np.column_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw]))

    return 0

#Output dvprof.dat
def sigv_out(cluster,fileout,do_beta=False):
    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    lrprofn,sigvprof=sigv_prof(cluster)

    if do_beta:
        lrprofn,betaprof=beta_prof(cluster)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for sig in sigvprof:
        fileout.write("%f " % sig)
    if do_beta:
        for beta in betaprof:
            fileout.write("%f " % beta)

    fileout.write("%f\n" % (cluster.rm))

#Output vprof.dat
def v_out(cluster,fileout,coord=None):
    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    lrprofn,vprof=v_prof(cluster,coord=coord)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for v in vprof:
        fileout.write("%f " % v)

    fileout.write("%f\n" % (cluster.rm))
