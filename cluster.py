import math
import numpy as np
from units import *

#Define the StarClustee class and add properties to the cluster
#With the exception of total cluster mass and core mass (if core radius is given), and other cluster properties
#
class StarCluster(object):
    def __init__(self,ntot=0,tphys=0.0,id=None,m=None,x=None,y=None,z=None,vx=None,vy=None,vz=None,units=None,origin=None,keyparams=True):
        #Total Number of Stars + Binaries in the cluster
        self.ntot=ntot
        #Age of cluster
        self.tphys=tphys
        
        #If array of stars and their properties are included when instance is defined:
        #Assume all coordinates are given in Msun, parsecs, km/s and clustercentric distance unless otherwise specified
        if id!= None:
            self.id=id
            self.m=m
            self.x=x
            self.y=y
            self.z=z
            self.vx=vx
            self.vy=vy
            self.vz=vz
        else:
            self.id=[]
            self.m=[]
            self.x=[]
            self.y=[]
            self.z=[]
            self.vx=[]
            self.vy=[]
            self.vz=[]
        
        self.units=units
        self.origin=origin
        self.keyparams=keyparams
        #Calculate key parameters
        if id!=None and self.keyparams:
            self.key_params()

    # Add key parameters from an Nbody6 output file
    def add_nbody6(self,nc=0,rc=0.0,rbar=1.0,rtide=0.0,xc=0.0,yc=0.0,zc=0.0,zmbar=1.0,vstar=1.0,rscale=1.0,ns=0.0,nb=0.0,np=0.0):
        #Number of stars in the core
        self.nc=nc
        #Core radius
        self.rc=rc
        #Distane scaling parameter
        self.rbar=rbar
        #Tidal limit from NBODY6 (not neccesarily a true tidal radius)
        self.rtide=rtide
        #Center of mass of cluster (x,yz)
        self.xc=xc
        self.yc=yc
        self.zc=zc
        #Mass scaling parameter
        self.zmbar=zmbar
        #Velocity scaling parameter
        self.vstar=vstar
        #Scale radius of cluster
        self.rscale=rscale
        #Number of single stars
        self.ns=ns
        #Number of binary stars
        self.nb=nb
        #Number of particles (from NBODY6 when tidal tail is being integrated)
        self.np=np
    
    #Add an array of stars with id's, masses, positions, and velocities
    def add_stars(self,id,m,x,y,z,vx,vy,vz):
        for i in range(0,len(id)):
            self.id.append(id[i])
            self.m.append(m[i])
            self.x.append(x[i])
            self.y.append(y[i])
            self.z.append(z[i])
            self.vx.append(vx[i])
            self.vy.append(vy[i])
            self.vz.append(vz[i])
        
        if len(self.id)!=self.ntot:
            print('Added %i stars to instance' % (len(self.id)-self.ntot))
            self.ntot=len(self.id)

        #Calculate key parameters
        if self.keyparams:
            self.key_params()

    def key_params(self):
        #Calculate Key Parameters
        self.v=[]
        self.rxy=[]
        self.r=[]
        self.mtot=0.0       
  
        for i in range(0,len(self.x)):
            self.v.append(math.sqrt(self.vx[i]**2.0+self.vy[i]**2.0+self.vz[i]**2.0))
            self.rxy.append(math.sqrt(self.x[i]**2.0+self.y[i]**2.0))
            self.r.append(math.sqrt(self.x[i]**2.0+self.y[i]**2.0+self.z[i]**2.0))
            self.mtot+=self.m[i]
        
        self.rmean=np.mean(self.r)

        if self.origin=='cluster':
            self.rmax=np.max(self.r)
            self.rm=self.rnfind(0.5)

    def rnfind(self,nfrac):
        
        if self.origin!='cluster':
            print('Cannot compute rm when origin = %s' % self.origin)
            return None
        else:
            #Radially order the stars to find half-mass radius
            self.rorder=sorted(range(0,self.ntot),key=lambda k:self.r[k])
            msum=0.0
            for i in range(0,self.ntot):
                msum+=self.m[self.rorder[i]]
                if msum >= nfrac*self.mtot:
                    rn=self.r[self.rorder[i]]
                    break
        return rn

    #Add stellar evolution parameters (type (NBODY6), luminosity, radius)
    #Units not specified
    def add_se(self,kw,logl,logr):
        self.kw=kw
        self.logl=logl
        self.logr=logr

    #Add binary evolution parameters (ids,type,eccentricity, semi major axis, masses, luminosities, and radii)
    #Units not specified
    def add_bse(self,id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2):
        self.id1=id1
        self.id2=id2
        self.kw1=kw1
        self.kw2=kw2
        self.kcm=kcm
        self.ecc=ecc
        self.pb=pb
        self.semi=semi
        self.m1=m1
        self.m2=m2
        self.logl1=logl1
        self.logl2=logl2
        self.logr1=logr1
        self.logr2=logr2

    #Add individual energies
    #Units not specified
    def add_energies(self,kin,pot,etot):
        self.kin=kin
        self.pot=pot
        self.etot=etot
    
    #Add cluster's orbital properties
    #Units not specified
    def add_orbit(self,xgc,ygc,zgc,vxgc,vygc,vzgc):
        self.xgc=xgc
        self.ygc=ygc
        self.zgc=zgc
        self.rgc=math.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
        self.vxgc=vxgc
        self.vygc=vygc
        self.vzgc=vzgc

#Routine for joining cluster timesteps

def join_clusters(clusters):

    base_cluster=clusters[0]

    id=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    
    for i in range(1,len(clusters)):
        for j in range(0,clusters[i].ntot):
            id.append(clusters[i].id[j])
            m.append(clusters[i].m[j])
            x.append(clusters[i].x[j])
            y.append(clusters[i].y[j])
            z.append(clusters[i].z[j])
            vx.append(clusters[i].vx[j])
            vy.append(clusters[i].vy[j])
            vz.append(clusters[i].vz[j])

    base_cluster.add_stars(id,m,x,y,z,vx,vy,vz)

    return base_cluster

def sub_cluster(cluster,rmin=None,rmax=None,mmin=None,mmax=None):

    id=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]

    for i in range(0,cluster.ntot):
        add_star=False
        if rmin!=None and mmin!= None:
            if cluster.r[i]>=rmin and cluster.r[i]<=rmax and cluster.m[i]>=mmin and cluster.m[i]<=mmax:
                add_star=True
        elif rmin!=None and cluster.r[i]>=rmin and cluster.r[i]<=rmax:
            add_star=True
        elif mmin!=None and cluster.m[i]>=mmin and cluster.m[i]<=mmax:
            add_star=True

        if add_star:
            id.append(clusters[i].id[i])
            m.append(clusters[i].m[i])
            x.append(clusters[i].x[i])
            y.append(clusters[i].y[i])
            z.append(clusters[i].z[i])
            vx.append(clusters[i].vx[i])
            vy.append(clusters[i].vy[i])
            vz.append(clusters[i].vz[i])

    return StarCluster(len(x),cluster.tphys,id,m,x,y,z,vx,vy,vz)





 
