import math
import numpy as np
from units import *

#Define the StarClustee class and add properties to the cluster
#With the exception of total cluster mass and core mass (if core radius is given), and other cluster properties
#
class StarCluster(object):
    def __init__(self,ntot=0,tphys=0.0,with_stars=False,id=None,m=None,x=None,y=None,z=None,vx=None,vy=None,vz=None,kw=None,units=None,origin=None,keyparams=True):
        #Set additional components as False initially
        self.nbody6=False
        self.se=False
        self.bse=False
        self.energies=False
        self.orbit=False
        self.rcom=False
        self.vcom=False
        #Total Number of Stars + Binaries in the cluster
        self.ntot=ntot
        #Age of cluster
        self.tphys=tphys
        
        #If array of stars and their properties are included when instance is defined:
        #Assume all coordinates are given in Msun, parsecs, km/s and clustercentric distance unless otherwise specified
        if with_stars:
            self.id=np.asarray(id)
            self.m=np.asarray(m)
            self.x=np.asarray(x)
            self.y=np.asarray(y)
            self.z=np.asarray(z)
            self.vx=np.asarray(vx)
            self.vy=np.asarray(vy)
            self.vz=np.asarray(vz)
            self.kw=np.asarray(kw)
        else:
            self.id=np.array([])
            self.m=np.array([])
            self.x=np.array([])
            self.y=np.array([])
            self.z=np.array([])
            self.vx=np.array([])
            self.vy=np.array([])
            self.vz=np.array([])
            self.kw=np.array([])

        self.units=units
        self.origin=origin
        self.keyparams=keyparams
        #Calculate key parameters
        if with_stars and self.keyparams:
            self.key_params()

    # Add key parameters from an Nbody6 output file
    def add_nbody6(self,nc=0,rc=0.0,rbar=1.0,rtide=0.0,xc=0.0,yc=0.0,zc=0.0,zmbar=1.0,vstar=1.0,rscale=1.0,ns=0.0,nb=0.0,np=0.0):
        self.nbody6=True

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
    def add_stars(self,id,m,x,y,z,vx,vy,vz,kw=None):
        self.id=np.append(self.id,np.asarray(id))
        self.m=np.append(self.m,np.asarray(m))
        self.x=np.append(self.x,np.asarray(x))
        self.y=np.append(self.y,np.asarray(y))
        self.z=np.append(self.z,np.asarray(z))
        self.vx=np.append(self.vx,np.asarray(vx))
        self.vy=np.append(self.vy,np.asarray(vy))
        self.vz=np.append(self.vz,np.asarray(vz))
        self.kw=np.append(self.kw,np.asarray(kw))
        

        if len(self.id)!=self.ntot:
            print('Added %i stars to instance' % (len(self.id)-self.ntot))
            self.ntot=len(self.id)

        #Calculate key parameters
        if self.keyparams:
            self.key_params()

    def key_params(self):
        #Calculate Key Parameters
        self.v=np.sqrt(self.vx**2.+self.vy**2.+self.vz**2.)
        self.r=np.sqrt(self.x**2.+self.y**2.+self.z**2.)
        self.rpro=np.sqrt(self.x**2.+self.y**2.)
        self.mtot=np.sum(self.m)
        self.lr=np.log10(self.r)
        self.lv=np.log10(self.v)
        self.rmean=np.mean(self.r)

        if self.origin!='cluster':
            print('Cannot compute rm when origin = %s' % self.origin)
        else:
            self.rmax=np.max(self.r)
            #Radially order the stars to find half-mass radius
            self.rorder=sorted(range(0,self.ntot),key=lambda k:self.r[k])
            msum=0.0
            for i in range(0,self.ntot):
                msum+=self.m[self.rorder[i]]
                if msum >= 0.5*self.mtot:
                    self.rm=self.r[self.rorder[i]]
                    break

            self.rproorder=sorted(range(0,self.ntot),key=lambda k:self.rpro[k])
            msum=0.0
            for i in range(0,self.ntot):
                msum+=self.m[self.rproorder[i]]
                if msum >= 0.5*self.mtot:
                    self.rmpro=self.rpro[self.rproorder[i]]
                    break
            if self.se:
                lsum=0.0
                for i in range(0,self.ntot):
                    lsum+=self.lum[self.rproorder[i]]
                    if lsum>=0.5*self.ltot:
                        self.rhpro=self.rpro[self.rproorder[i]]
                        break

    #Add stellar evolution parameters (type (NBODY6), luminosity, radius)
    #Units not specified
    def add_se(self,kw,logl,logr,ep,ospin):
        self.se=True
        self.kw=np.asarray(kw)
        self.logl=np.asarray(logl)
        self.logr=np.asarray(logr)
        self.lum=10.0**self.logl
        self.ltot=np.sum(self.lum)
        self.ep=np.asarray(ep)
        self.ospin=np.asarray(ospin)

    #Add binary evolution parameters (ids,type,eccentricity, semi major axis, masses, luminosities, and radii)
    #Units not specified
    def add_bse(self,id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2,ep1,ep2,ospin1,ospin2):
        self.bse=True
        self.id1=np.asarray(id1)
        self.id2=np.asarray(id2)
        self.kw1=np.asarray(kw1)
        self.kw2=np.asarray(kw2)
        self.kcm=np.asarray(kcm)
        self.ecc=np.asarray(ecc)
        self.pb=np.asarray(pb)
        self.semi=np.asarray(semi)
        self.m1=np.asarray(m1)
        self.m2=np.asarray(m2)
        self.logl1=np.asarray(logl1)
        self.logl2=np.asarray(logl2)
        self.logr1=np.asarray(logr1)
        self.logr2=np.asarray(logr2)
        self.ep1=np.asarray(ep1)
        self.ep2=np.asarray(ep2)
        self.ospin1=np.asarray(ospin1)
        self.ospin2=np.asarray(ospin2)

    #Add individual energies
    #Units not specified
    def add_energies(self,kin,pot,etot):
        self.energies=True
        self.kin=np.asarray(kin)
        self.pot=np.asarray(pot)
        self.etot=np.asarray(etot)

    #Add cluster's orbital properties
    #Units not specified
    def add_orbit(self,xgc,ygc,zgc,vxgc,vygc,vzgc):
        self.orbit=True
        self.xgc=xgc
        self.ygc=ygc
        self.zgc=zgc
        self.rgc=math.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
        self.vxgc=vxgc
        self.vygc=vygc
        self.vzgc=vzgc

   #Calculate location of center of mass:
    def r_com(self):
        self.rcom=True

        self.xcom=np.sum(self.m*self.x)/np.sum(self.m)
        self.ycom=np.sum(self.m*self.y)/np.sum(self.m)
        self.zcom=np.sum(self.m*self.z)/np.sum(self.m)

   #Calculate velocity of center of mass:
    def v_com(self):
        self.vcom=True
        self.xcom=np.sum(self.m*self.vx)/np.sum(self.m)
        self.ycom=np.sum(self.m*self.vy)/np.sum(self.m)
        self.zcom=np.sum(self.m*self.vz)/np.sum(self.m)

    def find_center(self,xgc,ygc,zgc):
        #Iterate to find center of cluster
        sigma_xgc=np.std(self.x)
        sigma_ygc=np.std(self.y)
        sigma_zgc=np.std(self.z)

        niterate=0

        while sigma_xgc>0.001 or sigma_ygc>0.001 or sigma_zgc>0.001:

            niterate+=1

            #Then find center of stars within 1sigma of median

            indx=(abs(self.x-xgc) < sigma_xgc) * (abs(self.y-ygc) < sigma_ygc) * (abs(self.z-zgc) < sigma_zgc)

            if len(self.x[indx])<=1:
                break

            xgc=np.mean(self.x[indx])
            ygc=np.mean(self.y[indx])
            zgc=np.mean(self.z[indx])
            vxgc=np.mean(self.vx[indx])
            vygc=np.mean(self.vy[indx])
            vzgc=np.mean(self.vz[indx])

            sigma_xgc=np.std(self.x[indx])
            sigma_ygc=np.std(self.y[indx])
            sigma_zgc=np.std(self.z[indx])

            if niterate>10:
                print('SEARCHING FOR CLUSTER CENTER....',niterate)

        return xgc,ygc,zgc,vxgc,vygc,vzgc


#Routine for joining cluster timesteps

def join_clusters(clusters):

    base_cluster=clusters[0]

    id=np.asarray([])
    m=np.asarray([])
    x=np.asarray([])
    y=np.asarray([])
    z=np.asarray([])
    vx=np.asarray([])
    vy=np.asarray([])
    vz=np.asarray([])
    kw=np.asarray([])
    
    for i in range(1,len(clusters)):
        id=np.append(id,clusters[i].id)
        m=np.append(m,clusters[i].m)
        x=np.append(x,clusters[i].x)
        y=np.append(y,clusters[i].y)
        z=np.append(z,clusters[i].z)
        x=np.append(vx,clusters[i].vx)
        y=np.append(vy,clusters[i].vy)
        z=np.append(vz,clusters[i].vz)
        kw=np.append(kw,clusters[i].kw)

    base_cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw)

    return base_cluster

# Extract a sub population of stars from a cluster
# Extraction criteria include radius, mass and stellar evolution (KW) type (default all stars)
def sub_cluster(cluster,rmin=None,rmax=None,mmin=None,mmax=None,se_min=0,se_max=100,e_min=None,e_max=None,projected=False):
 
    if not cluster.keyparams:
        cluster.key_params()
    
    if rmin==None and projected:
       rmin=np.min(cluster.rpro)
       rmax=np.max(cluster.rpro)
    elif rmin==None and not projected:
       rmin=np.min(cluster.r)
       rmax=np.max(cluster.r)

    if mmin==None:
        mmin=np.min(cluster.m)
    if mmax==None:
        mmax=np.max(cluster.m)

    if e_min==None and e_max!=None:
        eindx=cluster.etot<=e_max
    elif e_min!=None and e_max==None:
        eindx=cluster.etot>=e_min
    elif e_min!=None and e_max!=None:
        eindx=(cluster.etot<=e_max) * (cluster.etot>=e_min)
    else:
        eindx=cluster.id > -1

    indx=(cluster.r>=rmin) * (cluster.r<=rmax) * (cluster.m>=mmin) * (cluster.m<=mmax) * (cluster.kw>=se_min) * (cluster.kw<=se_max)


    indx=np.logical_and(indx,eindx)

    if len(cluster.id[indx])>0:
        subcluster=StarCluster(len(cluster.id[indx]),cluster.tphys,True,cluster.id[indx],cluster.m[indx],cluster.x[indx],cluster.y[indx],cluster.z[indx],cluster.vx[indx],cluster.vy[indx],cluster.vz[indx],cluster.kw[indx],units=cluster.units,origin=cluster.origin,keyparams=cluster.keyparams)
        if cluster.se:
            subcluster.add_se(cluster.kw[indx],cluster.logl[indx],cluster.logr[indx],cluster.ep[indx],cluster.ospin[indx])
        if cluster.bse:
            bindx=np.in1d(cluster.id1,cluster.id[indx])
            subcluster.add_bse(cluster.id1[bindx],cluster.id2[bindx],cluster.kw1[bindx],cluster.kw2[bindx],cluster.kcm[bindx],cluster.ecc[bindx],cluster.pb[bindx],cluster.semi[bindx],cluster.m1[bindx],cluster.m2[bindx],cluster.logl1[bindx],cluster.logl2[bindx],cluster.logr1[bindx],cluster.logr2[bindx],cluster.ep1[bindx],cluster.ep2[bindx],cluster.ospin1[bindx],cluster.ospin2[bindx])
        if cluster.energies:
            subcluster.add_energies(cluster.kin[indx],cluster.pot[indx],cluster.etot[indx])
    else:
        subcluster=StarCluster(0,cluster.tphys)

    return subcluster
