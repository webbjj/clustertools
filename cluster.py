import math
from units import *

#Define the StarCluster class and add properties to the cluster
class StarCluster(object):
    def __init__(self,ntot,tphys=0.0,nc=0,rc=0.0,rbar=1.0,rtide=0.0,xc=0.0,yc=0.0,zc=0.0,zmbar=1.0,vstar=1.0,rscale=1.0,ns=0.0,nb=0.0,np=0.0,units='nbody',origin='cluster'):
        #Total Number of Stars + Binaries
        self.ntot=ntot
        #Current Time
        self.tphys=tphys
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
        #Units (nbody/realpc/realkpc/galpy)
        self.units=units
        #Origin (cluster/galaxy)
        self.origin=origin
    
    #Add an array of stars with id's, masses, positions, and velocities
    def add_stars(self,id,m,x,y,z,vx,vy,vz):
        self.id=id
        self.m=m
        self.x=x
        self.y=y
        self.z=z
        self.vx=vx
        self.vy=vy
        self.vz=vz
    
        #Calculate Key Parameters
        self.v=[]
        self.rxy=[]
        self.r=[]
        self.mtot=0.0
        #Mass within rtide
        self.mgc=0.0
        #Mass within core
        self.mc=0.0
        #Mass beyond rtide
        self.mrt=0.0
        
        for i in range(0,len(x)):
            self.v.append(math.sqrt(self.vx[i]**2.0+self.vy[i]**2.0+self.vz[i]**2.0))
            self.rxy.append(math.sqrt(self.x[i]**2.0+self.y[i]**2.0))
            self.r.append(math.sqrt(self.x[i]**2.0+self.y[i]**2.0+self.z[i]**2.0))
            self.mtot+=self.m[i]
            if self.r[i]<self.rc:
                self.mc=self.mc+self.m[i]
            if self.r[i]<=self.rtide:
                self.mgc=self.mgc+self.m[i]
            else:
                self.mrt=self.mrt+self.m[i]

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

    #Calculate lagrange radii
    #Units will be whatever units the cluster is currently in
    def rlagrange(self,nlagrange=10):
        
        #Shift to clustercentric origin if not so already:
        origin0=self.origin
        units0=self.units
        if origin0 !='cluster':
            xvshift(self,self.xgc,self.ygc,self.zgc,self.vxgc,self.vygc,self.vzgc,'cluster')
        if units0 =='realkpc':
            kpctopc(self)
        
        #If rtide was not given:
        if self.mgc==0.0:
            self.mgc=self.mtot

        #Radially order the stars
        rorder=sorted(range(0,self.ntot),key=lambda k:self.r[k])
        msum=0.0
        nfrac=1
        self.rn=[]
        for i in range(0,self.ntot):
            mfrac=self.mgc*float(nfrac)/float(nlagrange)
            msum+=self.m[rorder[i]]
            if msum >= mfrac:
                self.rn.append(self.r[rorder[i]])
                nfrac+=1
            if nfrac>nlagrange:
                break

        #If rc was not given let rc=r10:
        if self.rc==0.0:
            for i in range(0,nlagrange):
                if float(i)/float(nlagrange) >= 0.1:
                    self.rc=self.rn[i]
                    self.mc=self.mgc*float(i)/float(nlagrange)
                    break

        #Shift back to original origin and units
        if units0 =='realkpc' and self.units=='realpc':
            pctokpc(self)
        if origin0 == 'galaxy' and self.origin=='cluster':
            xvshift(self,self.xgc,self.ygc,self.zgc,self.vxgc,self.vygc,self.vzgc,'galaxy')



 
