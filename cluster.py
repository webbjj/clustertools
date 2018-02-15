#Define the StarCluster class and add properties to the cluster
import math
from units import *

class StarCluster(object):
    def __init__(self,ntot,tphys=0.0,nc=0,rc=0.0,rbar=1.0,rtide=0.0,xc=0.0,yc=0.0,zc=0.0,zmbar=1.0,vstar=1.0,rscale=1.0,ns=0.0,nb=0.0,np=0.0,units='nbody',origin='cluster'):
        self.ntot=ntot
        self.tphys=tphys
        self.nc=nc
        self.rc=rc
        self.rbar=rbar
        self.rtide=rtide
        self.xc=xc
        self.yc=yc
        self.zc=zc
        self.zmbar=zmbar
        self.vstar=vstar
        self.rscale=rscale
        self.ns=ns
        self.nb=nb
        self.np=np
        self.units=units
        self.origin=origin
    
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
        self.mgc=0.0
        self.mc=0.0
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

    def add_se(self,kw,logl,logr):
        self.kw=kw
        self.logl=logl
        self.logr=logr

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

    def add_energies(self,kin,pot,etot):
        self.kin=kin
        self.pot=pot
        self.etot=etot
    
    def add_orbit(self,xgc,ygc,zgc,vxgc,vygc,vzgc):
        self.xgc=xgc
        self.ygc=ygc
        self.zgc=zgc
        self.rgc=math.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
        self.vxgc=vxgc
        self.vygc=vygc
        self.vzgc=vzgc

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



 
