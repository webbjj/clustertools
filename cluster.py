import numpy as np
from galpy.util import bovy_conversion


#Define the StarClustee class and add properties to the cluster
#With the exception of total cluster mass and core mass (if core radius is given), and other cluster properties
#
class StarCluster(object):
    def __init__(self,ntot=0,tphys=0.0,units=None,origin=None,center=False,**kwargs):
        
        #Total Number of Stars + Binaries in the cluster
        self.ntot=ntot
        self.nb=0
        #Age of cluster
        self.tphys=tphys
        
        if 'ctype' in kwargs:
            self.ctype=kwargs.get('ctype')
        if 'sfile' in kwargs:
            self.sfile=kwargs.get('sfile')
        if 'bfile' in kwargs:
            self.bfile=kwargs.get('bfile')
        if 'delimiter' in kwargs:
            self.delimiter=kwargs.get('delimiter')
        if 'nsnap' in kwargs:
             self.nsnap=kwargs.get('nsnap')
        if 'nzfill' in kwargs:
            self.nzfill=kwargs.get('nzfill')
        if 'wdir' in kwargs:
            self.wdir=kwargs.get('wdir')
 
        #Initial arrays
        self.id=np.array([])
        self.m=np.array([])
        self.x=np.array([])
        self.y=np.array([])
        self.z=np.array([])
        self.vx=np.array([])
        self.vy=np.array([])
        self.vz=np.array([])
        self.kw=np.array([])

        self.zmbar=1.
        self.rbar=1.
        self.vstar=1.
        self.tstar=1.

        self.xc=0.
        self.yc=0.
        self.zc=0.
        self.vxc=0.
        self.vyc=0.
        self.vzc=0.

        self.xgc=0.
        self.ygc=0.
        self.zgc=0.
        self.vxgc=0.
        self.vygc=0.
        self.vzgc=0.
        
        self.orbit=None

        #SE Arrays
        self.logl=np.asarray([])
        self.logr=np.asarray([])
        self.lum=np.asarray([])
        self.ep=np.asarray([])
        self.ospin=np.asarray([])
        
        #BSE Array
        self.id1=np.asarray([])
        self.id2=np.asarray([])
        self.kw1=np.asarray([])
        self.kw2=np.asarray([])
        self.kcm=np.asarray([])
        self.ecc=np.asarray([])
        self.pb=np.asarray([])
        self.semi=np.asarray([])
        self.m1=np.asarray([])
        self.m2=np.asarray([])
        self.logl1=np.asarray([])
        self.logl2=np.asarray([])
        self.logr1=np.asarray([])
        self.logr2=np.asarray([])
        self.ep1=np.asarray([])
        self.ep2=np.asarray([])
        self.ospin1=np.asarray([])
        self.ospin2=np.asarray([])


        #Energies
        self.kin=np.asarray([])
        self.pot=np.asarray([])
        self.etot=np.asarray([])
        self.bound=np.ones(self.ntot,dtype=bool)

        #Functions
        self.trh=None
        self.alpha=None
        self.ealpha=None
        self.dalpha=None
        self.edalpha=None
        self.rn=None
        self.eta=None
        self.eeta=None

        self.units=units
        self.origin=origin
        self.center=False

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


        self.key_params()

    #Add stellar evolution parameters (type (NBODY6), luminosity, radius)
    #Units not specified
    def add_se(self,kw,logl,logr,ep,ospin):
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
        
        self.kin=np.asarray(kin)
        self.pot=np.asarray(pot)
        self.etot=np.asarray(etot)

    #Add cluster's orbital properties
    #Units not specified
    def add_orbit(self,xgc,ygc,zgc,vxgc,vygc,vzgc):
        self.xgc=xgc
        self.ygc=ygc
        self.zgc=zgc
        self.rgc=np.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
        self.vxgc=vxgc
        self.vygc=vygc
        self.vzgc=vzgc

    def find_center(self,xgc,ygc,zgc,nsphere=100):
        #Iterate to find center of cluster
        
        x=self.x
        y=self.y
        z=self.z
        r=self.r
        id=self.id

        while len(r)>nsphere:
            sigma=np.std(r)
            indx=(r<sigma)
            
            if len(r[indx])>nsphere:
                id=id[indx]
                x=x[indx]-np.mean(x[indx])
                y=y[indx]-np.mean(y[indx])
                z=z[indx]-np.mean(z[indx])
                r=np.sqrt(x*x+y*y+z*z)
            else:
                break

        #Find center of mass and velocity of inner stars:
        indx=np.in1d(self.id,id)
        
        xgc=np.sum(self.m[indx]*self.x[indx])/np.sum(self.m[indx])
        ygc=np.sum(self.m[indx]*self.y[indx])/np.sum(self.m[indx])
        zgc=np.sum(self.m[indx]*self.z[indx])/np.sum(self.m[indx])

        vxgc=np.sum(self.m[indx]*self.vx[indx])/np.sum(self.m[indx])
        vygc=np.sum(self.m[indx]*self.vy[indx])/np.sum(self.m[indx])
        vzgc=np.sum(self.m[indx]*self.vz[indx])/np.sum(self.m[indx])
    
        return xgc,ygc,zgc,vxgc,vygc,vzgc

    def key_params(self):
        #Calculate Key Parameters
        self.v=np.sqrt(self.vx**2.+self.vy**2.+self.vz**2.)
        self.r=np.sqrt(self.x**2.+self.y**2.+self.z**2.)
        self.rpro=np.sqrt(self.x**2.+self.y**2.)
        self.mtot=np.sum(self.m)
        self.lr=np.log10(self.r)
        self.lv=np.log10(self.v)
        self.rmean=np.mean(self.r)
        self.rmax=np.max(self.r)

        if self.origin=='cluster':
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
            if len(self.logl)>0:
                lsum=0.0
                for i in range(0,self.ntot):
                    lsum+=self.lum[self.rproorder[i]]
                    if lsum>=0.5*self.ltot:
                        self.rhpro=self.rpro[self.rproorder[i]]
                        break

    #Change units
    def to_realpc(self):
        if self.units=='galpy': self.to_realkpc()

        if self.units=='nbody':
            self.m*=self.zmbar
            self.x*=self.rbar
            self.y*=self.rbar
            self.z*=self.rbar
            self.vx*=self.vstar
            self.vy*=self.vstar
            self.vz*=self.vstar
            
            self.xc*=self.rbar
            self.yc*=self.rbar
            self.zc*=self.rbar
            self.vxc*=self.vstar
            self.vyc*=self.vstar
            self.vzc*=self.vstar
            
            self.xgc*=self.rbar
            self.ygc*=self.rbar
            self.zgc*=self.rbar
            self.vxgc*=self.vstar
            self.vygc*=self.vstar
            self.vzgc*=self.vstar
            
            self.units='realpc'
            
            if self.nb>0:
                yrs = (self.rbar*1296000./(2.0*np.pi))**1.5/np.sqrt(self.zmbar)
                days = 365.25*yrs
                pctoau=206265.
                
                self.pb*=days
                self.semi*=self.rbar*pctoau
            
            self.key_params()
    
        elif self.units == 'realkpc':
            self.x*=1000.0
            self.y*=1000.0
            self.z*=1000.0
            
            self.xgc*=1000.0
            self.ygc*=1000.0
            self.zgc*=1000.0
            
            self.xc*=1000.0
            self.yc*=1000.0
            self.zc*=1000.0
            
            self.units='realpc'
            self.key_params()
        
        else:
            print('CLUSTER UNITS ALREADY EQAULS ',self.units)


    def to_realkpc(self,r0=8.,v0=220.):
        if self.units=='galpy':
            self.m*=bovy_conversion.mass_in_msol(ro=r0,vo=v0)
            self.x*=r0
            self.y*=r0
            self.z*=r0
            self.vx*=v0
            self.vy*=v0
            self.vz*=v0
            
            self.xc*=r0
            self.yc*=r0
            self.zc*=r0
            self.vxc*=v0
            self.vyc*=v0
            self.vzc*=v0
            
            self.xgc*=r0
            self.ygc*=r0
            self.zgc*=r0
            self.vxgc*=v0
            self.vygc*=v0
            self.vzgc*=v0

            self.units='realkpc'
            self.key_params()

        elif self.units=='nbody':
            self.m*=self.zmbar
            self.x*=(self.rbar/1000.0)
            self.y*=(self.rbar/1000.0)
            self.z*=(self.rbar/1000.0)
            self.vx*=self.vstar
            self.vy*=self.vstar
            self.vz*=self.vstar
            
            self.xc*=(self.rbar/1000.0)
            self.yc*=(self.rbar/1000.0)
            self.zc*=(self.rbar/1000.0)
            self.vxc*=self.vstar
            self.vyc*=self.vstar
            self.vzc*=self.vstar
            
            self.xgc*=(self.rbar/1000.0)
            self.ygc*=(self.rbar/1000.0)
            self.zgc*=(self.rbar/1000.0)
            self.vxgc*=self.vstar
            self.vygc*=self.vstar
            self.vzgc*=self.vstar
            
            self.units='realkpc'
            self.key_params()

        
        elif self.units=='realpc':
            self.x/=1000.0
            self.y/=1000.0
            self.z/=1000.0
            
            self.xgc/=1000.0
            self.ygc/=1000.0
            self.zgc/=1000.0
            
            self.xc/=1000.0
            self.yc/=1000.0
            self.zc/=1000.0
            
            self.units='realkpc'
            self.key_params()
        else:
            print('CLUSTER UNITS ALREADY EQAULS ',self.units)


    def to_nbody(self,r0=8.,v0=220.):
        
        if self.units=='galpy': self.to_realkpc()
        if self.units!='realpc': self.to_realpc()

        if self.units=='realpc':
            self.m/=self.zmbar
            self.x/=self.rbar
            self.y/=self.rbar
            self.z/=self.rbar
            self.vx/=self.vstar
            self.vy/=self.vstar
            self.vz/=self.vstar
            
            self.xc/=self.rbar
            self.yc/=self.rbar
            self.zc/=self.rbar
            self.vxc/=self.vstar
            self.vyc/=self.vstar
            self.vzc/=self.vstar
            
            self.xgc/=self.rbar
            self.ygc/=self.rbar
            self.zgc/=self.rbar
            self.vxgc/=self.vstar
            self.vygc/=self.vstar
            self.vzgc/=self.vstar
            
            self.units='nbody'
            self.key_params()

        else:
            print('CLUSTER UNITS ALREADY EQAULS ',self.units)


    def to_galpy(self,r0=8.,v0=220.):
        if self.units!='realkpc' and self.units!='galpy':
            self.to_realkpc()
        
        if self.units=='realkpc':
            self.m=self.m/bovy_conversion.mass_in_msol(ro=r0,vo=v0)
            self.x/=r0
            self.y/=r0
            self.z/=r0
            self.vx/=v0
            self.vy/=v0
            self.vz/=v0
            
            self.xc/=r0
            self.yc/=r0
            self.zc/=r0
            self.vxc/=v0
            self.vyc/=v0
            self.vzc/=v0
            
            self.xgc/=r0
            self.ygc/=r0
            self.zgc/=r0
            self.vxgc/=v0
            self.vygc/=v0
            self.vzgc/=v0
            
            self.units='galpy'
            self.key_params()
        else:
            print('CLUSTER UNITS ALREADY EQAULS ',self.units)

    #Change origin
    def to_center(self):
        if not self.center:

            if self.origin=='galaxy':
                self.to_cluster()
        
            self.x-=self.xc
            self.y-=self.yc
            self.z-=self.zc
            self.vx-=self.vxc
            self.vy-=self.vyc
            self.vz-=self.vzc

            self.key_params()
            
            self.center=True
        else:
            print('CLUSTER CENTER ALREADY EQAULS ',self.center)


    def from_center(self):
        
        if self.center and self.origin=='cluster':
            self.x+=self.xc
            self.y+=self.yc
            self.z+=self.zc
            self.vx+=self.vxc
            self.vy+=self.vyc
            self.vz+=self.vzc

            self.key_params()
            
            self.center=False
        else:
            print('CLUSTER CENTER ALREADY EQAULS ',self.center)



    def to_cluster(self,to_center=False):
        if self.origin!='cluster':
            self.x-=self.xgc
            self.y-=self.ygc
            self.z-=self.zgc
            self.vx-=self.vxgc
            self.vy-=self.vygc
            self.vz-=self.vzgc
            self.origin='cluster'

            if to_center:
                self.to_center()

            self.key_params()
        else:
            print('CLUSTER ORIGIN ALREADY EQAULS ',self.origin)


    def to_galaxy(self,from_center=False):

        if self.origin!='galaxy':
            if self.center:
                self.from_center()
            
            self.x+=self.xgc
            self.y+=self.ygc
            self.z+=self.zgc
            self.vx+=self.vxgc
            self.vy+=self.vygc
            self.vz+=self.vzgc
            
            self.key_params()


            self.origin='galaxy'
        else:
            print('CLUSTER ORIGIN ALREADY EQAULS ',self.origin)
