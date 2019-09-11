"""
CLUSTER

The StarCluster classes and key internal functions

"""
import numpy as np
from galpy.util import bovy_conversion,bovy_coords
from textwrap import dedent
from galpy.potential import MWPotential2014
from .orbit import rtidal,rlimiting,initialize_orbit,calc_actions
from .functions import *
from .profiles import *

class StarCluster(object):
    """
    Class representing a star cluster
    """    
    def __init__(self,ntot=0,tphys=0.0,units=None,origin=None,ctype='snapshot',**kwargs):
        """
        NAME:

           __init__

        PURPOSE:

           Initialize a Star Cluster
           Notes:
            - key stellar attributes are initialized so they can be added
            directly to the instance (self.x = x) or through one of the add
            functions below.

        INPUT:

           ntot - number of stars (default:0)

           tphys - age of cluster (default:0)

           units - units of input data (nbody/realpc/realkpc/galpy,default:None)

           origin - origin of input data (cluster/galaxy/default:None)

           ctype - code used to generate data (nbody6/nbody6se/gyrfalcon/...)


        KWARGS:

            sfile - name of file containing stellar data

            bfile - name of file contain binary star data

            ofilename - orbit filename if ofile is not given

            ounits - units of orbital information (else assumed equal to StarCluster.units)

            nsnap - if a specific snapshot is to be read in instead of starting from zero
        
            nzfill - value for zfill when reading and writing snapshots (Default: 5)
        
            delimiter - choice of delimiter when reading ascii/csv files (Default: ',')
        
            wdir - working directory of snapshots if not current directory

            intialize - initialize a galpy orbit after reading in orbital information (default: False)

            kwfile - open file containing stellar evolution type (kw) of individual stars (used if snapshot file does not contain kw and kw is needed)

            advance - is this a snapshot that has been advanced to from initial load_cluster?

            projected - calculate projected values as well as 3D values (Default: True)


        OUTPUT:

           instance

        HISTORY:

           2018 - Written - Webb (UofT)

        """        
        #Total Number of Stars + Binaries in the cluster
        self.ntot=ntot
        self.nb=0
        #Age of cluster
        self.tphys=tphys

        #Cluster Simulation Type
        self.ctype=ctype
        
        #Kwargs
        self.nsnap=int(kwargs.get('nsnap','0'))
        self.delimiter=kwargs.get('delimiter',None)
        self.wdir=kwargs.get('wdir','')
        self.nzfill=int(kwargs.get('nzfill','5'))
        self.snapbase=kwargs.get('snapbase','')
        self.snapend=kwargs.get('snapend','.dat')
        self.snapdir=kwargs.get('snapdir','')
        self.skiprows=kwargs.get('skiprows',0)
        self.sfile=kwargs.get('sfile','')
        self.bfile=kwargs.get('bfile','')
        self.projected=kwargs.get('projected',True)

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

        #Lagrange Radii,limiting radius, tidal radius, and virial radius
        self.rn=None
        self.rl=None
        self.rt=None
        self.rv=None
        self.rorder=None
        self.rproorder=None

        #Additional Parameters
        self.trh=None
        self.alpha=None
        self.dalpha=None
        self.eta=None

        self.units=units
        self.origin=origin

    def add_nbody6(self,nc=0,rc=0.0,rbar=1.0,rtide=0.0,xc=0.0,yc=0.0,zc=0.0,zmbar=1.0,vstar=1.0,rscale=1.0,ns=0.0,nb=0.0,np=0.0):
        """
        NAME:

           add_nbody6

        PURPOSE:

           For data generated using NBDOY6 (or one of its variants), add additional information from output files

        INPUT:

           nc - number of stars in core (default:0)

           rc - core radius (default:0)

           rbar - scaling factor between NBODY units and pc (default:1.)

           rtide - rtide set in the simulation (default:0)

           xc,yc,zc - position of cluster centre (default:0)

           zmbar - scaling factor between NBODY units and Msun (default:1.)

           vstar - scaling factor between NBODY units and km/s (default:1.)

           rscale - the scale radius of data (default:1.)

           ns - number of single stars (default:0)

           nb - number of binary stars (default:0)

           np - number of particles (default:0)

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """  
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
    
    def add_stars(self,id,m,x,y,z,vx,vy,vz,do_key_params=False,do_order=False):
        """
        NAME:

           add_stars

        PURPOSE:

           Add stars to the StarCluster instance

        INPUT:

           id - star id

           m - mass

           x,y,z - positions

           vx,vy,vz - velocities

           kw - stellar evolution tag (optional,for using with NBODY6) (default:0)

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """  
        self.id=np.append(self.id,np.asarray(id))
        self.m=np.append(self.m,np.asarray(m))
        self.x=np.append(self.x,np.asarray(x))
        self.y=np.append(self.y,np.asarray(y))
        self.z=np.append(self.z,np.asarray(z))
        self.vx=np.append(self.vx,np.asarray(vx))
        self.vy=np.append(self.vy,np.asarray(vy))
        self.vz=np.append(self.vz,np.asarray(vz))
        

        self.kw=np.append(self.kw,np.zeros(len(self.id)))

        if do_key_params:
            self.key_params(do_order=do_order)
        
        if len(self.id)!=self.ntot:
            print('Added %i stars to instance' % (len(self.id)-self.ntot))
            self.ntot=len(self.id)

    def add_se(self,kw,logl,logr,ep,ospin):
        """
        NAME:

           add_se

        PURPOSE:

           Add stellar evolution information to stars
           Notes:
            - parameters are common output variables in NBODY6
            - values are never adjusted during unit or coordinate changes

        INPUT:

           kw - stellar evolution tag (for using with NBODY6) 

           logl - log of luminosity

           logr - log of stellar radius

           ep - epoch

           ospin - ospin

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """  
        self.kw=np.asarray(kw)
        self.logl=np.asarray(logl)
        self.logr=np.asarray(logr)
        self.lum=10.0**self.logl
        self.ltot=np.sum(self.lum)
        self.ep=np.asarray(ep)
        self.ospin=np.asarray(ospin)

    def add_bse(self,id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2,ep1,ep2,ospin1,ospin2):
        """
        NAME:

           add_bse

        PURPOSE:

           Add binary star evolution information to stars
           Notes:
            - parameters are common output variables in NBODY6
            - values are never adjusted during unit or coordinate changes

        INPUT:

           id1/id2 - id of star1 and star2

           kw1/kw2 - stellar evolution tags (for using with NBODY6) 

           kcm - stellar evolution tag for binary star

           ecc - eccentricity of binary orbit

           pb - period of binary orbit

           semi - semi major axis of binary orbit

           m1/m2 - masses of binary stars

           logl1/logl2 - luminosities of binary stars

           logr1/logr2 - radii of binary stars

           ep1/ep2 - epochs of binary stars

           ospin1/ospin2 - opsin of binary stars

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)
        """
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

    def add_energies(self,kin,pot,etot):
        """
        NAME:

           add_energies

        PURPOSE:

           Add energy information to stars and calculate total energy and Q for the cluster
           Notes:
            - values are never adjusted during unit or coordinate changes

        INPUT:

           kin - kinetic energy 

           pot - potentail energy

           etot - total energy

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """ 

        self.kin=np.array(kin)
        self.pot=np.array(pot)
        self.etot=np.array(etot)
        self.ektot=np.sum(self.kin)
        self.ptot=np.sum(self.pot)/2.

        if self.ptot==0.:
            self.qvir=0.
        else:
            self.qvir=self.ektot/self.ptot

    def add_orbit(self,xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=None,initialize=False,r0=8.,v0=220.):
        """
        NAME:

           add_orbit

        PURPOSE:

           Add orbital information of the star cluster
           Notes:
            - the input units are assumed to be equal that of self.units

        INPUT:

           xgc,ygc,zgc - cluster's galactocentric position

           vxgc,vygc,vzgc - cluster's galactocentric velocity

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """ 
        if ounits!=None and ounits!=self.units:
              #First convert to realkpc
              if ounits!='realkpc':
                  if ounits=='nbody':
                      xgc*=self.rbar*1000.0
                      ygc*=self.rbar*1000.0
                      zgc*=self.rbar*1000.0
                      vxgc*=self.vstar
                      vygc*=self.vstar
                      vzgc*=self.vstar
                  elif ounits=='galpy':
                      xgc*=r0
                      ygc*=r0
                      zgc*=r0
                      vxgc*=v0
                      vygc*=v0
                      vzgc*=v0
                  elif ounits=='realpc':
                      xgc/=1000.
                      ygc/=1000.
                      zgc/=1000.

                  ounits='realkpc'

              if self.units=='realpc':
                  xgc*=1000.
                  ygc*=1000.
                  zgc*=1000.
              elif self.units=='nbody':
                  xgc*=1000./self.rbar
                  ygc*=1000./self.rbar
                  zgc*=1000./self.rbar
                  vxgc/=self.vstar
                  vygc/=self.vstar
                  vzgc/=self.vstar
              elif self.units=='galpy':
                  xgc/=r0
                  ygc/=r0
                  zgc/=r0
                  vxgc/=v0
                  vygc/=v0
                  vzgc/=v0


        self.xgc=xgc
        self.ygc=ygc
        self.zgc=zgc
        self.rgc=np.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
        self.vxgc=vxgc
        self.vygc=vygc
        self.vzgc=vzgc

        if initialize:
            initialize_orbit(self,from_centre=False)

    def find_centre_of_mass(self):
        """
        NAME:

           find_centre_of_mass

        PURPOSE:

           Find the centre of mass of the cluster

        INPUT:
           None

        OUTPUT:

           xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

        HISTORY:

           2018 - Written - Webb (UofT)

        """       
        xc=np.sum(self.m*self.x)/self.mtot
        yc=np.sum(self.m*self.y)/self.mtot
        zc=np.sum(self.m*self.z)/self.mtot

        vxc=np.sum(self.m*self.vx)/self.mtot
        vyc=np.sum(self.m*self.vy)/self.mtot
        vzc=np.sum(self.m*self.vz)/self.mtot
    
        return xc,yc,zc,vxc,vyc,vzc

    def find_centre_of_density(self,xstart=0.,ystart=0.,zstart=0.,vxstart=0.,vystart=0.,vzstart=0.,indx=None,rmin=0.1,nmax=100,r0=8.,v0=220.):
        """
        NAME:

           find_centre_of_density

        PURPOSE:

           Find the centre of density of the cluster
           Notes - the general framework for this code was taken from Yohai Meiron,
           who made a python adaptation of the centre of density finder in PhiGrape

        INPUT:
           xstart,ystart,zstart - starting position for centre
           vxstart,vystart,vzstart - starting velocity for centre
           rmin - minimum radius of sphere around which to estimate density centre (default: 0.1 pc)
           nmax - maximum number of iterations (default:100)

        OUTPUT:

           xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

        HISTORY:

           2019 - Written - Webb (UofT)
           - with kudos to Yohai Meiron (UofT)

        """  

        #Need to change rmin for pc to self.units:
        if self.units!='nbody':
            rmin/=self.rbar
        elif self.units=='realkpc':
            rmin/=1000.
        elif self.units=='galpy':
            rmin/=(1000.*r0)

        if indx is None:
            indx=np.ones(self.ntot,bool)

        m=self.m[indx]
        x=self.x[indx]-xstart
        y=self.y[indx]-ystart
        z=self.z[indx]-zstart
        vx=self.vx[indx]-vxstart
        vy=self.vy[indx]-vystart
        vz=self.vz[indx]-vzstart

        r=np.sqrt(x**2.+y**2.+z**2.)
        rlim=np.amax(r)

        xdc,ydc,zdc=xstart,ystart,zstart
        vxdc,vydc,vzdc=vxstart,vystart,vzstart

        n=0

        while (rlim > rmin) and (n < nmax):
            r2=x**2.+y**2.+z**2.
            indx = (r2 < rlim**2)
            nc = np.sum(indx)
            mc = np.sum(m[indx])

            if mc == 0:
                xc,yc,zc=0.0,0.0,0.0
                vxc,vyc,vzc=0.0,0.0,0.0
            else:

                xc=np.sum(m[indx]*x[indx])/mc
                yc=np.sum(m[indx]*y[indx])/mc
                zc=np.sum(m[indx]*z[indx])/mc

                vxc=np.sum(m[indx]*vx[indx])/mc
                vyc=np.sum(m[indx]*vy[indx])/mc
                vzc=np.sum(m[indx]*vz[indx])/mc       

            if((mc > 0) and (nc > 100)):
                x-=xc
                y-=yc
                z-=zc
                xdc+=xc
                ydc+=yc
                zdc+=zc

                vx-=vxc
                vy-=vyc
                vz-=vzc                
                vxdc+=vxc
                vydc+=vyc
                vzdc+=vzc

            else:
                break
            rlim *= 0.8
            n+=1
 
        self.xc,self.yc,self.zc=xdc,ydc,zdc
        self.vxc,self.vyc,self.vzc=vxdc,vydc,vzdc

        return xdc,ydc,zdc,vxdc,vydc,vzdc

    def find_centre(self,xstart=0.,ystart=0.,zstart=0.,vxstart=0.,vystart=0.,vzstart=0.,indx=None,nsigma=1.,nsphere=100,density=True,rmin=0.1,nmax=100,r0=8.,v0=220.):
        """
        NAME:

           find_centre

        PURPOSE:

           Find the centre of mass of the cluster
           Notes:
            - The routine first works to identify a sphere of nsphere stars around the centre in which
            to perform the C.O.M calculation. This step prevents long tidal tails from affecting the 
            calculation

        INPUT:
           xstart,ystart,zstart - starting position for centre
           vxstart,vystart,vzstart - starting velocity for centre
           indx - subset of stars to use when finding center
           nsigma - number of standard deviations to within which to keep stars
           nsphere - number of stars in centre sphere (default:100)
           density - use Yohai Meiron's centre of density calculator instead (Default: True)
           if density:
               - rmin - minimum radius to start looking for stars
               - nmax - maximum number of iterations to find centre
           r0,v0 - For converting to and from galpy units (Default: 8., 220.)

        OUTPUT:

           xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

        HISTORY:

           2019 - Written - Webb (UofT)

        """

        if indx is None:
            indx=np.ones(self.ntot,bool)
        elif np.sum(indx)==0.:
            print('NO SUBSET OF STARS GIVEN')
            return 0.,0.,0.,0.,0.,0.

        if density: 
           xc,yc,zc,vxc,vyc,vzc=self.find_centre_of_density(xstart=xstart,ystart=ystart,zstart=zstart,vxstart=vxstart,vystart=vystart,vzstart=vzstart,indx=indx,rmin=rmin,nmax=nmax,r0=r0,v0=v0)
        else:

            x=self.x[indx]-xstart
            y=self.y[indx]-ystart
            z=self.z[indx]-zstart
            r=np.sqrt(x**2.+y**2.+z**2.)
            i_d=self.id[indx]

            while len(r)>nsphere:
                sigma=nsigma*np.std(r)
                indx=(r<sigma)

                if len(r[indx])>nsphere:
                    i_d=i_d[indx]
                    x=x[indx]-np.mean(x[indx])
                    y=y[indx]-np.mean(y[indx])
                    z=z[indx]-np.mean(z[indx])
                    r=np.sqrt(x*x+y*y+z*z)
                else:
                    break

            #Find centre of mass and velocity of inner stars:
            indx=np.in1d(self.id,i_d)

            xc=np.sum(self.m[indx]*self.x[indx])/np.sum(self.m[indx])
            yc=np.sum(self.m[indx]*self.y[indx])/np.sum(self.m[indx])
            zc=np.sum(self.m[indx]*self.z[indx])/np.sum(self.m[indx])

            vxc=np.sum(self.m[indx]*self.vx[indx])/np.sum(self.m[indx])
            vyc=np.sum(self.m[indx]*self.vy[indx])/np.sum(self.m[indx])
            vzc=np.sum(self.m[indx]*self.vz[indx])/np.sum(self.m[indx])

        if self.origin=='galaxy':
            self.xgc,self.ygc,self.zgc=xc,yc,zc
            self.vxgc,self.vygc,self.vzgc=vxc,vyc,vzc

            self.xc,self.yc,self.zc=0.0,0.0,0.0
            self.vxc,self.vyc,self.vzc=0.0,0.0,0.0

        elif self.origin=='cluster':
            self.xc,self.yc,self.zc=xc,yc,zc
            self.vxc,self.vyc,self.vzc=vxc,vyc,vzc

        return xc,yc,zc,vxc,vyc,vzc

    def add_rotation(self,qrot):
        """
        NAME:

           add_rotation

        PURPOSE:

           Introduce a degree of rotation into the cluster

        INPUT:

           qrot - Degree of rotation, such that qrot% of stars with negative v_theta are switched positive

        OUTPUT:

            None

        HISTORY:

           2019 - Written - Webb (UofT)

        """   

        r,theta,z,vr,vtheta,vz=npy.cyl_coords(cluster)

        indx=(vtheta < 0.)
        rindx=(np.random.rand(cluster.ntot) < q)

        vtheta[indx*rindx]=np.sqrt(vtheta[indx*rindx]*vtheta[indx*rindx])
        cluster.x,cluster.y,cluster.z=bovy_coords.cyl_to_rect(r,theta,z)
        cluster.vx,cluster.vy,cluster.vz=bovy_coords.cyl_to_rect_vec(vr,vtheta,vz,theta)

    def add_actions(self,JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz):
        """
        NAME:

           add_actions

        PURPOSE:

           Add action angle values to the cluster instance

        INPUT:

           JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz

        OUTPUT:

            None

        HISTORY:

           2019 - Written - Webb (UofT)

        """  
        self.JR,self.Jphi,self.Jz=JR,Jphi,Jz
        self.OR,self.Ophi,self.Oz=OR,Ophi,Oz
        self.TR,self.Tphi,self.Tz=TR,Tphi,Tz


    def key_params(self,do_order=False):
        """
        NAME:

           key_params

        PURPOSE:

           Find key parameters of the cluster (mass,luminosity,r10,r50,rh10,rh50)

        INPUT:

           do_order - Perform the time consuming task of ordering stars based on radius to find r10,r50, etc. (default:False)

        OUTPUT:

            None

        HISTORY:

           2018 - Written - Webb (UofT)

        """  
        self.v=np.sqrt(self.vx**2.+self.vy**2.+self.vz**2.)
        self.vpro=np.sqrt(self.vx**2.+self.vy**2.)

        self.r=np.sqrt(self.x**2.+self.y**2.+self.z**2.)
        self.rpro=np.sqrt(self.x**2.+self.y**2.)
        self.mtot=np.sum(self.m)
        self.mmean=np.mean(self.m)
        self.rmean=np.mean(self.r)
        self.rmax=np.max(self.r)


        #Radially order the stars to find half-mass radius
        if do_order:
            self.rorder=np.argsort(self.r)
            if self.projected:
                self.rproorder=np.argsort(self.rpro)
 
        if self.rorder is not None:
            msum=np.cumsum(self.m[self.rorder])
            indx=(msum >= 0.5*self.mtot)
            self.rm=self.r[self.rorder[indx]][0]
            indx=(msum >= 0.1*self.mtot)
            self.r10=self.r[self.rorder[indx]][0]

        if self.projected and self.rproorder is not None:
            msum=np.cumsum(self.m[self.rproorder])
            indx=(msum >= 0.5*self.mtot)
            self.rmpro=self.r[self.rproorder[indx]][0]
            indx=(msum >= 0.1*self.mtot)
            self.r10pro=self.r[self.rproorder[indx]][0]
        else:
            self.rmpro=0.
            self.r10pro=0.

        if len(self.logl)>0 and self.rorder is not None:
            lsum=np.cumsum(self.lum[self.rorder])
            indx=(lsum >= 0.5*self.ltot)
            self.rh=self.r[self.rorder[indx]][0]
            indx=(lsum >= 0.1*self.ltot)
            self.rh10=self.r[self.rorder[indx]][0]

            if self.projected and self.rproorder is not None:
                lsum=np.cumsum(self.lum[self.rproorder])
                indx=(lsum >= 0.5*self.ltot)
                self.rhpro=self.rpro[self.rproorder[indx]][0]
                indx=(lsum >= 0.1*self.ltot)
                self.rh10pro=self.rpro[self.rproorder[indx]][0]
            else:
                self.rhpro=0.
                self.rh10pro=0.

    def to_realpc(self,do_key_params=True):
        """
        NAME:

           to_realpc

        PURPOSE:

           Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and km/s

        INPUT:

           None

        OUTPUT:

            None

        HISTORY:

           2018 - Written - Webb (UofT)

        """    
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

        if do_key_params:
            self.key_params()

    def to_realkpc(self,do_key_params=True,r0=8.,v0=220.):
        """
        NAME:

           to_realkpc

        PURPOSE:

           Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and km/s

        INPUT:

           None

        OUTPUT:

            None

        HISTORY:

           2018 - Written - Webb (UofT)

        """ 
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

        if do_key_params:
            self.key_params()

    def to_nbody(self,do_key_params=True,r0=8.,v0=220.):
        """
        NAME:

           to_nbody

        PURPOSE:

           Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Nbody units
           Notes:
            - requires that self.zmbar, self.rbar, self.vstar are set (defaults are 1 in add_nbody6)

        INPUT:

           None

        OUTPUT:

            None

        HISTORY:

           2018 - Written - Webb (UofT)

        """     
        if self.units=='galpy': self.to_realkpc(do_key_params=False)
        if self.units!='realpc': self.to_realpc(do_key_params=False)

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

        if do_key_params:
            self.key_params()

    def to_galpy(self,do_key_params=True,r0=8.,v0=220.):
        """
        NAME:

           to_galpy

        PURPOSE:

           Convert stellar positions/velocities, centre of mass, and orbital position and velocity to galpy units

        INPUT:

           None

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """ 
        if self.units!='realkpc' and self.units!='galpy':
            self.to_realkpc(do_key_params=False)
        
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

        if do_key_params:
            self.key_params()


    def to_units(self,units,do_key_params=True,r0=8.,v0=220.):
        """
        NAME:

           to_units

        PURPOSE:

           Convert stellar positions/velocities, centre of mass, and orbital position and velocity to user defined units

        INPUT:

           units - 'nbody','realpc','realkpc','galpy'

        OUTPUT:

            None

        HISTORY:

           2018 - Written - Webb (UofT)

        """    

        if units=='nbody':
            self.to_nbody(do_key_params=do_key_params)
        elif units=='galpy':
            self.to_galpy(do_key_params=do_key_params,r0=r0,v0=v0)
        elif units=='realpc':
            self.to_realpc(do_key_params=do_key_params)
        elif units=='realkpc':
            self.to_realkpc(do_key_params=do_key_params)


    def to_centre(self,do_order=False,do_key_params=False):
        """
        NAME:

           to_centre

        PURPOSE:

           Shift coordinates such that origin is the centre of mass

        INPUT:

           None

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)

        """ 
        if self.origin!='centre':

            if self.origin=='galaxy':
                self.to_cluster(do_key_params=False)
        
            self.x-=self.xc
            self.y-=self.yc
            self.z-=self.zc
            self.vx-=self.vxc
            self.vy-=self.vyc
            self.vz-=self.vzc
            
            self.origin='centre'

            if do_key_params:
                self.key_params(do_order=do_order)


    def to_cluster(self,do_order=False,do_key_params=False):
        """
        NAME:

           to_cluster

        PURPOSE:

           Shift coordinates to clustercentric reference frame

        INPUT:

           to_centre - shift to centre of mass as well (default: False)

        OUTPUT:

            None

        HISTORY:

           2018 - Written - Webb (UofT)
        """
        if self.origin!='cluster':

            if self.origin=='centre':
              self.x+=self.xc
              self.y+=self.yc
              self.z+=self.zc
              self.vx+=self.vxc
              self.vy+=self.vyc
              self.vz+=self.vzc
            elif self.origin=='galaxy':
              self.x-=self.xgc
              self.y-=self.ygc
              self.z-=self.zgc
              self.vx-=self.vxgc
              self.vy-=self.vygc
              self.vz-=self.vzgc
            self.origin='cluster'

            if do_key_params:
                self.key_params(do_order=do_order)

    def to_galaxy(self,do_order=False,do_key_params=False):
        """
        NAME:

           to_galaxy

        PURPOSE:

           Shift coordinates to galactocentric reference frame

        INPUT:

           None

        OUTPUT:

           None

        HISTORY:

           2018 - Written - Webb (UofT)
        """
        if self.origin!='galaxy':
            if self.origin=='centre':
                self.to_cluster(do_key_params=False)
            
            self.x+=self.xgc
            self.y+=self.ygc
            self.z+=self.zgc
            self.vx+=self.vxgc
            self.vy+=self.vygc
            self.vz+=self.vzgc
            

            self.origin='galaxy'
        
            if do_key_params:
                self.key_params(do_order=do_order)

    def to_sky(self):
        """
        NAME:

           to_sky

        PURPOSE:

           Calculate on-sky position, proper motion, and radial velocity of cluster

        INPUT:

           None

        OUTPUT:

            None

        HISTORY:

           2019 - Written - Webb (UofT)

        """    
        units0,origin0=self.units,self.origin

        self.to_galaxy()
        self.to_realkpc()

        x0,y0,z0=bovy_coords.galcenrect_to_XYZ(self.x,self.y,self.z,Xsun=8.,Zsun=0.025).T
        vx0,vy0,vz0=bovy_coords.galcenrect_to_vxvyvz(self.vx,self.vy,self.vz,Xsun=8.,Zsun=0.025,vsun=[-11.1,244.,7.25]).T

        self.vlos=(vx0*x0+vy0*y0+vz0*z0)/np.sqrt(x0**2.+y0**2.+z0**2.)

        l0,b0,self.dist=bovy_coords.XYZ_to_lbd(x0,y0,z0,degree=True).T
        self.ra,self.dec=bovy_coords.lb_to_radec(l0,b0,degree=True).T

        vr0,pmll0,pmbb0=bovy_coords.vxvyvz_to_vrpmllpmbb(vx0,vy0,vz0,l0,b0,self.dist,degree=True).T
        self.pmra,self.pmdec=bovy_coords.pmllpmbb_to_pmrapmdec(pmll0,pmbb0,l0,b0,degree=True).T

        x0,y0,z0=bovy_coords.galcenrect_to_XYZ(self.xgc,self.ygc,self.zgc,Xsun=8.,Zsun=0.025)
        vx0,vy0,vz0=bovy_coords.galcenrect_to_vxvyvz(self.vxgc,self.vygc,self.vzgc,Xsun=8.,Zsun=0.025,vsun=[-11.1,244.,7.25])

        self.vlosgc=(vx0*x0+vy0*y0+vz0*z0)/np.sqrt(x0**2.+y0**2.+z0**2.)

        l0,b0,self.distgc=bovy_coords.XYZ_to_lbd(x0,y0,z0,degree=True)
        self.ragc,self.decgc=bovy_coords.lb_to_radec(l0,b0,degree=True)

        vr0,pmll0,pmbb0=bovy_coords.vxvyvz_to_vrpmllpmbb(vx0,vy0,vz0,l0,b0,self.distgc,degree=True)
        self.pmragc,self.pmdecgc=bovy_coords.pmllpmbb_to_pmrapmdec(pmll0,pmbb0,l0,b0,degree=True)

       
        if origin0=='centre':
            self.to_centre()
        elif origin0=='cluster':
            self.to_cluster()

        if units0=='realpc':
            self.to_realpc()
        elif units0=='nbody':
            self.to_nbody()
        elif units0=='galpy':
            self.to_galpy()

    def to_origin(self,origin,do_order=False):
        """
        NAME:

           to_origin

        PURPOSE:

           Shift cluster to origin as defined by keyword

        INPUT:

           origin - accepts 'cluster','centre','galaxy'

        OUTPUT:

            None

        HISTORY:

           2019 - Written - Webb (UofT)

        """    

        if origin=='centre':
            self.to_centre(do_order=do_order,do_key_params=False)
        elif origin=='cluster':
            self.to_cluster(do_order=do_order,do_key_params=False)
        elif origin=='galaxy':
            self.to_galaxy(do_order=do_order,do_key_params=False)
        elif origin=='to_sky':
            self.to_sky()

    #Directly call from functions.py and profiles.py (see respective files for documenation):
    def energies(self,specific=True,i_d=None,full=True,parallel=False):
        ek,pot,etot=energies(self,specific=specific,i_d=i_d,full=full,parallel=parallel)
        self.add_energies(ek,pot,etot)

    def rlagrange(self,nlagrange=10,projected=False):
        self.rn=rlagrange(self,nlagrange=nlagrange,projected=projected)
        
    def rvirial(self,H=70.0,Om=0.3,overdens=200.,nrad=20,projected=False):
        self.rv=rvirial(self,H=H,Om=Om,overdens=overdens,nrad=nrad,projected=projected)

    def virial_radius(self,projected=False):
        self.rv=virial_radius(self,projected=projected)      

    def rtidal(self,pot=MWPotential2014 ,rtiterate=0,rgc=None,r0=8.,v0=220.):
        self.rt=rtidal(self,pot=pot,rtiterate=rtiterate,rgc=rgc,r0=r0,v0=v0)

    def rlimiting(self,pot=MWPotential2014 ,rgc=None,r0=8.,v0=220.,nrad=20,projected=False,obs_cut=False,plot=False,**kwargs):
        self.rl=rlimiting(self,pot=pot,rgc=rgc,r0=r0,v0=v0,nrad=nrad,projected=projected,obs_cut=obs_cut,plot=plot,**kwargs)

    def mass_function(self,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None,plot=False,**kwargs):
        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(self,mmin=mmin,mmax=mmax,nmass=nmass,rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut,plot=plot,title='GLOBAL',**kwargs)
        self.alpha=alpha 
        self.rn=rlagrange(self,nlagrange=10,projected=projected)
        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=mass_function(self,mmin=mmin,mmax=mmax,nmass=nmass,rmin=self.rn[4],rmax=self.rn[6],vmin=vmin,vmax=vmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut,plot=plot,title='AT R50',**kwargs)
        self.alpha50=alpha 

    def eta_function(self,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None,plot=False,**kwargs):
        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(self,mmin=mmin,mmax=mmax,nmass=nmass,rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut,plot=plot,**kwargs)
        self.eta=eta

    def alpha_prof(self,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,projected=False,obs_cut=None,plot=False,**kwargs):
        lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha=alpha_prof(self,mmin=mmin,mmax=mmax,nmass=nmass,rmin=rmin,rmax=rmax,nrad=nrad,vmin=vmin,vmax=vmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut,plot=plot,**kwargs)
        self.dalpha=dalpha

    def calc_actions(self,pot=MWPotential2014,r0=8.,v0=220.,**kwargs):
        JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz=calc_actions(self,pot=pot,r0=r0,v0=v0,**kwargs)
        self.add_actions(JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz)

    def vcirc_prof(self,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
        rprof,vcprof,rvmax,vmax=vcirc_prof(self,mmin=mmin,mmax=mmax,rmin=rmin,rmax=rmax,nrad=nrad,vmin=vmin,vmax=vmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,indx=indx,projected=projected,obs_cut=obs_cut,plot=plot,**kwargs)
        self.rvmax=rvmax
        self.vmax=vmax

def sub_cluster(cluster,rmin=None,rmax=None,mmin=None,mmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,indx=[None],projected=False,reset_centre=False,reset_nbody_scale=False,reset_nbody_mass=False,reset_nbody_radii=False):
    """
    NAME:

       sub_cluster

    PURPOSE:

       Extract a sub population of stars from a StarCluster
       Notes:
        -- automatically moves cluster to centre of mass, so all constraints are in clustercentric coordinates and current cluster.units

    INPUT:

       rmin/rmax - minimum and maximum stellar radii

       mmin/mmax - minimum and maximum stellar mass

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (default:False)

       reset_centre - re-calculate cluster centre after extraction (default:False)

    OUTPUT:

       instance

    HISTORY:

       2018 - Written - Webb (UofT)

    """     

    units0,origin0=cluster.units,cluster.origin
    cluster.to_centre()

    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if emin==None and emax!=None:
        eindx=(cluster.etot<=emax)
    elif emin!=None and emax==None:
        eindx=(cluster.etot>=emin)
    elif emin!=None and emax!=None:
        eindx=(cluster.etot<=emax) * (cluster.etot>=emin)
    else:
        eindx=(cluster.id > -1)

    if None in indx:
        indx=(cluster.id > -1)

    indx*=(r>=rmin) * (r<=rmax) * (cluster.m>=mmin) * (cluster.m<=mmax) * (v>=vmin) * (v<=vmax) * (cluster.kw>=kwmin) * (cluster.kw<=kwmax) * eindx

    if np.sum(indx)>0:
        subcluster=StarCluster(len(cluster.id[indx]),cluster.tphys,units=cluster.units,origin=cluster.origin)
        subcluster.add_stars(cluster.id[indx],cluster.m[indx],cluster.x[indx],cluster.y[indx],cluster.z[indx],cluster.vx[indx],cluster.vy[indx],cluster.vz[indx])
        subcluster.kw=cluster.kw[indx]

        subcluster.zmbar=cluster.zmbar
        subcluster.rbar=cluster.rbar
        subcluster.vstar=cluster.vstar
        subcluster.tstar=cluster.tstar
        subcluster.projected=cluster.projected

        if len(cluster.logl)>0:
            subcluster.add_se(cluster.kw[indx],cluster.logl[indx],cluster.logr[indx],cluster.ep[indx],cluster.ospin[indx])
        if len(cluster.id2)>0:
            bindx=np.in1d(cluster.id1,cluster.id[indx])
            subcluster.add_bse(cluster.id1[bindx],cluster.id2[bindx],cluster.kw1[bindx],cluster.kw2[bindx],cluster.kcm[bindx],cluster.ecc[bindx],cluster.pb[bindx],cluster.semi[bindx],cluster.m1[bindx],cluster.m2[bindx],cluster.logl1[bindx],cluster.logl2[bindx],cluster.logr1[bindx],cluster.logr2[bindx],cluster.ep1[bindx],cluster.ep2[bindx],cluster.ospin1[bindx],cluster.ospin2[bindx])
        if len(cluster.etot)>0:
            subcluster.add_energies(cluster.kin[indx],cluster.pot[indx],cluster.etot[indx])

        if reset_centre:
            subcluster.add_orbit(cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc,cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc)
            subcluster.xc,subcluster.yc,subcluster.zc=0.0,0.0,0.0
            subcluster.vxc,subcluster.vyc,subcluster.vzc=0.0,0.0,0.0
            subcluster.xc,subcluster.yc,subcluster.zc,subcluster.vxc,subcluster.vyc,subcluster.vzc=subcluster.find_centre(0.0,0.0,0.0)

        else:
            subcluster.add_orbit(cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc)
            subcluster.xc,subcluster.yc,subcluster.zc=cluster.xc,cluster.yc,cluster.zc
            subcluster.vxc,subcluster.vyc,subcluster.vzc=cluster.vxc,cluster.vyc,cluster.vzc

        if reset_nbody_scale or reset_nbody_mass or reset_nbody_radii:
            subcluster.to_realpc()
            subcluster.key_params(do_order=True)

            if reset_nbody_scale or reset_nbody_mass:
                subcluster.zmbar=subcluster.mtot
            if reset_nbody_scale or reset_nbody_radii:
                subcluster.rbar=4.0*subcluster.rm/3.0

            subcluster.vstar=0.06557*np.sqrt(subcluster.zmbar/subcluster.rbar)
            subcluster.tstar=subcluster.rbar/subcluster.vstar
        else:

            subcluster.key_params(do_order=True)

    else:
        subcluster=StarCluster(0,cluster.tphys)

    cluster.to_origin(origin0)
    cluster.to_units(units0)

    if subcluster.ntot > 0:
        subcluster.to_origin(origin0)
        subcluster.to_units(units0)

    return subcluster

def kwtypes():
    """
    NAME:

       kwtypes

    PURPOSE:

       Print legend for converting kwtype (from NBODY6) to stellar evolution type

    INPUT:

       None

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)

    """   

    print (dedent("""\
        *       Stellar evolution types
        *       ***********************
        *
        *       ---------------------------------------------------------------------
        *       0       Low main sequence (M < 0.7).
        *       1       Main sequence.
        *       2       Hertzsprung gap (HG).
        *       3       Red giant.
        *       4       Core Helium burning.
        *       5       First AGB.
        *       6       Second AGB.
        *       7       Helium main sequence.
        *       8       Helium HG.
        *       9       Helium GB.
        *      10       Helium white dwarf.
        *      11       Carbon-Oxygen white dwarf.
        *      12       Oxygen-Neon white dwarf.
        *      13       Neutron star.
        *      14       Black hole.
        *      15       Massless supernova remnant.
        *       ---------------------------------------------------------------------
        *
        *       Binary types
        *       ************
        *
        *       ---------------------------------------------------------------------
        *       0       Standard case.
        *      -1       Chaotic (option 27 = 2).
        *      -2       Continuous circularizing (option 27 = 2).
        *       9       Sequential circularization (option 27 = 1).
        *      10       Circularized.
        *      11       First Roche stage (option 34 = 1/2).
        *      12       End of first Roche stage.
        *      13       Start of second Roche stage.
        *      xx       Further Roche stages.
        *       ---------------------------------------------------------------------
        *



        """))
