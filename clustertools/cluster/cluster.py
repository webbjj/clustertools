""" The StarCluster class and key internal functions

"""

__author__ = "Jeremy J Webb"

__all__ = [
    "StarCluster",
    "sub_cluster",
    "overlap_cluster",
]

import numpy as np
try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion
from textwrap import dedent
from ..analysis.orbits import *
from ..analysis.functions import *
from .operations import *
from ..tidaltail.tails import *
from copy import copy
from galpy.orbit import Orbit
from ..util.constants import *

try:
    from amuse.lab import *
    import amuse.units.units as u
    from amuse.datamodel import Particles
    from ..util.units import _convert_amuse,_convert_length,_convert_velocity
    _hasamuse=True
except:
    _hasamuse=False

class StarCluster(object):
    """ A class that represents a star cluster population that ooperations and functions can be performed on
    
    Parameters
    ----------
    tphys : float
        Time (units not necessary) associated with the population (default: 0)
    units : str
        Units of stellar positions and velocties. Options include 'pckms','pcmyr',
        'kpckms','kpcgyr','radec','nbody',and 'galpy'. For 'pckms' and 'kpckms', 
        stellar velocities are assumed to be km/s. (default: None)
    origin : str
        Origin of coordinate system within which stellar positions and velocities are defined. 
        Options include 'centre', 'cluster', 'galaxy' and 'sky'. Note that 'centre' corresponds
        to the systems centre of density or mass (as calculated by clustertools) while 'cluster' corresponds
        to the orbital position of the cluster's initial centre of mass given 'tphys'. In some cases
        these two coordinate systems may be the same. (default: None)
    ctype : str
        Code type used to generate the cluster population. Current options include 'snapshot',
        'nbody6','nbody6se','nemo','gyrfalcon','snaptrim', amuse','clustertools','snapauto'. The parameter
        informs clustertools how to load the stellar popuation and advance to the next snapshot.
        (default: 'snapshot')
    projected : bool
        return projected values instead of 3D value. (default: False)

    Other Parameters
    ----------------
    sfile : str
        file containing single star data
    bfile : str
        file contain binary star data
    ssefile : str
        file containing single star evolution data
    bsefile : str
        file contain binary star evolution data
    tfile : str
        name of file contain tail star data
    ofile : str
        orbit file
    ofilename : str
        orbit filename if ofile is not given
    orbit : str
        galpy orbit instance
    ounits : str
        {'pckms','pcmyr','kpckms','kpcgyr','radec','nbody','galpy'} units of orbital information (else assumed equal to StarCluster.units)
    ro : float
        distance to the Galactic centre (Default: solar_ro)
    vo : float
        circular velocity at ro (Default: solar_vo)
    zo : float
        Sun's distance above the Galactic plane (default: solar_zo)
    solarmotion : float
        array representing U,V,W of Sun (default: solar_motion)
    nsnap : int
        if a specific snapshot is to be read in instead of starting from zero
    nzfill : int
        value for zfill when reading and writing snapshots (default: 5)
    snapbase : str
        string of characters in filename before nsnap (default: '')
    snapend : str
        string of character in filename after nsnap (default: '.dat')
    snapdir : str
        string name of directory of snapshots if different than wdir (Default: '')
    delimiter : str
        choice of delimiter when reading ascii/csv files (default: ',')
    wdir : str
        working directory of snapshots if not current directory
    intialize : bool
        initialize a galpy orbit after reading in orbital information (default: False)
    advance : bool
        set to True if this a snapshot that has been advanced to from an initial one? (default: False)
    centre_method: str
        {None,'orthographic','VandeVen'} method to convert to clustercentric coordinates when units are in radec (Default: None)
    give : str
        set what parameters are read in from nemo/gyrfalcon (default: 'mxv')
        Currently only accepts 'mxvpqael' as an alternative.
    History
    -------
    2018 - Written - Webb (UofT)

    """

    def __init__(
        self,
        tphys=0.0,
        units=None,
        origin=None,
        ctype="snapshot",
        projected=False,
        **kwargs,
    ):

        # Age of cluster
        self.tphys = tphys
        self.dt = None

        # Units and origin
        self.units = units
        self.bunits = units
        self.origin = origin
        self.units_init=units
        self.origin_init=origin

        # Cluster Simulation Type
        self.ctype = ctype

        # Return projected values only
        self.projected = projected

        # Kwargs
        self.nsnap = int(kwargs.get("nsnap", 0))
        self.delimiter = kwargs.get("delimiter", None)
        self.wdir = kwargs.get("wdir", "./")
        self.nzfill = int(kwargs.get("nzfill", 5))
        self.snapbase = kwargs.get("snapbase", "")
        self.snapend = kwargs.get("snapend", ".dat")
        self.snapdir = kwargs.get("snapdir", "")
        self.skiprows = kwargs.get("skiprows", 0)
        self.sfile = kwargs.get("sfile", None)
        self.bfile = kwargs.get("bfile", None)
        self.ssefile = kwargs.get("ssefile", None)
        self.bsefile = kwargs.get("bsefile", None)
        self.ofile = kwargs.get("ofile", None)
        self.ofilename = kwargs.get("ofilename", None)
        self.orbit = kwargs.get("orbit", None)


        self._ro=kwargs.get('ro',solar_ro)
        self._vo=kwargs.get('vo',solar_vo)
        self._zo=kwargs.get('zo',solar_zo)
        self._solarmotion=kwargs.get('solarmotion',solar_motion)

        self.give=kwargs.get('give','mxv')

        self.centre_method = kwargs.get("centre_method", None)

        # Total Number of Stars + Binaries in the cluster
        self.ntot = 0
        self.nb = 0

        # variables for add_stars
        self.id = np.array([])
        self.m = np.array([])
        self.x = np.array([])
        self.y = np.array([])
        self.z = np.array([])
        self.vx = np.array([])
        self.vy = np.array([])
        self.vz = np.array([])
        self.m0 = np.array([])

        # variables for add_binary_stars
        self.mb1 = np.array([])
        self.xb1 = np.array([])
        self.yb1 = np.array([])
        self.zb1 = np.array([])
        self.vxb1 = np.array([])
        self.vyb1 = np.array([])
        self.vzb1 = np.array([])
        self.mb2 = np.array([])
        self.xb2 = np.array([])
        self.yb2 = np.array([])
        self.zb2 = np.array([])
        self.vxb2 = np.array([])
        self.vyb2 = np.array([])
        self.vzb2 = np.array([])

        #If using Nbodypp
        self.rhos=np.array([])


        # variables for add_nbody
        self.zmbar = 1.0
        self.rbar = 1.0
        self.vbar = 1.0
        self.tbar = 1.0


        # variables for centre of cluster
        self.xc = 0.0
        self.yc = 0.0
        self.zc = 0.0
        self.vxc = 0.0
        self.vyc = 0.0
        self.vzc = 0.0

        # variable for galpy orbits
        self.orbits= None

        # variables for orbital position and kinematics
        if self.orbit is None:
            self.xgc = 0.0
            self.ygc = 0.0
            self.zgc = 0.0
            self.vxgc = 0.0
            self.vygc = 0.0
            self.vzgc = 0.0
        elif units=='radec':
            self.xgc = self.orbit.ra()
            self.ygc = self.orbit.dec()
            self.zgc = self.orbit.dist()
            self.vxgc = self.orbit.pmra()
            self.vygc = self.orbit.pmdec()
            self.vzgc = self.orbit.vlos()   
        else:
            self.xgc = self.orbit.x()
            self.ygc = self.orbit.y()
            self.zgc = self.orbit.z()
            self.vxgc = self.orbit.vx()
            self.vygc = self.orbit.vy()
            self.vzgc = self.orbit.vz()

        # variable for cluster's on-sky coordinates
        self.ra = np.array([])
        self.dec = np.array([])
        self.dist = np.array([])
        self.pmra = np.array([])
        self.pmdec = np.array([])
        self.vlos = np.array([])

        if self.orbit is None:
            self.ra_gc = 0.0
            self.dec_gc = 0.0
            self.dist_gc = 0.0
            self.pmra_gc = 0.0
            self.pmdec_gc = 0.0
            self.vlos_gc = 0.0
        else:
            self.ra_gc = self.orbit.ra()
            self.dec_gc = self.orbit.dec()
            self.dist_gc = self.orbit.dist()
            self.pmra_gc = self.orbit.pmra()
            self.pmdec_gc = self.orbit.pmdec()
            self.vlos_gc = self.orbit.vlos()

        self.ra_c = 0.0
        self.dec_c = 0.0
        self.dist_c = 0.0
        self.pmra_c = 0.0
        self.pmdec_c = 0.0
        self.vlos_c = 0.0   

        # variables for add_nbody6
        # Number of stars in the core
        self.nc = 0
        # Core radius
        self.rc = 0
        # Distance scaling parameter
        self.rbar = 1.
        self.rbar_su=1.
        self.rbar_au=1.
        # Tidal limit from NBODY6 (not neccesarily a true tidal radius)
        self.rtide = 0.
        # Center of mass of cluster (x,yz)
        self.xc = 0.
        self.yc = 0.
        self.zc = 0.
        self.xcn = None
        self.ycn = None
        self.zcn = None
        # Mass scaling parameter
        self.zmbar = 1.
        # Velocity scaling parameter
        self.vbar = 1.
        # Time scaling parameter
        self.tbar = 1.
        self.tbar_days=1.
        # Scale radius of cluster
        self.rscale = 1.
        # Number of single stars
        self.ns = 0
        # Number of binary stars
        self.nb = 0
        # Number of particles (from NBODY6 when tidal tail is being integrated)
        self.n_p = 0

        # variables for add_sse (stellar evolution information)
        self.kw = np.array([])
        self.logl = np.array([])
        self.logr = np.array([])
        self.ep = np.array([])
        self.ospin = np.array([])
        self.lum = np.array([])

        # variables for add_bse (binary star evolution information)
        self.id1 = np.array([])
        self.id2 = np.array([])
        self.kw1 = np.array([])
        self.kw2 = np.array([])
        self.kcm = np.array([])
        self.ecc = np.array([])
        self.pb = np.array([])
        self.semi = np.array([])
        self.m1 = np.array([])
        self.m2 = np.array([])
        self.m01 = np.array([])
        self.m02 = np.array([])
        self.logl1 = np.array([])
        self.logl2 = np.array([])
        self.logr1 = np.array([])
        self.logr2 = np.array([])
        self.ep1 = np.array([])
        self.ep2 = np.array([])
        self.ospin1 = np.array([])
        self.ospin2 = np.array([])
        self.npop1 = np.array([])
        self.npop2 = np.array([])

        # variables of energies
        self.kin = np.array([])
        self.pot = np.array([])
        self.etot = np.array([])

        # Lagrange Radii,10% lagrage radius, half-mass radius, limiting radius, tidal radius, and virial radius
        self.rn = None
        self.r10 = None
        self.r10pro=None
        self.rm = None
        self.rmpro = None
        self.rh = None
        self.rhpro = None
        self.rl = None
        self.rt = None
        self.rv = None

        #3D and projected order of stars with respect to origin
        self.rorder = None
        self.rorder_origin = None
        self.rproorder = None

        # Additional variables for operation and function calls
        self.trelax = None
        self.trh = None
        self.trc = None
        self.qv = None
        self.alpha = None
        self.eta = None
        self.rvmax = None
        self.vmax = None

        #For use with multiple populations
        self.npop = np.array([])

        #For use with extended nemo/gyrfalcon output
        if self.give == 'mxvpqael':
            self.gyrpot=np.array([])
            self.gyrq=np.array([])
            self.gyracc=np.array([])
            self.eps=np.array([])
            self.gyrlev=np.array([])
        elif self.give =='mxve':
            self.eps=np.array([])

        #For use with HDF5
        self.hdf5=False
        self.ngroups=0
        self.ngroup=0

    def add_stars(
        self, x, y, z, vx, vy, vz,m=None,id=None,m0=None,npop=None,nb=0,sortstars=False,analyze=True
    ):
        """Add stars to StarCluster.

        Parameters
        ----------
        x,y,z: float
            stellar positions. Input is assumed to be in cartesian coordinates unless self.units=='radec' 
            and self.origin=='sky', then positions are assumed to be ra,dec,dist (degrees, degrees, kpc)
        vx,vy,vz: float
            atellar velocities. Input is assumed to be in cartesian coordinates unless self.units=='radec' 
            and self.origin=='sky', then positions are assumed to be pmra,pmdec,vlos (mas/yr, mas/yr, km/s)
        m: float int
            stellar mass
        id: int 
            star id
        m0: float int
            initial stellar mass
        npop : int
            population number, for use with multiple populations
        nb: int
            number of binary stars, which are assumed to be the first 2*nb stars in the array (default:0)
        sortstars: bool
            order stars by radius (default: False)
        analyze : bool
            perform analysis on cluster after stars are added

        Notes
        -----
        History:

            - 2018 - Written - Webb (UofT)

        """

        #Check for float entries
        params=([x,y,z,vx,vy,vz,m,id,m0,npop])
        nfloat=np.zeros(len(params),dtype=bool)
        npmax=0
        for i,p in enumerate(params):
            if p is not None:
                if _hasamuse:
                    if isinstance(p,ScalarQuantity):
                        nfloat[i]=True
                        npmax=int(np.maximum(npmax,1))    

                if isinstance(p,float) or isinstance(p,int):
                    nfloat[i]=True
                    npmax=int(np.maximum(npmax,1))
                else:
                    npmax=int(np.maximum(npmax,len(p)))

        if nfloat[0]: x=np.ones(npmax)*x
        if nfloat[1]: y=np.ones(npmax)*y
        if nfloat[2]: z=np.ones(npmax)*z
        if nfloat[3]: vx=np.ones(npmax)*vx
        if nfloat[4]: vy=np.ones(npmax)*vy
        if nfloat[5]: vz=np.ones(npmax)*vz
        if nfloat[6]: m=np.ones(npmax)*m
        if nfloat[7]: id=np.ones(npmax)*id
        if nfloat[8]: m0=np.ones(npmax)*m0
        if nfloat[9]: npop=np.ones(npmax)*npop

        #Check for AMUSE units:
        isamuse=False

        if _hasamuse and (self.units=='amuse' or (self.units is None and self.ctype=='amuse')):
            if isinstance(x,VectorQuantity):
                stars=Particles(len(x))
                stars.x,stars.y,stars.z=x,y,z
                stars.vx,stars.vy,stars.vz=vx,vy,vz
                stars.mass=m
                stars.key=id

                if isinstance(self.tphys,ScalarQuantity):
                    self.tphys=self.tphys.value_in(u.Myr)

                self.units='pckms'                
                x,y,z,vx,vy,vz,m,id=_convert_amuse(stars,self)
                isamuse=True


        #Check for binaries
        if nb>0:
            if nb==1:
                arg1,arg2=0,1
            else:
                arg1=np.arange(0,2*nb-1,2)
                arg2=arg1+1

            if m is None:
                xcom,ycom,zcom,vxcom,vycom,vzcom,mcom=self.add_binary_stars(x[arg1],y[arg1],z[arg1],vx[arg1],vy[arg1],vz[arg1],x[arg2],y[arg2],z[arg2],vx[arg2],vy[arg2],vz[arg2],return_com=True)
            else:
                xcom,ycom,zcom,vxcom,vycom,vzcom,mcom=self.add_binary_stars(x[arg1],y[arg1],z[arg1],vx[arg1],vy[arg1],vz[arg1],x[arg2],y[arg2],z[arg2],vx[arg2],vy[arg2],vz[arg2],m[arg1],m[arg2],return_com=True)


            self.x = np.append(xcom,self.x)
            self.y = np.append(ycom,self.y)
            self.z = np.append(zcom,self.z)
            self.vx = np.append(vxcom,self.vx)
            self.vy = np.append(vycom,self.vy)
            self.vz = np.append(vzcom,self.vz) 
            self.m = np.append(mcom,self.m) 

            if id is None:
                if len(self.id)!=0:
                    idstart=np.amax(self.id)+1
                else:
                    idstart=0

                ids=idstart+np.arange(0,2*nb-1,2)
            else:
                ids=id[arg1]

            self.id=np.append(ids,self.id)

            if m0 is not None:
                self.m0=np.append(m0[arg1]+m0[arg2],self.m0)
            else:
                self.m0=np.append(np.zeros(nb),self.m0)

            if npop is None:
                npopb=np.ones(nb,int)
            else:
                npopb=npop[arg1]

            self.npop=np.append(npopb,self.npop).astype(int)

            args=2*nb
        else:
            args=0


        #Add single stars
        self.x = np.append(self.x, x[args:])
        self.y = np.append(self.y, y[args:])
        self.z = np.append(self.z, z[args:])
        self.vx = np.append(self.vx, vx[args:])
        self.vy = np.append(self.vy, vy[args:])
        self.vz = np.append(self.vz, vz[args:])

        if m is None:
            ms = np.ones(len(x[args:]),float)
        else:
            ms=m[args:]

        self.m = np.append(self.m, ms)

        if m0 is not None:
            self.m0=np.append(self.m0,m0[args:])
        else:
            self.m0=np.append(self.m0,np.zeros(len(x[args:])))

        if npop is not None:
            self.npop=np.append(self.npop,npop[args:]).astype(int)
        else:
            self.npop=np.append(self.npop,np.ones(len(x[args:]))).astype(int)


        if id is None:

            if len(self.id)!=0:
                idstart=np.amax(self.id)+1
                if self.nb>0:
                    idstart+=1
            else:
                idstart=0

            ids=idstart+np.arange(0, len(x[args:]), dtype=int)
        else:
            ids=id[args:]

        self.id = np.append(self.id, ids)

        self.id = self.id.astype(int)

        # Check lengths

        length_error=False
        nmax = np.amax(
            [
                len(self.id),
                len(self.m),
                len(self.x),
                len(self.y),
                len(self.z),
                len(self.vx),
                len(self.vy),
                len(self.vz),
                len(self.m0),
                len(self.npop),
            ]
        )
        nmin = np.amin(
            [
                len(self.id),
                len(self.m),
                len(self.x),
                len(self.y),
                len(self.z),
                len(self.vx),
                len(self.vy),
                len(self.vz),
                len(self.m0),
                len(self.npop),
            ]
        )

        if nmax != nmin:
            if len(self.id) == 1:
                self.id = np.linspace(0, nmax - 1, nmax, dtype=int)
            elif len(self.id)<nmax:
                length_error=True

            if len(self.m) == 1:
                self.m = np.ones(nmax) * self.m[0]
            elif len(self.m) <nmax:
                length_error=True

            if len(self.x) == 1:
                self.x = np.ones(nmax) * self.x[0]
            elif len(self.x) <nmax:
                length_error=True

            if len(self.y) == 1:
                self.y = np.ones(nmax) * self.y[0]
            elif len(self.y) <nmax:
                length_error=True

            if len(self.z) == 1:
                self.z = np.ones(nmax) * self.z[0]
            elif len(self.z) <nmax:
                length_error=True

            if len(self.vx) == 1:
                self.vx = np.ones(nmax) * self.vx[0]
            elif len(self.vx) <nmax:
                length_error=True

            if len(self.vy) == 1:
                self.vy = np.ones(nmax) * self.vy[0]
            elif len(self.vy) <nmax:
                length_error=True

            if len(self.vz) == 1:
                self.vz = np.ones(nmax) * self.vz[0]
            elif len(self.vz) <nmax:
                length_error=True

            if len(self.m0) == 1:   
                self.m0 = np.ones(nmax) * self.m0[0]
            elif len(self.m0) <nmax:
                length_error=True

            if len(self.npop) == 1:   
                self.npop = np.ones(nmax) * self.npop[0]
            elif len(self.npop) <nmax:
                length_error=True

        if length_error:
            print('ONE OR MORE INPUT ARRAY HAS INCORRECT LENGTH: ',nmin,nmax)
            print(len(self.id),len(self.m),len(self.x),len(self.y),len(self.z),len(self.x),len(self.y),len(self.z),len(self.m0),len(self.npop))

        if self.units == "radec" and self.origin == "sky":
            self.ra = np.append(self.ra, x)
            self.dec = np.append(self.dec, y)
            self.dist = np.append(self.dist, z)
            self.pmra = np.append(self.pmra, vx)
            self.pmdec = np.append(self.pmdec, vy)
            self.vlos = np.append(self.vlos, vz)


        if isamuse: 
            self.to_amuse()
            self.units_init='amuse'
        if analyze: self.analyze(sortstars=sortstars)


        self.ntot = len(self.x)
        self.ns=self.ntot-self.nb

    def add_binary_stars(
        self, xb1, yb1, zb1, vxb1, vyb1, vzb1, xb2, yb2, zb2, vxb2, vyb2, vzb2,mb1=None,mb2=None,id1=None,id2=None,m01=None,m02=None,npop1=None,npop2=None,set_com=True,return_com=False,
    ):
        """Individually add binary stars to StarCluster.

        Parameters
        ----------
        xb1,yb1,zb1: float
            stellar positions of primary. Input is assumed to be in cartesian coordinates unless self.units=='radec' 
            and self.origin=='sky', then positions are assumed to be ra,dec,dist (degrees, degrees, kpc)
        vxb1,vyb1,vzb1: float
            stellar velocities of primary. Input is assumed to be in cartesian coordinates unless self.units=='radec' 
            and self.origin=='sky', then positions are assumed to be pmra,pmdec,vlos (mas/yr, mas/yr, km/s)
        xb2,yb2,zb2: float
            stellar positions of secondary. Input is assumed to be in cartesian coordinates unless self.units=='radec' 
            and self.origin=='sky', then positions are assumed to be ra,dec,dist (degrees, degrees, kpc)
        vxb2,vyb2,vzb2: float
            stellar velocities of secondary. Input is assumed to be in cartesian coordinates unless self.units=='radec' 
            and self.origin=='sky', then positions are assumed to be pmra,pmdec,vlos (mas/yr, mas/yr, km/s)
        mb1. mb2 : float
            mass of primary and secondary (default: None)
        id1,id2: int 
            primary and secondary star ids (default: None)
        m01,m02: float
            initial primary and secondary stellar masses (default: None)
        npop1,npop2 : int
            population number of primary and secondary star, for use with multiple populations (default: None)
        set_com : bool
            set coordinates of binary's centre of mass in main arrays (default: False)
        return_com : bool
            reuturn coordinates of binary's centre of mass (default: False)

        Notes
        -----
        History:

        - 2022 - Written - Webb (UofT)

        """

        if isinstance(xb1,float):
            self.nb+=1
        else:
            self.nb+=len(xb1)

        self.xb1=np.append(self.xb1,xb1)
        self.yb1=np.append(self.yb1,yb1)
        self.zb1=np.append(self.zb1,zb1)
        self.vxb1=np.append(self.vxb1,vxb1)
        self.vyb1=np.append(self.vyb1,vyb1)
        self.vzb1=np.append(self.vzb1,vzb1)

        self.xb2=np.append(self.xb2,xb2)
        self.yb2=np.append(self.yb2,yb2)
        self.zb2=np.append(self.zb2,zb2)
        self.vxb2=np.append(self.vxb2,vxb2)
        self.vyb2=np.append(self.vyb2,vyb2)
        self.vzb2=np.append(self.vzb2,vzb2)

        if mb1 is not None:
            self.mb1=mb1
        else:
            self.mb1=np.ones(len(xb1))
        if mb2 is not None:
            self.mb2=mb2
        else:
            self.mb2=np.ones(len(xb2))

        if set_com and not return_com:
            mcom=self.mb1+self.mb2
            xcom=(xb1*self.mb1+xb2*self.mb2)/mcom
            ycom=(yb1*self.mb1+yb2*self.mb2)/mcom
            zcom=(zb1*self.mb1+zb2*self.mb2)/mcom
            vxcom=(vxb1*self.mb1+vxb2*self.mb2)/mcom
            vycom=(vyb1*self.mb1+vyb2*self.mb2)/mcom
            vzcom=(vzb1*self.mb1+vzb2*self.mb2)/mcom

            self.x = np.append(xcom,self.x)
            self.y = np.append(ycom,self.y)
            self.z = np.append(zcom,self.z)
            self.vx = np.append(vxcom,self.vx)
            self.vy = np.append(vycom,self.vy)
            self.vz = np.append(vzcom,self.vz) 
            self.m = np.append(mcom,self.m) 

            if id1 is None:
                if len(self.id)!=0:
                    idstart=np.amax(self.id)+1
                else:
                    idstart=0

                ids=idstart+np.linspace(0,len(xcom)-1,len(xcom))*2
            else:
                ids=id1
                self.id1=id1

            self.id=np.append(ids,self.id)

            if id2 is not None:
                self.id2=id2
            else:
                self.id2=self.id1+1

            if m01 is not None and m02 is not None:
                self.m0=np.append(m01+m02,self.m0)
            else:
                self.m0=np.append(np.zeros(len(xcom)),self.m0)

            if m01 is not None:
                self.m01=np.append(m01,self.m01)
            if m02 is not None:
                self.m02=np.append(m02,self.m02)

            if npop1 is not None and npop2 is not None:
                npopb=np.maximum(npop1,npop2)
                self.npop=np.append(npopb,self.npop).astype(int)

            if npop1 is None:
                self.npop1=np.append(np.ones(len(xcom)),self.npop1)
            elif isinstance(npop1,float):
                self.npop1=np.append(np.ones(len(xcom))*npop1,self.npop1)
            else:
                self.npop1=np.append(npop1,self.npop1)

            if npop2 is None:
                self.npop2=np.append(np.ones(len(xcom)),self.npop2)
            elif isinstance(npop1,float):
                self.npop2=np.append(np.ones(len(xcom))*npop2,self.npop2)
            else:
                self.npop2=np.append(npop2,self.npop2)

            self.ntot = len(self.x)
            self.ns=self.ntot-self.nb

        elif return_com:
            mcom=self.mb1+self.mb2
            xcom=(xb1*self.mb1+xb2*self.mb2)/mcom
            ycom=(yb1*self.mb1+yb2*self.mb2)/mcom
            zcom=(zb1*self.mb1+zb2*self.mb2)/mcom
            vxcom=(vxb1*self.mb1+vxb2*self.mb2)/mcom
            vycom=(vyb1*self.mb1+vyb2*self.mb2)/mcom
            vzcom=(vzb1*self.mb1+vzb2*self.mb2)/mcom

            return xcom,ycom,zcom,vxcom,vycom,vzcom,mcom

    def add_orbit(
        self,
        xgc,
        ygc,
        zgc,
        vxgc,
        vygc,
        vzgc,
        ounits=None,
        initialize=False,
        from_centre=False,
        tphys=None,
    ):
        """ Add orbit properties to StarCluster

        Parameters
        ----------
        xgc,ygc,zgc: float
            cluster's galactocentric position
        vxgc,vygc,vzgc: float
            cluster's galactocentric velocity
        ounits: str
            units of position and velocity. Options include 'pckms','pcmyr'
            'kpckms','kpcgyr','radec','nbody',and 'galpy'. Values will be converted 
            to match self.units
        initialize: bool
            Initialize a galpy orbit for self.orbit (default: False)
        from_centre : bool
            genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
        tphys : float
            physical time as per the orbit file

        Returns
        ----------

        None

        History:
        ----------

        2018 - Written - Webb (UofT)

        """

        ro,vo,zo,solarmotion=self._ro,self._vo,self._zo,self._solarmotion

        if tphys is not None:
            otime=True
        else:
            tphys=0.
            otime=False

        if ounits != None and ounits != self.units:
            # First convert to kpckms
            if ounits != "kpckms":
                if ounits == "nbody":
                    xgc *= self.rbar / 1000.0
                    ygc *= self.rbar / 1000.0
                    zgc *= self.rbar / 1000.0
                    vxgc *= self.vbar
                    vygc *= self.vbar
                    vzgc *= self.vbar
                    tphys *= self.tbar
                elif ounits == "galpy":
                    xgc *= ro
                    ygc *= ro
                    zgc *= ro
                    vxgc *= vo
                    vygc *= vo
                    vzgc *= vo
                    tphys *= conversion.time_in_Gyr(ro=ro,vo=vo)
                elif ounits == "pckms":
                    xgc /= 1000.0
                    ygc /= 1000.0
                    zgc /= 1000.0
                    tphys /= 1000.0
                elif ounits == "pcmyr":
                    xgc /= 1000.0
                    ygc /= 1000.0
                    zgc /= 1000.0
                    vxgc/=1.022712165045695
                    vygc/=1.022712165045695
                    vzgc/=1.02271216504569
                    tphys /= 1000.0
                elif ounits == 'kpcgyr':
                    vxgc/=1.022712165045695
                    vygc/=1.022712165045695
                    vzgc/=1.022712165045695


                elif ounits == 'radec':
                    o=Orbit([xgc,ygc,zgc,vxgc,vygc,vzgc],radec=True,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
                    xgc=o.x()
                    ygc=o.y()
                    zgc=o.z()
                    vxgc=o.vx()
                    vygc=o.vy()
                    vzgc=o.vz()

                ounits = "kpckms"

            if self.units == "pckms":
                xgc *= 1000.0
                ygc *= 1000.0
                zgc *= 1000.0
                tphys *= 1000.0
            elif self.units == "nbody":
                xgc *= 1000.0 / self.rbar
                ygc *= 1000.0 / self.rbar
                zgc *= 1000.0 / self.rbar
                vxgc /= self.vbar
                vygc /= self.vbar
                vzgc /= self.vbar
                tphys *= 1000.0/self.tbar
            elif self.units == "galpy":
                xgc /= ro
                ygc /= ro
                zgc /= ro
                vxgc /= vo
                vygc /= vo
                vzgc /= vo
                tphys /= conversion.time_in_Gyr(ro=ro,vo=vo)


        self.xgc = xgc
        self.ygc = ygc
        self.zgc = zgc
        self.rgc = np.sqrt(xgc ** 2.0 + ygc ** 2.0 + zgc ** 2.0)
        self.vxgc = vxgc
        self.vygc = vygc
        self.vzgc = vzgc

        if otime:
            self.tphys=tphys

        if self.units == "radec":
            self.ra_gc = xgc
            self.dec_gc = ygc
            self.dist_gc = zgc
            self.pmra_gc = vxgc
            self.pmdec_gc = vygc
            self.vlos_gc = vzgc


            #Add on orbital parameters for missing data
            if self.origin == 'sky':
                if np.array_equal(self.ra,np.zeros(len(self.ra))):
                    self.ra += self.ra_gc
                if np.array_equal(self.dec,np.zeros(len(self.dec))):
                    self.dec += self.dec_gc
                if np.array_equal(self.dist, np.zeros(len(self.dist))):
                    self.dist += self.dist_gc
                if np.array_equal(self.pmra,np.zeros(len(self.pmra))):
                    self.pmra += self.pmra_gc
                if np.array_equal(self.pmdec,np.zeros(len(self.pmdec))):
                    self.pmdec += self.pmdec_gc
                if np.array_equal(self.vlos,np.zeros(len(self.vlos))):
                    self.vlos += self.vlos_gc

        if initialize:
            if self.origin=='galaxy' or self.origin=='sky':
                self.initialize_orbit(from_centre=from_centre,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
            else:
                origin0=self.origin
                self.to_galaxy()
                self.initialize_orbit(from_centre=from_centre,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
                self.to_origin(origin0)


    def add_nbody6(
        self,
        nc=0,
        rc=0.0,
        rbar=1.0,
        rtide=0.0,
        xc=0.0,
        yc=0.0,
        zc=0.0,
        zmbar=1.0,
        vbar=1.0,
        tbar=1.0,
        rscale=1.0,
        ns=0.0,
        nb=0.0,
        n_p=0.0,
    ):
        """ Add additional information to StarCluster

        - parameters are common output variables in NBODY6
        - values are never adjusted during unit or coordinate changes

        Parameters
        ----------
        nc : int
            number of stars in core (default:0)
        rc : int
            core radius (default:0)
        rbar : float
            scaling factor between NBODY units and pc (default:1.)
        rtide :
            rtide set in the simulation (default:0)
        xc,yc,zc : float
            position of cluster centre (default:0)
        zmbar : float
            scaling factor between NBODY units and Msun (default:1.)
        vbar : float
            scaling factor between NBODY units and km/s (default:1.)
        tbar : float
            scaling factor between NBODY units and Myr (default:1.)
        rscale : float
            the scale radius of data (default:1.)
        ns : int
            number of single stars (default:0)
        nb : int
            number of binary stars (default:0)
        np : int
            number of particles (default:0)

        Returns
        ----------
        None

        History:
        ----------
        2018 - Written - Webb (UofT)
        """

        if isinstance(nc,list):
            nc,rc,rbar,rtide,xc,yc,zc,zmbar,vbar,tbar,rscale,ns,nb,n_p=nc

        # Number of stars in the core
        self.nc = nc
        # Core radius
        self.rc = rc
        # Distane scaling parameter
        self.rbar = rbar
        # Tidal limit from NBODY6 (not neccesarily a true tidal radius)
        self.rtide = rtide
        # Center of mass of cluster (x,yz)
        self.xc = xc
        self.yc = yc
        self.zc = zc
        self.xcn = xc
        self.ycn = yc
        self.zcn = zc
        # Mass scaling parameter
        self.zmbar = zmbar
        # Velocity scaling parameter
        self.vbar = vbar
        # Time scaling parameter
        self.tbar=tbar

        # Scale radius of cluster
        self.rscale = rscale
        # Number of single stars
        self.ns = ns
        # Number of binary stars
        self.nb = nb
        # Number of particles (from NBODY6 when tidal tail is being integrated)
        self.n_p = n_p

        au_to_cm = 1.49597870700e13
        pc_to_cm = 1296000.0/(2.0*np.pi)*au_to_cm

        nbody_to_years = (self.rbar*1296000.0/(2.0*np.pi))**1.5/np.sqrt(self.zmbar)
        self.tbar_days = 365.25*nbody_to_years
        rsun_to_cm = 6.957e10
        self.rbar_su = pc_to_cm/rsun_to_cm*self.rbar
        self.rbar_au= pc_to_cm/au_to_cm*self.rbar

    def add_sse(self, kw, logl, logr, ep = None, ospin = None, arg = None):
        """Add stellar evolution information to stars
        
        - parameters are common output variables in NBODY6
        - values are never adjusted during unit or coordinate changes

        Parameters
        ----------
        kw : int
            stellar evolution type (based on NBODY6) 
        logl : float
            log of luminosity
        logr : float
            log of stellar radius
        ep : float
            epoch
        ospin : float
            ospin
        arg : int
            array address arguments if SSE parameters do not necessarily match position and velocity arrays

        Returns
        ----------
        None

        History:
        ----------
        2018 - Written - Webb (UofT)

        """


        if arg is not None:

            if len(self.kw) == 0:

                nstart=0

                self.kw=np.zeros(len(kw))
                self.logl=np.zeros(len(kw))
                self.logr=np.zeros(len(kw))
                self.lum=np.zeros(len(kw))
                self.ep=np.zeros(len(kw))
                self.ospin=np.zeros(len(kw))
            else:

                nstart=len(self.kw)

                self.kw=np.append(self.kw,np.zeros(len(kw)))
                self.logl=np.append(self.logl,np.zeros(len(kw)))
                self.logr=np.append(self.logr,np.zeros(len(kw)))
                self.lum=np.append(self.lum,np.zeros(len(kw)))
                self.ep=np.append( self.ep,np.zeros(len(kw)))
                self.ospin=np.append(self.ospin,np.zeros(len(kw)))

            self.kw[nstart+arg.astype(int)] = np.asarray(kw)
            self.logl[nstart+arg.astype(int)] = np.asarray(logl)
            self.logr[nstart+arg.astype(int)] = np.asarray(logr)
            self.ep[nstart+arg.astype(int)]= np.array(ep)
            self.ospin[nstart+arg.astype(int)] = np.array(ospin)

        else:
            self.kw = np.append(self.kw,kw)
            self.logl = np.append(self.logl,logl)
            self.logr = np.append(self.logr,logr)

        self.lum = 10.0 ** self.logl
        self.ltot = np.sum(self.lum)

        if ep is not None: 
            self.ep= np.append(self.ep,np.array(ep))
        else:
            self.ep= np.append(self.ep,np.zeros(len(kw)))

        if ospin is not None: 
            self.ospin = np.append(self.ospin,np.array(ospin))
        else:
            self.ospin = np.append(self.ospin,np.zeros(len(kw)))

    def add_bse(
        self,
        id1,
        id2,
        kw1,
        kw2,
        kcm,
        ecc,
        pb,
        semi,
        m1,
        m2,
        logl1,
        logl2,
        logr1,
        logr2,
        ep1=None,
        ep2=None,
        ospin1=None,
        ospin2=None,
    ):
        """Add binary star evolution information to stars

        - parameters are common output variables in NBODY6
        - values are never adjusted during unit or coordinate changes

        Parameters
        ----------

        id1/id2 : int
            id of star1 and star2
        kw1/kw2 : int
            stellar evolution tags (for using with NBODY6) 
        kcm : int
            stellar evolution tag for binary star
        ecc : float 
            eccentricity of binary orbit
        pb: float
            period of binary orbit
        semi : float 
            semi major axis of binary orbit
        m1/m2 : float
            masses of binary stars
        logl1/logl2 : float
            luminosities of binary stars
        logr1/logr2 : float
            radii of binary stars
        ep1/ep2 : float
            epochs of binary stars
        ospin1/ospin2 : float
            opsin of binary stars

        Returns
        ----------

        None

        History:
        ----------

        2018 - Written - Webb (UofT)
        """

        self.id1 = np.append(self.id1,np.array(id1))
        self.id2 = np.append(self.id2,np.array(id2))
        self.kw1 = np.append(self.kw1,np.array(kw1))
        self.kw2 = np.append(self.kw2,np.array(kw2))
        self.kcm = np.append(self.kcm,np.array(kcm))
        self.ecc = np.append(self.ecc,np.array(ecc))
        self.pb = np.append(self.pb,np.array(pb))
        self.semi = np.append(self.semi,np.array(semi))
        self.m1 = np.append(self.m1,np.array(m1))
        self.m2 = np.append(self.m2,np.array(m2))
        self.logl1 = np.append(self.logl1,np.array(logl1))
        self.logl2 = np.append(self.logl2,np.array(logl2))
        self.logr1 = np.append(self.logr1,np.array(logr1))
        self.logr2 = np.append(self.logr2,np.array(logr2))

        if ep1 is not None: self.ep1 = np.append(self.ep1,np.array(ep1))
        if ep2 is not None: self.ep2 = np.append(self.ep2,np.array(ep2))
        if ospin1 is not None: self.ospin1 = np.append(self.ospin1,np.array(ospin1))
        if ospin2 is not None: self.ospin2 = np.append(self.ospin2,np.array(ospin2))

        self.eb=0.5*self.m1*self.m2/self.semi

        self.nb=len(self.id1)

    def add_energies(self, kin, pot, etot=None):
        """ Add energy information for stars 

        - total energy and Q for the cluster are also calculated
        - values are never adjusted during unit or coordinate changes

        Parameters
        ----------
        kin : float
            kinetic energy 
        pot : float
            potentail energy
        etot : float
            total energy - calculated as kin+pot if not given

        Returns
        ----------
        None

        History
        ----------

        2018 - Written - Webb (UofT)

        """

        if _hasamuse and (self.units=='amuse' or (self.units is None and self.ctype=='amuse')):
            if isinstance(kin,VectorQuantity) or isinstance(kin,ScalarQuantity):
                self.kin=kin
                self.pot=pot
            else:
                self.kin = np.array(kin)
                self.pot = np.array(pot)               
        else:
            self.kin = np.array(kin)
            self.pot = np.array(pot)

        if etot is None:
            self.etot=self.kin+self.pot 
        else:
            self.etot = np.array(etot)

        self.ektot = np.sum(np.nan_to_num(self.kin))
        self.ptot = np.sum(np.nan_to_num(self.pot)) / 2.0

        try:
            self.qvir = self.ektot / self.ptot
        except:
            self.qvir = 0.0

    def add_action(self, JR, Jphi, Jz, OR=None, Ophi=None, Oz=None, TR=None, Tphi=None, Tz=None):
        """ Add action values to the cluster instance

        Parameters
        ----------
        JR,Jphi,Jz : float
            orbit actions
        OR,Ophi,Oz : float
            orbit frequencies
        TR,Tphi,Tz : float
            orbit periods

        Returns
        -------
        None

        History
        --------
        2019 - Written - Webb (UofT)
        """
        self.JR, self.Jphi, self.Jz = JR, Jphi, Jz
        self.OR, self.Ophi, self.Oz = OR, Ophi, Oz
        self.TR, self.Tphi, self.Tz = TR, Tphi, Tz

    def add_actions(self, JR, Jphi, Jz, OR=None, Ophi=None, Oz=None, TR=None, Tphi=None, Tz=None):
        """ Add action values to the cluster instance

        Parameters
        ----------
        JR,Jphi,Jz : float
            orbit actions
        OR,Ophi,Oz : float
            orbit frequencies
        TR,Tphi,Tz : float
            orbit periods

        Returns
        -------
        None

        History
        --------
        2019 - Written - Webb (UofT)
        """
        self.JRs, self.Jphis, self.Jzs = JR, Jphi, Jz
        self.ORs, self.Ophis, self.Ozs = OR, Ophi, Oz
        self.TRs, self.Tphis, self.Tzs = TR, Tphi, Tz

    def analyze(self, sortstars = True, projected = True):
        """ Calculate properties related to mass, radius, and velocity

        Parameters
        ----------
        sortstars : bool
            sort star by radius after coordinate change (default: True)
        projected : bool
            sort projected radii as well, but do not change self.projected (default: True) 

        Returns
        ----------
        None

        History
        ----------

        2019 - Written - Webb (UofT)
        """

        self.analyze_units=self.units
        self.analyze_origin=self.origin

        self.r = np.sqrt(self.x ** 2.0 + self.y ** 2.0 + self.z ** 2.0)
        self.rpro = np.sqrt(self.x ** 2.0 + self.y ** 2.0)
        self.v = np.sqrt(self.vx ** 2.0 + self.vy ** 2.0 + self.vz ** 2.0)
        self.vpro = np.sqrt(self.vx ** 2.0 + self.vy ** 2.0)

        self.mtot = np.sum(self.m)
        self.mmean = np.mean(self.m)

        self.rmean = np.mean(self.r)
        self.rmax = np.max(self.r)
        self.rmeanpro = np.mean(self.rpro)
        self.rmaxpro = np.max(self.rpro)

        if sortstars: self.sortstars(projected=projected)

        # Find half-mass radius

        if self.rorder is not None:
            msum = np.cumsum(self.m[self.rorder])
            indx = msum >= 0.5 * self.mtot
            self.rm = self.r[self.rorder[indx][0]]  
            indx = msum >= 0.1 * self.mtot
            self.r10 = self.r[self.rorder[indx][0]]

        if self.rproorder is not None:
            msum = np.cumsum(self.m[self.rproorder])
            indx = msum >= 0.5 * self.mtot
            self.rmpro = self.rpro[self.rproorder[indx][0]]
            indx = msum >= 0.1 * self.mtot
            self.r10pro = self.rpro[self.rproorder[indx][0]]


        if len(self.lum) > 0 and self.rorder is not None:
            lsum = np.cumsum(self.lum[self.rorder])
            indx = lsum >= 0.5 * self.ltot
            self.rh = self.r[self.rorder[indx][0]]
            indx = lsum >= 0.1 * self.ltot
            self.rh10 = self.r[self.rorder[indx][0]]

            if self.rproorder is not None:
                lsum = np.cumsum(self.lum[self.rproorder])
                indx = lsum >= 0.5 * self.ltot
                self.rhpro = self.rpro[self.rproorder[indx][0]]
                indx = lsum >= 0.1 * self.ltot
                self.rh10pro = self.rpro[self.rproorder[indx][0]]
            else:
                self.rhpro = 0.0
                self.rh10pro = 0.0

    def analyse(self, sortstars = True, projected=True):
        """Call analyze with alternative spelling

        Parameters
        ----------
        sortstars : bool
            sort star by radius after coordinate change (default: True)
        projected : bool
            sort projected radii as well, but do not change self.projected (default: True) 

        Returns
        ----------

        None

        History
        ----------

        2020 - Written - Webb (UofT)
        """
        analyze(self, sortstars = sortstars, projected=projected)

    def key_params(self, do_order=True, projected=True):
        """Call analyze with key_params for backwards compatibility

        Parameters
        ----------
        do_order : bool
            sort star by radius after coordinate change (default: True)
        projected : bool
            sort projected radii as well, but do not change self.projected (default: True) 

        Returns
        ----------

        None

        History
        ----------
 
        2020 - Written - Webb (UofT)
        """

        analyze(self,sortstars=do_order, projected=projected)


    def sortstars(self, projected=True):
        """sort stars based on radius 

        Parameters
        ----------
        projected : bool
            sort projected radii as well, but do not change self.projected (default: True) 

        Returns
        ----------

        None

        History
        ----------

        2018 - Written - Webb (UofT)

        """
        if self. rorder is None or self.rorder_origin!=self.origin or len(self.r)!=len(self.rorder):

            self.rorder = np.argsort(self.r)
            self.rorder_origin=self.origin

        if (self. rproorder is None or self.rorder_origin!=self.origin or len(self.rpro)!=len(self.rproorder)) and (projected or self.projected):
            self.rproorder = np.argsort(self.rpro)

    def _analysis_check(self,sortstars=True,projected=False):
        if self.units!=self.analyze_units or self.analyze_units!=self.units:
            self.analyze(sortstars=sortstars,projected=projected)

    def subset(self,**kwargs):
        """Generate a boolean array that corresponds to subset of star cluster members that meet a certain criteria

        Parameters
        ----------
        rmin/rmax : float
            minimum and maximum stellar radii
        mmin/mmax : float
            minimum and maximum stellar mass
        vmin/vmax : float
            minimum and maximum stellar velocity
        emin/emax : float
            minimum and maximum stellar energy
        kwmin/kwmax : int
            minimum and maximum stellar type (kw)
        npop : int
            population number
        indx : bool
            user defined boolean array from which to extract the subset
        projected : bool 
            use projected values and constraints (default:False)

        Returns
        -------
        indx : bool
            boolean array that is True for stars that meet the criteria

        History
        -------
        2022 - Written - Webb (UofT)

        """ 
        self.indx=_get_subset(self,**kwargs)
        return self.indx

    # Directly call from operations.py (see operations.py files for documenation):

    def to_pckms(self,analyze=True):
        """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and km/s

        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History
        -------
        2018 - Written - Webb (UofT)

        """
        to_pckms(self,analyze=analyze)

    def to_kpckms(self,analyze=True):
        """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and km/s

        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_kpckms(self,analyze=analyze)

    def to_pcmyr(self,analyze=True):
        """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and pc/Myr
           
        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2022 - Written - Webb (UofT)

        """
        to_pcmyr(self,analyze=analyze)

    def to_kpcgyr(self,analyze=True):
        """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and kpc/Gyr
           
        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2022 - Written - Webb (UofT)

        """
        to_kpcgyr(self,analyze=analyze)


    def to_nbody(self,analyze=True):
        """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Nbody units
       
        - requires that cluster.zmbar, cluster.rbar, cluster.vbar are set (defaults are 1)

        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_nbody(self,analyze=analyze)

    def to_radec(self, sortstars=True,centre_method=None,analyze=True):
        """Convert to on-sky position, proper motion, and radial velocity of cluster
        
        Parameters
        ----------
        cluster : class
            StarCluster
        sortstars : bool
            sort star by radius after coordinate change (default: False)
        centre_method : str
            method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)
        """
        to_radec(self, sortstars=sortstars,centre_method=centre_method,analyze=analyze)

    def to_galpy(self,analyze=True):
        """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to galpy units
        
        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_galpy(self,analyze=analyze)

    def to_WDunits(self,analyze=True):
        """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Walter Dehnen Units

        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2022 - Written - Webb (UofT)

        """
        to_WDunits(self,analyze=analyze)

    def to_amuse(self,analyze=True):
        """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to AMUSE Vector and Scalar Quantities

        Parameters
        ----------
        cluster : class
            StarCluster
        analyze : bool
            run analysis function (default: True)

        Returns
        -------
        None

        History:
        -------
        2022 - Written - Webb (UofT)

        """

        to_amuse(self,analyze=analyze)

    def to_units(self, units):
        """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to user defined units

        Parameters
        ----------
        cluster : class
            StarCluster
        units : str
            units to be converted to (nbody,galpy,pckms,kpckms,radec)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_units(self, units=units)

    def to_sudays(self):
        """ Convert binary star semi-major axis and periods to solar radii and days
            Note: Masses will be converted to solar masses

        Parameters
        ----------
        cluster : class
            StarCluster

        Returns
        -------
        None

        History
        -------
        2022 - Written - Webb (UofT)

        """
        to_sudays(self)

    def to_audays(self):
        """ Convert binary star semi-major axis and periods to solar radii and days

        Parameters
        ----------
        cluster : class
            StarCluster

        Returns
        -------
        None

        History
        -------
        2022 - Written - Webb (UofT)

        """
        to_audays(self)

    def to_centre(self, sortstars=True, centre_method=None):
        """Shift coordinates such that origin is the centre of mass

        Parameters
        ----------
        cluster : class
            StarCluster
        sortstars : bool
            sort star by radius after coordinate change (default: False)
        centre_method : str
            method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_centre(self, sortstars=sortstars, centre_method=centre_method)

    def to_center(self, sortstars=True, centre_method=None):
        """Shift coordinates such that origin is the centre of mass

        - this is the same to to_centre, but allows for alternate spelling

        Parameters
        ----------
        cluster : class
            StarCluster
        sortstars : bool
            sort star by radius after coordinate change (default: False)
        centre_method : str
            method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_centre(self, sortstars=sortstars, centre_method=centre_method)

    def to_cluster(self, sortstars=True, centre_method=None):
        """Shift coordinates to clustercentric reference frame

        - When units='radec' and origin='sky', the cluster will be shifted to clustercentric coordinates using either:
        --centre_method=None: angular distance between each star's RA/DEC and the RA/DEC of the cluster's centre with proper motions directly subtracted off
        --centre_method='orthographic' , positions and velocities changed to orthnormal coordinates (Helmi et al. 2018)
        -- centre_method='VandeVen' , positions and velocities changed to clustercentric coordinates using method outlined by Van de Ven et al. 2005

        - Note the the line of site positions and velocities will just have the galactocentric coordinates of the cluster
        subtracted off. Be sure to set projected=True when making any calculations to use only x and y coordinates

        Parameters
        ----------
        cluster : class
            StarCluster
        sortstars : bool
            sort star by radius after coordinate change (default: False)
        centre_method : str
            method for shifting coordinates to clustercentric coordinates. (default: None)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_cluster(self, sortstars=sortstars, centre_method=centre_method)

    def to_galaxy(self, sortstars=True):
        """Shift coordinates to galactocentric reference frame

        Parameters
        ----------
        cluster : class
            StarCluster
        sortstars : bool
            sort star by radius after coordinate change (default: False)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_galaxy(self, sortstars=sortstars)

    def to_sky(self, sortstars=True, centre_method=None):
        """Calculate on-sky position, proper motion, and radial velocity of cluster
            
        - Also changes units to radec

        Parameters
        ----------
        cluster : class
            StarCluster
        sortstars : bool
            sort star by radius after coordinate change (default: False)
        centre_method : str
            method for shifting coordinates to clustercentric coordinates. (default: None)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)
        """
        to_sky(self, sortstars=sortstars,centre_method=centre_method)

    def to_origin(self, origin, sortstars=True, centre_method=None):
        """Shift cluster to origin as defined by keyword

        Parameters
        ----------
        cluster : class
            StarCluster
        origin : str
            origin of coordinate system to be shifted to (centre, cluster, galaxy, sky)
        sortstars : bool
            sort star by radius after coordinate change (default: False)
        centre_method : str
            method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        to_origin(self, origin, sortstars=sortstars, centre_method=centre_method)

    def save_cluster(self):
        """Save cluster's units and origin

        Parameters
        ----------
        cluster : class
            StarCluster

        Returns
        -------
        cluster.units, cluster.origin : str
            units and origin of StarCluster

        History:
        -------
        2018 - Written - Webb (UofT)
        """
        self.units0,self.origin0, self.rorder0, self.rorder_origin0=save_cluster(self)

    def return_cluster(self, units0=None, origin0=None, rorder0=None, rorder_origin0=None, sortstars=False):
        """ return cluster to a specific combination of units and origin

        Parameters
        ----------

        cluster : class
            StarCluster
        units0 : str
            units that StarCluster will be changed to
        origin0 : str
            origin that StarCluster will be changed to
        sortstars : bool
            sort star by radius after coordinate change (default: False)

        Returns
        -------
        None

        History:
        -------
        2018 - Written - Webb (UofT)
        """
        if units0 is None: units0=self.units0
        if origin0 is None: origin0=self.origin0
        if rorder0 is None: rorder0=self.rorder0
        if rorder_origin0 is None: rorder_origin0=self.rorder_origin0

        return_cluster(self, units0, origin0, rorder0, rorder_origin0, sortstars=sortstars)

    def reset_nbody_scale(self, mass=True, radii=True, rvirial=True, projected=None, **kwargs):
        """ Assign new conversions for real mass, size, and velocity to Nbody units

        - kwargs are passed to the virial_radius function. See the virial_radius documenation in functions.py

        Parameters
        ----------

        cluster : class
            StarCluster instance
        mass : bool
            find new mass conversion (default: True)
        radii : bool
            find new radius conversion (default: True)
        rvirial : bool (default: True)
            use virial radius to set conversion rate for radii as opposed to the approximation that rbar=4/3 rm
        projected : bool
            use projected values to calculate radius conversion (default: False)

        Returns
        -------

        zmbar : float
            mass conversion
        rbar : float
            radius conversion
        vbar : float
            velocity conversion
        tbar : float
            time conversion

        History:

        2018 - Written - Webb (UofT)
        """
        if projected==None:
            projected=self.projected

        self.zmbar,self.rbar,self.vbar,self.tbar=reset_nbody_scale(self, mass=mass, radii=radii, rvirial=rvirial,projected=projected,**kwargs)

    def add_rotation(self, qrot):
        """Add a degree of rotation to an already generated StarCluster

        Parameters
        ----------
        cluster : class
            StarCluster
        qrot : float
            fraction of stars with v_phi < 0 that are switched to vphi > 0

        Returns
        -------
        x,y,z,vx,vy,vz : float
            stellar positions and velocities (now with rotation)

        History
        -------
        2019 - Written - Webb (UofT)
        """
        self.x,self.y,self.z,self.vx,self.vy,self.vz=add_rotation(self, qrot)
        self.qrot=qrot

    def virialize(self, qvir=0.5,specific=True, full=True, projected=None, softening=0.0):
        """ Adjust stellar velocities so cluster is in virial equilibrium

        Parameters
        ----------
        cluster : class
            StarCluster
        qvir : float
            value you wish to virial parameter to be (default: 0.5)
        specific : bool
            find specific energies (default: True)
        full: bool
            do full array of stars at once with numba (default: full_default)
        projected : bool
            use projected values when calculating energies (default: False)
        softening : float
          Plummer softening length in cluster.units (default: 0.0)
        Returns
        -------
        qv : float
            scaling factor used to adjust velocities

        History
        -------
        2019 - Written - Webb (UofT)
        """
        if projected==None:
            projected=self.projected

        self.qv=virialize(self, qvir=qvir, specific=True, full=True, projected=projected,softening=softening)

        self.save_cluster()
        units0,origin0, rorder0, rorder_origin0 = self.units0,self.origin0, self.rorder0, self.rorder_origin0

        if self.origin0 != 'cluster' and self.origin0 != 'centre':
            self.to_centre(sortstars=False)

        self.vx *= self.qv
        self.vy *= self.qv
        self.vz *= self.qv

        self.return_cluster(units0,origin0, rorder0, rorder_origin0)

    # Directly call from functions.py (see functions.py files for documenation):

    def find_centre(
        self,
        xstart=0.0,
        ystart=0.0,
        zstart=0.0,
        vxstart=0.0,
        vystart=0.0,
        vzstart=0.0,
        indx=None,
        nsigma=1.0,
        nsphere=100,
        density=True,
        rmin=0.1,
        rmax=None,
        nmax=100,
        method='harfst',
        nneighbour=6,
        reset_centre=False
    ):
        """Find the cluster's centre

        - The default assumes the cluster's centre is the centre of density, calculated via the find_centre_of_density function.
        - For density=False, the routine first works to identify a sphere of nsphere stars around the centre in which to perform a centre of mass calculation (similar to NBODY6). Stars beyond nsigma standard deviations are removed from the calculation until only nsphere stars remain. This step prevents long tidal tails from affecting the calculation

        Parameters
        ----------
        cluster : class
            StarCluster
        xstart,ystart,zstart : float
            starting position for centre (default: 0,0,0)
        vxstart,vystart,vzstart : float
            starting velocity for centre (default: 0,0,0)
        indx: bool
            subset of stars to perform centre of density calculation on (default: None)
        nsphere : int
            number of stars in centre sphere (default:100)
        rmin : float
            minimum radius of sphere around which to estimate density centre (default: 0.1 cluster.units)
        rmax : float
            maxmimum radius of sphere around which to estimate density centre (default: None cluster.units, uses maximum r)
        nmax : float
            maximum number of iterations (default:100)

        method : str
            method for finding the centre of density ('harfst' (default), 'casertano')

        if method=='casertano'
            nneighbour : int
                number of neighbours for calculation local densities


        reset_centre : bool
            forcibly reset cluster's centre of mass (default: False)

        Returns
        -------
        xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

        History
        -------
        2019 - Written - Webb (UofT)
        """

        xc,yc,zc,vxc,vyc,vzc=find_centre(self,xstart=xstart,
            ystart=ystart,zstart=zstart,vxstart=vxstart,vystart=vystart,vzstart=vzstart,indx=indx,
            nsigma=nsigma,nsphere=nsphere,density=density,
            rmin=rmin,rmax=rmax,nmax=nmax,method=method,nneighbour=nneighbour)

        if self.origin=='centre':

            warning=False
            if xc!=0.0 or yc!=0.0 or zc!=0.0:
                warning=True
            if vxc!=0.0 or vyc!=0.0 or vzc!=0.0:
                warning=True

            if warning:
                print('Centre is not at origin')

            return xc,yc,zc,vxc,vyc,vzc

        elif self.origin=='cluster':
            self.xc, self.yc, self.zc = xc,yc,zc
            self.vxc, self.vyc, self.vzc = vxc,vyc,vzc
            
            return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc


        elif self.origin == "galaxy":

            if self.units!='amuse':
                if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0.,0.,0.,0.,0.,0.) or reset_centre:
                    self.xgc, self.ygc, self.zgc = xc,yc,zc
                    self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                    self.xc, self.yc, self.zc = 0.0, 0.0, 0.0
                    self.vxc, self.vyc, self.vzc = 0.0, 0.0, 0.0

                    return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc
                else:
                    self.xc,self.yc,self.zc=xc-self.xgc,yc-self.ygc,zc-self.zgc
                    self.vxc,self.vyc,self.vzc=vxc-self.vxgc,vyc-self.vygc,vzc-self.vzgc

            elif self.units=='amuse':
                if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0. | u.pc ,0. | u.pc,0. | u.pc,0. | u.kms,0. | u.kms ,0. | u.kms) or reset_centre:
                    self.xgc, self.ygc, self.zgc = xc,yc,zc
                    self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                    self.xc, self.yc, self.zc = 0.0 | u.pc, 0.0 | u.pc, 0.0 | u.pc
                    self.vxc, self.vyc, self.vzc = 0.0 | u.kms, 0.0 | u.kms , 0.0 | u.kms

                    return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc
     
                else:
                    self.xc,self.yc,self.zc=xc-self.xgc,yc-self.ygc,zc-self.zgc
                    self.vxc,self.vyc,self.vzc=vxc-self.vxgc,vyc-self.vygc,vzc-self.vzgc


                return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc

        elif self.origin=='sky':
            if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0.,0.,0.,0.,0.,0.) or reset_centre:
                self.xgc, self.ygc, self.zgc = xc,yc,zc
                self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                self.xc, self.yc, self.zc = 0.0, 0.0, 0.0
                self.vxc, self.vyc, self.vzc = 0.0, 0.0, 0.0

                return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc

            else:
                print('No Cluster Variables Set')
                return xc,yc,zc,vxc,vyc,vzc
        elif self.origin is None:
            self.xc, self.yc, self.zc = xc,yc,zc
            self.vxc, self.vyc, self.vzc = vxc,vyc,vzc

            return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc

        else:
            print('No Cluster Variables Set')
            return xc,yc,zc,vxc,vyc,vzc

    def find_centre_of_density(
        self,
        xstart=0.0,
        ystart=0.0,
        zstart=0.0,
        vxstart=0.0,
        vystart=0.0,
        vzstart=0.0,
        indx=None,
        nsphere=100,
        rmin=0.1,
        rmax=None,
        nmax=100,
        method='harfst',
        nneighbour=6,
        reset_centre=False,
    ):
        """Find cluster's centre of density

        - Find cluster's centre of density:
            if method=='harfst' use (Harfst, S., Gualandris, A., Merritt, D., et al. 2007, NewA, 12, 357) courtesy of Yohai Meiron
            if method=='casertano' use (Casertano, S., Hut, P. 1985, ApJ, 298, 80)

        Parameters
        ----------
        cluster : class
            StarCluster
        xstart,ystart,zstart : float
            starting position for centre (default: 0,0,0)
        vxstart,vystart,vzstart : float
            starting velocity for centre (default: 0,0,0)
        indx: bool
            subset of stars to perform centre of density calculation on (default: None)
        nsphere : int
            number of stars in centre sphere (default:100)
        rmin : float
            minimum radius of sphere around which to estimate density centre (default: 0.1 cluster.units)
        rmax : float
            maxmimum radius of sphere around which to estimate density centre (default: None cluster.units, uses maximum r)
        nmax : float
            maximum number of iterations (default:100)

        method : str
            method for finding the centre of density ('harfst' (default), 'casertano')

        if method=='casertano'
            nneighbour : int
                number of neighbours for calculation local densities

        reset_centre : bool
            forcibly reset cluster's centre of mass (default: False)


        Returns
        -------
        xc,yc,zc,vxc,vyc,vzc : float
            coordinates of centre of mass

        HISTORY
        -------
        2019 - Written - Webb (UofT) with Yohai Meiron (UofT)
        2022 - Written - Webb (UofT) - add method=='casertano'
        """

        xc,yc,zc,vxc,vyc,vzc=find_centre_of_density(self,xstart=xstart,
            ystart=ystart,zstart=zstart,vxstart=vxstart,vystart=vystart,vzstart=vzstart,indx=indx,
            nsphere=nsphere,rmin=rmin,rmax=rmax,nmax=nmax)

        if self.origin=='cluster':
            self.xc, self.yc, self.zc = xc,yc,zc
            self.vxc, self.vyc, self.vzc = vxc,vyc,vzc

            return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc


        elif self.origin == "galaxy":

            if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0.,0.,0.,0.,0.,0.) or reset_centre:
                self.xgc, self.ygc, self.zgc = xc,yc,zc
                self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                self.xc, self.yc, self.zc = 0.0, 0.0, 0.0
                self.vxc, self.vyc, self.vzc = 0.0, 0.0, 0.0

                return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc
     
            else:
                self.xc,self.yc,self.zc=xc-self.xgc,yc-self.ygc,zc-self.zgc
                self.vxc,self.vyc,self.vzc=vxc-self.vxgc,vyc-self.vygc,vzc-self.vzgc


                return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc

        elif self.origin=='sky':
            if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0.,0.,0.,0.,0.,0.) or reset_centre:
                self.xgc, self.ygc, self.zgc = xc,yc,zc
                self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                self.xc, self.yc, self.zc = 0.0, 0.0, 0.0
                self.vxc, self.vyc, self.vzc = 0.0, 0.0, 0.0

                return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc
            else:
                print('No Cluster Variables Set')
                return xc,yc,zc,vxc,vyc,vzc

        else:
            print('No Cluster Variables Set')
            return xc,yc,zc,vxc,vyc,vzc

    def find_centre_of_mass(self,reset_centre=False):
        """ Find the centre of mass of the cluster

        Parameters
        ----------
        cluster : class
            StarCluster

        Returns
        -------
        xc,yc,zc,vxc,vyc,vzc : float
            coordinates of centre of mass

        HISTORY
        -------
        2018 - Written - Webb (UofT)
        """
        xc,yc,zc,vxc,vyc,vzc=find_centre_of_mass(self)

        if self.origin=='cluster':
            self.xc, self.yc, self.zc = xc,yc,zc
            self.vxc, self.vyc, self.vzc = vxc,vyc,vzc

            return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc


        elif self.origin == "galaxy":

            if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0.,0.,0.,0.,0.,0.) or reset_centre:
                self.xgc, self.ygc, self.zgc = xc,yc,zc
                self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                self.xc, self.yc, self.zc = 0.0, 0.0, 0.0
                self.vxc, self.vyc, self.vzc = 0.0, 0.0, 0.0

                return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc
     
            else:
                self.xc,self.yc,self.zc=xc-self.xgc,yc-self.ygc,zc-self.zgc
                self.vxc,self.vyc,self.vzc=vxc-self.vxgc,vyc-self.vygc,vzc-self.vzgc


                return self.xc, self.yc, self.zc,self.vxc, self.vyc, self.vzc

        elif self.origin=='sky':
            if (self.xgc, self.ygc, self.zgc, self.vxgc, self.vygc, self.vzgc)==(0.,0.,0.,0.,0.,0.) or reset_centre:
                self.xgc, self.ygc, self.zgc = xc,yc,zc
                self.vxgc, self.vygc, self.vzgc = vxc, vyc, vzc
                self.xc, self.yc, self.zc = 0.0, 0.0, 0.0
                self.vxc, self.vyc, self.vzc = 0.0, 0.0, 0.0

                return self.xgc, self.ygc, self.zgc,self.vxgc, self.vygc, self.vzgc
            else:
                print('No Cluster Variables Set')
                return xc,yc,zc,vxc,vyc,vzc
        else:
            print('No Cluster Variables Set')
            return xc,yc,zc,vxc,vyc,vzc

    def relaxation_time(self, rad=None, coulomb=0.4, projected=None,method='spitzer'):
        """Calculate the relaxation time (Spitzer & Hart 1971) within a given radius of the cluster

        - Spitzer, L. Jr, Hart, M.H. 1971, ApJ, 164, 399 (Equation 5)
        - Need to adjust amplitude for different input units
        Parameters
        ----------
        cluster : class
          StarCluster
        rad : float
          radius within which to calculate the relaxation time (defult: cluster.rm)
        coulomb : float
          Coulomb parameter (default: 0.4)
        projected : bool
          use projected values (default: False)
        method : str
          choose between Spitzer & Hart 1971 and other methods (in development)

        Returns
        -------
           trelax : float
              relaxation time within radius rad

        History
        -------
        2020 - Written - Webb (UofT)

        """
        if projected==None:
            projected=self.projected

        self.trelax=relaxation_time(self, rad=rad, coulomb=coulomb, projected=projected,method='spitzer')

        return self.trelax

    def half_mass_relaxation_time(self, coulomb=0.4, projected=None):
        """ Calculate the half-mass relaxation time (Spitzer 1987) of the cluster
        - Spitzer, L. 1987, Dynamical evolution of globular clusters
        - Need to adjust amplitude for different input units

        Parameters
        ----------
        cluster : class
          StarCluster

        coulomb : float
          Coulomb parameter (default: 0.4)

        projected : bool
          use projected values (default: False)

        Returns
        -------
           trh : float
              half-mass relaxation time within radius rad

        History
        -------
           2019 - Written - Webb (UofT)

        """
        if projected==None:
            projected=self.projected
        self.trh=half_mass_relaxation_time(self, coulomb=coulomb, projected=projected)

        return self.trh


    def core_relaxation_time(self, coulomb=0.4, projected=None):
        """ Calculate the core relaxation time (Stone & Ostriker 2015) of the cluster
        
        - Stone, N.C. & Ostriker, J.P. 2015, ApJ, 806, 28

        Parameters

        cluster : class
          StarCluster

        coulomb : float
          Coulomb parameter (default: 0.4)

        projected : bool
          use projected values (default: False)

        method : str
          choose between Stone & Ostriker 2015 and other methods (in development)

        Returns

         trc

        History
        -------
        2019 - Written - Webb (UofT)

        """
        if projected==None:
            projected=self.projected
        self.trc=core_relaxation_time(self, coulomb=coulomb, projected=projected)

        return self.trc

    def energies(self, specific=True, ids=None,  projected=None,softening=0.0, full=True, parallel=False, **kwargs):
        """Calculate kinetic and potential energy of every star

        Parameters
        ----------
        cluster : class
          StarCluster instance
        specific : bool
          find specific energies (default: True)
        ids: boolean array or integer array
          if given, find the energues of a subset of stars defined either by an array of
          star ids, or a boolean array that can be used to slice the cluster (default: None)
        projected : bool
          use projected values (default: False)
        softening : float
          Plummer softening length in cluster.units (default: 0.0)
        full : bool
          calculate distance of full array of stars at once with numbra (default: True)
        parallel : bool
          calculate distances in parallel (default: False)

        Returns
        -------
        kin,pot : float
          kinetic and potential energy of every star if the i_d argument is not used. If i_d
          argument is used, return an arrays with potential and kinetic energy in the same shape
          of i_d
        History
        -------
           2019 - Written - Webb (UofT)
           2022 - Updated with support for multiple ids or an idexing array - Gillis (UofT)

        """

        if ids is None:
            ids=kwargs.get('i_d',None)

        if projected == None:
            projected=self.projected

        self.save_cluster()
        units0,origin0, rorder0, rorder_origin0 = self.units0,self.origin0, self.rorder0, self.rorder_origin0

        if self.origin0 != 'cluster' and self.origin0 != 'centre':
            self.to_cluster(sortstars=False)
        if self.units=='amuse':
            self.to_pckms()

        kin, pot = energies(self, specific=specific, ids=ids, softening=softening, full=full, 
                            projected=projected, parallel=parallel)
        
        if ids is not None:

            # Convert ids to boolean array if given as an array of ids
            if isinstance(ids,int) or isinstance(ids,float) or isinstance(ids,np.int64) or isinstance(ids,np.float64):
                ids = self.id == ids
            elif isinstance(ids[0],int) or isinstance(ids[0],float) or isinstance(ids[0],np.int64) or isinstance(ids[0],np.float64):
                ids = np.in1d(self.id, ids)

            kin_full, pot_full = np.zeros(self.ntot), np.zeros(self.ntot)
            kin_full[ids], pot_full[ids] = kin, pot
            

            self.return_cluster(units0,origin0, rorder0, rorder_origin0)

            if self.units=='amuse' and specific:
                kin_full=kin_full | (u.kms*u.kms)
                pot_full=pot_full | (u.kms*u.kms)
                kin=kin | (u.kms*u.kms)
                pot=pot | (u.kms*u.kms)
            elif self.units=='amuse' and not specific:
                kin_full=kin_full | (u.MSun*u.kms*u.kms)
                pot_full=pot_full | (u.MSun*u.kms*u.kms)
                kin=kin | (u.MSun*u.kms*u.kms)
                pot=pot | (u.MSun*u.kms*u.kms)

            self.add_energies(kin_full, pot_full)
            
            return kin, pot

        else:
            
            self.return_cluster(units0,origin0, rorder0, rorder_origin0)

            if self.units=='amuse' and specific:
                kin=kin | (u.kms*u.kms)
                pot=pot | (u.kms*u.kms)
            elif self.units=='amuse' and not specific:
                kin=kin | (u.MSun*u.kms*u.kms)
                pot=pot | (u.MSun*u.kms*u.kms)

            self.add_energies(kin, pot)

            return self.kin, self.pot

    def closest_star(self, projected=None, argmument=False):
        """Find distance to closest star for each star
        - uses numba

        Parameters
        ----------
        cluster : class
            positions of stars within the StarCluster
        projected : bool
          use projected values (default: False)
        argument : bool
          return argument of closest star as well (default: False)

        Returns
        -------
            minimum_distance : float
                distance to closest star for each star

            if argument:
                arg : int
                    argument of closest star for each star            

        History
        -------
           2019 - Written - Webb (UofT)
        """
        if projected==None:
            projected=self.projected

        if argument:
            self.dclosest,self.argclosest=closest_star(self, projected=projected,argument=argument)
            return self.dclosest,self.argclosest
        else:
            self.dclosest=closest_star(self, projected=projected)
            return self.dclosest

    def rlagrange(self, nlagrange=10, projected=None):
        """Calculate lagrange radii of the cluster by mass

        Parameters
        ----------
        cluster : class
            StarCluster
        nlagrange : int 
            number of lagrange radii bins (default: 10)
        mfrac : float
            Exact masss fraction to calculate radius. Will supersede nlagrange if not None (default : None)
        projected : bool
            calculate projected lagrange radii (default: False)

        Returns
        -------
        rn : float
            lagrange radii

        History
        -------
           2019 - Written - Webb (UofT)
        """
        if projected==None:
            projected=self.projected
        self.rn = rlagrange(self, nlagrange=nlagrange, projected=projected)

        return self.rn

    def virial_radius(self,method='inverse_distance', full=True, H=70.0, Om=0.3, overdens=200.0,
        nrad=20, projected=None, plot=False, **kwargs):
        """Calculate virial radius of the cluster
        - Virial radius is calculated using either:
        -- the average inverse distance between particles, weighted by their masses (default)
        -- the radius at which the density is equal to the critical density of the Universe at the redshift of the system, multiplied by an overdensity constant
        Parameters
        ----------
        cluster : class
            StarCluster
        method : str
            method for calculating virial radius (default: 'inverse_distance', alternative: 'critical_density')

        Returns
        -------
        rv : float
            virial radius

        Other Parameters
        ----------------
        full : bool
            Use Numba to calculate average inverse distance between stars (default:True)
        H : float
            Hubble constant
        Om : float
            density of matter
        overdens : float
            overdensity constant
        projected : bool
            calculate projected virial radius (default: False)
        plot : bool
            plot cluster density profile and illustrate virial radius calculation
        kwargs : str
            key word arguments for plotting function

        History
        -------
           2019 - Written - Webb (UofT)
        """
        if projected==None:
            projected=self.projected

        self.rv = virial_radius(self,method=method, full=full, H=H, Om=Om, overdens=overdens,
        nrad=nrad, projected=projected, plot=plot, **kwargs)

        return self.rv

    def mass_function(
        self,
        mmin=None,
        mmax=None,
        nmass=10,
        rmin=None,
        rmax=None,
        vmin=None,
        vmax=None,
        emin=None,
        emax=None,
        kwmin=0,
        kwmax=1,
        indx=None,
        mcorr=None,
        projected=False,
        plot=False,
        **kwargs
    ):
        """Find mass function over a given mass range

        - mass bins are set up so that there are an equal number of stars in each bin

        Parameters
        ----------
        cluster : class
            StarCluster instance
        mmin/mmax : float
            specific mass range
        nmass : 
            number of mass bins used to calculate alpha
        rmin/rmax : 
            specific radial range
        vmin/vmax : float
            specific velocity range
        emin/emax : float
            specific energy range
        kwmin/kwmax : int
            specific stellar evolution type range
        npop : int
            population number
        indx : bool 
            specific subset of stars
        projected : bool 
            use projected values (default: False)
        mcorr : float
            completeness correction for masses
        plot : bool 
            plot the mass function

        Returns
        -------
        m_mean : float
            mean mass in each bin
        m_hist : float
            number of stars in each bin
        dm : float
            dN/dm of each bin
        alpha : float
            power-law slope of the mass function (dN/dm ~ m^alpha)
        ealpha : float
            error in alpha
        yalpha : float
            y-intercept of fit to log(dN/dm) vs log(m)
        eyalpha : float
            error in yalpha

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2018 - Written - Webb (UofT)
        """
        self.mmin=mmin
        self.mmax=mmax

        if projected==None:
            projected=self.projected

        m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(
            self,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rmin,
            rmax=rmax,
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            indx=indx,
            mcorr=mcorr,
            projected=projected,
            plot=plot,
            title="GLOBAL",
            **kwargs
        )
        self.alpha = alpha

        return self.alpha

    def tapered_mass_function(
        cluster,
        mmin=None,
        mmax=None,
        nmass=10,
        rmin=None,
        rmax=None,
        vmin=None,
        vmax=None,
        emin=None,
        emax=None,
        kwmin=0,
        kwmax=1,
        npop=None,
        indx=None,
        projected=False,
        mcorr=None,
        plot=False,
        **kwargs
    ):
        """Find a tapered mass function over a given mass range

        - mass bins are set up so that there are an equal number of stars in each bin
        - functional form of the tapered mass function is taken from De Marchi, Paresce & Portegies Zwart 2010
        Parameters
        ----------
        cluster : class
            StarCluster instance
        mmin/mmax : float
            specific mass range
        nmass : 
            number of mass bins used to calculate alpha
        rmin/rmax : 
            specific radial range
        vmin/vmax : float
            specific velocity range
        emin/emax : float
            specific energy range
        kwmin/kwmax : int
            specific stellar evolution type range
        npop : int
            population number
        indx : bool 
            specific subset of stars
        projected : bool 
            use projected values (default: False)
        mcorr : float
            completeness correction for masses
        plot : bool 
            plot the mass function

        Returns
        -------
        m_mean : float
            mean mass in each bin
        m_hist : float
            number of stars in each bin
        dm : float
            dN/dm of each bin
        alpha : float
            power-law slope of the mass function (dN/dm ~ m^alpha)
        ealpha : float
            error in alpha
        yalpha : float
            y-intercept of fit to log(dN/dm) vs log(m)
        eyalpha : float
            error in yalpha

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2018 - Written - Webb (UofT)
        """

        self.mmin=mmin
        self.mmax=mmax

        if projected==None:
            projected=self.projected

        m_mean, m_hist, dm, A, eA, alpha, ealpha, mc, emc, beta, ebeta= tapered_mass_function(
            cluster,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rmin,
            rmax=rmax,
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            npop=npop,
            indx=indx,
            projected=projected,
            mcorr=mcorr,
            plot=plot,
            **kwargs
        )

        self.alpha=alpha

        return self.alpha

    def eta_function(
        self,
        mmin=None,
        mmax=None,
        nmass=10,
        rmin=None,
        rmax=None,
        vmin=None,
        vmax=None,
        emin=None,
        emax=None,
        kwmin=0,
        kwmax=1,
        indx=None,
        projected=False,
        plot=False,
        **kwargs
    ):
        """
        NAME: Find power slope of velocity dispersion versus mass
        
        - mass bins are set up so that there are an equal number of stars in each bin

        Parameters
        ----------
        cluster : class
            StarCluster instance
        mmin/mmax : float
            specific mass range
        nmass : 
            number of mass bins used to calculate alpha
        rmin/rmax : 
            specific radial range
        vmin/vmax : float
            specific velocity range
        emin/emax : float
            specific energy range
        kwmin/kwmax : int
            specific stellar evolution type range
        npop : int
            population number
        indx : bool 
            specific subset of stars
        projected : bool 
            use projected values (default: False)
        plot : bool 
            plot the mass function

        Returns
        -------
        m_mean : float
            mean mass in each bin
        sigvm : float
            velocity dispersion of stars in each bin
        eta : float
            power-law slope of (sigvm ~ m^eta)
        eeta : float
            error in eta
        yeta : float
            y-intercept of fit to log(sigvm) vs log(m)
        eeta : float
            error in yeta

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2018 - Written - Webb (UofT)
        """
        self.mmin=mmin
        self.mmax=mmax

        m_mean, sigvm, eta, eeta, yeta, eyeta = eta_function(
            self,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rmin,
            rmax=rmax,
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            indx=indx,
            projected=projected,
            plot=plot,
            **kwargs
        )

        self.eta = eta

        return self.eta

    def meq_function(
        self,
        mmin=None,
        mmax=None,
        nmass=10,
        rmin=None,
        rmax=None,
        vmin=None,
        vmax=None,
        emin=None,
        emax=None,
        kwmin=0,
        kwmax=1,
        indx=None,
        projected=False,
        plot=False,
        **kwargs
    ):

        """
        NAME: Find meq from velocity dispersion versus mass
        
        - mass bins are set up so that there are an equal number of stars in each bin
        - As per Bianchini, P. et al. 2016, MNRAS, 458, 3644, velocity dispersion 
          versus mass is fit with the following:
          sigma(m)= sigma e^(-1/2 m/meq) if m<= meq
                  = sigma0 e^(-1/2) (m/meq)^-1/2 if m > meq

        Parameters
        ----------
        cluster : class
            StarCluster instance
        mmin/mmax : float
            specific mass range
        nmass : 
            number of mass bins used to calculate alpha
        rmin/rmax : 
            specific radial range
        vmin/vmax : float
            specific velocity range
        emin/emax : float
            specific energy range
        kwmin/kwmax : int
            specific stellar evolution type range
        npop : int
            population number
        indx : bool 
            specific subset of stars
        projected : bool 
            use projected values (default: False)
        plot : bool 
            plot the mass function

        Returns
        -------
        m_mean : float
            mean mass in each bin
        sigvm : float
            velocity dispersion of stars in each bin
        meq : float
            Bianchini fit to sigvm vs m
        emeq : float
            error in Bianchini fit to sigvm vs m
        sigma0 : float
            Bianchini fit to sigvm vs m
        esigma0 : float
            error in Bianchini fit to sigvm vs m

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2020
        """

        self.mmin=mmin
        self.mmax=mmax

        m_mean, sigvm, meq, emeq, sigma0, esigma0 = meq_function(
            self,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rmin,
            rmax=rmax,
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            indx=indx,
            projected=projected,
            plot=plot,
            **kwargs
        )
        self.meq = meq

        return self.meq

    def ckin(
        self,
        mmin=None,
        mmax=None,
        nmass=10,
        rmin=None,
        rmax=None,
        vmin=None,
        vmax=None,
        emin=None,
        emax=None,
        kwmin=0,
        kwmax=1,
        indx=None,
        projected=False,
    ):
        """
        NAME: Find the kinematic concentration parameter ck
        
        - see Bianchini et al. 2018, MNRAS, 475, 96

        Parameters
        ----------
        cluster : class
            StarCluster instance
        mmin/mmax : float
            specific mass range
        nmass : 
            number of mass bins used to calculate alpha
        rmin/rmax : 
            specific radial range
        vmin/vmax : float
            specific velocity range
        emin/emax : float
            specific energy range
        kwmin/kwmax : int
            specific stellar evolution type range
        npop : int
            population number
        indx : bool 
            specific subset of stars
        projected : bool 
            use projected values (default: False)

        Returns
        -------
        ck : float
            kinematic concentration

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2020
        """
        c_k = ckin(
            self,
            mmin=mmin,
            mmax=mmax,
            nmass=nmass,
            rmin=rmin,
            rmax=rmax,
            vmin=vmin,
            vmax=vmax,
            emin=emin,
            emax=emax,
            kwmin=kwmin,
            kwmax=kwmax,
            indx=indx,
            projected=projected,
        )

        self.ck=c_k

        return self.ck

    def rcore(
        self,
        method='casertano',
        nneighbour=6,
        mfrac=0.1,
        projected=False,
        plot=False,
        **kwargs
    ):
        """Calculate core radius of the cluster
        -- The default method (method='casertano') follows Casertano, S., Hut, P. 1985, ApJ, 298, 80 to find the core
        -- An alternative metrhod (method=='isothermal') assumes the cluster is an isothermal sphere the core radius is where density drops to 1/3 central value
        --- For projected core radius, the core radius is where the surface density profile drops to 1/2 the central value
        --- Note that the inner mass fraction of stars used to calculate central density is set by mfrac (default 0.1 = 10%)

        Parameters
        ----------
        cluster : class
            StarCluster instance
        method : string
            method of calculating the core radius of a star cluster (default 'casertano')
        if method =='casertano':
            nneighbour : int
                number of neighbours for calculation local densities
        if method=='isothermal':
            mfrac : float
                inner mass fraction to be used to establish the central density
        projected : bool
            use projected values (default: False)
        plot : bool
            plot the density profile and mark the core radius of the cluster (default: False)
        Returns
        -------
        rc : float
            core radius

        Other Parameters
        ----------------
        None

        History
        -------
        2021 - Written - Webb (UofT)
        """


        self.rc = rcore(
            self,
            method=method,
            nneighbour=nneighbour,
            mfrac=mfrac,
            projected=projected,
            plot=plot,
            **kwargs
        )

        return self.rc

    def rtidal(self, pot=None, rtiterate=0, rtconverge=0.9, indx=None, rgc=None, from_centre=False, plot=False, verbose=False, **kwargs):
        """Calculate tidal radius of the cluster
        - The calculation uses Galpy (Bovy 2015_, which takes the formalism of Bertin & Varri 2008 to calculate the tidal radius
        -- Bertin, G. & Varri, A.L. 2008, ApJ, 689, 1005
        -- Bovy J., 2015, ApJS, 216, 29
        - riterate = 0 corresponds to a single calculation of the tidal radius based on the cluster's mass (np.sum(cluster.m))
        -- Additional iterations take the mass within the previous iteration's calculation of the tidal radius and calculates the tidal
           radius again using the new mass until the change is less than 90%
        - for cases where the cluster's orbital parameters are not set, it is possible to manually set rgc which is assumed to be in kpc.

        Parameters
        ----------
        cluster : class
            StarCluster instance
        pot : class 
            GALPY potential used to calculate tidal radius (default: None)
        rtiterate : int
            how many times to iterate on the calculation of r_t (default: 0)
        rtconverge : float
            criteria for tidal radius convergence within iterations (default 0.9)
        indx : bool
            subset of stars to use when calculate the tidal radius (default: None)
        rgc : float
            Manually set galactocentric distance in kpc at which the tidal radius is to be evaluated (default: None)
        zgc : float
            For non-spherically symmetric potentials, manually set distance in kpc above disk at which the tidal radius is to be evaluated. When set, rgc becomes radius in cylindrical coordinates (default: None)
        from_centre : bool
            calculate tidal radius based on location of cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
        plot : bool
            plot the x and y coordinates of stars and mark the tidal radius of the cluster (default: False)
        verbose : bool
            Print information about iterative calculation of rt

        Returns
        -------
        rt : float
            tidal radius

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2019 - Written - Webb (UofT)
        """

        self.rt = rtidal(self, pot=pot, rtiterate=rtiterate,rtconverge=rtconverge, indx=indx, rgc=rgc, from_centre=from_centre, plot=plot, verbose=verbose, **kwargs)

        return self.rt

    def rlimiting(
        self,
        pot=None,
        rgc=None,
        nrad=20,
        projected=False,
        plot=False,
        from_centre=False,
        verbose=False,
        **kwargs
    ):
        """Calculate limiting radius of the cluster
           
        - The limiting radius is defined to be where the cluster's density reaches the local background density of the host galaxy
        - for cases where the cluster's orbital parameters are not set, it is possible to manually set rgc which is assumed to be in kpc.

        Parameters
        ----------

        cluster : class
            StarCluster
        pot : class 
            GALPY potential used to calculate actions
        rgc : 
            Manually set galactocentric distance in kpc at which the tidal radius is to be evaluated (default: None)
        zgc : float
            For non-spherically symmetric potentials, manually set distance in kpc above disk at which the tidal radius is to be evaluated. When set, rgc becomes radius in cylindrical coordinates (default: None)
        nrad : int
            number of radial bins used to calculate density profile (Default: 20)
        projected : bool
            use projected values (default: False)
        plot : bool
            plot the density profile and mark the limiting radius of the cluster (default: False)
        from_centre : bool
            calculate tidal radius based on location of cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
        verbose : bool
        Returns
        -------
            rl : float
                limiting radius

        Other Parameters
        ----------------
        kwargs : str
            key words for plotting

        History
        -------
        2019 - Written - Webb (UofT)
        """
        self.rl = rlimiting(
            self,
            pot=pot,
            rgc=rgc,
            nrad=nrad,
            projected=projected,
            plot=plot,
            from_centre=from_centre,
            verbose=verbose,
            **kwargs
        )

        return self.rl

    def initialize_orbit(self, from_centre=False, ro=None, vo=None,zo=None, solarmotion=None):
        """ Initialize a galpy orbit instance for the cluster

        Parameters
        ----------
        cluster : class
            StarCluster
        from_centre : bool
            intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
        ro : float
            galpy distance scale (default: None)
        vo : float
            galpy velocity scale (default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)

        Returns
        -------
        orbit : class
            GALPY orbit

        History
        -------
        2018 - Written - Webb (UofT)
        """
        self.orbit=initialize_orbit(self, from_centre=from_centre, ro=ro, vo=vo, zo=zo, solarmotion=solarmotion)
        return self.orbit

    def initialize_orbits(self,ro=None, vo=None, zo=None,solarmotion=None):
        """Initialize a galpy orbit for every star in the cluster

        Parameters
        ----------
        cluster : class
            StarCluster
        ro : float
            galpy distance scale (default: None)
        vo : float
            galpy velocity scale (default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)

        Returns
        -------
        orbit : class
            GALPY orbit

        History
        -------
        2018 - Written - Webb (UofT)
        """
        self.orbits=initialize_orbits(self,ro=ro, vo=vo, zo=zo, solarmotion=solarmotion)
        return self.orbits

    def interpolate_orbit(self,pot=None,tfinal=None,nt=1000,from_centre=False, ro=None,vo=None,zo=None,solarmotion=None):
        """
        Interpolate past or future position of cluster and escaped stars

        Parameters
        ----------
        cluster : class
            StarCluster
        cluster_pot : class
            Galpy potential for host cluster that orbit is to be integrated in
            if None, assume a Plumme Potential
        pot : class
            galpy Potential that orbit is to be integrate in (default: None)
        tfinal : float
            final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
        nt : int
            number of timesteps
        from_centre : bool
            intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)

        Returns
        -------
        x,y,z : float
            interpolated positions of each star
        vx,vy,vz : float
            interpolated velocities of each star    

        History
        -------
        2021 - Written - Webb (UofT)
        """
        xgc,ygc,zgc,vxgc,vygc,vzgc=interpolate_orbit(self,pot=pot,tfinal=tfinal,nt=nt,from_centre=from_centre, ro=ro,vo=vo, zo=zo, solarmotion=solarmotion)

        origin0=self.origin
        self.to_cluster()

        self.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

        if self.units=='radec' and self.origin=='sky':
            self.ra_gc,self.dec_gc,self.dist_gc=xgc,ygc,zgc
            self.pmra_gc,self.pmdec_gc,self.vlos_gc=vxgc,vygc,vzgc

        self.to_origin(origin0)

        return self.xgc,self.ygc,self.zgc,self.vxgc,self.vygc,self.vzgc

    def orbit_interpolate(self,pot=None,tfinal=None,nt=1000,from_centre=False,ro=None,vo=None,zo=None,solarmotion=None):
        """
        Interpolate past or future position of cluster and escaped stars

        - same as interpolate_orbit, but included for legacy purposes

        Parameters
        ----------
        cluster : class
            StarCluster
        cluster_pot : class
            Galpy potential for host cluster that orbit is to be integrated in
            if None, assume a Plumme Potential
        pot : class
            galpy Potential that orbit is to be integrate in (default: None)
        tfinal : float
            final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
        nt : int
            number of timesteps
        from_centre : bool
            intialize orbits from cluster's exact centre instead of cluster's position in galaxy (default :False)
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)

        Returns
        -------
        x,y,z : float
            interpolated positions of each star
        vx,vy,vz : float
            interpolated velocities of each star    

        History
        -------
        2021 - Written - Webb (UofT)
        """
        self.interpolate_orbit(self,pot=pot,tfinal=tfinal,nt=nt,from_centre=from_centre,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)


    def interpolate_orbits(self,pot=None,tfinal=None,nt=1000,ro=None,vo=None,zo=None,solarmotion=None):
        """
        Interpolate past or future position of stars within the cluster

        Parameters
        ----------
        cluster : class
            StarCluster
        pot : class
            Galpy potential for host cluster that orbit is to be integrated in
            if None, assume a Plumme Potential
        tfinal : float
            final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
        nt : int
            number of timesteps
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)

        Returns
        -------
        x,y,z : float
            interpolated positions of each star
        vx,vy,vz : float
            interpolated velocities of each star    

        History
        -------
        2021 - Written - Webb (UofT)
        """
        x,y,z,vx,vy,vz=interpolate_orbits(self,pot=pot,tfinal=tfinal,nt=nt,ro=ro,vo=vo, zo=zo, solarmotion=solarmotion)

        if (self.origin=='centre' or self.origin=='cluster') and self.units!='radec':
            self.x,self.y,self.z=x,y,z
            self.vx,self.vy,self.vz=vx,vy,vz
        elif self.origin=='sky' and self.units=='radec':
            self.to_sky()
            self.ra,self.dec,self.dist=x,y,z
            self.pmra,self.pmdec,self.vlos=vx,vy,vz
            self.x,self.y,self.z=x,y,z
            self.vx,self.vy,self.vz=vx,vy,vz
        elif self.units=='radec' and (self.origin=='centre' or self.origin=='cluster'):
            print('CANT INTEGRATE ORBITS WITH FROM_CENTRE OR FROM_CLUSTER AND RETURN IN SKY COORDINATES')
        else:
            self.x,self.y,self.z=x,y,z
            self.vx,self.vy,self.vz=vx,vy,vz


        if self.origin!='centre' and self.origin!='cluster' :

            if pot is None:
                xgc,ygc,zgc,vxgc,vygc,vzgc=interpolate_orbit(self,pot=None,tfinal=tfinal,nt=nt, ro=ro,vo=vo,zo=zo, solarmotion=solarmotion)
            else:
                xgc,ygc,zgc,vxgc,vygc,vzgc=interpolate_orbit(self,pot=pot,tfinal=tfinal,nt=nt, ro=ro,vo=vo,zo=zo, solarmotion=solarmotion)

            self.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

            if self.units=='radec' and self.origin=='sky':
                self.ra_gc,self.dec_gc,self.dist_gc=xgc,ygc,zgc
                self.pmra_gc,self.pmdec_gc,self.vlos_gc=vxgc,vygc,vzgc

        return self.x,self.y,self.z,self.vx,self.vy,self.vz


    def orbits_interpolate(self,pot=None,tfinal=None,nt=1000,ro=None,vo=None,zo=None,solarmotion=None):
        """
        Interpolate past or future position of stars within the cluster

        - same as interpolate_orbits, but kept for legacy purposes

        Parameters
        ----------
        cluster : class
            StarCluster
        pot : class
            Galpy potential for host cluster that orbit is to be integrated in
            if None, assume a Plumme Potential
        tfinal : float
            final time (in cluster.units) to integrate orbit to (default: 12 Gyr)
        nt : int
            number of timesteps
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)

        Returns
        -------
        x,y,z : float
            interpolated positions of each star
        vx,vy,vz : float
            interpolated velocities of each star    

        History
        -------
        2021 - Written - Webb (UofT)
        """
        self.interpolate_orbits(self,pot=pot,tfinal=tfinal,nt=nt,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)

    def orbital_path(self,tfinal=0.1,nt=1000,pot=None,from_centre=False,
        skypath=False,initialize=False,ro=None,vo=None,zo=None,solarmotion=None, plot=False,**kwargs):
        """Calculate the cluster's orbital path

        Parameters
        ----------
        cluster : class
            StarCluster
        tfinal : float
            final time (in cluster.units) to integrate orbit to (default: 0.1 Gyr)
        nt : int
            number of timesteps
        pot : class
            galpy Potential that orbit is to be integrate in (default: None)
        from_centre : bool
            genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
        skypath : bool
            return sky coordinates instead of cartesian coordinates (default: False)
        initialize : bool
            Initialize and return Orbit (default: False)
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)
        plot : bool
            plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)

        Returns
        -------
        t : float
            times for which path is provided
        x,y,z : float
            orbit positions
        vx,vy,vz : float
            orbit velocity
        o : class
            galpy orbit (if initialize==True)
        History
        -------
        2018 - Written - Webb (UofT)
        """
        if initialize:
            self.tpath,self.xpath,self.ypath,self.zpath,self.vxpath,self.vypath,self.vzpath,self.orbit=orbital_path(self,tfinal=tfinal,nt=nt,pot=pot,from_centre=from_centre,skypath=skypath,initialize=initialize,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion,plot=plot,**kwargs)
        else:
            self.tpath,self.xpath,self.ypath,self.zpath,self.vxpath,self.vypath,self.vzpath=orbital_path(self,tfinal=tfinal,nt=nt,pot=pot,from_centre=from_centre,skypath=skypath,initialize=initialize,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion,plot=plot,**kwargs)

        return  self.tpath,self.xpath,self.ypath,self.zpath,self.vxpath,self.vypath,self.vzpath

    def orbital_path_match(self,tfinal=0.1,nt=1000,pot=None,path=None,from_centre=False,
        skypath=False,to_path=False,do_full=False,ro=None,vo=None,zo=None,solarmotion=None,plot=False,projected=False,**kwargs):

        """Match stars to a position along the orbital path of the cluster

        Parameters
        ----------
        cluster : class
            StarCluster
        tfinal : float
            final time (in cluster.units) to integrate orbit to (default: 0.1 Gyr)
        nt : int
            number of timesteps
        pot : class
            galpy Potential that orbit is to be integrate in (default: None)
        path : array
            array of (t,x,y,x,vx,vy,vz) corresponding to the tail path. If none path is calculated (default: None)
        from_centre : bool
            genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
        skypath : bool
            return sky coordinates instead of cartesian coordinates (default: False)
            if True, projected is set to True
        to_path : bool
            measure distance to the path itself instead of distance to central point along the path (default: False)
        do_full : bool
            calculate dpath all at once in a single numpy array (can be memory intensive) (default:False)
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)
        plot : bool
            plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)
        projected : bool
            match to projected orbital path, which means matching just x and y coordinates or Ra and Dec coordinates (not z, or dist) (default:False)

        Returns
        -------
        tstar : float
            orbital time associated with star
        dprog : float
            distance along the orbit to the progenitor
        dpath : 
            distance to centre of the orbital path bin (Default) or the orbit path (to_path = True)

        History
        -------
        2018 - Written - Webb (UofT)
        """

        self.tstar,self.dprog,self.dpath=orbital_path_match(self,tfinal=tfinal,nt=nt,pot=pot,path=path,
        from_centre=from_centre,skypath=skypath,to_path=to_path,do_full=do_full,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion,plot=plot,projected=projected,**kwargs)

        return self.tstar,self.dprog,self.dpath

    def to_tail(self):
        """Calculate positions and velocities of stars when rotated such that clusters velocity vector
           points along x-axis

        - no change to coordinates in StarCluster

        Parameters
        ----------
        cluster : class
            StarCluster

        Returns
        -------
        x_tail,y_tail,z_tail,vx_tail,vy_tail,vz_tail : float
            rotated coordinates with cluster's velocity vector point along x-axis

        History:
        -------
        2018 - Written - Webb (UofT)

        """
        self.x_tail,self.y_tail,self.z_tail,self.vx_tail,self.vy_tail,self.vz_tail=to_tail(self)
        self.r_tail = np.sqrt(self.x_tail ** 2.0 + self.y_tail ** 2.0 + self.z_tail ** 2.0)
        self.v_tail = np.sqrt(self.vx_tail ** 2.0 + self.vy_tail ** 2.0 + self.vz_tail ** 2.0)

    def tail_path(self,dt=0.1,no=1000,nt=100,ntail=100,pot=None,dmax=None,bintype='fix',from_centre=False,skypath=False,
        to_path=False,
        do_full=False,
        ro=None,vo=None,zo=None,solarmotion=None,plot=False,**kwargs):

        """Calculate tail path +/- dt Gyr around the cluster

                Parameters
            ----------
            cluster : class
                StarCluster
            dt : float
                timestep that StarCluster is to be moved to
            no : int
                number of timesteps for orbit integration (default:1000)
            nt : int
                number of points along the tail to set the tail spacing (default: 100)
            ntail : int
                number of points along the tail with roaming average (default: 1000)
            pot : class
                galpy Potential that orbit is to be integrate in (default: None)
            dmax : float
                maximum distance (assumed to be same units as cluster) from orbital path to be included in generating tail path (default: None)
            bintype : str
                type of binning for tail stars (default : 'fix')
            from_centre : bool
                genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
            skypath : bool
                return sky coordinates instead of cartesian coordinates (default: False)
            to_path : bool
                measure distance to the path itself instead of distance to central point along the path (default: False)
            do_full : bool
                calculate dpath all at once in a single numpy array (can be memory intensive) (default:False)
            ro :float 
                galpy distance scale (Default: None)
            vo : float
                galpy velocity scale (Default: None)
            zo : float
                Sun's distance above the Galactic plane (default: None)
            solarmotion : float
                array representing U,V,W of Sun (default: None)
            plot : bool
                plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)
            projected : bool
                match to projected orbital path, which means matching just x and y coordinates or Ra and Dec coordinates (not z, or dist) (default:False)


            Returns
            -------
            t : float
                times for which path is provided
            x,y,z : float
                tail path positions
            vx,vy,vz : float
                tail path velocities
            History
            -------
            2018 - Written - Webb (UofT)
            2019 - Implemented numpy array preallocation to minimize runtime - Nathaniel Starkman (UofT)
            """

        self.tpath,self.xpath,self.ypath,self.zpath,self.vxpath,self.vypath,self.vzpath=tail_path(self,dt=dt,no=no,nt=nt,ntail=ntail,pot=pot,dmax=dmax,bintype=bintype,from_centre=from_centre,skypath=skypath,to_path=to_path,do_full=do_full,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion,plot=plot,**kwargs)

        return self.tpath,self.xpath,self.ypath,self.zpath,self.vxpath,self.vypath,self.vzpath

    def tail_path_match(self,dt=0.1,no=1000,nt=100,ntail=100, pot=None,dmax=None,from_centre=False,
        to_path=False,do_full=False,ro=None,vo=None,zo=None,solarmotion=None,plot=False,**kwargs,):

        """Match stars to a position along the tail path of the cluster

        Parameters
        ----------
        cluster : class
            StarCluster
        dt : float
            timestep that StarCluster is to be moved to
        no : int
            number of timesteps for orbit integration (default:1000)
        nt : int
            number of points along the tail to set the tail spacing (default: 100)
        ntail : int
            number of points along the tail with roaming average (default: 1000)
        pot : class
            galpy Potential that orbit is to be integrate in (default: None)
        path : array
            array of (t,x,y,x,vx,vy,vz) corresponding to the tail path. If none path is calculated (default: None)
        from_centre : bool
            genrate orbit from cluster's exact centre instead of its assigned galactocentric coordinates (default: False)
        skypath : bool
            return sky coordinates instead of cartesian coordinates (default: False)
            if True, projected is set to True
        to_path : bool
            measure distance to the path itself instead of distance to central point along the path (default: False)
        do_full : bool
            calculate dpath all at once in a single numpy array (can be memory intensive) (default:False)
        ro :float 
            galpy distance scale (Default: None)
        vo : float
            galpy velocity scale (Default: None)
        zo : float
            Sun's distance above the Galactic plane (default: None)
        solarmotion : float
            array representing U,V,W of Sun (default: None)
        plot : bool
            plot a snapshot of the cluster in galactocentric coordinates with the orbital path (defualt: False)
        projected : bool
            match to projected orbital path, which means matching just x and y coordinates or Ra and Dec coordinates (not z, or dist) (default:False)

        Returns
        -------
        tstar : float
            orbital time associated with star
        dprog : float
            distance along the path to the progenitor
        dpath : 
            distance to centre of the tail path bin (default) or the tail path (to_path = True)

        History
        -------
        2018 - Written - Webb (UofT)
        """

        self.tstar,self.dprog,self.dpath=tail_path_match(self,dt=dt,no=no,nt=nt,ntail=ntail,pot=pot,dmax=dmax,
        from_centre=from_centre,to_path=to_path,do_full=do_full,ro=ro,vo=vo,zo=zo,solarmotion=solarmotion,plot=plot,**kwargs)

        return self.tstar,self.dprog,self.dpath

def _get_subset(
    cluster,
    rmin=None,
    rmax=None,
    mmin=None,
    mmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=15,
    npop=None,
    indx=[None],
    projected=False,
    **kwargs,
):

    """Generate a boolean array that corresponds to subset of star cluster members that meet a certain criteria

    Parameters
    ----------
    rmin/rmax : float
        minimum and maximum stellar radii
    mmin/mmax : float
        minimum and maximum stellar mass
    vmin/vmax : float
        minimum and maximum stellar velocity
    emin/emax : float
        minimum and maximum stellar energy
    kwmin/kwmax : int
        minimum and maximum stellar type (kw)
    npop : int
        population number
    indx : bool
        user defined boolean array from which to extract the subset
    projected : bool 
        use projected values and constraints (default:False)

    Returns
    -------
    indx : bool
        boolean array that is True for stars that meet the criteria

    History
    -------
    2022 - Written - Webb (UofT)

    """    
    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    if rmin == None:
        rmin = np.amin(r)
    if rmax == None:
        rmax = np.amax(r)
    if vmin == None:
        vmin = np.amin(v)
    if vmax == None:
        vmax = np.amax(v)
    if mmin == None:
        mmin = np.amin(cluster.m)
    if mmax == None:
        mmax = np.amax(cluster.m)

    if npop == None:
        npopindx = np.ones(cluster.ntot,dtype=bool)
    else:
        npopindx=(cluster.npop == npop)

    if emin == None and emax != None:
        eindx = cluster.etot <= emax
    elif emin != None and emax == None:
        eindx = cluster.etot >= emin
    elif emin != None and emax != None:
        eindx = (cluster.etot <= emax) * (cluster.etot >= emin)
    else:
        eindx = np.ones(cluster.ntot,dtype=bool)

    if len(cluster.kw) > 0:
        kwindx=((cluster.kw >= kwmin) * (cluster.kw <= kwmax))
    else:
        kwindx = np.ones(cluster.ntot,dtype=bool)

    if indx is None:
        indx = np.ones(cluster.ntot,dtype=bool)
    elif None in indx:
        indx = np.ones(cluster.ntot,dtype=bool)

    # Build subcluster containing only stars in the full radial and mass range:

    try:

        indx *= (
            (r >= rmin)
            * (r <= rmax)
            * (cluster.m >= mmin)
            * (cluster.m <= mmax)
            * (v >= vmin)
            * (v <= vmax)
            * npopindx
            * kwindx
            * eindx
        )

    except:
        print('SUBSET ERROR: ',rmin,rmax,mmin,mmax,vmin,vmax,np.sum(kwindx),np.sum(eindx))
        indx=-1

    return indx

def sub_cluster(
    cluster,
    rmin=None,
    rmax=None,
    mmin=None,
    mmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=15,
    npop=None,
    indx=[None],
    projected=False,
    sortstars=True,
    reset_centre=False,
    reset_nbody=False,
    reset_nbody_mass=False,
    reset_nbody_radii=False,
    reset_rvirial=False,
    reset_projected=False,
    **kwargs
):
    """Extract a sub population of stars from a StarCluster

    - automatically moves cluster to centre of mass, so all constraints are in clustercentric coordinates and current StarCluster.units

    Parameters
    ----------
    rmin/rmax : float
        minimum and maximum stellar radii
    mmin/mmax : float
        minimum and maximum stellar mass
    vmin/vmax : float
        minimum and maximum stellar velocity
    emin/emax : float
        minimum and maximum stellar energy
    kwmin/kwmax : int
        minimum and maximum stellar type (kw)
    npop : int
        population number
    indx : bool
        user defined boolean array from which to extract the subset
    projected : bool 
        use projected values and constraints (default:False)
    sortstars: bool
        order stars by radius (default: True)
    reset_centre : bool
        re-calculate cluster centre after extraction (default:False)
    reset_nbody : bool
        reset nbody scaling factors (default:False)
    reset_nbody_mass : bool 
        find new nbody scaling mass (default:False)
    reset_nbody_radii : bool
        find new nbody scaling radius (default:False)
    reset_rvirial : bool
        use virial radius to find nbody scaling radius (default: True)
    reset_projected : bool 
        use projected radii to find nbody scaling radius (default: False)

    Returns
    -------
    StarCluster

    History
    -------
    2018 - Written - Webb (UofT)

    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0


    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    """
    if rmin == None:
        rmin = np.amin(r)
    if rmax == None:
        rmax = np.amax(r)
    if vmin == None:
        vmin = np.amin(v)
    if vmax == None:
        vmax = np.amax(v)
    if mmin == None:
        mmin = np.amin(cluster.m)
    if mmax == None:
        mmax = np.amax(cluster.m)

    if emin == None and emax != None:
        eindx = cluster.etot <= emax
    elif emin != None and emax == None:
        eindx = cluster.etot >= emin
    elif emin != None and emax != None:
        eindx = (cluster.etot <= emax) * (cluster.etot >= emin)
    else:
        eindx = cluster.id > -1

    if None in indx:
        indx = cluster.id > -1

    indx *= (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m <= mmax)
        * (v >= vmin)
        * (v <= vmax)
        * eindx
    )

    if len(cluster.kw) > 0:
        indx*=((cluster.kw >= kwmin) * (cluster.kw <= kwmax))
    """

    indx=cluster.subset(rmin=rmin,rmax=rmax,vmin=vmin,vmax=vmax,mmin=mmin,mmax=mmax,emin=emin,emax=emax,kwmin=kwmin,kwmax=kwmax,npop=npop,indx=indx,projected=projected)


    if np.sum(indx) > 0:


        subcluster = StarCluster(
            cluster.tphys,
            units=cluster.units,
            origin=cluster.origin,
            ctype=cluster.ctype,
            ro=cluster._ro,
            vo=cluster._vo,
            zo=cluster._zo,
            solarmotion=cluster._solarmotion,
        )

        subcluster.add_stars(
            cluster.x[indx],
            cluster.y[indx],
            cluster.z[indx],
            cluster.vx[indx],
            cluster.vy[indx],
            cluster.vz[indx],
            cluster.m[indx],
            cluster.id[indx],
            cluster.m0[indx],
            cluster.npop[indx],
            sortstars=sortstars,
        )

        if len(cluster.ra)==len(cluster.x):
            subcluster.ra, subcluster.dec, subcluster.dist = (
                cluster.ra[indx],
                cluster.dec[indx],
                cluster.dist[indx],
            )
            subcluster.pmra, subcluster.pmdec, subcluster.vlos = (
                cluster.pmra[indx],
                cluster.pmdec[indx],
                cluster.vlos[indx],
            )

        subcluster.add_nbody6(cluster.nc,cluster.rc,cluster.rbar,
            cluster.rtide,cluster.xc,cluster.yc,cluster.zc,
            cluster.zmbar,cluster.vbar,cluster.tbar,cluster.rscale,
            cluster.ns,cluster.nb,cluster.n_p)

        subcluster.projected = cluster.projected
        subcluster.centre_method = cluster.centre_method

        if len(cluster.logl) > 0:
            if len(cluster.ep) !=0 and len(cluster.ospin) != 0:
                subcluster.add_sse(
                    cluster.kw[indx],
                    cluster.logl[indx],
                    cluster.logr[indx],
                    cluster.ep[indx],
                    cluster.ospin[indx],
                )
            else:
                subcluster.add_sse(
                        cluster.kw[indx],
                        cluster.logl[indx],
                        cluster.logr[indx],
                    )
        elif len(cluster.kw) > 0:
            subcluster.kw = cluster.kw[indx]


        if len(cluster.id2) > 0:
            bindx1 = np.in1d(cluster.id1, cluster.id[indx])
            bindx2 = np.in1d(cluster.id2, cluster.id[indx])
            bindx=np.logical_or(bindx1,bindx2)


            if len(cluster.ep1) !=0 and len(cluster.ospin1) != 0:

                subcluster.add_bse(
                    cluster.id1[bindx],
                    cluster.id2[bindx],
                    cluster.kw1[bindx],
                    cluster.kw2[bindx],
                    cluster.kcm[bindx],
                    cluster.ecc[bindx],
                    cluster.pb[bindx],
                    cluster.semi[bindx],
                    cluster.m1[bindx],
                    cluster.m2[bindx],
                    cluster.logl1[bindx],
                    cluster.logl2[bindx],
                    cluster.logr1[bindx],
                    cluster.logr2[bindx],
                    cluster.ep1[bindx],
                    cluster.ep2[bindx],
                    cluster.ospin1[bindx],
                    cluster.ospin2[bindx],
                )
            else:
                subcluster.add_bse(
                    cluster.id1[bindx],
                    cluster.id2[bindx],
                    cluster.kw1[bindx],
                    cluster.kw2[bindx],
                    cluster.kcm[bindx],
                    cluster.ecc[bindx],
                    cluster.pb[bindx],
                    cluster.semi[bindx],
                    cluster.m1[bindx],
                    cluster.m2[bindx],
                    cluster.logl1[bindx],
                    cluster.logl2[bindx],
                    cluster.logr1[bindx],
                    cluster.logr2[bindx],
                )

        if len(cluster.etot) > 0:
            subcluster.add_energies(
                cluster.kin[indx], cluster.pot[indx],
            )

        if cluster.give == 'mxvpqael':
            subcluster.give=cluster.give
            subcluster.gyrpot=cluster.gyrpot[indx]
            subcluster.gyrq=cluster.gyrq[indx]
            subcluster.gyracc=cluster.gyracc[indx]
            subcluster.eps=cluster.eps[indx]
            subcluster.gyrlev=cluster.gyrlev[indx]
        elif cluster.give =='mxve':
            subcluster.give=cluster.give
            subcluster.eps=cluster.eps[indx]


        if reset_centre:
            subcluster.add_orbit(
                cluster.xgc,
                cluster.ygc,
                cluster.zgc,
                cluster.vxgc,
                cluster.vygc,
                cluster.vzgc,
            )

            if cluster.origin=='centre' or cluster.origin=='cluster':
                subcluster.find_centre(0.0, 0.0, 0.0, reset_centre=reset_centre)
            else:
                subcluster.find_centre(reset_centre=reset_centre)

        else:
            subcluster.add_orbit(
                cluster.xgc,
                cluster.ygc,
                cluster.zgc,
                cluster.vxgc,
                cluster.vygc,
                cluster.vzgc,
            )
            subcluster.xc, subcluster.yc, subcluster.zc = (
                cluster.xc,
                cluster.yc,
                cluster.zc,
            )
            subcluster.vxc, subcluster.vyc, subcluster.vzc = (
                cluster.vxc,
                cluster.vyc,
                cluster.vzc,
            )

            subcluster.ra_gc, subcluster.dec_gc, subcluster.dist_gc = cluster.ra_gc, cluster.dec_gc, cluster.dist_gc
            subcluster.pmra_gc, subcluster.pmdec_gc, subcluster.vlos_gc = (
                cluster.pmra_gc,
                cluster.pmdec_gc,
                cluster.vlos_gc,
            )

        if reset_nbody:
            subcluster.to_pckms()
            subcluster.analyze()
            subcluster.reset_nbody_scale(mass=True,radius=True,rvirial=reset_rvirial,projected=reset_projected,**kwargs)
        elif reset_nbody_mass or reset_nbody_radii:
            subcluster.to_pckms()
            subcluster.analyze()
            subcluster.reset_nbody_scale(mass=reset_nbody_mass,radius=reset_nbody_radii,rvirial=reset_rvirial,projected=reset_projected,**kwargs)

    else:
        subcluster = StarCluster(cluster.tphys)

    if subcluster.ntot > 0:
        if subcluster.units!=units0: subcluster.to_units(units0)
        if subcluster.origin!=origin0: subcluster.to_origin(origin0)
        subcluster.analyze(sortstars=sortstars)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return subcluster

def overlap_cluster(cluster1,cluster2,tol=0.1,projected=False,return_cluster=True):
    """Extract a sub population of stars from cluster1 that spatially overlaps with cluster2

    Parameters
    ----------
    cluster1 : StarCluster
        cluster from which stars are to be extracted
    cluster2 : StarCluster
        cluster from which overlap region is defined
    tol: float
        tolerance parameter for how much regions need to overlap (default: 0.1)
    projected : bool 
        use projected values and constraints (default:False)
    return_cluster: bool
        returns a sub_cluster if True, otherwise return the boolean array (default:True)

    Returns
    -------
    StarCluster

    History
    -------
    2021 - Written - Webb (UofT)

    """
    indx=np.zeros(cluster1.ntot,dtype=bool)
    drmin=np.zeros(cluster1.ntot)

    for i in range(0,cluster1.ntot):
        dx=cluster1.x[i]-cluster2.x
        dy=cluster1.y[i]-cluster2.y

        if not projected:
            dz=cluster1.z[i]-cluster2.z
            dr=np.sqrt(dx**2.+dy**2.+dz**2.)
        else:
            dr=np.sqrt(dx**2.+dy**2.)

        drmin[i]=np.amin(dr)
        if np.amin(dr) < tol:
            indx[i]=True

    if return_cluster:
        return sub_cluster(cluster1,indx=indx)
    else:
        return indx







