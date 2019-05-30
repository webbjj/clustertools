import numpy as np
from galpy.util import bovy_conversion,bovy_coords
from textwrap import dedent
from galpy.potential import MWPotential2014

from ..nbodypy.orbit import rtidal,rlimiting,initialize_orbit
from ..nbodypy.functions import *
from ..nbodypy.profiles import *

class StarClusterSimulation(object):
    """
    Class representing a full star cluster simulation
    """    
    def __init__(self,cluster=None,exfilename=None,t=None,m=None,n=None,rm=None,rc=None):
        """
        NAME:

           __init__

        PURPOSE:

           Initialize a Star Cluster Simulation, which keeps track of the time evolution of key cluster parameters

        INPUT:

           cluster - initial cluster instance (optional)
           
           exfilename - if files are to be read from a saved file (format set to extrcy.npy in custom_output.npy) (optional)

           t,m,n,rm,rc - set time,mass,n,half-mass radius,core radius manually (optional)

        OUTPUT:

           instance

        HISTORY:

           2019 - Written - Webb (UofT)

        """
        if exfilename!=None:
            self.load_extrct(exfilename)
        elif cluster!=None:
            self.cluster=np.array([cluster])
            self.t=np.array([cluster.tphys])
            self.m=np.array([cluster.mtot])
            self.n=np.array([cluster.ntot])
            self.rm=np.array([cluster.rm])
            self.rc=np.array([cluster.r10])
        elif isinstance(t,float):
            self.cluster=None
            self.t=np.array([t])
            self.m=np.array([m])
            self.n=np.array([n])
            self.rm=np.array([rm])
            self.rc=np.array([rc])
        else:
            self.cluster=None
            self.t=np.array(t)
            self.m=np.array(m)
            self.n=np.array(n)
            self.rm=np.array(rm)
            self.rc=np.array(rc)

    def add_snapshot(self,cluster=None,t=None,m=None,n=None,rm=None,rc=None):
        if cluster!=None:
            self.cluster=np.append(self.cluster,cluster)
            self.t=np.append(self.t,cluster.tphys)
            self.m=np.append(self.m,cluster.mtot)
            self.n=np.append(self.n,cluster.ntot)
            self.rm=np.append(self.rm,cluster.rm)
            self.rc=np.append(self.rc,cluster.r10)

        else:
            self.cluster=np.append(self.cluster,None)
            self.t=np.append(self.t,t)
            self.m=np.append(self.m,m)
            self.n=np.append(self.n,n)
            self.rm=np.append(self.rm,rm)
            self.rc=np.append(self.rc,rc)

    def load_extrct(self,exfilename='extrct.npy'):
        exdata=np.loadtxt(exfilename)
        self.n=exdata[:,0]        
        self.t=exdata[:,2]
        self.m=exdata[:,4]
        self.rm=exdata[:,9]
        self.rc=exdata[:,5]

    def load_alpha_prof(self,apfilename='alpha_prof.npy',nmass=10,nrad=20):
        data=np.loadtxt(apfilename)

        self.taprof=data[:,0]
        self.alpha=data[:,1]
        self.ealpha=data[:,2]
        self.yalpha=data[:,3]
        self.eyalpha=data[:,4]

        nstart,nstop=5,5+nmass
        self.mmean=data[:,nstart:nstop]
      
        nstart=nstop
        nstop=nstart+nmass
        self.dm=data[:,nstart:nstop]

        nstart=nstop
        nstop=nstart+nrad
        self.lrprof=data[:,nstart:nstop]

        nstart=nstop
        nstop=nstart+nrad
        self.aprof=data[:,nstart:nstop]

        nstart=nstop
        self.dalpha=data[:,nstart]
        self.edalpha=data[:,nstart+1]
        self.ydalpha=data[:,nstart+2]
        self.eydalpha=data[:,nstart+3]

#    def load_dalpha(cluster=None,dafilename='dalpha_prof.npy',nmass=3,**kwargs):
#        dafile=np.loadtxt(dafilename)
#        self.alpha=dafile(:,2:nmass)
#        self.dalpha=dafile(:,2+nmass-2+2*nmass)
