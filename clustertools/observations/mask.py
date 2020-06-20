from ..analysis.cluster import sub_cluster
from ..util.recipes import *
import numpy as np
from scipy.spatial import ConvexHull

# Place a mask over the simulated data in order to replicate an observational dataset


class ObservationalMask(object):
    def __init__(self, name):

        self.name = name
        self.rnorm=False

        self.rmlower=None
        self.rmupper=None

    def set_orbit(self,ra,dec,dist,pmra,pmdec,vlos):

        self.ra = ra
        self.dec = dec
        self.dist = dist
        self.pmra = pmra
        self.pmdec = pmdec
        self.vlos = vlos

    def set_key_params(self,**kwargs):

        if 'rm' in kwargs:
            self.rm = kwargs.get("rm")
        if 'rt' in kwargs:
            self.rt = kwargs.get("rt")
        if 'rmlower' in kwargs:
            self.rmlower = kwargs.get("rmlower")
        if 'rmupper' in kwargs:
            self.rmupper = kwargs.get("rmupper")
 
    def set_mcorr(self, mcorr):
        # Se mass correction for entire dataset
        self.mcorr = mcorr

    def set_rlims(self,rmin,rmax):
        self.rmin=rmin
        self.rmax=rmax

    def set_rbins(self,r_lower,r_mean,r_upper):
        self.r_lower=r_lower
        self.r_mean=r_mean
        self.r_upper=r_upper

    def calc_rbins(self,cluster,nsplit,nrbin,rmin,rmax,indx=None,projected=True):

        if indx is None:
            indx=np.ones(cluster.ntot,bool)

        if projected:
            r=cluster.rpro
        else:
            r=cluster.r

        if nsplit==0:
            if not isinstance(rmin,float): rmin=rmin[0]
            if not isinstance(rmax,float): rmax=rmax[0]
            if not isinstance(nrbin,int): rmax=int(nrbin[0])

            rindx=(r>=rmin) * (r<=rmax) * indx

            self.r_lower,self.r_mean,self.r_upper,r_hist=npy.nbinmaker(r[rindx],nrbin)
        else:
            self.r_lower=np.array([])
            self.r_mean=np.array([])
            self.r_upper=np.array([])

            for i in range(0,nsplit):
                rindx=(r>=rmin[i]) * (r<=rmax[i])
                r_lower,r_mean,r_upper,r_hist=nbinmaker(r[rindx],nrbin[i])
                self.r_lower=np.append(self.r_lower,r_lower[:])
                self.r_mean=np.append(self.r_mean,r_mean[:])
                self.r_upper=np.append(self.r_upper,r_upper[:])

    def set_mbins(self,m_lower,m_mean,m_upper):
        self.m_lower=m_lower
        self.m_mean=m_mean
        self.m_upper=m_upper

    def calc_mbins(self,cluster,mmin=0.1,mmax=0.8):
        mindx=(cluster.m >= mmin) * (cluster.m<=mmax)
        self.m_lower,self.m_mean,self.m_upper=nbinmaker(cluster.m[mindx])

    def apply_orbit(self,cluster):
        # Apply orbit to self.cluster based on mask
        cluster.add_orbit(
            self.ra,
            self.dec,
            self.dist,
            self.pmra,
            self.pmdec,
            self.vlos,
        )

    def apply_key_params(self,cluster):
        # Force key parameters to be specific values - useful if only part of the data is provided
        cluster.rm = self.rm
        cluster.rt = self.rt

    def apply_rnorm(self):
        self.rnorm=True

        self.r_lower/=self.rm
        self.r_mean/=self.rm
        self.r_upper/=self.rm

    def undo_rnorm(self):
        self.rnorm=False

        self.r_lower*=self.rm
        self.r_mean*=self.rm
        self.r_upper*=self.rm
