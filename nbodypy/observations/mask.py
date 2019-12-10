from ..main.cluster import sub_cluster

import numpy as np
from scipy.spatial import ConvexHull

# Place a mask over the simulated data in order to replicate an observational dataset


class ObservationalMask(object):
    def __init__(self, name):

        self.name = name
        self.rnorm=False

    def set_orbit(self,ra,dec,dist,pmra,pmdec,vlos):

        self.ra = ra
        self.dec = dec
        self.dist = dist
        self.pmra = pmra
        self.pmdec = pmdec
        self.vlos = vlos

    def set_key_params(self,**kwargs):

        self.rm = kwargs.get("rm", 0.)
        self.rt = kwargs.get("rt", 0.)

    def set_mcorr(self, mcorr):
        # Se mass correction for entire dataset
        self.mcorr = mcorr

    def set_cluster(self, cluster):
        # Set base cluster to use values of rm for
        self.cluster = cluster

    def set_rlims(self,rmin,rmax):
        self.rmin=rmin
        self.rmax=rmax

    def set_rbins(self,r_lower,r_mean,r_upper):
        self.r_lower=r_lower
        self.r_mean=r_mean
        self.r_upper=r_upper

    def set_mbins(self,m_lower,m_mean,m_upper):
        self.m_lower=m_lower
        self.m_mean=m_mean
        self.m_upper=m_upper

    def apply_mask(self):
        # Mask simulated data to reflect fields of view of observational

        self.apply_orbit()
        self.apply_key_params()

    def apply_orbit(self):
        # Apply orbit to self.cluster based on mask
        self.cluster.add_orbit(
            self.ra,
            self.dec,
            self.dist,
            self.pmra,
            self.pmdec,
            self.vlos,
        )

    def apply_key_params(self):
        # Force key parameters to be specific values - useful if only part of the data is provided
        self.cluster.rm = self.rm
        self.cluster.rt = self.rt

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
