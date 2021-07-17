from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,KingPotential,MovingObjectPotential
from galpy.util import bovy_conversion,conversion,bovy_plot
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation

from ..util.recipes import *
from ..util.coordinates import *
from ..analysis.cluster import StarCluster

#from streamtools.df import streamspraydf
import streamspraydf_Webb


class tailspray(object):

	""" A class that is responsible for periodcially ejecting stars from a cluster via tidal stripping
	- Note this code is effectively a wrapper around Jo Bovy's streamtools python package that puts the 
	- output into a StarCluster

	Parameters
	----------


	History
	-------
	2021 - Written - El-Falou/Webb (UofT)

	"""

	def __init__(self,gcname,mgc=None,rgc=None,pot=MWPotential2014,rtpot=MWPotential2014,ro=8.,vo=220.):
		
		#mu0 - average 1D velocity in core (km/s)
		#sig0 - 1D velocity dispersion in core (km/s)
		#vesc0 - core velocity dispersion (km/s)
		#rho0 - core density (Msun/pc^3)
		#mgc - mass of GC
		#rgc - size of gc (assumed to be rh if Plummer or rt if King)
		#W0 - central potential if King model
		#mmin,mmax - minimum and maximum mass in the core (MSun)
		#alpha - slope of core mass function

		self.gcname=gcname
		self.init_orbit(self.gcname)

		self.mgc=mgc
		self.rgc=rgc

		self.ro,self.vo=ro,vo
		self.to=conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.
		self.mo=conversion.mass_in_msol(ro=self.ro,vo=self.vo)

		self.pot=pot
		self.rtpot=rtpot

	def sample(self,tstart= -1000., tend=0., nstar=100, nper=None, integ=True, verbose=False):
		"""
		Makes stream using streamspraydf

		Returns globular cluster orbit integratrated forward to tdisrupt, 
		globular cluster orbit integratrated backward to -tdisrupt, 6D numpy 
		array of Orbit's vxvv of each star in leading stream at t=0, array 
		of Orbit's vxvv of each star in trailing stream at t=0, array of 
		times when stars in leading stream escaped, and array of times 
		when stars in trailing stream escaped.

		If not integ, RvR and RvRt are Orbit's vxvv at t=-tdisrupt
		"""

		if nper is not None:

			ts=np.linspace(0.,1.,1000)
			self.o.integrate(ts,self.pot)

			if nper < 0:
				self.tstart=nper*self.o.Tp()*1000.0
				self.tend=tend
			else:
				self.tstart=tstart
				self.tend=nper*self.o.Tp()*1000.0

			self.tdisrupt=self.tend-self.tstart

			self.init_orbit(self.gcname)

		else:
			self.tstart=tstart
			self.tend=tend
			self.tdisrupt=tend-tstart

		spdf= streamspraydf_Webb.streamspraydf(self.mgc,
                           progenitor=self.o,
                           pot=self.pot,
                           tdisrupt=self.tdisrupt, 
                           rtpot=self.rtpot)
		spdft= streamspraydf_Webb.streamspraydf(self.mgc,
                           progenitor=self.o,
                           pot=self.pot,
                           tdisrupt=self.tdisrupt, 
                           rtpot=self.rtpot,
                           leading=False)

		RvR,dt= spdf.sample(n=nstar,returndt=True,integrate=integ)
		RvRt,dtt= spdft.sample(n=nstar,returndt=True,integrate=integ)

		#vxvv=np.column_stack([RvR[0]/self.ro,RvR[1]/self.vo,RvR[2]/self.vo,RvR[3]/self.ro,RvR[4]/self.vo,RvR[5]])
		#vxvvt=np.column_stack([RvRt[0]/self.ro,RvRt[1]/self.vo,RvRt[2]/self.vo,RvRt[3]/self.ro,RvRt[4]/self.vo,RvRt[5]])
		vxvv=np.column_stack([RvR[0],RvR[1],RvR[2],RvR[3],RvR[4],RvR[5]])
		vxvvt=np.column_stack([RvRt[0],RvRt[1],RvRt[2],RvRt[3],RvRt[4],RvRt[5]])


		osleading=Orbit(vxvv,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
		ostrailing=Orbit(vxvvt,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

		self.init_orbit(self.gcname,use_physical=True)

		cluster=StarCluster(units='kpckms',origin='galaxy')
		cluster.add_stars(osleading.x(),osleading.y(),osleading.z(),osleading.vx(),osleading.vy(),osleading.vz())
		cluster.add_stars(ostrailing.x(),ostrailing.y(),ostrailing.z(),ostrailing.vx(),ostrailing.vy(),ostrailing.vz())
		cluster.add_orbit(self.o.x(),self.o.y(),self.o.z(),self.o.vx(),self.o.vy(),self.o.vz())

		self.dt=np.append(dt,dtt)

		return cluster,RvR,RvRt

	def init_orbit(self,gcname,use_physical=False):

		if use_physical:
			self.o=Orbit.from_name(gcname,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
		else:
			temp = Orbit.from_name(gcname)
			self.o = Orbit([temp.R(use_physical=False), 
			       temp.vR(use_physical=False), 
			       temp.vT(use_physical=False), 
			       temp.z(use_physical=False), 
			       temp.vz(use_physical=False), 
			       temp.phi()])

	def _init_fig(self,xlim=(-20,20),ylim=(-20,20)):
	    self.fig = plt.figure()
	    self.ax = plt.axes(xlim=xlim, ylim=ylim)
	    self.ax.set_xlabel('X (kpc)')
	    self.ax.set_ylabel('Y (kpc)')
	    self.txt_title=self.ax.set_title('')
	    self.line, = self.ax.plot([], [], lw=2)
	    self.pt, = self.ax.plot([],[],'.')

	def _set_data(self,gcdata,sdata):
	    self.gcdata = gcdata
	    self.sdata=sdata

	def _ani_init(self):
	    self.line.set_data([], [])
	    self.pt.set_data([],[])
	    return self.line,self.pt

	def _ani_update(self, i):

		if i < 5:
		    x = self.gcdata[0:i,0]
		    y = self.gcdata[0:i,1]
		else:
			x = self.gcdata[i-5:i,0]
			y = self.gcdata[i-5:i,1]	
		self.line.set_data(x, y)

		escindx=self.tesc/self.to <= self.ts[i]


		if np.sum(escindx)>0:
			self.pt.set_data(self.sdata[i][0][escindx],self.sdata[i][1][escindx])
		else:
			self.pt.set_data([],[])

		self.txt_title.set_text('%s' % str (self.ts[i]*self.to))

		return self.line,self.pt


	def animate(self,frames=100,interval=50,xlim=(-20,20),ylim=(-20,20)):

		self._init_fig(xlim,ylim)

		self.ts=np.linspace(self.tstart/self.to,self.tend/self.to,frames)
		self.oe.integrate(np.flip(self.ts),self.pot)

		gcdata=np.zeros(shape=(frames,2))

		for i in range(0,frames):
		    gcdata[i]=[self.o.x(self.ts[i]),self.o.y(self.ts[i])]

		sdata=np.zeros(shape=(frames,2,self.nstar))

		for i in range(0,frames):
			sdata[i]=[self.oe.x(self.ts[i]),self.oe.y(self.ts[i])]

		self._set_data(gcdata,sdata)

		self.anim = animation.FuncAnimation(self.fig, self._ani_update, init_func=self._ani_init, frames=frames, interval=interval, blit=False)

	def snapout(self,filename='corespray.dat'):
		R=np.append(self.o.R(0.),self.oe.R(0.))
		vR=np.append(self.o.vR(0.),self.oe.vR(0.))
		vT=np.append(self.o.vT(0.),self.oe.vT(0.))
		z=np.append(self.o.z(0.),self.oe.z(0.))
		vz=np.append(self.o.vz(0.),self.oe.vz(0.))
		phi=np.append(self.o.phi(0.),self.oe.phi(0.))

		vesc=np.append(0.,self.vesc)
		tesc=np.append(0.,self.tesc)

		np.savetxt(filename,np.column_stack([R,vR,vT,z,vz,phi,vesc,tesc]))