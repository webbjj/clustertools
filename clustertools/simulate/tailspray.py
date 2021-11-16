""" The tailspray d class
  
"""

__author__ = "Nada El-Falou & Jeremy J Webb"

__all__ = [
    "tailspray",
]

from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,KingPotential,MovingObjectPotential
from galpy.util import bovy_conversion,conversion,bovy_plot
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation

from ..util.recipes import *
from ..util.coordinates import *
from ..cluster.cluster import StarCluster
from ..tidaltail.tails import *
from ..analysis.orbits import orbital_path

from streamtools.df import streamspraydf

from astropy import units

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

	def __init__(self,gcorbit,mgc=None,pot=MWPotential2014,rtpot=MWPotential2014,ro=8.,vo=220.):

		# - name of or orbit instance corresponding to GC
		#mgc - mass of GC
		#pot - potential for integrating cluster
		#rtpot - potential for calcuating the tidal radius


		self.ro,self.vo=ro,vo
		self.to=conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.
		self.mo=conversion.mass_in_msol(ro=self.ro,vo=self.vo)

		if isinstance(gcorbit,str):
			self.gcname=gcorbit
			self.o=Orbit.from_name(self.gcname, ro=self.ro, vo=self.vo, solarmotion=[-11.1, 24.0, 7.25])
		else:
			self.gcname='unknown'
			self.o=gcorbit
		
		self.o.turn_physical_off()

		self.mgc=mgc

		self.pot=pot
		self.rtpot=rtpot

	def sample(self,tdisrupt,nstar=100,integ=True, verbos=False):
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

		self.tdisrupt=tdisrupt
		self.nstar=nstar

		spdf= streamspraydf(self.mgc/self.mo,
                           progenitor=self.o,
                           pot=self.pot,
                           tdisrupt=self.tdisrupt/self.to, 
                           rtpot=self.rtpot)
		spdft= streamspraydf(self.mgc/self.mo,
                           progenitor=self.o,
                           pot=self.pot,
                           tdisrupt=self.tdisrupt/self.to, 
                           rtpot=self.rtpot,
                           leading=False)

		RvR,dt= spdf.sample(n=self.nstar,returndt=True,integrate=integ)
		RvRt,dtt= spdft.sample(n=self.nstar,returndt=True,integrate=integ)

		vxvv=np.column_stack([RvR[0],RvR[1],RvR[2],RvR[3],RvR[4],RvR[5]])
		vxvvt=np.column_stack([RvRt[0],RvRt[1],RvRt[2],RvRt[3],RvRt[4],RvRt[5]])

		ol=Orbit(vxvv,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
		ot=Orbit(vxvvt,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])


		self.cluster=StarCluster(origin='galaxy',units='kpckms')
		self.cluster.add_stars(ol.x(),ol.y(),ol.z(),ol.vx(),ol.vy(),ol.vz())
		self.cluster.add_stars(ot.x(),ot.y(),ot.z(),ot.vx(),ot.vy(),ot.vz())
		self.cluster.add_orbit(self.o.x()*self.ro,self.o.y()*self.ro,self.o.z()*self.ro,self.o.vx()*self.vo,self.o.vy()*self.vo,self.o.vz()*self.vo)

		self.tesc=-1.*np.append(dt*self.to,dtt*self.to)

		return self.cluster

	def orbital_path(self, 
		dt=0.1, 
		nt=100, 
		pot=None, 
		from_centre=False, 
		skypath=False, 
		initialize=False,
		ro=None, 
		vo=None,
		plot=False,
		**kwargs):

			if pot is None: pot=self.pot
			if ro is None: ro=self.ro
			if vo is None: vo=self.vo

			self.torbit, self.xorbit, self.yorbit, self.zorbit, self.vxorbit, self.vyorbit, self.vzorbit=orbital_path(self.cluster,dt=dt,nt=nt,pot=pot,from_centre=from_centre,skypath=skypath,initialize=initialize,ro=ro,vo=vo,plot=plot,**kwargs)

	def tail_path(self, 
		dt=0.1, 
		nt=100, 
		pot=None, 
		from_centre=False, 
		skypath=False, 
		ro=None, 
		vo=None,
		plot=False,
		**kwargs):

			if pot is None: pot=self.pot
			if ro is None: ro=self.ro
			if vo is None: vo=self.vo

			self.ttail, self.xtail, self.ytail, self.ztail, self.vxtail, self.vytail, self.vztail=tail_path(self.cluster,dt=dt,nt=nt,pot=pot,from_centre=from_centre,skypath=skypath,ro=ro,vo=vo,plot=plot,**kwargs)

	def tail_path_match(self,
	    dt=0.1,
	    nt=100,
	    pot=None,
	    path=None,
	    from_centre=False,
	    skypath=False,
	    to_path=False,
	    do_full=False,
	    ro=None,
	    vo=None,
	    plot=False,
	    projected=False,
	    **kwargs,
		):

			if pot is None: pot=self.pot
			if ro is None: ro=self.ro
			if vo is None: vo=self.vo

			self.tstar,self.dprog,self.dpath=tail_path_match(self.cluster,dt=dt,nt=nt,pot=pot,path=path,from_centre=from_centre,skypath=skypath,to_path=to_path,do_full=do_full,ro=ro,vo=vo,plot=plot,projected=projected,**kwargs)

	def _init_fig(self,xlim=(-20,20),ylim=(-20,20)):
	    self.fig = plt.figure()
	    self.ax = plt.axes(xlim=xlim, ylim=ylim)
	    self.ax.set_xlabel('X (kpc)')
	    self.ax.set_ylabel('Y (kpc)')
	    self.txt_title=self.ax.set_title('')
	    self.pt, = self.ax.plot([],[],'.')
	    self.line, = self.ax.plot([], [], lw=2)

	def _set_data(self,gcdata,sdata):
	    self.gcdata = gcdata
	    self.sdata=sdata

	def _ani_init(self):
	    self.line.set_data([], [])
	    self.pt.set_data([],[])
	    return self.line,self.pt

	def _ani_update(self, i):

		if i < 5:
		    x = self.gcdata[0:i+1,0]
		    y = self.gcdata[0:i+1,1]
		else:
			x = self.gcdata[i-5:i+1,0]
			y = self.gcdata[i-5:i+1,1]

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

		self.ts=np.linspace(-1.*self.tdisrupt/self.to,0.,frames)

		tsint=np.linspace(0,-1.*self.tdisrupt/self.to,1000)

		self.o.integrate(tsint,self.pot)

		gcdata=np.zeros(shape=(frames,2))

		for i in range(0,frames):
		    gcdata[i]=[self.o.x(self.ts[i])*self.ro,self.o.y(self.ts[i])*self.ro]

		sdata=np.zeros(shape=(frames,2,int(2*self.nstar)))

		self.oe=self.cluster.initialize_orbits()
		self.oe.integrate(tsint,self.pot)

		for i in range(0,frames):
			sdata[i]=[self.oe.x(self.ts[i]),self.oe.y(self.ts[i])]

		self._set_data(gcdata,sdata)

		self.anim = animation.FuncAnimation(self.fig, self._ani_update, init_func=self._ani_init, frames=frames, interval=interval, blit=False)

	def snapout(self,filename='corespray.dat'):
		#Need positions, velocities, and escape times!
		self.oe=self.cluster.initialize_orbits()

		R=np.append(self.o.R(0.),self.oe.R(0.))
		vR=np.append(self.o.vR(0.),self.oe.vR(0.))
		vT=np.append(self.o.vT(0.),self.oe.vT(0.))
		z=np.append(self.o.z(0.),self.oe.z(0.))
		vz=np.append(self.o.vz(0.),self.oe.vz(0.))
		phi=np.append(self.o.phi(0.),self.oe.phi(0.))

		tesc=np.append(0.,self.tesc)

		np.savetxt(filename,np.column_stack([R,vR,vT,z,vz,phi,tesc]))
