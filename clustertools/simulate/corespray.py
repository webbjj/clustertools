from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,KingPotential,MovingObjectPotential
from galpy.util import bovy_conversion,conversion,bovy_plot
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import animation

from ..util.recipes import *
from ..util.coordinates import *
from ..analysis.cluster import StarCluster

class corespray(object):

	""" A class that is responsible for periodcially ejecting stars from a cluster's core

	Parameters
	----------


	History
	-------
	2021 - Written - Grandin/Webb (UofT)

	"""

	def __init__(self,gcorbit,pot=MWPotential2014,mu0=0.,sig0=10.0,vesc0=10.0,rho0=1.,mgc=None,rgc=None,W0=None,mmin=0.1,mmax=1.4,alpha=-1.35,ro=8.,vo=220.,q=-2):
		
		#gc orbit - name of or orbit instance corresponding to GC
		#mu0 - average 1D velocity in core (km/s)
		#sig0 - 1D velocity dispersion in core (km/s)
		#vesc0 - core velocity dispersion (km/s)
		#rho0 - core density (Msun/pc^3)
		#mgc - mass of GC
		#rgc - size of gc (assumed to be rh if Plummer or rt if King)
		#W0 - central potential if King model
		#mmin,mmax - minimum and maximum mass in the core (MSun)
		#alpha - slope of core mass function

		if isinstance(gcorbit,str):
			self.gcname=gcorbit
			self.o = Orbit.from_name(self.gcname, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])
		else:
			self.gcname='unknown'
			self.o=gcorbit

		self.mu0,self.sig0,self.vesc0,self.rho0=mu0,sig0,vesc0,rho0

		self.mmin,self.mmax,self.alpha=mmin,mmax,alpha

		#Mean separation of star's in the core equal to twice the radius of a sphere that contains one star
		#Assume all stars in the core have mass equal to the mean mass
		masses=power_law_distribution_function(1000, self.alpha, self.mmin, self.mmax)
		self.mbar=np.mean(masses)
		self.rsep=2.*((self.mbar/self.rho0)/(4.*np.pi/3.))**(1./3.)


		self.ro,self.vo=ro,vo
		self.to=conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.
		self.mo=bovy_conversion.mass_in_msol(ro=self.ro,vo=self.vo)


		self.q=q

		self.mwpot=pot

		if mgc is None:
			self.gcpot=None
		else:
			if W0 is None:
				ra=rgc/1.3
				self.gcpot=PlummerPotential(mgc/self.mo,ra/self.ro,ro=self.ro,vo=self.vo)
			else:
				self.gcpot=KingPotential(W0,mgc/self.mo,rgc/self.ro,ro=self.ro,vo=self.vo)

	def sample(self,tdisrupt=1000.,rate=1.,nstar=None,npeak=5.,verbose=False):
		#tdisrupt - time over which sampling begins (Myr)
		#rate - ejection rate (default 1 per Myr)
		#nstar - if set, nstar stars will be ejected randomly between start and tend
		#npeak - when sampling vs, sampling range will be from 0 to npeak*vpeak, where vpeak is the peak in the distribution 

		grav=4.302e-3 #pc/Msun (km/s)^2

		self.tdisrupt=tdisrupt

		#Select escape times
		#If nstar is not None, randomly select escapers between tstart and tend
		if nstar is not None:
			self.nstar=nstar
			self.rate=nstar/self.tdisrupt
		else:
			self.rate=rate
			self.nstar=self.tdisrupt*rate
			
		self.tesc=-1.*self.tdisrupt*np.random.rand(self.nstar)

		ts=np.linspace(0.,-1.*self.tdisrupt/self.to,1000)
		self.o.integrate(ts,self.mwpot)

		moving_pot=MovingObjectPotential(self.o,self.gcpot,ro=self.ro,vo=self.vo)
		self.pot=[self.mwpot,moving_pot]
	    
	    #Generate positions and new velocities for escaped stars with velocity dispersions of 100/root(3), based on
	    #gcname's position 1 Gyr ago

		vxesc=np.zeros(self.nstar)
		vyesc=np.zeros(self.nstar)
		vzesc=np.zeros(self.nstar)
	    
		self.vesc=np.array([])

		nescape=0


		self.mstar=np.zeros(self.nstar)
		self.mb1=np.zeros(self.nstar)
		self.mb2=np.zeros(self.nstar)
		self.eb=np.zeros(self.nstar)

		while nescape < self.nstar:
			masses=power_law_distribution_function(3, self.alpha, self.mmin, self.mmax)

			mindx=(masses==np.amin(masses))
			ms=masses[mindx][0]
			m_a,m_b=masses[np.invert(mindx)]
			        
			mb=m_a+m_b
			M=ms+mb

			prob=self.prob_three_body_escape(ms,m_a,m_b,self.q)

			if np.random.rand() < prob:
			    
				vxs,vys,vzs=np.random.normal(self.mu0,self.sig0,3)
				vstar=np.sqrt(vxs**2.+vys**2.+vzs**2.)
				vxb,vyb,vzb=np.random.normal(self.mu0,self.sig0,3)
				rdot=np.sqrt((vxs-vxb)**2.+(vys-vyb)**2.+(vzs-vzb)**2.)

				ebin,semi=self.sample_binding_energy(m_a,m_b,alpha=-1,xmin=10.0**36,xmax=10.0**40.0)

				e0=0.5*(mb*ms/M)*(rdot**2.)-grav*ms*mb/self.rsep + ebin

				vs=self.sample_escape_velocity(e0,ms,mb,npeak)

				if vs >self.vesc0:  
				    self.vesc=np.append(self.vesc,vs)

				    vxesc[nescape]=vs*(vxs/vstar)
				    vyesc[nescape]=vs*(vys/vstar)
				    vzesc[nescape]=vs*(vzs/vstar)

				    self.mstar[nescape]=ms
				    self.mb1[nescape]=m_a
				    self.mb2[nescape]=m_b
				    self.eb[nescape]=ebin

				    nescape+=1

				if verbose: print('DEBUG: ',nescape,prob,vs,self.vesc0)
		

		Re0, phie0, ze0, vRe0, vTe0, vze0=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])

		self.o0=[]

		for i in range(0,self.nstar):
			xe,ye,ze=self.o.x(self.tesc[i]/self.to),self.o.y(self.tesc[i]/self.to),self.o.z(self.tesc[i]/self.to)
			vxe=vxesc[i]+self.o.vx(self.tesc[i]/self.to)
			vye=vyesc[i]+self.o.vy(self.tesc[i]/self.to)
			vze=vzesc[i]+self.o.vy(self.tesc[i]/self.to)

			if verbose: print(i,self.tesc[i],xe,ye,ze,vxe,vye,vze)

			Re, phie, ze, vRe, vTe, vze=cart_to_cyl(xe,ye,ze,vxe,vye,vze)
			oe=Orbit([Re/self.ro, vRe/self.vo, vTe/self.vo, ze/self.ro, vze/self.vo, phie],ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

			self.o0.append(oe)

			ts=np.linspace(self.tesc[i]/self.to,0.,1000)
			oe.integrate(ts,self.pot)


			Re0=np.append(Re0,oe.R(0.))
			phie0=np.append(phie0,oe.phi(0.))
			ze0=np.append(ze0,oe.z(0.))
			vRe0=np.append(vRe0,oe.vR(0.))
			vTe0=np.append(vTe0,oe.vT(0.))
			vze0=np.append(vze0,oe.vz(0.))


		vxvve=np.column_stack([Re0/self.ro,vRe0/self.vo,vTe0/self.vo,ze0/self.ro,vze0/self.vo,phie0])
		self.oe=Orbit(vxvve,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
		

		cluster=StarCluster(units='kpckms',origin='galaxy')
		cluster.add_stars(self.oe.x(),self.oe.y(),self.oe.z(),self.oe.vx(),self.oe.vy(),self.oe.vz(),m=self.mstar)
		cluster.add_orbit(self.o.x(),self.o.y(),self.o.z(),self.o.vx(),self.o.vy(),self.o.vz())

		return cluster

	def prob_three_body_escape(self,ms,m_a,m_b,q):
		#Equation 7.23
		prob=(ms**q)/(ms**q+m_a**q+m_b**q)
		return prob

	def sample_binding_energy(self,mb1,mb2,alpha=-1,xmin=10.0**36,xmax=10.0**40.0):
	    #Opik's Law
		#Default binding energy distribution is:
		# power law of slope -1 
		# Between 10.0**36 and 10.0**40 J

		grav=4.302e-3 #pc/Msun (km/s)^2


		if isinstance(mb1,float):
			n=1
		else:
			n=len(mb1)

		ebin_si=self.xfunc(n,alpha,xmin,xmax) #Joules = kg (m/s)^2
		ebin=ebin_si/1.9891e30 # Msun (m/s)^2
		ebin/=(1000.0*1000.0) #Msun (km/s)^2


		#Log normal a:
		semi_pc=(0.5*grav*mb1*mb2)/ebin
		semi_au=semi_pc/4.84814e-6    

		semi=semi_pc

		return ebin,semi

	def sample_escape_velocity(self,e0,ms,mb,npeak=5):
		#randomly sample between npeak*vs_peak

		vs_peak=self.escape_velocity_distribution_peak(e0,ms,mb)
		match=False

		while not match:
			vstemp=np.random.rand()*npeak*vs_peak
			amptemp=np.random.rand()*vs_peak

			if amptemp < self.escape_velocity_distribution(vstemp,e0,ms,mb):
				vs=vstemp
				match=True

		return vs


	def escape_velocity_distribution(self,vs,e0,ms,mb):
		#Equation 7.19
		M=ms+mb
		fv=(3.5*(np.fabs(e0)**(7./2.))*ms*M/mb)*vs/((np.fabs(e0)+0.5*(ms*M/mb)*(vs**2.))**(9./2.))
		return fv

	def escape_velocity_distribution_peak(self,e0,ms,mb):
		M=ms+mb
		vs_peak=0.5*np.sqrt((M-ms)/(ms*M))*np.sqrt(np.fabs(e0))

		return vs_peak

	def xfunc(self,n,alpha,xmin,xmax):

	    eta=alpha+1.
	    
	    if xmin==xmax:
	        x=xmin
	    elif alpha==0:
	        x=xmin+np.random.random(n)*(xmax-xmin)
	    elif alpha>0:
	        x=xmin+np.random.power(eta,n)*(xmax-xmin)
	    elif alpha<0 and alpha!=-1.:
	        x=(xmin**eta + (xmax**eta - xmin**eta)*np.random.rand(n))**(1./eta)
	    elif alpha==-1:
	        x=np.log10(xmin)+np.random.random(n)*(np.log10(xmax)-np.log10(xmin))
	        x=10.0**x
	        
	    if n==1:
	        return x
	    else:      
	        return np.array(x)

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
		tsint=np.linspace(0.,-1.*self.tdisrupt/self.to,1000)
		self.oe.integrate(tsint,self.pot)

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