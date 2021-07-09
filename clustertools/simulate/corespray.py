from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.util import bovy_conversion,bovy_plot
import numpy as np

import matplotlib.pyplot as plt

from ..util.recipes import *
from ..util.coordinates import *

class corespray(object):

	""" A class that is responsible for periodcially ejecting stars from a cluster's core

	Parameters
	----------


	History
	-------
	2021 - Written - Grandin/Webb (UofT)

	"""

	def __init__(self,gcname,pot=MWPotential2014,mu0=0.,sig0=10.0,vesc0=10.0,rho0=1.,mmin=0.1,mmax=1.4,alpha=-1.35,ro=8.,vo=220.):
		
		#mu0 - average 1D velocity in core (km/s)
		#sig0 - 1D velocity dispersion in core (km/s)
		#vesc0 - core velocity dispersion (km/s)
		#rho0 - core density (Msun/pc^3)
		#mmin,mmax - minimum and maximum mass in the core (MSun)
		#alpha - slope of core mass function

		to=bovy_conversion.time_in_Gyr(ro=ro,vo=vo)*1000.

		self.gcname=gcname
		self.o = Orbit.from_name(gcname, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])		

		self.pot=pot

		self.mu0,self.sig0,self.vesc0,self.rho0=mu0,sig0,vesc0,rho0

		self.mmin,self.mmax,self.alpha=mmin,mmax,alpha

		#Mean separation of star's in the core equal to twice the radius of a sphere that contains one star
		#Assume all stars in the core have mass equal to the mean mass
		masses=power_law_distribution_function(1000, self.alpha, self.mmin, self.mmax)
		self.mbar=np.mean(masses)
		self.rsep=2.*((self.mbar/self.rho0)/(4.*np.pi/3.))**(1./3.)


		self.ro,self.vo=ro,vo


	def sample(self,tstart=-1000.,tend=0.,rate=1.,nstar=None,npeak=5.,verbose=False):
		#tstart - initial time when sampling begins (Myr)
		#tend - final time when sampling ends (Myr)
		#rate - ejection rate (default 1 per Myr)
		#nstar - if set, nstar stars will be ejected randomly between start and tend
		#npeak - when sampling vs, sampling range will be from 0 to npeak*vpeak, where vpeak is the peak in the distribution 

		grav=4.302e-3 #pc/Msun (km/s)^2

		to=bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.0

		self.tstart=tstart
		self.tend=tend
		self.dt=tend-tstart

		#Select escape times
		#If nstar is not None, randomly select escapers between tstart and tend
		if nstar is not None:
			self.nstar=nstar
			self.rate=self.dt/nstar
			self.tesc=tstart+np.random.rand(nstar)*self.dt
		else:
			self.rate=rate
			self.nstar=(tend-tstart)/rate
			self.tesc=tstart+np.random.rand(nstar)*(tend-tstart)

		if tstart<tend:
			ts=np.linspace(self.tend/to,self.tstart/to,1000)
			self.o.integrate(ts,self.pot)
		else:
			print('FORWARD INTEGRATION HAS NOT YET BEEN IMPLEMENTED')
			return -1
	    
	    #Generate positions and new velocities for escaped stars with velocity dispersions of 100/root(3), based on
	    #gcname's position 1 Gyr ago
		self.xe=self.o.x(self.tesc/to)
		self.ye=self.o.y(self.tesc/to)
		self.ze=self.o.z(self.tesc/to)

		self.vxe=self.o.vx(self.tesc/to)
		self.vye=self.o.vy(self.tesc/to)
		self.vze=self.o.vz(self.tesc/to)
	    
		self.vescape=np.array([])

		nescape=0

		while nescape < self.nstar:
			masses=power_law_distribution_function(3, self.alpha, self.mmin, self.mmax)

			mindx=(masses==np.amin(masses))
			ms=masses[mindx][0]
			m_a,m_b=masses[np.invert(mindx)]
			        
			mb=m_a+m_b
			M=ms+mb

			prob=self.prob_three_body_escape(ms,m_a,m_b)

			if np.random.rand() < prob:
			    
				vxs,vys,vzs=np.random.normal(self.mu0,self.sig0,3)
				vstar=np.sqrt(vxs**2.+vys**2.+vzs**2.)
				vxb,vyb,vzb=np.random.normal(self.mu0,self.sig0,3)
				rdot=np.sqrt((vxs-vxb)**2.+(vys-vyb)**2.+(vzs-vzb)**2.)

				ebin,semi=self.sample_binding_energy(m_a,m_b,alpha=-1,xmin=10.0**36,xmax=10.0**40.0)

				e0=0.5*(mb*ms/M)*(rdot**2.)-grav*ms*mb/self.rsep + ebin

				vs=self.sample_escape_velocity(e0,ms,mb,npeak)

				if vs >self.vesc0:  
				    self.vescape=np.append(self.vescape,vs)

				    self.vxe[nescape]+=vs*(vxs/vstar)
				    self.vye[nescape]+=vs*(vys/vstar)
				    self.vze[nescape]+=vs*(vzs/vstar)

				    nescape+=1

				if verbose: print('DEBUG: ',nescape,prob,vs,self.vesc0)
		

		#Need to do an orbit for each star
		Re, phie, ze, vRe, vTe, vze=cart_to_cyl(self.xe,self.ye,self.ze,self.vxe,self.vye,self.vze)

		#Make array of star orbits
		#Note positions and velocities are scaled by ro and vo respectively to be on Galpy units
		vxvv=np.column_stack([Re/self.ro,vRe/self.vo,vTe/self.vo,ze/self.ro,vze/self.vo,phie])
		os=Orbit(vxvv,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])

		#Integrate all stars forward by dt
		ts=np.linspace(0,self.dt/to,1000)
		os.integrate(ts,self.pot)

		#Get each stars orbit at time=tend
		Re0, phie0, ze0, vRe0, vTe0, vze0=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
		for i in range(0,len(os)):
			t0=(tend-self.tesc[i])/to

			Re0=np.append(Re0,os[i].R(t0))
			phie0=np.append(phie0,os[i].phi(t0))
			ze0=np.append(ze0,os[i].z(t0))
			vRe0=np.append(vRe0,os[i].vR(t0))
			vTe0=np.append(vTe0,os[i].vT(t0))
			vze0=np.append(vze0,os[i].vz(t0))

		vxvve=np.column_stack([Re0/self.ro,vRe0/self.vo,vTe0/self.vo,ze0/self.ro,vze0/self.vo,phie0])
		self.oe=Orbit(vxvve,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
		ts=np.linspace(self.tend/to,self.tstart/to,1000)
		self.oe.integrate(ts,self.pot)

		return self.o,self.oe,self.tesc,self.vescape

	def prob_three_body_escape(self,ms,m_a,m_b,q=-2):
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
		vs_peak=0.5*np.sqrt((M-ms)/(ms*M))*np.sqrt(e0)

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

	def animate(self,nsnap=100,xlim=10,ylim=10):

		to=bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)*1000.0
		ts=np.linspace(self.tend/to,self.tstart/to,nsnap)


		for i in range(0,len(ts)):
			plt.plot(self.o.x(ts[i]),o.y(ts[i]),'o')

			escindx=self.tesc/to<ts[i]

			for j in range(0,np.sum(escindx)):
				plt.plot(self.os[escindx][j].x(ts[i]/to),self.os[escindx][j].y(ts[i]/to),'k.')

			plt.xlabel('X (kpc)')
			plt.ylabel('Y (kpc)')
			plt.xlim(-xlim,xlim)
			plt.ylim(-ylim,ylim)

			plt.savefig('%s.png' % str(i).zfill(5))
			plt.close()






