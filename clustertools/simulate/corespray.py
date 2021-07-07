from galpy.orbit import Orbit
from galpy.potential import MWPotential2014
from galpy.util import bovy_conversion,bovy_plot
import numpy as np

from ..util.recipes import *

class corespray(object):

    """ A class that is responsible for periodcially ejecting stars from a cluster's core
    
    Parameters
    ----------


    History
    -------
    2021 - Written - Grandin/Webb (UofT)

    """

	def __init__(self,gcname,pot=MWPotential2014,mu0=0.,sig0=10.0,vesc0=10.0,mmin=0.1,mmax=1.4,alpha=-1.35,ro=8.,vo=220.):
		
		to=bovy_conversion.time_in_Gyr(ro=ro,vo=vo)

		self.gcname=gcname
	    self.o = Orbit.from_name(gcname, ro=ro, vo=vo, solarmotion=[-11.1, 24.0, 7.25])		

		self.pot=pot

		self.mu0,self.sig0,self.vesc0=mu0,sig0,vesc0
		self.mmin,self.mmax,self.alpha=mmin,mmax,alpha

		self.ro,self.vo=ro,vo


	def sample(self,tstart=-1000.,tend=0.,rate=1.,nstar=None):

		to=bovy_conversion.time_in_Gyr(ro=self.ro,vo=self.vo)

		self.tstart=tstart
		self.tend=tend
		self.dt=tend-tstart

		#Select escape times
		#If nstar is not None, randomly select escapers between tstart and tend
		if nstar is not None:
			self.nstar=nstar
			self.rate=(tend-tstart)/nstar
			self.tesc=tstart+np.random.rand(nstar)*(tend-tstart)
		else:
			self.rate=rate
			self.nstar=(tend-tstart)/rate
			self.tesc=tstart+np.random.rand(nstar)*(tend-tstart)

	    if tstart<tend:
	    	ts=np.linspace(tend,tstart/(to*1000.0),1000)
		    self.o.integrate(ts,self.pot)
    	else:
	    	ts=np.linspace(tstart,tend/(to*1000.0),1000)
		    self.o.integrate(ts,MWPotential2014)
	    
	    #Generate positions and new velocities for escaped stars with velocity dispersions of 100/root(3), based on
	    #gcname's position 1 Gyr ago
	    self.xe=self.o.x(self.tesc/(to*1000.0))
	    self.ye=self.o.y(self.tesc/(to*1000.0))
	    self.ze=self.o.z(self.tesc/(to*1000.0))
	    
	    self.vxe=self.o.vx(self.tesc/(to*1000.0))
	    self.vye=self.o.vy(self.tesc/(to*1000.0))
	    self.vze=self.o.vz(self.tesc/(to*1000.0))
	    
	    self.vescape=np.array([])

	    nescape=0

		while nescape < self.nstar:
			masses=power_law_distribution_function(3, self.alpha, self.mmin, self.mmax)

			mindx=(masses==np.amin(masses))
			ms=masses[mindx][0]
			m_a,m_b=masses[np.invert(mindx)]
			        
			mb=m_a+m_b
			M=ms+mb

			prob=three_body_escape(ms,m_a,m_b)

			if np.random.rand() < prob:
			    
			    vxs,vys,vzs=np.random.normal(0.,self.sig0,3)
			    vxb,vyb,vzb=np.random.normal(0.,self.sig0,3)
			    dv=np.sqrt((vxs-vxb)**2.+(vys-vyb)**2.+(vzs-vzb)**2.)
			    
			    vs=mb/M*dv
			    
			    if vs >self.vesc0:  
			        self.vescape=np.append(self.vescape,vs)
			        r_hat=1.
			        phihat=2.0*np.pi*np.random.rand()
			        thetahat=np.arccos(1.0-2.0*np.random.rand())

			        Rhat=r_hat*np.sin(thetahat)
			        xhat=Rhat*np.cos(phihat)
			        yhat=Rhat*np.sin(phihat)
			        zhat=r_hat*np.cos(thetahat)

			        self.vxe[nescape]+=vs*xhat
			        self.vye[nescape]+=vs*yhat
			        self.vze[nescape]+=vs*zhat

			        nescape+=1
		

		#Need to do an orbit for each star

	    Re, phie, ze, vRe, vTe, vze= ctools.cart_to_cyl(self.xe,self.ye,self.ze,self.vxe,self.vye,self.vze)
	    
	    #Make array of star orbits
	    #Note positions and velocities are scaled by ro and vo respectively to be on Galpy units
	    vxvv=np.column_stack([Re/self.ro,vRe/self.vo,vTe/self.vo,ze/self.ro,vze/self.vo,phie])
	    os=Orbit(vxvv,ro=self.ro,vo=self.vo,solarmotion=[-11.1, 24.0, 7.25])
	    ts=np.linspace(0,tesc/to,1000)
	    os.integrate(ts,pot)

	    return o,os,ts,vescape

	    def three_body_escape(self,ms,m_a,m_b,q=-3):
			prob=(ms**q)/(ms**q+m_a**q+m_b**q)
			return prob

		def sample_binding_energy(mb1,mb2,alpha=-1,xmin=10.0**36,xmax=10.0**40.0):
	        
			#Default binding energy distribution is:
			# power law of slope -1 
			# Between 10.0**36 and 10.0**40 J

			grav=4.302e-3 #pc/Msun (km/s)^2


	        n=len(mb1)

	        ebin_si=xfunc(n,alpha,xmin,xmax) #Joules = kg (m/s)^2
	        ebin=ebin_si/1.9891e30 # Msun (m/s)^2
	        ebin/=(1000.0*1000.0) #Msun (km/s)^2
	        
	    
	        #Log normal a:
	        semi_pc=(0.5*grav*mb1*mb2)/ebin
	        semi_au=semi_pc/4.84814e-6    

	        self.ebin=ebin
	        self.semi=semi_au

		def xfunc(n,alpha,xmin,xmax):

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





