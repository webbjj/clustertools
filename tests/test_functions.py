import clustertools as ctools
from clustertools.analysis.functions import tpl_func
from scipy.optimize import curve_fit

import numpy as np
from galpy.orbit import Orbit
from galpy.potential import NFWPotential,MWPotential2014,rtide,evaluateDensities
from galpy.df import isotropicNFWdf
from galpy.util import bovy_conversion

import limepy
from limepy import limepy

try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion

try:
    import amuse.units.units as u
except:
		pass

import time as time

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]
solar_zo=0.0208 #Bennett & Bovy 2019
to=conversion.time_in_Gyr(ro=solar_ro,vo=solar_vo)

def test_find_centre_of_density(tol=0.01):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_density(cluster)

	assert np.fabs(1.0-xc/o.x()) <= tol
	assert np.fabs(1.0-yc/o.y()) <= tol
	assert np.fabs(1.0-zc/o.z()) <= tol


	assert np.fabs(1.0-vxc/o.vx()) <= tol
	assert np.fabs(1.0-vyc/o.vy()) <= tol
	assert np.fabs(1.0-vzc/o.vz()) <= tol

def test_find_centre_of_density_casertano(tol=0.01):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_density(cluster,method='casertano')

	assert np.fabs(1.0-xc/o.x()) <= tol
	assert np.fabs(1.0-yc/o.y()) <= tol
	assert np.fabs(1.0-zc/o.z()) <= tol


	assert np.fabs(1.0-vxc/o.vx()) <= tol
	assert np.fabs(1.0-vyc/o.vy()) <= tol
	assert np.fabs(1.0-vzc/o.vz()) <= tol


def test_find_centre_of_mass(tol=0.01):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_mass(cluster)

	assert np.fabs(1.0-xc/o.x()) <= tol
	assert np.fabs(1.0-yc/o.y()) <= tol
	assert np.fabs(1.0-zc/o.z()) <= tol


	assert np.fabs(1.0-vxc/o.vx()) <= tol
	assert np.fabs(1.0-vyc/o.vy()) <= tol
	assert np.fabs(1.0-vzc/o.vz()) <= tol

def test_find_centre(tol=0.01):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster.find_centre()

	assert np.fabs(1.0-cluster.xgc/o.x()) <= tol
	assert np.fabs(1.0-cluster.ygc/o.y()) <= tol
	assert np.fabs(1.0-cluster.zgc/o.z()) <= tol


	assert np.fabs(1.0-cluster.vxgc/o.vx()) <= tol
	assert np.fabs(1.0-cluster.vygc/o.vy()) <= tol
	assert np.fabs(1.0-cluster.vzgc/o.vz()) <= tol

def test_relaxation_time(tol=0.01):
	cluster=ctools.StarCluster(units='pckms',origin='centre')

	#Setup cluster with rad=1, mean velocity=1, mass=ntot
	ntot=1000
	grav=4.302e-3
	x=np.append(np.ones(int(ntot/2)),np.ones(int(ntot/2))*2.)
	y,z=np.zeros(ntot),np.zeros(ntot)

	m=np.ones(ntot)
	vx=np.ones(ntot)
	vy,vz=np.zeros(ntot),np.zeros(ntot)

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.analyze()

	coulomb=0.4
	rad=1.
	mbar=1.
	vol=4.0*np.pi/3.0
	rho=(ntot/2.)/vol
	v2=1.
	lnlambda=np.log(coulomb*ntot)

	trelax=1./(15.4*grav**2.*rho*lnlambda)

	# Units of Myr
	trelax*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

	assert np.fabs(1.0-trelax/cluster.relaxation_time()) <= tol

	#Test calculation in Gyr
	cluster.to_kpckms()
	assert np.fabs(1.0-(trelax/1000.0)/cluster.relaxation_time()) <= tol

	#Test calculation in galpy
	cluster.to_galpy()
	assert np.fabs(1.0-(trelax/1000.0/to)/cluster.relaxation_time()) <= tol

	#Test calculation in amuse
	cluster.to_amuse()

	assert np.fabs(1.0-(trelax | u.Myr)/cluster.relaxation_time()) <= tol

def test_projected_relaxation_time(tol=0.01):
	cluster=ctools.StarCluster(units='pckms',origin='centre')

	#Setup cluster with rad=1, mean velocity=1, mass=ntot
	ntot=1000
	grav=4.302e-3
	x=np.append(np.ones(int(ntot/2)),np.ones(int(ntot/2))*2.)
	y,z=np.zeros(ntot),np.random.rand(ntot)

	m=np.ones(ntot)
	vx=np.ones(ntot)
	vy,vz=np.zeros(ntot),np.random.rand(ntot)

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.analyze()

	coulomb=0.4
	rad=1.
	mbar=1.
	vol=4.0*np.pi/3.0
	rho=(ntot/2.)/vol
	v2=1.
	lnlambda=np.log(coulomb*ntot)

	trelax=1./(15.4*grav**2.*rho*lnlambda)

	# Units of Myr
	trelax*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

	assert np.fabs(1.0-trelax/cluster.relaxation_time(projected=True)) <= tol


	#Test calculation in Gyr
	cluster.to_kpckms()
	assert np.fabs(1.0-(trelax/1000.0)/cluster.relaxation_time(projected=True)) <= tol

def test_half_mass_relaxation_time(tol=0.01):
	cluster=ctools.StarCluster(units='pckms',origin='centre')

	#Setup cluster with rad=1, mean velocity=1, mass=ntot
	ntot=1000
	grav=4.302e-3
	x=np.append(np.ones(int(ntot/2)),np.ones(int(ntot/2))*2.)
	y,z=np.zeros(ntot),np.zeros(ntot)

	m=np.ones(ntot)
	vx=np.ones(ntot)
	vy,vz=np.zeros(ntot),np.zeros(ntot)

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.analyze()

	coulomb=0.4
	rad=1.
	mbar=1.
	vol=4.0*np.pi/3.0
	rho=(ntot/2.)/vol
	v2=1.
	lnlambda=np.log(coulomb*ntot)

	trelax=0.138*((mbar+ntot)**0.5)/(np.sqrt(grav)*lnlambda)

	# Units of Myr
	trelax*= 3.086e13 / (3600.0 * 24.0 * 365.0 * 1000000.0)

	assert np.fabs(1.0-trelax/cluster.half_mass_relaxation_time()) <= tol

def test_core_relaxation_time(tol=0.01):
	cluster=ctools.StarCluster(units='pckms',origin='centre')

	#Setup cluster with rad=1, mean velocity=1, mass=ntot
	ntot=1000
	grav=4.302e-3
	x=np.append(np.random.rand(int(ntot/2)),1.0+np.random.rand(int(ntot/2)))
	y,z=np.zeros(ntot),np.zeros(ntot)

	m=np.ones(ntot)
	vx=np.ones(ntot)
	vy,vz=np.zeros(ntot),np.zeros(ntot)

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.analyze()

	coulomb=0.4
	rh=1.
	rc=cluster.rcore()
	mbar=1.
	vol=4.0*np.pi/3.0
	rho=(ntot/2.)/vol
	v2=1.
	lnlambda=np.log(coulomb*ntot)

	trelax=(0.39/lnlambda)*np.sqrt(rc**3./(grav*(ntot)))*(ntot)*np.sqrt(rc*rh)/(rc+rh)

	assert np.fabs(1.0-trelax/cluster.core_relaxation_time()) <= tol

	lnlambda=np.log(0.2*ntot)
	trelax=(0.39/lnlambda)*np.sqrt(rc**3./(grav*(ntot)))*(ntot)*np.sqrt(rc*rh)/(rc+rh)
	assert np.fabs(1.0-trelax/cluster.core_relaxation_time(coulomb=0.2)) <= tol

def test_energies(tol=0.01):

	cluster=ctools.StarCluster(units='nbody',origin='centre')

	m=np.array([1.,1.])
	x=np.array([-1.,1.])
	vx=np.array([-2.,2.])

	y,z=np.zeros(2)
	vy,vz=np.zeros(2)

	dx=np.array([2.,2.])

	kin=0.5*vx**2.
	pot=-1./dx
	etot=kin+pot
	ektot=np.sum(kin)
	ptot=np.sum(pot)/2.

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.analyze()
	cluster.energies()

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

	cluster.energies(parallel=True)

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

	kid,pid=cluster.energies(ids=cluster.id[0])
	eid=kid+pid


	assert(eid[0]==etot[0])

	kid,pid=cluster.energies(ids=[cluster.id[0],cluster.id[1]])
	eid=kid+pid

	assert(eid[0]==etot[0])
	assert(eid[1]==etot[1])

	ids=np.ones(cluster.ntot,dtype=bool)

	kid,pid=cluster.energies(ids=ids)
	eid=kid+pid
	assert(eid[0]==etot[0])
	assert(eid[1]==etot[1])

	ids[0]=False
	kid,pid=cluster.energies(ids=ids)
	eid=kid+pid
	assert(eid[0]==etot[1])

	cluster.z=x
	cluster.vz=vx
	cluster.analyze()
	cluster.energies(projected=True)

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

	cluster.m=np.array([0.1,0.1])
	cluster.analyze()
	cluster.energies(projected=True,specific=False)

	assert np.fabs(1.0-0.1*ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-0.1*0.1*ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(0.1*kin+0.1*0.1*pot,cluster.etot)

def test_softened_energies(tol=0.01):

	softening=0.1

	cluster=ctools.StarCluster(units='nbody',origin='centre')

	m=np.array([1.,1.])
	x=np.array([-1.,1.])
	vx=np.array([-2.,2.])

	y,z=np.zeros(2)
	vy,vz=np.zeros(2)

	dx=np.array([2.,2.])

	kin=0.5*vx**2.
	pot=-1./np.sqrt(dx**2.+softening**2.)
	etot=kin+pot
	ektot=np.sum(kin)
	ptot=np.sum(pot)/2.

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.analyze()
	cluster.energies(softening=softening)

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

	cluster.energies(parallel=True,softening=softening)

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

	kid,pid=cluster.energies(ids=cluster.id[0],softening=softening)
	eid=kid+pid


	assert(eid[0]==etot[0])

	kid,pid=cluster.energies(ids=[cluster.id[0],cluster.id[1]],softening=softening)
	eid=kid+pid

	assert(eid[0]==etot[0])
	assert(eid[1]==etot[1])

	ids=np.ones(cluster.ntot,dtype=bool)

	kid,pid=cluster.energies(ids=ids,softening=softening)
	eid=kid+pid
	assert(eid[0]==etot[0])
	assert(eid[1]==etot[1])

	ids[0]=False
	kid,pid=cluster.energies(ids=ids,softening=softening)
	eid=kid+pid
	assert(eid[0]==etot[1])

	cluster.z=x
	cluster.vz=vx
	cluster.analyze()
	cluster.energies(projected=True,softening=softening)

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

	cluster.m=np.array([0.1,0.1])
	cluster.analyze()
	cluster.energies(projected=True,specific=False,softening=softening)

	assert np.fabs(1.0-0.1*ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-0.1*0.1*ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(0.1*kin+0.1*0.1*pot,cluster.etot)

def test_closest_star(tol=0.01):

	x,y,z=np.linspace(0,4,5),np.zeros(5),np.zeros(5)
	vx,vy,vz=np.zeros(5),np.zeros(5),np.zeros(5)
	cluster=ctools.StarCluster(units='nbody',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.analyze()

	mindx=ctools.closest_star(cluster)
	np.testing.assert_array_equal(mindx,np.ones(5))

	dist,arg=ctools.closest_star(cluster,argument=True)
	np.testing.assert_array_equal(mindx,dist)

	cluster.z=cluster.x
	mindx=ctools.closest_star(cluster,projected=True)
	np.testing.assert_array_equal(mindx,np.ones(5))
	
	dist,arg=ctools.closest_star(cluster,argument=True,projected=True)
	np.testing.assert_array_equal(mindx,dist)


def test_rlagrange(tol=0.01):
	x,y,z=np.linspace(0,10,21),np.zeros(21),np.zeros(21)
	vx,vy,vz=np.zeros(21),np.zeros(21),np.zeros(21)
	cluster=ctools.StarCluster(units='nbody',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.analyze()

	rn=ctools.rlagrange(cluster)
	np.testing.assert_array_equal(rn,np.linspace(1,10,10))

	cluster.z=np.linspace(0,10,21)
	cluster.analyze(sortstars=True)

	rn=ctools.rlagrange(cluster,projected=True)
	np.testing.assert_array_equal(rn,np.linspace(1,10,10))

	r10=ctools.rlagrange(cluster,projected=True,mfrac=0.1)
	assert rn[0]==r10

def test_virial_radius(tol=0.01):
	cluster=ctools.load_cluster(ctype='limepy',g=1.,phi0=1.,rv=1.,m=1.,N=10000)

	rv=ctools.virial_radius(cluster,method='inverse_distance')
	assert np.fabs(rv-1.)<=tol

	np.random.seed(5)

	#Generate GALPY NFW Halo and calculate virial radius via critical density method
	mhalo=2562109.8143978487
	chalo=24.31377803
	mo=bovy_conversion.mass_in_msol(ro=solar_ro,vo=solar_vo)

	nfw=NFWPotential(mvir=mhalo/1e12,conc=chalo,wrtcrit=True,ro=solar_ro,vo=solar_vo)
	rvnfw=nfw.rvir(wrtcrit=True,ro=solar_ro,vo=solar_vo)
	ndf= isotropicNFWdf(pot=nfw,rmax=rvnfw/8.,ro=solar_ro,vo=solar_vo) 
	os= ndf.sample(n=100000)
	msub=nfw.mass(rvnfw/8.)/100000
	cluster=ctools.StarCluster(units='kpckms',origin='cluster')
	cluster.add_stars(os.x(),os.y(),os.z(),os.vx(),os.vy(),os.vz(),m=msub)
	cluster.analyze()
	rv=ctools.virial_radius(cluster,method='critical_density')

	assert np.fabs(rv-rvnfw) <= tol

def test_mass_function(tol=0.01):

	alphatest=-1
	m=ctools.power_law_distribution_function(10000,alphatest,0.1,1.)

	x,y,z=np.zeros(len(m)),np.zeros(len(m)),np.zeros(len(m))
	vx,vy,vz=np.zeros(len(m)),np.zeros(len(m)),np.zeros(len(m))

	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m)
	cluster.analyze()
	m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha=ctools.mass_function(cluster)

	assert np.fabs(alpha-alphatest) <= tol

def test_tapered_mass_function(tol=3.):

	msample=np.linspace(0.1,1.,100)
	deltam=(msample[1]-msample[0])/2.
	A=10000
	alpha=-1.
	mc=0.5
	beta=2.

	m=np.array([])
	for i in range(0,len(msample)-1):
	    
	    dm=tpl_func(msample[i]+deltam,A,alpha,mc,beta)*deltam
	    m=np.append(m,np.random.uniform(msample[i],msample[i+1],int(np.ceil(dm))))
 	    

	x,y,z=np.zeros(len(m)),np.zeros(len(m)),np.zeros(len(m))
	vx,vy,vz=np.zeros(len(m)),np.zeros(len(m)),np.zeros(len(m))

	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)

	lower_bounds=[5000.,-2.,np.amin(cluster.m),1.]
	upper_bounds=[15000.,0.,np.amax(cluster.m),3.]

	m_mean, m_hist, dm, Afit, eA, alphafit, ealpha, mcfit, emc, betafit, ebeta=ctools.tapered_mass_function(cluster,lower_bounds=lower_bounds,upper_bounds=upper_bounds)

	assert np.fabs(alphafit-alpha) <= ealpha*tol
	assert np.fabs(mcfit-mc) <= emc*tol
	assert np.fabs(betafit-beta) <= ebeta*tol


def test_eta_function(tol=0.1):
	alphatest=-1
	m=ctools.power_law_distribution_function(10000,alphatest,0.1,1.)
	x,y,z=np.zeros(10000),np.zeros(10000),np.zeros(10000)
	vx=np.array([])

	etatest=-0.25
	for i in range(0,len(x)):
		sig=m[i]**etatest
		vx=np.append(vx,np.random.normal(0.,sig))

	vy,vz=np.zeros(10000),np.zeros(10000)
	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m)
	cluster.analyze()
	m_mean, sigvm, eta, eeta, yeta, eyeta=ctools.eta_function(cluster)

	assert np.fabs(eta-etatest) <= tol

def test_meq_function(tol=0.01):
	"""
	    - As per Bianchini, P. et al. 2016, MNRAS, 458, 3644, velocity dispersion 
      versus mass is fit with the following:
      sigma(m)= sigma e^(-1/2 m/meq) if m<= meq
              = sigma0 e^(-1/2) (m/meq)^-1/2 if m > meq
    """

	alphatest=-1
	m=ctools.power_law_distribution_function(10000,alphatest,0.1,1.)
	x,y,z=np.zeros(10000),np.zeros(10000),np.zeros(10000)
	vx=np.array([])

	meqtest=0.5
	sigma0test=1.

	for i in range(0,len(x)):
		if m[i]<=meqtest:
			sig=sigma0test*np.exp(-0.5*m[i]/meqtest)
		else:
			sig=sigma0test*np.exp(-0.5)*((m[i]/meqtest)**(-0.5))
       
		vxtemp=np.random.normal(0.,sig)
		while vxtemp < 0:
			vxtemp=np.random.normal(0.,sig)

		vx=np.append(vx,vxtemp)

	vy,vz=np.zeros(10000),np.zeros(10000)

	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m)
	cluster.analyze()
	m_mean, sigvm, meq, emq, sigma0, esigma0=ctools.meq_function(cluster)

	assert np.fabs(meq-meqtest) <= tol

def test_ckin(tol=0.1):
	alphatest=-1
	m=ctools.power_law_distribution_function(10000,alphatest,0.1,1.)
	x,y,z=np.zeros(10000),np.zeros(10000),np.zeros(10000)
	vx=np.array([])

	meqtest=0.5
	sigma0test=1.

	for i in range(0,len(x)):
		if m[i]<=meqtest:
			sig=sigma0test*np.exp(-0.5*m[i]/meqtest)
		else:
			sig=sigma0test*np.exp(-0.5)*((m[i]/meqtest)**(-0.5))
	   
		vxtemp=np.random.normal(0.,sig)
		while vxtemp < 0:
			vxtemp=np.random.normal(0.,sig)

		vx=np.append(vx,vxtemp)

	vy,vz=np.zeros(10000),np.zeros(10000)

	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m)
	cluster.analyze()

	ck=ctools.ckin(cluster)

	assert np.fabs(ck-1.)<=tol

def test_rcore(tol=0.1):

	wdir='../docs/source/notebooks/'
	rc0=0.4465

	#Plummer sphere with rm=1.
	a0=0.76628
	rc0=0.64*a0

	m,x,y,z,vx,vy,vz=np.loadtxt(wdir+'N10k.dat.10',unpack=True)
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m)

	rc=ctools.rcore(cluster)

	print(rc,rc0)

	assert np.fabs(rc-rc0) <= tol

	m,x,y,z,vx,vy,vz=np.loadtxt(wdir+'N1k.dat.10',unpack=True)
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	rc=ctools.rcore(cluster,method='isothermal')
	assert np.fabs(rc-rc0) <= tol

def test_rtide(tol=0.1):
	#test rtide is the same as galpy calculates
	mo=bovy_conversion.mass_in_msol(ro=solar_ro,vo=solar_vo)

	cluster=ctools.load_cluster('limepy',gcname='NGC5466')
	m=cluster.mtot
	Rgc=np.sqrt(cluster.xgc**2.+cluster.ygc**2.)/1000.
	zgc=cluster.zgc/1000.

	rtgalpy=rtide(MWPotential2014,Rgc/solar_ro,zgc/solar_ro,M=m/mo,ro=solar_ro,vo=solar_vo)*1000.0
	rtctools=ctools.rtidal(cluster,MWPotential2014)

	assert np.fabs(rtgalpy-rtctools) < tol

def test_rlimiting(tol=0.1):
	#test the densities are equal
	ro,vo=8.,220.
	mo=bovy_conversion.mass_in_msol(ro=ro,vo=vo)

	cluster=ctools.load_cluster('limepy',gcname='NGC5466')

	m=cluster.mtot
	Rgc=np.sqrt(cluster.xgc**2.+cluster.ygc**2.)/1000.
	zgc=cluster.zgc/1000.
	pot=MWPotential2014

	rholocal=evaluateDensities(pot,Rgc/ro,zgc/ro,ro=ro,vo=vo)

	rl=ctools.rlimiting(cluster,MWPotential2014)

	rprof, pprof, nprof = ctools.rho_prof(cluster, nrad=20, projected=False)

	rlarg=np.argmin(np.fabs(rprof-rl))

	if rprof[rlarg] > rl:
		assert pprof[rlarg-1] > rholocal
		assert pprof[rlarg] < rholocal
	else:
		assert pprof[rlarg] > rholocal
		assert pprof[rlarg+1] < rholocal

	rl=ctools.rlimiting(cluster,MWPotential2014,projected=True)

	rprof, pprof, nprof = ctools.rho_prof(cluster, nrad=20, projected=True)

	rlarg=np.argmin(np.fabs(rprof-rl))

	#Approximate projected local density across area of cluster
	rholocal_pro=rholocal*(4.*rprof[-1]/3)

	if rprof[rlarg] > rl:
		assert pprof[rlarg-1] > rholocal_pro
		assert pprof[rlarg] < rholocal_pro
	else:
		assert pprof[rlarg] > rholocal_pro
		assert pprof[rlarg+1] < rholocal_pro
