import clustertools as ctools
from clustertools.analysis.functions import tpl_func
from scipy.optimize import curve_fit

import pytest

import numpy as np
from galpy.orbit import Orbit
from galpy.potential import NFWPotential,MWPotential2014,rtide,evaluateDensities
from galpy.df import isotropicNFWdf
from galpy.util import bovy_conversion

try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion

try:
	import amuse.units.units as u
	noamuse=False
except:
	noamuse=True

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]
solar_zo=0.0208 #Bennett & Bovy 2019
to=conversion.time_in_Gyr(ro=solar_ro,vo=solar_vo)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_find_centre_of_density(tol=0.01):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster.to_amuse()

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_density(cluster)

	assert np.fabs(1.0-xc/(o.x() | u.kpc)) <= tol
	assert np.fabs(1.0-yc/(o.y() | u.kpc)) <= tol
	assert np.fabs(1.0-zc/(o.z() | u.kpc)) <= tol


	assert np.fabs(1.0-vxc/(o.vx() | u.kms)) <= tol
	assert np.fabs(1.0-vyc/(o.vy() | u.kms)) <= tol
	assert np.fabs(1.0-vzc/(o.vz() | u.kms)) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_find_centre_of_density_casertano(tol=0.01):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster.to_amuse()

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_density(cluster,method='casertano')

	assert np.fabs(1.0-xc/(o.x() | u.kpc)) <= tol
	assert np.fabs(1.0-yc/(o.y() | u.kpc)) <= tol
	assert np.fabs(1.0-zc/(o.z() | u.kpc)) <= tol


	assert np.fabs(1.0-vxc/(o.vx() | u.kms)) <= tol
	assert np.fabs(1.0-vyc/(o.vy() | u.kms)) <= tol
	assert np.fabs(1.0-vzc/(o.vz() | u.kms)) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_find_centre_of_mass(tol=0.01):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster.to_amuse()

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_mass(cluster)

	assert np.fabs(1.0-xc/(o.x() | u.kpc)) <= tol
	assert np.fabs(1.0-yc/(o.y() | u.kpc)) <= tol
	assert np.fabs(1.0-zc/(o.z() | u.kpc)) <= tol


	assert np.fabs(1.0-vxc/(o.vx() | u.kms)) <= tol
	assert np.fabs(1.0-vyc/(o.vy() | u.kms)) <= tol
	assert np.fabs(1.0-vzc/(o.vz() | u.kms)) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_find_centre(tol=0.01):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster.to_amuse()

	print('DEBUG TEST1',cluster.xgc,cluster.vxgc,cluster.xc,cluster.vxc)


	cluster.find_centre()

	print('DEBUG TEST2',cluster.xgc,cluster.vxgc,cluster.xc,cluster.vxc)

	assert np.fabs(1.0-cluster.xgc/(o.x() | u.kpc)) <= tol
	assert np.fabs(1.0-cluster.ygc/(o.y() | u.kpc)) <= tol
	assert np.fabs(1.0-cluster.zgc/(o.z() | u.kpc)) <= tol


	assert np.fabs(1.0-cluster.vxgc/(o.vx() | u.kms)) <= tol
	assert np.fabs(1.0-cluster.vygc/(o.vy() | u.kms)) <= tol
	assert np.fabs(1.0-cluster.vzgc/(o.vz() | u.kms)) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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

	#Test calculation in amuse
	cluster.to_amuse()

	assert np.fabs(1.0-(trelax | u.Myr)/cluster.relaxation_time()) <= tol
	assert np.fabs(1.0-(trelax | u.Myr)/cluster.trelax) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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

	#Test calculation in amuse
	cluster.to_amuse()

	assert np.fabs(1.0-(trelax | u.Myr)/cluster.relaxation_time(projected=True)) <= tol
	assert np.fabs(1.0-(trelax | u.Myr)/cluster.trelax) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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

	#Test calculation in amuse
	cluster.to_amuse()

	assert np.fabs(1.0-(trelax | u.Myr)/cluster.half_mass_relaxation_time()) <= tol
	assert np.fabs(1.0-(trelax | u.Myr)/cluster.trh) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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

	#Test calculation in amuse
	cluster.to_amuse()

	assert np.fabs(1.0-(trelax | u.Myr)/cluster.core_relaxation_time()) <= tol
	assert np.fabs(1.0-(trelax | u.Myr)/cluster.trc) <= tol

	lnlambda=np.log(0.2*ntot)
	trelax=(0.39/lnlambda)*np.sqrt(rc**3./(grav*(ntot)))*(ntot)*np.sqrt(rc*rh)/(rc+rh)

	assert np.fabs(1.0-(trelax | u.Myr)/cluster.core_relaxation_time(coulomb=0.2)) <= tol
	assert np.fabs(1.0-(trelax | u.Myr)/cluster.trc) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_energies(tol=0.01):

	cluster=ctools.StarCluster(units='pckms',origin='centre')

	m=np.array([1.,1.])
	x=np.array([-1.,1.])
	vx=np.array([-2.,2.])

	y,z=np.zeros(2)
	vy,vz=np.zeros(2)

	dx=np.array([2.,2.])

	kin=0.5*vx**2.
	pot=-4.302e-3/dx
	etot=kin+pot
	ektot=np.sum(kin)
	ptot=np.sum(pot)/2.

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)

	cluster.to_amuse()

	cluster.analyze()
	cluster.energies()

	assert np.fabs(1.0-(ektot | (u.kms*u.kms)) /cluster.ektot) <= tol
	assert np.fabs(1.0-(ptot | (u.kms*u.kms)) /cluster.ptot) <= tol
	np.testing.assert_array_equal(etot | (u.kms*u.kms) ,cluster.etot)

	cluster.energies(parallel=True)

	assert np.fabs(1.0-(ektot | (u.kms*u.kms)) /cluster.ektot) <= tol
	assert np.fabs(1.0-(ptot | (u.kms*u.kms)) /cluster.ptot) <= tol
	np.testing.assert_array_equal(etot | (u.kms*u.kms) ,cluster.etot)

	kid,pid=cluster.energies(ids=cluster.id[0])
	eid=kid+pid


	assert(eid[0]==etot[0] | (u.kms*u.kms))

	kid,pid=cluster.energies(ids=[cluster.id[0],cluster.id[1]])
	eid=kid+pid

	assert(eid[0]==etot[0] | (u.kms*u.kms))
	assert(eid[1]==etot[1] | (u.kms*u.kms))

	ids=np.ones(cluster.ntot,dtype=bool)

	kid,pid=cluster.energies(ids=ids)
	eid=kid+pid
	assert(eid[0]==etot[0] | (u.kms*u.kms))
	assert(eid[1]==etot[1] | (u.kms*u.kms))

	ids[0]=False
	kid,pid=cluster.energies(ids=ids)
	eid=kid+pid
	assert(eid[0]==etot[1] | (u.kms*u.kms))

	cluster.z=x | u.parsec
	cluster.vz=vx | u.kms
	cluster.analyze()
	cluster.energies(projected=True)

	assert np.fabs(1.0-(ektot | (u.kms*u.kms)) /cluster.ektot) <= tol
	assert np.fabs(1.0-(ptot | (u.kms*u.kms)) /cluster.ptot) <= tol
	np.testing.assert_array_equal(etot | (u.kms*u.kms),cluster.etot)

	cluster.m=np.array([0.1,0.1]) | u.MSun
	cluster.analyze()
	cluster.energies(projected=True,specific=False)

	assert np.fabs(1.0-0.1*(ektot | (u.MSun*u.kms*u.kms)) /cluster.ektot) <= tol
	assert np.fabs(1.0-0.1*0.1*(ptot | (u.MSun*u.kms*u.kms)) /cluster.ptot) <= tol
	np.testing.assert_array_equal(0.1*(kin| (u.MSun*u.kms*u.kms))+0.1*0.1*(pot| (u.MSun*u.kms*u.kms)),cluster.etot)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_closest_star(tol=0.01):

	x,y,z=np.linspace(0,4,5),np.zeros(5),np.zeros(5)
	vx,vy,vz=np.zeros(5),np.zeros(5),np.zeros(5)
	cluster=ctools.StarCluster(units='nbody',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.analyze()

	mindx=ctools.closest_star(cluster)
	np.testing.assert_array_equal(mindx,np.ones(5))

	cluster.z=cluster.x
	mindx=ctools.closest_star(cluster,projected=True)
	np.testing.assert_array_equal(mindx,np.ones(5))

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_rlagrange(tol=0.01):
	x,y,z=np.linspace(0,10,21),np.zeros(21),np.zeros(21)
	vx,vy,vz=np.zeros(21),np.zeros(21),np.zeros(21)
	cluster=ctools.StarCluster(units='nbody',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.analyze()
	cluster.to_amuse()

	rn=ctools.rlagrange(cluster)
	np.testing.assert_array_equal(rn.value_in(u.parsec),np.linspace(1,10,10))

	cluster.z=np.linspace(0,10,21) | u.parsec
	cluster.analyze(sortstars=True)

	rn=ctools.rlagrange(cluster,projected=True)
	np.testing.assert_array_equal(rn.value_in(u.parsec),np.linspace(1,10,10))

	r10=ctools.rlagrange(cluster,projected=True,mfrac=0.1)
	assert rn[0]==r10

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_virial_radius(tol=0.01):
	cluster=ctools.load_cluster(ctype='snapshot',filename='g1phi1rv1m1N10000.dat',units='nbody',origin='cluster')
	cluster.zmbar=10000
	cluster.reset_nbody_scale(mass=False,radii=False)
	cluster.to_amuse()
	rv=ctools.virial_radius(cluster,method='inverse_distance')
	assert np.fabs(rv.value_in(u.parsec)-1.)<=tol

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
	cluster.to_amuse()

	rv=ctools.virial_radius(cluster,method='critical_density')

	assert np.fabs(rv.value_in(u.kpc)-rvnfw) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_mass_function(tol=0.01):

	alphatest=-1
	m=ctools.power_law_distribution_function(10000,alphatest,0.1,1.)

	x,y,z=np.zeros(len(m)),np.zeros(len(m)),np.zeros(len(m))
	vx,vy,vz=np.zeros(len(m)),np.zeros(len(m)),np.zeros(len(m))

	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m)
	cluster.analyze()
	cluster.to_amuse()

	m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha=ctools.mass_function(cluster)

	assert np.fabs(alpha-alphatest) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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
	cluster.to_amuse()

	lower_bounds=[5000.,-2.,np.amin(cluster.m.value_in(u.MSun)),1.]
	upper_bounds=[15000.,0.,np.amax(cluster.m.value_in(u.MSun)),3.]

	m_mean, m_hist, dm, Afit, eA, alphafit, ealpha, mcfit, emc, betafit, ebeta=ctools.tapered_mass_function(cluster,lower_bounds=lower_bounds,upper_bounds=upper_bounds)

	assert np.fabs(alphafit-alpha) <= ealpha*tol
	assert np.fabs(mcfit-mc) <= emc*tol
	assert np.fabs(betafit-beta) <= ebeta*tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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
	cluster.to_amuse()
	m_mean, sigvm, eta, eeta, yeta, eyeta=ctools.eta_function(cluster)

	assert np.fabs(eta-etatest) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_meq_function(tol=0.1):
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
	cluster.to_amuse()
	m_mean, sigvm, meq, emq, sigma0, esigma0=ctools.meq_function(cluster)

	assert np.fabs(meq-meqtest) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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
	cluster.to_amuse()

	ck=ctools.ckin(cluster)

	assert np.fabs(ck-1.)<=tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_rcore(tol=0.1):

	wdir='../docs/source/notebooks/'
	rc0=0.4465

	#Plummer sphere with rm=1.
	a0=0.76628
	rc0=0.64*a0

	m,x,y,z,vx,vy,vz=np.loadtxt(wdir+'N10k.dat.10',unpack=True)
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.to_amuse()

	rc=ctools.rcore(cluster)

	print(rc,rc0)

	assert np.fabs(rc.value_in(u.parsec)-rc0) <= tol

	m,x,y,z,vx,vy,vz=np.loadtxt(wdir+'N1k.dat.10',unpack=True)
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.to_amuse()
	rc=ctools.rcore(cluster,method='isothermal')
	assert np.fabs(rc.value_in(u.parsec)-rc0) <= tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_rtide(tol=0.1):
	#test rtide is the same as galpy calculates
	mo=bovy_conversion.mass_in_msol(ro=solar_ro,vo=solar_vo)

	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	m=cluster.mtot
	Rgc=np.sqrt(cluster.xgc**2.+cluster.ygc**2.)/1000.
	zgc=cluster.zgc/1000.

	rtgalpy=rtide(MWPotential2014,Rgc/solar_ro,zgc/solar_ro,M=m/mo,ro=solar_ro,vo=solar_vo)*1000.0
	
	cluster.to_amuse()

	rtctools=ctools.rtidal(cluster,MWPotential2014)

	assert np.fabs(rtgalpy-rtctools.value_in(u.parsec)) < tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_rlimiting(tol=0.1):
	#test the densities are equal

	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)

	ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion
	mo=bovy_conversion.mass_in_msol(ro=ro,vo=vo)

	m=cluster.mtot
	Rgc=np.sqrt(cluster.xgc**2.+cluster.ygc**2.)/1000.
	zgc=cluster.zgc/1000.
	pot=MWPotential2014

	rholocal=evaluateDensities(pot,Rgc/ro,zgc/ro,ro=ro,vo=vo)

	cluster.to_amuse()
	rl=ctools.rlimiting(cluster,MWPotential2014)

	rprof, pprof, nprof = ctools.rho_prof(cluster, nrad=20, projected=False)

	rlarg=np.argmin(np.fabs(rprof.value_in(u.parsec)-rl.value_in(u.parsec)))

	if rprof[rlarg] > rl:
		assert pprof[rlarg-1].value_in(u.MSun/(u.parsec**3.)) > rholocal
		assert pprof[rlarg].value_in(u.MSun/(u.parsec**3.)) < rholocal
	else:
		assert pprof[rlarg].value_in(u.MSun/(u.parsec**3.)) > rholocal
		assert pprof[rlarg+1].value_in(u.MSun/(u.parsec**3.)) < rholocal

	rl=ctools.rlimiting(cluster,MWPotential2014,projected=True)

	rprof, pprof, nprof = ctools.rho_prof(cluster, nrad=20, projected=True)

	rlarg=np.argmin(np.fabs(rprof.value_in(u.parsec)-rl.value_in(u.parsec)))

	#Approximate projected local density across area of cluster
	rholocal_pro=rholocal*(4.*rprof[-1].value_in(u.parsec)/3)

	if rprof[rlarg] > rl:
		assert pprof[rlarg-1].value_in(u.MSun/(u.parsec**2.)) > rholocal_pro
		assert pprof[rlarg].value_in(u.MSun/(u.parsec**2.)) < rholocal_pro
	else:
		assert pprof[rlarg].value_in(u.MSun/(u.parsec**2.)) > rholocal_pro
		assert pprof[rlarg+1].value_in(u.MSun/(u.parsec**2.)) < rholocal_pro
