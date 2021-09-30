import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit

def test_find_centre_of_density(tol=0.01):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_density(cluster)

	assert np.fabs(1.0-xc/o.x()) <= tol
	assert np.fabs(1.0-yc/o.y()) <= tol
	assert np.fabs(1.0-zc/o.z()) <= tol


	assert np.fabs(1.0-vxc/o.vx()) <= tol
	assert np.fabs(1.0-vyc/o.vy()) <= tol
	assert np.fabs(1.0-vzc/o.vz()) <= tol

def test_find_centre_of_mass(tol=0.01):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])

	xc, yc, zc, vxc, vyc, vzc=ctools.find_centre_of_mass(cluster)

	assert np.fabs(1.0-xc/o.x()) <= tol
	assert np.fabs(1.0-yc/o.y()) <= tol
	assert np.fabs(1.0-zc/o.z()) <= tol


	assert np.fabs(1.0-vxc/o.vx()) <= tol
	assert np.fabs(1.0-vyc/o.vy()) <= tol
	assert np.fabs(1.0-vzc/o.vz()) <= tol

def test_find_centre(tol=0.01):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])

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
	x=np.append(np.ones(int(ntot/2)),np.ones(int(ntot/2))*2.)
	y,z=np.zeros(ntot),np.zeros(ntot)

	m=np.ones(ntot)
	vx=np.ones(ntot)
	vy,vz=np.zeros(ntot),np.zeros(ntot)

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)

	coulomb=0.4
	rc,rh=1.,1.
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
	vx=np.array([-1.,1.])

	y,z=np.zeros(2)
	vy,vz=np.zeros(2)

	dx=np.array([2.,2.])

	kin=0.5*vx**2.
	pot=-1./dx
	etot=kin+pot
	ektot=np.sum(kin)
	ptot=np.sum(pot)/2.

	cluster.add_stars(x,y,z,vx,vy,vz,m=m)
	cluster.energies()

	print(ektot,cluster.ektot)

	assert np.fabs(1.0-ektot/cluster.ektot) <= tol
	assert np.fabs(1.0-ptot/cluster.ptot) <= tol
	np.testing.assert_array_equal(etot,cluster.etot)

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

def test_closest_star(tol=0.01):

	x,y,z=np.linspace(0,4,5),np.zeros(5),np.zeros(5)
	vx,vy,vz=np.zeros(5),np.zeros(5),np.zeros(5)
	cluster=ctools.StarCluster(units='nbody',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz)

	mindx=ctools.closest_star(cluster)
	np.testing.assert_array_equal(mindx,np.ones(5))

	cluster.z=cluster.x
	mindx=ctools.closest_star(cluster,projected=True)
	np.testing.assert_array_equal(mindx,np.ones(5))


def test_rlagrange(tol=0.01):
	x,y,z=np.linspace(0,10,21),np.zeros(21),np.zeros(21)
	vx,vy,vz=np.zeros(21),np.zeros(21),np.zeros(21)
	cluster=ctools.StarCluster(units='nbody',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz)

	rn=ctools.rlagrange(cluster)
	np.testing.assert_array_equal(rn,np.linspace(1,10,10))

	cluster.z=np.linspace(0,10,21)
	cluster.analyze(sortstars=True)

	rn=ctools.rlagrange(cluster,projected=True)
	np.testing.assert_array_equal(rn,np.linspace(1,10,10))

#def test_virial_radius():


