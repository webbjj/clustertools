import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014

def test_initialize_orbit(tol=0.001):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])

	ocluster=cluster.initialize_orbit()

	assert np.fabs(o.x()-ocluster.x()) <= tol
	assert np.fabs(o.y()-ocluster.y()) <= tol
	assert np.fabs(o.z()-ocluster.z()) <= tol
	assert np.fabs(o.vx()-ocluster.vx()) <= tol
	assert np.fabs(o.vy()-ocluster.vy()) <= tol
	assert np.fabs(o.vz()-ocluster.vz()) <= tol

	ocluster=cluster.initialize_orbit(from_centre=True)

	assert np.fabs(o.x()-ocluster.x()+cluster.xc) <= tol
	assert np.fabs(o.y()-ocluster.y()+cluster.yc) <= tol
	assert np.fabs(o.z()-ocluster.z()+cluster.zc) <= tol
	assert np.fabs(o.vx()-ocluster.vx()+cluster.vxc) <= tol
	assert np.fabs(o.vy()-ocluster.vy()+cluster.vyc) <= tol
	assert np.fabs(o.vz()-ocluster.vz()+cluster.vzc) <= tol

def test_initialize_orbits(tol=0.1):
	
	x=np.random.rand(1000)
	y=np.random.rand(1000)
	z=np.random.rand(1000)
	vx=np.random.rand(1000)
	vy=np.random.rand(1000)
	vz=np.random.rand(1000)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=1.,sortstars=True,analyze=True)

	cluster.add_orbit(8000.,0.,0.,0.,220.,0.)
	cluster.to_galaxy()
	cluster.to_kpckms()

	ocluster=cluster.initialize_orbits()

	np.testing.assert_allclose(8.0+x/1000.,ocluster.x(),rtol=0.0001)
	np.testing.assert_allclose(y/1000.,ocluster.y(),rtol=0.0001)
	np.testing.assert_allclose(z/1000.,ocluster.z(),rtol=0.0001)

	np.testing.assert_allclose(vx,ocluster.vx(),rtol=0.0001)
	np.testing.assert_allclose(220.+vy,ocluster.vy(),rtol=0.0001)
	np.testing.assert_allclose(vz,ocluster.vz(),rtol=0.0001)

	ocluster=cluster.initialize_orbits(from_centre=True)

	np.testing.assert_allclose(8.0+x/1000.,ocluster.x()+cluster.xc,rtol=0.0001)
	np.testing.assert_allclose(y/1000.,ocluster.y()+cluster.yc,rtol=0.0001)
	np.testing.assert_allclose(z/1000.,ocluster.z()+cluster.zc,rtol=0.0001)

	np.testing.assert_allclose(vx,ocluster.vx()+cluster.vxc,rtol=0.0001)
	np.testing.assert_allclose(220.+vy,ocluster.vy()+cluster.vyc,rtol=0.0001)
	np.testing.assert_allclose(vz,ocluster.vz()+cluster.vzc,rtol=0.0001)

def test_integrate_orbit(tol=0.1):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_kpckms()
	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])

	ts,ocluster=cluster.integrate_orbit(nt=1000)

	assert len(ts)==1000
	assert cluster.xgc == o.x()
	assert cluster.xgc = ocluster.x(ts[0])

	o.integrate(ts,MWPotential2014)

	assert ocluster.x(ts[-1])==o.x(ts[-1])

def test_integrate_orbits(tol=0.1):
	pass

def test_orbit_interpolate(tol=0.1):
	pass

def test_orbital_path(tol=0.1):
	pass

def test_orbital_path_match(tol=0.1):
	pass

def test_calc_actions(tol=0.1):
	pass
	
def test_ttensor(tol=0.1):
	pass