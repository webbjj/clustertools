import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit


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
	pass

def test_integrate_orbit(tol=0.1):
	pass

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