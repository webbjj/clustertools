import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit

def test_find_centre_of_density(tol=0.01):
	cluster=ctools.setup_cluster(ctype='NGC6101')
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
	cluster=ctools.setup_cluster(ctype='NGC6101')
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
	cluster=ctools.setup_cluster(ctype='NGC6101')
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