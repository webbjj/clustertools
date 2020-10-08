import clustertools as ctools


def test_init_default():
	cluster=ctools.StarCluster()

	assert cluster.ntot == 0
	assert cluster.tphys == 0.0
	assert cluster.units is None
	assert cluster.origin is None
	assert cluster.ctype == 'snapshot'
	assert cluster.projected == False

	return None

def test_init_custom():
	cluster=ctools.StarCluster(1000,10.0,'cluster','nbody','nbody6',True)

	assert cluster.ntot == 1000
	assert cluster.tphys == 10.0
	assert cluster.units == 'cluster'
	assert cluster.origin == 'nbody'
	assert cluster.ctype == 'nbody6'
	assert cluster.projected == True

	return None