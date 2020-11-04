import clustertools as ctools
import numpy as np

units=['pckms','kpckms','nbody','galpy']

def test_init_default():
	cluster=ctools.StarCluster()

	assert cluster.tphys == 0.0
	assert cluster.units is None
	assert cluster.origin is None
	assert cluster.ctype == 'snapshot'
	assert cluster.projected == False

	return None

def test_init_custom():

	for u in units:

		cluster=ctools.StarCluster(10.0,'cluster',u,'nbody6',True)

		assert cluster.tphys == 10.0
		assert cluster.units == 'cluster'
		assert cluster.origin == u
		assert cluster.ctype == 'nbody6'
		assert cluster.projected == True

	return None

def test_add_stars_default():
	nstar=100
	cluster=ctools.StarCluster()
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)

	np.testing.assert_array_equal(x,cluster.x)
	np.testing.assert_array_equal(y,cluster.y)
	np.testing.assert_array_equal(z,cluster.z)
	np.testing.assert_array_equal(vx,cluster.vx)
	np.testing.assert_array_equal(vy,cluster.vy)
	np.testing.assert_array_equal(vz,cluster.vz)

	np.testing.assert_array_equal(cluster.m,np.ones(len(x)))
	np.testing.assert_array_equal(cluster.id, np.linspace(0, len(x) - 1, len(x), dtype=int))

	assert cluster.ntot == nstar

	np.testing.assert_array_equal(cluster.kw,np.zeros(nstar))


def test_add_stars_custom():
	nstar=100
	cluster=ctools.StarCluster()
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	m=np.ones(nstar)*0.1
	id=np.linspace(0,nstar,nstar)
	cluster.add_stars(x,y,z,vx,vy,vz,m,id)

	np.testing.assert_array_equal(m,cluster.m)
	np.testing.assert_array_equal(id,cluster.id)

def test_add_stars_lenghtdif():
	nstar=100
	cluster=ctools.StarCluster()
	x,y=np.ones(nstar),np.ones(nstar)
	vx,vy=np.ones(nstar),np.ones(nstar)
	z,vz,m=1.,1.,1.
	cluster.add_stars(x,y,z,vx,vy,vz,m)

	np.testing.assert_array_equal(cluster.m,np.ones(nstar)*m)
	np.testing.assert_array_equal(cluster.z,np.ones(nstar)*z)
	np.testing.assert_array_equal(cluster.vz,np.ones(nstar)*vz)

def test_add_stars_radec():
	nstar=100
	cluster=ctools.StarCluster(units='radec',origin='sky')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)

	np.testing.assert_array_equal(x,cluster.ra)
	np.testing.assert_array_equal(y,cluster.dec)
	np.testing.assert_array_equal(z,cluster.dist)
	np.testing.assert_array_equal(vx,cluster.pmra)
	np.testing.assert_array_equal(vy,cluster.pmdec)
	np.testing.assert_array_equal(vz,cluster.vlos)

def test_add_orbit():
	nstar=100
	cluster=ctools.StarCluster(units='kpckms',origin='cluster')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(8.,0.,0.,220.,0.,0.,initialize=True)

	assert cluster.xgc == 8.
	assert cluster.ygc == 0.
	assert cluster.zgc == 0.
	assert cluster.rgc == 8.
	assert cluster.vxgc == 220.
	assert cluster.vygc == 0.
	assert cluster.vzgc == 0.
	assert cluster.orbit.r() == cluster.rgc

def test_add_orbit_from_centre():
	nstar=100
	cluster=ctools.StarCluster(units='kpckms',origin='cluster')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.xc,cluster.yc,cluster.zc=1.,0.0,0.
	cluster.vxc,cluster.vyc,cluster.vzc=10.,0.0,0.

	cluster.add_orbit(8.,0.,0.,220.,0.,0.,initialize=True,from_centre=True)

	assert cluster.orbit.x() == 9.
	assert cluster.orbit.vx() == 230.

def test_add_orbit_units():
	for u in units:
		nstar=100
		cluster=ctools.StarCluster(units='pckms',origin='cluster')
		x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
		vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
		cluster.add_stars(x,y,z,vx,vy,vz)

		if u == 'pckms':
			cluster.add_orbit(8.,0.,0.,220.,0.,0.,ounits=u)
			assert cluster.xgc == 8.
			assert cluster.vxgc == 220.
		elif u=='kpckms':
			cluster.add_orbit(8.,0.,0.,220.,0.,0.,ounits=u)
			assert cluster.xgc == 8000.
			assert cluster.vxgc == 220.
		elif u=='nbody':
			cluster.add_orbit(8.,0.,0.,220.,0.,0.,ounits=u)
			assert cluster.xgc == 8.
			assert cluster.vxgc == 220.
		elif u=='galpy':
			cluster.add_orbit(1.,0.,0.,1.,0.,0.,ounits=u)

			assert cluster.xgc == 8000.
			assert cluster.vxgc == 220.

def test_add_orbit_radec():
	nstar=100
	cluster=ctools.StarCluster(units='radec',origin='sky')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(8.,0.,0.,220.,0.,0.,initialize=True)

	assert cluster.ra_gc == 8.
	assert cluster.dec_gc == 0.
	assert cluster.dist_gc == 0.
	assert cluster.pmra_gc == 220.
	assert cluster.pmdec_gc == 0.
	assert cluster.vlos_gc == 0.

def test_add_nbody6():
	cluster=ctools.StarCluster()
	cluster.add_nbody6(1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.)

	assert cluster.nc == 1.
	assert cluster.rc == 1.
	assert cluster.rbar == 1.
	assert cluster.rtide == 1.
	assert cluster.xc == 1.
	assert cluster.yc == 1.
	assert cluster.zc == 1.
	assert cluster.zmbar == 1.
	assert cluster.vbar == 1.
	assert cluster.rscale == 1.
	assert cluster.ns == 1.
	assert cluster.nb == 1.
	assert cluster.np == 1.


def test_add_nbody6():
	nstar = 100
	cluster=ctools.StarCluster()
	kw = np.ones(nstar)
	logl = np.ones(nstar)
	logr = np.ones(nstar)
	ep = np.ones(nstar)
	ospin = np.ones(nstar)

	cluster.add_sse(kw, logl, logr, ep, ospin)

	np.testing.assert_array_equal(kw,cluster.kw)
	np.testing.assert_array_equal(logl,cluster.logl)
	np.testing.assert_array_equal(logr,cluster.logr)
	np.testing.assert_array_equal(10.0 ** logl,cluster.lum)
	assert cluster.ltot == np.sum(10.0 ** logl)
	np.testing.assert_array_equal(ep,cluster.ep)
	np.testing.assert_array_equal(ospin,cluster.ospin)

def test_add_energies():
	cluster=ctools.StarCluster()
	kin=np.ones(100)
	pot=np.ones(100)
	cluster.add_energies(kin,pot)
	np.testing.assert_array_equal(kin,cluster.kin)
	np.testing.assert_array_equal(pot,cluster.pot)
	np.testing.assert_array_equal(kin+pot,cluster.etot)
	assert cluster.ektot == 100.
	assert cluster.ptot == 50.
	assert cluster.qvir == 2.

def test_add_energies_etot():
	cluster=ctools.StarCluster()
	kin=np.ones(100)
	pot=np.ones(100)
	etot=np.ones(100)
	cluster.add_energies(kin,pot,etot)
	np.testing.assert_array_equal(etot,cluster.etot)

def test_add_actions():
	cluster=ctools.StarCluster()
	JR,Jphi,Jz=np.ones(100),np.ones(100),np.ones(100)
	OR,Ophi,Oz=np.ones(100),np.ones(100),np.ones(100)
	TR,Tphi,Tz=np.ones(100),np.ones(100),np.ones(100)

	cluster.add_actions(JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz)
	np.testing.assert_array_equal(JR,cluster.JR)
	np.testing.assert_array_equal(Jphi,cluster.Jphi)
	np.testing.assert_array_equal(Jz,cluster.Jz)
	np.testing.assert_array_equal(OR,cluster.OR)
	np.testing.assert_array_equal(Ophi,cluster.Ophi)
	np.testing.assert_array_equal(Oz,cluster.Oz)
	np.testing.assert_array_equal(TR,cluster.TR)
	np.testing.assert_array_equal(Tphi,cluster.Tphi)
	np.testing.assert_array_equal(Tz,cluster.Tz)

def test_analyze(tol=0.01):
	cluster=ctools.setup_cluster(ctype='king',phi0=5.,rh=3.,M=10000,N=10000)
	assert cluster.ctype == 'king'
	assert float(np.fabs(cluster.ntot-10000)/10000) <= tol

	assert np.fabs(cluster.rm-3.)/3. <= tol
	assert np.fabs(cluster.mtot-10000.0)/10000.0 <= tol

	r = np.sqrt(cluster.x ** 2.0 + cluster.y ** 2.0 + cluster.z ** 2.0)
	rpro = np.sqrt(cluster.x ** 2.0 + cluster.y ** 2.0)
	v = np.sqrt(cluster.vx ** 2.0 + cluster.vy ** 2.0 + cluster.vz ** 2.0)
	vpro = np.sqrt(cluster.vx ** 2.0 + cluster.vy ** 2.0)


	np.testing.assert_array_equal(r,cluster.r)
	np.testing.assert_array_equal(rpro,cluster.rpro)
	np.testing.assert_array_equal(v,cluster.v)
	np.testing.assert_array_equal(vpro,cluster.vpro)

	assert np.mean(r) == cluster.rmean
	assert np.amax(r) == cluster.rmax
	assert np.mean(rpro) == cluster.rpromean
	assert np.amax(rpro) == cluster.rpromax
	np.testing.assert_array_equal(np.argsort(r),cluster.rorder)
	np.testing.assert_array_equal(np.argsort(rpro),cluster.rproorder)

	assert np.fabs(cluster.rm-r[np.argsort(r)][5000])/r[np.argsort(r)][5000] <= tol
	assert np.fabs(cluster.r10-r[np.argsort(r)][1000])/r[np.argsort(r)][5000] <= tol
	assert np.fabs(cluster.rmpro-rpro[np.argsort(rpro)][5000])/rpro[np.argsort(rpro)][5000] <= tol
	assert np.fabs(cluster.r10pro-rpro[np.argsort(rpro)][1000])/rpro[np.argsort(rpro)][5000] <= tol

	kw = np.ones(10000)
	logl = np.ones(10000)
	logr = np.ones(10000)
	ep = np.ones(10000)
	ospin = np.ones(10000)
	cluster.add_sse(kw, logl, logr, ep, ospin)
	cluster.analyze()

	assert np.fabs(cluster.rh-r[np.argsort(r)][5000])/r[np.argsort(r)][5000] <= tol
	assert np.fabs(cluster.rh10-r[np.argsort(r)][1000])/r[np.argsort(r)][5000] <= tol
	assert np.fabs(cluster.rhpro-rpro[np.argsort(rpro)][5000])/rpro[np.argsort(rpro)][5000] <= tol
	assert np.fabs(cluster.rh10pro-rpro[np.argsort(rpro)][1000])/rpro[np.argsort(rpro)][5000] <= tol

def test_sortstars():
	cluster=ctools.setup_cluster(ctype='king',phi0=5.,rh=3.,M=10000,N=10000)
	np.testing.assert_array_equal(np.argsort(cluster.r),cluster.rorder)
	np.testing.assert_array_equal(np.argsort(cluster.rpro),cluster.rproorder)

def test_subcluster():
	cluster=ctools.setup_cluster(ctype='NGC6101')
	cluster.find_centre()
	cluster.reset_nbody_scale()

	subcluster=ctools.sub_cluster(cluster,rmax=cluster.rm)

	#Assert meta data is correct
	assert subcluster.tphys == cluster.tphys
	assert subcluster.units == cluster.units
	assert subcluster.origin == cluster.origin
	assert subcluster.ctype == cluster.ctype

	#Assert orbit is correct
	assert subcluster.xgc == cluster.xgc
	assert subcluster.ygc == cluster.ygc
	assert subcluster.zgc == cluster.zgc
	assert subcluster.vxgc == cluster.vxgc
	assert subcluster.vygc == cluster.vygc
	assert subcluster.vzgc == cluster.vzgc

	#Assert centre of mass is unchanged
	assert subcluster.xc == cluster.xc
	assert subcluster.yc == cluster.yc
	assert subcluster.zc == cluster.zc
	assert subcluster.vxc == cluster.vxc
	assert subcluster.vyc == cluster.vyc
	assert subcluster.vzc == cluster.vzc

	#Assert subset of stars is correct
	assert subcluster.ntot == cluster.ntot/2
	assert np.amin(subcluster.r) == np.amin(cluster.r)
	assert subcluster.r[subcluster.rorder[0]] == cluster.r[cluster.rorder[0]]
	assert subcluster.rpro[subcluster.rproorder[0]] == cluster.rpro[cluster.rproorder[0]]

	kin=np.random.rand(cluster.ntot)
	pot=np.random.rand(cluster.ntot)
	cluster.add_energies(kin,pot)
	cluster.m=np.random.rand(cluster.ntot)
	cluster.kw=np.random.rand(cluster.ntot)

	#Assert different cuts are working
	subcluster=ctools.sub_cluster(cluster,rmin=cluster.rm)
	assert np.amin(subcluster.r) == np.amin(cluster.r[cluster.r>=cluster.rm])

	subcluster=ctools.sub_cluster(cluster,vmin=np.mean(cluster.v))
	assert subcluster.ntot == np.sum(cluster.v >= np.mean(cluster.v))

	subcluster=ctools.sub_cluster(cluster,vmax=np.mean(cluster.v))
	assert subcluster.ntot == np.sum(cluster.v <= np.mean(cluster.v))

	subcluster=ctools.sub_cluster(cluster,emin=np.mean(cluster.etot))
	assert subcluster.ntot == np.sum(cluster.etot >= np.mean(cluster.etot))

	subcluster=ctools.sub_cluster(cluster,emax=np.mean(cluster.etot))
	assert subcluster.ntot == np.sum(cluster.etot <= np.mean(cluster.etot))

	subcluster=ctools.sub_cluster(cluster,mmin=np.mean(cluster.m))
	assert subcluster.ntot == np.sum(cluster.m >= np.mean(cluster.m))

	subcluster=ctools.sub_cluster(cluster,mmax=np.mean(cluster.m))
	assert subcluster.ntot == np.sum(cluster.m <= np.mean(cluster.m))

	subcluster=ctools.sub_cluster(cluster,kwmin=np.mean(cluster.kw))
	assert subcluster.ntot == np.sum(cluster.kw >= np.mean(cluster.kw))

	subcluster=ctools.sub_cluster(cluster,kwmax=np.mean(cluster.kw))
	assert subcluster.ntot == np.sum(cluster.kw <= np.mean(cluster.kw))

	indx=np.append(np.ones(int(cluster.ntot/2),dtype=bool),np.zeros(int(cluster.ntot/2),dtype=bool))
	subcluster=ctools.sub_cluster(cluster,indx=indx)
	assert subcluster.ntot == cluster.ntot/2

	#Assert projected cuts are working
	subcluster=ctools.sub_cluster(cluster,rmin=cluster.rm, projected=True)
	assert np.amin(subcluster.r) == np.amin(cluster.r[cluster.rpro>=cluster.rm])

	subcluster=ctools.sub_cluster(cluster,vmin=np.mean(cluster.v),projected=True)
	assert subcluster.ntot == np.sum(cluster.vpro >= np.mean(cluster.v))

	#Assert that Nbody conversion can be rescaled and centre recalculated
	subcluster=ctools.sub_cluster(cluster,rmin=cluster.rm,reset_nbody_scale=True,reset_centre=True)

	print(np.sum(subcluster.m),subcluster.mtot,np.sum(cluster.m[cluster.r>=cluster.rm]))

	assert subcluster.zmbar == np.sum(subcluster.m)
	assert subcluster.rbar != cluster.rbar
	assert subcluster.vbar != cluster.vbar
	assert subcluster.tbar != cluster.tbar

	assert subcluster.xc != cluster.xc
	assert subcluster.yc != cluster.yc
	assert subcluster.zc != cluster.zc
	assert subcluster.vxc != cluster.vxc
	assert subcluster.vyc != cluster.vyc
	assert subcluster.vzc != cluster.vzc

	# test ra_dec, parameters in add_nbody6 with non-default
	#test add_sse, add_bse, add_energies, analyze

