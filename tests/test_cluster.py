import clustertools as ctools
import numpy as np

import pytest

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]

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

def test_add_stars_custom():
	nstar=100
	cluster=ctools.StarCluster()
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	m=np.ones(nstar)*0.1
	id=np.linspace(0,nstar,nstar+1)
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

def test_add_stars_binaries():
	nstar=100
	nb=10
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	x[:2*nb]*=2

	cluster.add_stars(x,y,z,vx,vy,vz,nb=nb)

	assert(len(cluster.x)==nstar-nb)
	assert(len(cluster.y)==nstar-nb)
	assert(len(cluster.z)==nstar-nb)
	assert(len(cluster.vx)==nstar-nb)
	assert(len(cluster.vy)==nstar-nb)
	assert(len(cluster.vz)==nstar-nb)
	assert(len(cluster.m)==nstar-nb)

	assert(len(cluster.xb1)==nb)
	assert(len(cluster.yb1)==nb)
	assert(len(cluster.zb1)==nb)
	assert(len(cluster.vxb1)==nb)
	assert(len(cluster.vyb1)==nb)
	assert(len(cluster.vzb1)==nb)
	assert(len(cluster.xb2)==nb)
	assert(len(cluster.yb2)==nb)
	assert(len(cluster.zb2)==nb)
	assert(len(cluster.vxb2)==nb)
	assert(len(cluster.vyb2)==nb)
	assert(len(cluster.vzb2)==nb)
	assert(len(cluster.mb1)==nb)
	assert(len(cluster.mb2)==nb)

	arg1=np.arange(0,2*nb,2)
	arg2=arg1+1

	xcom=(x[arg1]+x[arg2])/2.
	ycom=(y[arg1]+y[arg2])/2.
	zcom=(z[arg1]+z[arg2])/2.
	vxcom=(vx[arg1]+vx[arg2])/2.
	vycom=(vy[arg1]+vy[arg2])/2.
	vzcom=(vz[arg1]+vz[arg2])/2.

	np.testing.assert_array_equal(xcom,cluster.x[:nb])
	np.testing.assert_array_equal(np.ones(nb)*2,cluster.x[:nb])
	np.testing.assert_array_equal(np.ones(nstar-2*nb),cluster.x[nb:])

	np.testing.assert_array_equal(ycom,cluster.y[:nb])
	np.testing.assert_array_equal(zcom,cluster.z[:nb])

	np.testing.assert_array_equal(vxcom,cluster.vx[:nb])
	np.testing.assert_array_equal(vycom,cluster.vy[:nb])
	np.testing.assert_array_equal(vzcom,cluster.vz[:nb])

	assert(len(np.unique(cluster.id))==nstar-nb)
	assert np.amax(cluster.id)==nstar-1


	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	x[:2*nb]*=2

	arg1=np.arange(0,2*nb-1,2)
	arg2=arg1+1
	args=arg2[-1]+1

	cluster.add_stars(x[args:],y[args:],z[args:],vx[args:],vy[args:],vz[args:])
	cluster.add_binary_stars(x[arg1],y[arg1],z[arg1],vx[arg1],vy[arg1],vz[arg1],x[arg2],y[arg2],z[arg2],vx[arg2],vy[arg2],vz[arg2])

	assert(len(cluster.x)==nstar-nb)
	assert(len(cluster.y)==nstar-nb)
	assert(len(cluster.z)==nstar-nb)
	assert(len(cluster.vx)==nstar-nb)
	assert(len(cluster.vy)==nstar-nb)
	assert(len(cluster.vz)==nstar-nb)
	assert(len(cluster.m)==nstar-nb)

	assert(len(cluster.xb1)==nb)
	assert(len(cluster.yb1)==nb)
	assert(len(cluster.zb1)==nb)
	assert(len(cluster.vxb1)==nb)
	assert(len(cluster.vyb1)==nb)
	assert(len(cluster.vzb1)==nb)
	assert(len(cluster.xb2)==nb)
	assert(len(cluster.yb2)==nb)
	assert(len(cluster.zb2)==nb)
	assert(len(cluster.vxb2)==nb)
	assert(len(cluster.vyb2)==nb)
	assert(len(cluster.vzb2)==nb)
	assert(len(cluster.mb1)==nb)
	assert(len(cluster.mb2)==nb)

	arg1=np.arange(0,2*nb,2)
	arg2=arg1+1

	xcom=(x[arg1]+x[arg2])/2.
	ycom=(y[arg1]+y[arg2])/2.
	zcom=(z[arg1]+z[arg2])/2.
	vxcom=(vx[arg1]+vx[arg2])/2.
	vycom=(vy[arg1]+vy[arg2])/2.
	vzcom=(vz[arg1]+vz[arg2])/2.

	np.testing.assert_array_equal(xcom,cluster.x[:nb])
	np.testing.assert_array_equal(np.ones(nb)*2,cluster.x[:nb])
	np.testing.assert_array_equal(np.ones(nstar-2*nb),cluster.x[nb:])

	np.testing.assert_array_equal(ycom,cluster.y[:nb])
	np.testing.assert_array_equal(zcom,cluster.z[:nb])

	np.testing.assert_array_equal(vxcom,cluster.vx[:nb])
	np.testing.assert_array_equal(vycom,cluster.vy[:nb])
	np.testing.assert_array_equal(vzcom,cluster.vz[:nb])

	assert(len(np.unique(cluster.id))==nstar-nb)

	print(cluster.id,cluster.ntot,cluster.nb)

	assert np.amax(cluster.id)==nstar-2

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	x[:2*nb]*=2

	arg1=np.arange(0,2*nb-1,2)
	arg2=arg1+1
	args=arg2[-1]+1

	cluster.add_binary_stars(x[arg1],y[arg1],z[arg1],vx[arg1],vy[arg1],vz[arg1],x[arg2],y[arg2],z[arg2],vx[arg2],vy[arg2],vz[arg2])
	cluster.add_stars(x[args:],y[args:],z[args:],vx[args:],vy[args:],vz[args:])

	assert(len(cluster.x)==nstar-nb)
	assert(len(cluster.y)==nstar-nb)
	assert(len(cluster.z)==nstar-nb)
	assert(len(cluster.vx)==nstar-nb)
	assert(len(cluster.vy)==nstar-nb)
	assert(len(cluster.vz)==nstar-nb)
	assert(len(cluster.m)==nstar-nb)

	assert(len(cluster.xb1)==nb)
	assert(len(cluster.yb1)==nb)
	assert(len(cluster.zb1)==nb)
	assert(len(cluster.vxb1)==nb)
	assert(len(cluster.vyb1)==nb)
	assert(len(cluster.vzb1)==nb)
	assert(len(cluster.xb2)==nb)
	assert(len(cluster.yb2)==nb)
	assert(len(cluster.zb2)==nb)
	assert(len(cluster.vxb2)==nb)
	assert(len(cluster.vyb2)==nb)
	assert(len(cluster.vzb2)==nb)
	assert(len(cluster.mb1)==nb)
	assert(len(cluster.mb2)==nb)

	arg1=np.arange(0,2*nb,2)
	arg2=arg1+1

	xcom=(x[arg1]+x[arg2])/2.
	ycom=(y[arg1]+y[arg2])/2.
	zcom=(z[arg1]+z[arg2])/2.
	vxcom=(vx[arg1]+vx[arg2])/2.
	vycom=(vy[arg1]+vy[arg2])/2.
	vzcom=(vz[arg1]+vz[arg2])/2.

	np.testing.assert_array_equal(xcom,cluster.x[:nb])
	np.testing.assert_array_equal(np.ones(nb)*2,cluster.x[:nb])
	np.testing.assert_array_equal(np.ones(nstar-2*nb),cluster.x[nb:])

	np.testing.assert_array_equal(ycom,cluster.y[:nb])
	np.testing.assert_array_equal(zcom,cluster.z[:nb])

	np.testing.assert_array_equal(vxcom,cluster.vx[:nb])
	np.testing.assert_array_equal(vycom,cluster.vy[:nb])
	np.testing.assert_array_equal(vzcom,cluster.vz[:nb])

	assert(len(np.unique(cluster.id))==nstar-nb)

	print(cluster.id,cluster.ntot,cluster.nb)

	assert np.amax(cluster.id)==nstar-1

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

def test_add_orbit_from_centre(tol=0.01):
	nstar=100
	cluster=ctools.StarCluster(units='kpckms',origin='cluster')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.xc,cluster.yc,cluster.zc=1.,0.0,0.
	cluster.vxc,cluster.vyc,cluster.vzc=10.,0.0,0.

	cluster.add_orbit(solar_ro,0.,0.,solar_vo,0.,0.,initialize=True,from_centre=True)

	assert np.fabs(cluster.orbit.x() / (solar_ro+1)) - 1. <=tol
	assert cluster.orbit.vx() == solar_vo+10

def test_add_orbit_units():
	for u in units:
		nstar=100
		cluster=ctools.StarCluster(units='pckms',origin='cluster')
		x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
		vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
		cluster.add_stars(x,y,z,vx,vy,vz)

		if u == 'pckms':
			cluster.add_orbit(solar_ro,0.,0.,solar_vo,0.,0.,ounits=u)
			assert cluster.xgc == solar_ro
			assert cluster.vxgc == solar_vo
		elif u=='kpckms':
			cluster.add_orbit(solar_ro,0.,0.,solar_vo,0.,0.,ounits=u)
			assert cluster.xgc == solar_ro*1000.0
			assert cluster.vxgc == solar_vo
		elif u=='nbody':
			cluster.add_orbit(solar_ro,0.,0.,solar_vo,0.,0.,ounits=u)
			assert cluster.xgc == solar_ro
			assert cluster.vxgc == solar_vo
		elif u=='galpy':
			cluster.add_orbit(1.,0.,0.,1.,0.,0.,ounits=u)

			assert cluster.xgc == solar_ro*1000.0
			assert cluster.vxgc == solar_vo

def test_add_orbit_radec():
	nstar=100
	cluster=ctools.StarCluster(units='radec',origin='sky')
	x,y,z=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	vx,vy,vz=np.ones(nstar),np.ones(nstar),np.ones(nstar)
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(solar_ro,0.,0.,solar_vo,0.,0.,initialize=True)

	assert cluster.ra_gc == solar_ro
	assert cluster.dec_gc == 0.
	assert cluster.dist_gc == 0.
	assert cluster.pmra_gc == solar_vo
	assert cluster.pmdec_gc == 0.
	assert cluster.vlos_gc == 0.

def test_add_nbody6():
	cluster=ctools.StarCluster()
	cluster.add_nbody6(1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.)

	assert cluster.nc == 1.
	assert cluster.rc == 1.
	assert cluster.rbar == 1.
	assert cluster.rtide == 1.
	assert cluster.xc == 1.
	assert cluster.yc == 1.
	assert cluster.zc == 1.
	assert cluster.zmbar == 1.
	assert cluster.vbar == 1.
	assert cluster.tbar == 1.
	assert cluster.rscale == 1.
	assert cluster.ns == 1.
	assert cluster.nb == 1.
	assert cluster.n_p == 1.


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

def test_add_action():
	cluster=ctools.StarCluster()
	JR,Jphi,Jz=1.,1.,1.
	OR,Ophi,Oz=1.,1.,1.
	TR,Tphi,Tz=1.,1.,1.

	cluster.add_action(JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz)
	np.testing.assert_array_equal(JR,cluster.JR)
	np.testing.assert_array_equal(Jphi,cluster.Jphi)
	np.testing.assert_array_equal(Jz,cluster.Jz)
	np.testing.assert_array_equal(OR,cluster.OR)
	np.testing.assert_array_equal(Ophi,cluster.Ophi)
	np.testing.assert_array_equal(Oz,cluster.Oz)
	np.testing.assert_array_equal(TR,cluster.TR)
	np.testing.assert_array_equal(Tphi,cluster.Tphi)
	np.testing.assert_array_equal(Tz,cluster.Tz)

def test_add_actions():
	cluster=ctools.StarCluster()
	JR,Jphi,Jz=np.ones(100),np.ones(100),np.ones(100)
	OR,Ophi,Oz=np.ones(100),np.ones(100),np.ones(100)
	TR,Tphi,Tz=np.ones(100),np.ones(100),np.ones(100)

	cluster.add_actions(JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz)
	np.testing.assert_array_equal(JR,cluster.JRs)
	np.testing.assert_array_equal(Jphi,cluster.Jphis)
	np.testing.assert_array_equal(Jz,cluster.Jzs)
	np.testing.assert_array_equal(OR,cluster.ORs)
	np.testing.assert_array_equal(Ophi,cluster.Ophis)
	np.testing.assert_array_equal(Oz,cluster.Ozs)
	np.testing.assert_array_equal(TR,cluster.TRs)
	np.testing.assert_array_equal(Tphi,cluster.Tphis)
	np.testing.assert_array_equal(Tz,cluster.Tzs)

def test_analyze(tol=0.01):

	cluster=ctools.load_cluster('snapshot',filename='g1phi5rh3m10000.dat',origin='cluster')

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
	assert np.mean(rpro) == cluster.rmeanpro
	assert np.amax(rpro) == cluster.rmaxpro
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
	cluster=ctools.load_cluster('snapshot',filename='g1phi5rh3m10000.dat',origin='cluster')
	np.testing.assert_array_equal(np.argsort(cluster.r),cluster.rorder)
	np.testing.assert_array_equal(np.argsort(cluster.rpro),cluster.rproorder)

def test_subset(tol=0.001):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)

	cluster.to_cluster()

	kin=np.random.rand(cluster.ntot)
	pot=np.random.rand(cluster.ntot)
	cluster.add_energies(kin,pot)
	cluster.m=np.random.rand(cluster.ntot)
	cluster.kw=np.random.rand(cluster.ntot)
	cluster.npop=(np.ones(cluster.ntot)*2).astype(int)
	cluster.analyze()

	#Assert different cuts are working

	indx=cluster.subset(rmin=cluster.rm)
	assert np.amin(cluster.r[indx]) == np.amin(cluster.r[cluster.r>=cluster.rm])

	indx=cluster.subset(vmin=np.mean(cluster.v))
	assert np.sum(indx) == np.sum(cluster.v >= np.mean(cluster.v))

	indx=cluster.subset(vmax=np.mean(cluster.v))
	assert np.sum(indx) == np.sum(cluster.v <= np.mean(cluster.v))

	indx=cluster.subset(emin=np.mean(cluster.etot))
	assert np.sum(indx) == np.sum(cluster.etot >= np.mean(cluster.etot))

	indx=cluster.subset(emax=np.mean(cluster.etot))
	assert np.sum(indx) == np.sum(cluster.etot <= np.mean(cluster.etot))

	indx=cluster.subset(mmin=np.mean(cluster.m))
	assert np.sum(indx) == np.sum(cluster.m >= np.mean(cluster.m))

	indx=cluster.subset(mmax=np.mean(cluster.m))
	assert np.sum(indx) == np.sum(cluster.m <= np.mean(cluster.m))

	indx=cluster.subset(kwmin=np.mean(cluster.kw))
	assert np.sum(indx) == np.sum(cluster.kw >= np.mean(cluster.kw))

	indx=cluster.subset(kwmax=np.mean(cluster.kw))
	assert np.sum(indx) == np.sum(cluster.kw <= np.mean(cluster.kw))

	indx=np.append(np.ones(int(cluster.ntot/2),dtype=bool),np.zeros(int(cluster.ntot/2),dtype=bool))
	newindx=cluster.subset(indx=indx)
	assert np.sum(newindx) == cluster.ntot/2

	indx=cluster.subset(npop=1)
	assert np.sum(indx)==0
	indx=cluster.subset(npop=2)
	assert np.sum(indx)==cluster.ntot

	#Assert projected cuts are working
	indx=cluster.subset(rmin=cluster.rm, projected=True)
	assert np.amin(cluster.r[indx]) == np.amin(cluster.r[cluster.rpro>=cluster.rm])

	indx=cluster.subset(vmin=np.mean(cluster.v),projected=True)
	assert np.sum(indx) == np.sum(cluster.vpro >= np.mean(cluster.v))
	
def test_subcluster():
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.to_radec()

	cluster.to_cluster(sortstars=True)

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
	assert subcluster.ntot == np.sum(cluster.r <= cluster.rm)
	assert np.amin(subcluster.r) == np.amin(cluster.r)
	assert subcluster.r[subcluster.rorder[0]] == cluster.r[cluster.rorder[0]]

	kin=np.random.rand(cluster.ntot)
	pot=np.random.rand(cluster.ntot)
	cluster.add_energies(kin,pot)
	cluster.m=np.random.rand(cluster.ntot)
	cluster.kw=np.random.rand(cluster.ntot)
	cluster.analyze()

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
	subcluster=ctools.sub_cluster(cluster,rmin=cluster.rm,reset_nbody=True,reset_centre=True)

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

	#Test radec are properly transferred
	indx=cluster.r < cluster.rm
	subcluster=ctools.sub_cluster(cluster,indx=indx)

	np.testing.assert_array_equal(subcluster.ra,cluster.ra[indx])
	np.testing.assert_array_equal(subcluster.dec,cluster.dec[indx])
	np.testing.assert_array_equal(subcluster.dist,cluster.dist[indx])
	np.testing.assert_array_equal(subcluster.pmra,cluster.pmra[indx])
	np.testing.assert_array_equal(subcluster.pmdec,cluster.pmdec[indx])
	np.testing.assert_array_equal(subcluster.vlos,cluster.vlos[indx])

	assert subcluster.ra_gc==cluster.ra_gc
	assert subcluster.dec_gc==cluster.dec_gc
	assert subcluster.dist_gc==cluster.dist_gc
	assert subcluster.pmra_gc==cluster.pmra_gc
	assert subcluster.pmdec_gc==cluster.pmdec_gc
	assert subcluster.vlos_gc==cluster.vlos_gc

	cluster.add_nbody6(nc=100,rc=cluster.r10,rbar=2.,rtide=100.,xc=cluster.xc, yc=cluster.yc,zc=cluster.zc,zmbar=cluster.mtot,vbar=5.,tbar=2.,rscale=2.,ns=cluster.ntot,nb=2,n_p=1.)
	subcluster=ctools.sub_cluster(cluster,indx=indx)

	assert subcluster.nc==cluster.nc
	assert subcluster.rc==cluster.rc
	assert subcluster.rbar==cluster.rbar
	assert subcluster.rtide==cluster.rtide
	assert subcluster.xc==cluster.xc
	assert subcluster.yc==cluster.yc
	assert subcluster.zc==cluster. zc
	assert subcluster.zmbar==cluster.zmbar
	assert subcluster.vbar==cluster.vbar
	assert subcluster.tbar==cluster.tbar
	assert subcluster.rscale==cluster.rscale
	assert subcluster.ns==cluster.ns
	assert subcluster.nb==cluster.nb
	assert subcluster.n_p==cluster.n_p

	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_pckms_cluster.dat',units='pckms',origin='cluster')
	cluster.add_orbit(-2005.2100994789871, -9348.1814843660959, -3945.4681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.to_cluster(sortstars=True)

	kw=np.random.randint(0,10,cluster.ntot)
	logl=np.random.rand(cluster.ntot)
	logr=np.random.rand(cluster.ntot)
	ep=np.random.rand(cluster.ntot)
	ospin=np.random.rand(cluster.ntot)

	cluster.add_sse(kw,logl,logr,ep,ospin)
	btot=100
	#Manually set ids
	id1=2*np.linspace(1,btot,btot,dtype=int)-1
	id2=2*np.linspace(1,btot,btot,dtype=int)

	cluster.id=np.append(id1,np.linspace(2*btot+1,len(kw),len(kw)-btot,dtype=int))

	kw1=np.random.randint(0,10,btot)
	kw2=np.random.randint(0,10,btot)
	kcm=np.random.randint(0,10,btot)
	ecc=np.random.rand(btot)
	pb=np.random.rand(btot)
	semi=np.random.rand(btot)
	m1=np.random.rand(btot)
	m2=np.random.rand(btot)
	logl1=np.random.rand(btot)
	logl2=np.random.rand(btot)
	logr1=np.random.rand(btot)
	logr2=np.random.rand(btot)
	ep1=np.random.rand(btot)
	ep2=np.random.rand(btot)
	ospin1=np.random.rand(btot)
	ospin2=np.random.rand(btot)

	cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,
		logr1,logr2,ep1,ep2,ospin1,ospin2)

	indx=cluster.r<=cluster.rm
	bindx=cluster.r[0:btot] <=cluster.rm 

	subcluster=ctools.sub_cluster(cluster,rmax=cluster.rm)

	print(np.sum(indx),np.sum(bindx),subcluster.ntot,len(subcluster.id1))

	np.testing.assert_array_equal(subcluster.kw,cluster.kw[indx])
	np.testing.assert_array_equal(subcluster.logl,cluster.logl[indx])
	np.testing.assert_array_equal(subcluster.logr,cluster.logr[indx])
	np.testing.assert_array_equal(subcluster.ep,cluster.ep[indx])
	np.testing.assert_array_equal(subcluster.ospin,cluster.ospin[indx])
	np.testing.assert_array_equal(subcluster.id1,cluster.id1[bindx])
	np.testing.assert_array_equal(subcluster.id2,cluster.id2[bindx])
	np.testing.assert_array_equal(subcluster.kw1,cluster.kw1[bindx])
	np.testing.assert_array_equal(subcluster.kw2,cluster.kw2[bindx])
	np.testing.assert_array_equal(subcluster.kcm,cluster.kcm[bindx])
	np.testing.assert_array_equal(subcluster.ecc,cluster.ecc[bindx])
	np.testing.assert_array_equal(subcluster.pb,cluster.pb[bindx])
	np.testing.assert_array_equal(subcluster.semi,cluster.semi[bindx])
	np.testing.assert_array_equal(subcluster.m1,cluster.m1[bindx])
	np.testing.assert_array_equal(subcluster.m2,cluster.m2[bindx])
	np.testing.assert_array_equal(subcluster.logl1,cluster.logl1[bindx])
	np.testing.assert_array_equal(subcluster.logl2,cluster.logl2[bindx])
	np.testing.assert_array_equal(subcluster.logr1,cluster.logr1[bindx])
	np.testing.assert_array_equal(subcluster.logr2,cluster.logr2[bindx])
	np.testing.assert_array_equal(subcluster.ep1,cluster.ep1[bindx])
	np.testing.assert_array_equal(subcluster.ep2,cluster.ep2[bindx])
	np.testing.assert_array_equal(subcluster.ospin1,cluster.ospin1[bindx])
	np.testing.assert_array_equal(subcluster.ospin2,cluster.ospin2[bindx])


