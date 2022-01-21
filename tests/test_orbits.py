import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential
try:
	from galpy.util import conversion
except:
	import galpy.util.bovy_conversion as conversion

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

def test_initialize_orbits(tol=0.0001):
	
	x=np.random.rand(1000)
	y=np.random.rand(1000)
	z=np.random.rand(1000)
	vx=np.random.rand(1000)
	vy=np.random.rand(1000)
	vz=np.random.rand(1000)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=1.,sortstars=True,analyze=True)

	cluster.add_orbit(8000.,0.,0.,0.,220.,0.)
	cluster.find_centre()
	cluster.to_galaxy()
	cluster.to_kpckms()

	ocluster=cluster.initialize_orbits()

	np.testing.assert_allclose(8.0+x/1000.,ocluster.x(),rtol=tol)
	np.testing.assert_allclose(y/1000.,ocluster.y(),rtol=tol)
	np.testing.assert_allclose(z/1000.,ocluster.z(),rtol=tol)

	np.testing.assert_allclose(vx,ocluster.vx(),rtol=tol)
	np.testing.assert_allclose(220.+vy,ocluster.vy(),rtol=tol)
	np.testing.assert_allclose(vz,ocluster.vz(),rtol=tol)

	cluster.to_centre()
	print(cluster.xc,cluster.yc,cluster.zc)
	ocluster=cluster.initialize_orbits()

	np.testing.assert_allclose(x/1000.,ocluster.x()+cluster.xc,rtol=tol)
	np.testing.assert_allclose(y/1000.,ocluster.y()+cluster.yc,rtol=tol)
	np.testing.assert_allclose(z/1000.,ocluster.z()+cluster.zc,rtol=tol)

	np.testing.assert_allclose(vx,ocluster.vx()+cluster.vxc,rtol=tol)
	np.testing.assert_allclose(vy,ocluster.vy()+cluster.vyc,rtol=tol)
	np.testing.assert_allclose(vz,ocluster.vz()+cluster.vzc,rtol=tol)

	cluster.to_cluster()
	ocluster=cluster.initialize_orbits()

	np.testing.assert_allclose(x/1000.,ocluster.x(),rtol=tol)
	np.testing.assert_allclose(y/1000.,ocluster.y(),rtol=tol)
	np.testing.assert_allclose(z/1000.,ocluster.z(),rtol=tol)

	np.testing.assert_allclose(vx,ocluster.vx(),rtol=tol)
	np.testing.assert_allclose(vy,ocluster.vy(),rtol=tol)
	np.testing.assert_allclose(vz,ocluster.vz(),rtol=tol)

def test_integrate_orbit(tol=0.1,ro=8.,vo=220.):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy')
	ocluster=cluster.integrate_orbit(tfinal=1.,nt=1000)

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo))

	assert cluster.xgc == o.x()
	assert cluster.xgc == ocluster.x(ts[0])

	o.integrate(ts,MWPotential2014)

	assert np.fabs(ocluster.x(ts[-1])-o.x(ts[-1])) <= tol
	assert np.fabs(ocluster.y(ts[-1])-o.y(ts[-1])) <= tol
	assert np.fabs(ocluster.z(ts[-1])-o.z(ts[-1])) <= tol
	assert np.fabs(ocluster.vx(ts[-1])-o.vx(ts[-1])) <= tol
	assert np.fabs(ocluster.vy(ts[-1])-o.vy(ts[-1])) <= tol
	assert np.fabs(ocluster.vz(ts[-1])-o.vz(ts[-1])) <= tol

def test_integrate_orbits(tol=0.0001,ro=8,vo=220.):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',mbar=10)
	cluster.to_galaxy()
	cluster.to_kpckms()

	ocluster=cluster.integrate_orbits(tfinal=1,nt=1000)

	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)
	vxvv=np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])

	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),1000)
	o.integrate(ts,MWPotential2014)

	np.testing.assert_allclose(ocluster.x(ts[-1]),o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.y(ts[-1]),o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.z(ts[-1]),o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.vx(ts[-1]),o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.vy(ts[-1]),o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.vz(ts[-1]),o.vz(ts[-1]),rtol=tol)

def test_interpolate_orbit(tol=0.1,ro=8.,vo=220.):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy')
	x,y,z,vx,vy,vz=cluster.interpolate_orbit(tfinal=1.,nt=1000)

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo))
	o.integrate(ts,MWPotential2014)


	assert np.fabs(cluster.xgc-o.x(ts[-1])) <= tol
	assert np.fabs(cluster.ygc-o.y(ts[-1])) <= tol
	assert np.fabs(cluster.zgc-o.z(ts[-1])) <= tol

	assert np.fabs(cluster.vxgc-o.vx(ts[-1])) <= tol
	assert np.fabs(cluster.vygc-o.vy(ts[-1])) <= tol
	assert np.fabs(cluster.vzgc-o.vz(ts[-1])) <= tol

def test_interpolate_orbits(tol=0.1,ro=8.,vo=220.):

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy')
	xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.interpolate_orbit(tfinal=1.,nt=1000)

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy',mbar=10)
	
	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)

	x,y,z,vx,vy,vz=cluster.interpolate_orbits(tfinal=1.,nt=1000)
	np.testing.assert_allclose(cluster.x,x,rtol=tol)
	np.testing.assert_allclose(cluster.y,y,rtol=tol)
	np.testing.assert_allclose(cluster.z,z,rtol=tol)
	np.testing.assert_allclose(cluster.vx,vx,rtol=tol)
	np.testing.assert_allclose(cluster.vy,vy,rtol=tol)
	np.testing.assert_allclose(cluster.vz,vz,rtol=tol)

	vxvv=np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo))
	o.integrate(ts,MWPotential2014)

	np.testing.assert_allclose(cluster.x,o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.y,o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.z,o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vx,o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vy,o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vz,o.vz(ts[-1]),rtol=tol)

	assert np.fabs(xgc-cluster.xgc) <= tol
	assert np.fabs(ygc-cluster.ygc) <= tol
	assert np.fabs(zgc-cluster.zgc) <= tol
	assert np.fabs(vxgc-cluster.vxgc) <= tol
	assert np.fabs(vygc-cluster.vygc) <= tol
	assert np.fabs(vzgc-cluster.vzgc) <= tol

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='cluster',mbar=10)
	xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc
	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)

	cluster.to_galpy()
	pot=PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
	cluster.to_kpckms()

	vxvv=np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=8.,vo=220.)
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),1000)
	o.integrate(ts,pot)

	cluster.interpolate_orbits(tfinal=1,nt=1000)

	np.testing.assert_allclose(cluster.x,o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.y,o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.z,o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vx,o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vy,o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vz,o.vz(ts[-1]),rtol=tol)

	assert np.fabs(xgc-cluster.xgc) <= tol
	assert np.fabs(ygc-cluster.ygc) <= tol
	assert np.fabs(zgc-cluster.zgc) <= tol
	assert np.fabs(vxgc-cluster.vxgc) <= tol
	assert np.fabs(vygc-cluster.vygc) <= tol
	assert np.fabs(vzgc-cluster.vzgc) <= tol	


	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='pckms',origin='centre',mbar=10)
	xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc

	print(cluster.xgc,cluster.ygc,cluster.zgc)
	print(cluster.xc,cluster.yc,cluster.zc)

	cluster.to_kpckms()
	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)
	cluster.to_pckms()

	cluster.to_galpy()
	print(cluster.mtot,cluster.rm)
	pot=PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
	cluster.to_pckms()

	vxvv=np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=8.,vo=220.)
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),1000)
	o.integrate(ts,pot)

	cluster.interpolate_orbits(tfinal=1000,nt=1000)

	np.testing.assert_allclose(cluster.x/1000,o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.y/1000,o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.z/1000,o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vx,o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vy,o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vz,o.vz(ts[-1]),rtol=tol)

	assert np.fabs(xgc-cluster.xgc) <= tol
	assert np.fabs(ygc-cluster.ygc) <= tol
	assert np.fabs(zgc-cluster.zgc) <= tol
	assert np.fabs(vxgc-cluster.vxgc) <= tol
	assert np.fabs(vygc-cluster.vygc) <= tol
	assert np.fabs(vzgc-cluster.vzgc) <= tol

def test_orbital_path(tol=0.1):
	pass

def test_orbital_path_match(tol=0.1):
	pass

def test_calc_actions(tol=0.1):
	#Test internal and external
	pass
	
def test_ttensor(tol=0.1):
	#Test internal and external
	pass