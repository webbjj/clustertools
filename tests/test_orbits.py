import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,IsochronePotential,KeplerPotential
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

def test_orbital_path(tol=0.1,ro=8.,vo=220.):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy',mbar=10)

	tfinal=0.1

	t, x, y, z, vx, vy, vz=cluster.orbital_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=False,
	    skypath=False,
	    initialize=True,
	    ro=8.0,
	    vo=220.0,
	    solarmotion=[-11.1, 24.0, 7.25],
	    plot=False,
	)

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts = np.linspace(0, -1.0 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[0]-o.x(ts[-1])) <= tol
	assert np.fabs(y[0]-o.y(ts[-1])) <= tol
	assert np.fabs(z[0]-o.z(ts[-1])) <= tol
	assert np.fabs(vx[0]-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[0]-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[0]-o.vz(ts[-1])) <= tol

	assert np.fabs(cluster.xpath[0]-o.x(ts[-1])) <= tol
	assert np.fabs(cluster.ypath[0]-o.y(ts[-1])) <= tol
	assert np.fabs(cluster.zpath[0]-o.z(ts[-1])) <= tol
	assert np.fabs(cluster.vxpath[0]-o.vx(ts[-1])) <= tol
	assert np.fabs(cluster.vypath[0]-o.vy(ts[-1])) <= tol
	assert np.fabs(cluster.vzpath[0]-o.vz(ts[-1])) <= tol

	assert np.fabs(cluster.orbit.x(ts[-1])-o.x(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.y(ts[-1])-o.y(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.z(ts[-1])-o.z(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.vx(ts[-1])-o.vx(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.vy(ts[-1])-o.vy(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.vz(ts[-1])-o.vz(ts[-1])) <= tol

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts = np.linspace(0, 1.0 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[-1]-o.x(ts[-1])) <= tol
	assert np.fabs(y[-1]-o.y(ts[-1])) <= tol
	assert np.fabs(z[-1]-o.z(ts[-1])) <= tol
	assert np.fabs(vx[-1]-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[-1]-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[-1]-o.vz(ts[-1])) <= tol

	#Check skypath
	t, x, y, z, vx, vy, vz=cluster.orbital_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=False,
	    skypath=True,
	    initialize=True,
	    ro=8.0,
	    vo=220.0,
	    solarmotion=[-11.1, 24.0, 7.25],
	    plot=False,
	)

	assert np.fabs(x[-1]-o.ra(ts[-1])) <= tol
	assert np.fabs(y[-1]-o.dec(ts[-1])) <= tol
	assert np.fabs(z[-1]-o.dist(ts[-1])) <= tol
	assert np.fabs(vx[-1]-o.pmra(ts[-1])) <= tol
	assert np.fabs(vy[-1]-o.pmdec(ts[-1])) <= tol
	assert np.fabs(vz[-1]-o.vlos(ts[-1])) <= tol

	#Check from_centre:

	t, x, y, z, vx, vy, vz=cluster.orbital_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=True,
	    skypath=False,
	    initialize=False,
	    ro=8.0,
	    vo=220.0,
	    solarmotion=[-11.1, 24.0, 7.25],
	    plot=False,
	)

	xgc,ygc,zgc=cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc
	vxgc,vygc,vzgc=cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc

	rad,phi,zed,vR,vT,vzed=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi]

	o=Orbit(vxvv,ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts = np.linspace(0, -1.0 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[0]-o.x(ts[-1])) <= tol
	assert np.fabs(y[0]-o.y(ts[-1])) <= tol
	assert np.fabs(z[0]-o.z(ts[-1])) <= tol
	assert np.fabs(vx[0]-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[0]-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[0]-o.vz(ts[-1])) <= tol

def test_orbital_path_match(tol=0.1,ro=8.,vo=220.):

	tfinal=0.1
	nt=1000

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts = np.linspace(0, 0.5 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), nt)
	o.integrate(ts,MWPotential2014)

	x,y,z=o.x(ts),o.y(ts),o.z(ts)
	vx,vy,vz=o.vx(ts),o.vy(ts),o.vz(ts)

	cluster=ctools.StarCluster(units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,sortstars=True,analyze=True)
	cluster.add_orbit(o.x(),o.y(),o.z(),o.vx(),o.vy(),o.vz(),ounits='kpckms')

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts = np.linspace(0, -0.5 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), nt)
	o.integrate(ts,MWPotential2014)

	x,y,z=o.x(ts),o.y(ts),o.z(ts)
	vx,vy,vz=o.vx(ts),o.vy(ts),o.vz(ts)
	cluster.add_stars(x,y,z,vx,vy,vz,sortstars=True,analyze=True)

	cluster.xc,cluster.yc,cluster.zc=0.,0.,0.
	cluster.vxc,cluster.vyc,cluster.vzc=0.,0.,0.

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=True,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=False,
    skypath=True,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=True,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	n=1000
	x=np.linspace(-10,10,n)
	y,z=np.zeros(n),np.zeros(n)
	vx,vy,vz=np.zeros(n),np.zeros(n),np.zeros(n)
	m=np.ones(n)

	cluster=ctools.StarCluster(units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(0.,0.,0.,0.,0.,0.)
	cluster.xc,cluster.yc,cluster.zc=0.,0.,0.
	cluster.vxc,cluster.vyc,cluster.vzc=0.,0.,0.

	tpath=np.linspace(-1,1,1000)
	xpath=np.linspace(-20,20,1000)
	ypath,zpath=np.zeros(n),np.zeros(n)
	vxpath,vypath,vzpath=np.zeros(n),np.zeros(n),np.zeros(n)

	path=[tpath,xpath,ypath,zpath,vxpath,vypath,vzpath]

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=path,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=path,
    from_centre=False,
    skypath=False,
    to_path=True,
    do_full=False,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.orbital_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=path,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=True,
    ro=8.0,
    vo=220.0,
    solarmotion=[-11.1, 24.0, 7.25],
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)

def test_calc_action(tol=0.1,ro=8.,vo=220.):

	mo=conversion.mass_in_msol(ro=ro,vo=vo)

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',mbar=100)
	cluster.to_galaxy()
	cluster.to_kpckms()

	o=Orbit.from_name('NGC6101',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	atype = "staeckel"
	delta = 0.45
	c = True
	pot=MWPotential2014

	J_R = o.jr(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo)
	J_phi = o.jp(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo)
	J_z = o.jz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo)


	jr,jp,jz=cluster.calc_action(pot=pot, ro=8.0, vo=220.0,solarmotion=[-11.1, 24.0, 7.25],full=False)

	assert(np.fabs(jr-J_R)<tol)
	assert(np.fabs(jp-J_phi)<tol)
	assert(np.fabs(jz-J_z)<tol)

	ppot=PlummerPotential(1e13/mo, 1.0/ro, ro=ro,vo=vo)
	o=Orbit.from_name('NGC6101',ro=ro,vo=vo,solarmotion=[-11.1, 24.0, 7.25])
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),100)
	o.integrate(ts,ppot)

	jr,jp,jz=cluster.calc_action(pot=ppot, ro=ro, vo=vo,solarmotion=[-11.1, 24.0, 7.25],full=False)

	assert(np.fabs(jr-o.jr())<tol)
	assert(np.fabs(jp-o.jp())<tol)
	assert(np.fabs(jz-o.jz())<tol)

	jr,jp,jz,OR, Ophi, Oz, TR, Tphi, Tz=cluster.calc_action(pot=ppot, ro=ro, vo=vo,solarmotion=[-11.1, 24.0, 7.25],full=True)

	assert(np.fabs(OR-o.Or())<tol)
	assert(np.fabs(Ophi-o.Op())<tol)
	assert(np.fabs(Oz-o.Oz())<tol)
	assert(np.fabs(TR-o.Tr())<tol)
	assert(np.fabs(Tphi-o.Tp())<tol)
	assert(np.fabs(Tz-o.Tz())<tol)

def test_calc_actions(tol=0.1,ro=8.,vo=220.):

	mo=conversion.mass_in_msol(ro=ro,vo=vo)
	pot=MWPotential2014

	o=Orbit.from_name('MWglobularclusters',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	atype = "staeckel"
	delta = 0.45
	c = True

	J_R = o.jr(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo)
	J_phi = o.jp(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo)
	J_z = o.jz(pot=pot, type=atype, delta=delta, c=c, ro=ro, vo=vo)

	cluster=ctools.load_cluster(ctype='galpy',particles=o,units='kpckms',origin='galaxy')
	jr,jp,jz=cluster.calc_actions(pot=pot, ro=ro, vo=vo,solarmotion=[-11.1, 24.0, 7.25],full=False)


	np.testing.assert_allclose(jr,J_R,rtol=tol,atol=1.)
	np.testing.assert_allclose(jp,J_phi,rtol=tol,atol=1.)
	np.testing.assert_allclose(jz,J_z,rtol=tol,atol=1.)


	ppot=PlummerPotential(1e13/mo, 1.0/ro, ro=ro,vo=vo)
	o=Orbit.from_name('MWglobularclusters',ro=8.,vo=220.,solarmotion=[-11.1, 24.0, 7.25])
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),100)
	o.integrate(ts,ppot)
	cluster=ctools.load_cluster(ctype='galpy',particles=o,units='kpckms',origin='galaxy')
	jr,jp,jz,OR, Ophi, Oz, TR, Tphi, Tz=cluster.calc_actions(pot=ppot, ro=ro, vo=vo,solarmotion=[-11.1, 24.0, 7.25],full=True)

	np.testing.assert_allclose(jr,o.jr(),rtol=tol,atol=1.)
	np.testing.assert_allclose(jp,o.jp(),rtol=tol,atol=1.)
	np.testing.assert_allclose(jz,o.jz(),rtol=tol,atol=1.)
	np.testing.assert_allclose(OR,o.Or(),rtol=tol,atol=1.)
	np.testing.assert_allclose(Ophi,o.Op(),rtol=tol,atol=1.)
	np.testing.assert_allclose(Oz,o.Oz(),rtol=tol,atol=1.)
	np.testing.assert_allclose(TR,o.Tr(),rtol=tol,atol=1.)
	np.testing.assert_allclose(Tphi,o.Tp(),rtol=tol,atol=1.)
	np.testing.assert_allclose(Tz,o.Tz(),rtol=tol,atol=1.)

def test_ttensor(tol=0.1,ro=8.,vo=220.):
	
	pmass= KeplerPotential(normalize=1.)
	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(8.,0.,0.,0.,220.,0.,ounits='kpckms')

	tij=cluster.ttensor(pot=pmass)

	assert np.all(np.fabs(tij-np.diag([2,-1,-1])) < 1e-10)

	tij=cluster.ttensor(pot=pmass,eigenval=True)

	assert np.all(np.fabs(tij-np.array([2,-1,-1])) < 1e-10)

def test_ttensors(tol=0.1,ro=8.,vo=220.):

	pmass= KeplerPotential(normalize=1.)
	x,y,z=np.ones(100)*8.,np.zeros(100),np.zeros(100)
	vx,vy,vz=np.zeros(100),np.ones(100)*220.,np.zeros(100)

	cluster=ctools.StarCluster(units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,analyze=True)
	cluster.add_orbit(8.,0.,0.,0.,220.,0.,ounits='kpckms')

	tij=cluster.ttensors(pot=pmass,eigenval=False)

	for i in range(0,cluster.ntot):
		assert np.all(np.fabs(tij[:,:,i]-np.diag([2,-1,-1])) < 1e-10)



