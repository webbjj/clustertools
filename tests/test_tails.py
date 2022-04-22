import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,IsochronePotential,KeplerPotential
try:
	from galpy.util import conversion
except:
	import galpy.util.bovy_conversion as conversion

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]

def test_to_tail(tol=0.1):
	
	x=np.arange(-5,5,0.1)
	y=np.zeros(len(x))
	z=np.zeros(len(x))
	vx=np.ones(len(x))
	vy=np.zeros(len(x))
	vz=np.zeros(len(x))

	cluster=ctools.StarCluster(origin='galaxy',units='kpckms')

	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(0.,0.,0.,1.0,0.,0.)

	cluster.to_tail()

	np.testing.assert_allclose(x,cluster.x_tail,rtol=tol,atol=1.)

	y=np.arange(-5,5,0.1)
	x=np.zeros(len(x))
	z=np.zeros(len(x))
	vy=np.ones(len(x))
	vx=np.zeros(len(x))
	vz=np.zeros(len(x))

	cluster=ctools.StarCluster(origin='galaxy',units='kpckms')

	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(0.,0.,0.,0.,1.,0.)

	cluster.to_tail()

	np.testing.assert_allclose(y,-cluster.x_tail,rtol=tol,atol=1.)

	z=np.arange(-5,5,0.1)
	x=np.zeros(len(x))
	y=np.zeros(len(x))
	vz=np.ones(len(x))
	vx=np.zeros(len(x))
	vy=np.zeros(len(x))

	cluster=ctools.StarCluster(origin='galaxy',units='kpckms')

	cluster.add_stars(x,y,z,vx,vy,vz)
	cluster.add_orbit(0.,0.,0.,0.,0.,1.)

	cluster.to_tail()

	np.testing.assert_allclose(z,-cluster.x_tail,rtol=tol,atol=1.)


def test_tail_path(tol=0.1,ro=solar_ro,vo=solar_vo):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy',mbar=10)

	tfinal=0.1

	t, x, y, z, vx, vy, vz=cluster.orbital_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=False,
	    skypath=False,
	    initialize=True,
	    ro=solar_ro,
	    vo=solar_vo,
	    solarmotion=solar_motion,
	    plot=False,
	)

	cluster=ctools.StarCluster(units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz)

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, -1.0 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)
	cluster.add_orbit(o.x(),o.y(),o.z(),o.vx(),o.vy(),o.vz())

	t, x, y, z, vx, vy, vz=cluster.tail_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=False,
	    skypath=False,
	    initialize=True,
	    ro=solar_ro,
	    vo=solar_vo,
	    solarmotion=solar_motion,
	    plot=False,
	)


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


	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, 1.0 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[-1]-o.x(ts[-1])) <= tol
	assert np.fabs(y[-1]-o.y(ts[-1])) <= tol
	assert np.fabs(z[-1]-o.z(ts[-1])) <= tol
	assert np.fabs(vx[-1]-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[-1]-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[-1]-o.vz(ts[-1])) <= tol

	#Check skypath

	t, x, y, z, vx, vy, vz=cluster.tail_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=False,
	    skypath=True,
	    initialize=True,
	    ro=solar_ro,
	    vo=solar_vo,
	    solarmotion=solar_motion,
	    plot=False,
	)

	print(x[-5:],o.ra(ts)[-5:])
	print(y[-5:],o.dec(ts)[-5:])
	print(z[-5:],o.dist(ts)[-5:])
	print(vx[-5:],o.pmra(ts)[-5:])
	print(vy[-5:],o.pmdec(ts)[-5:])
	print(vz[-5:],o.vlos(ts)[-5:])


	assert np.fabs(x[-1]-o.ra(ts[-1])) <= tol
	assert np.fabs(y[-1]-o.dec(ts[-1])) <= tol
	assert np.fabs(z[-1]-o.dist(ts[-1])) <= tol
	assert np.fabs(vx[-1]-o.pmra(ts[-1])) <= tol
	assert np.fabs(vy[-1]-o.pmdec(ts[-1])) <= tol
	assert np.fabs(vz[-1]-o.vlos(ts[-1])) <= tol

	#Check from_centre:

	cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=0.,0.,0.,0.,0.,0.


	t, x, y, z, vx, vy, vz=cluster.tail_path(
	    tfinal=tfinal,
	    nt=1000,
	    pot=MWPotential2014,
	    from_centre=True,
	    skypath=False,
	    initialize=False,
	    ro=solar_ro,
	    vo=solar_vo,
	    solarmotion=solar_motion,
	    plot=False,
	)

	xgc,ygc,zgc=cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc
	vxgc,vygc,vzgc=cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc

	rad,phi,zed,vR,vT,vzed=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi]

	o=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, -1.0 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[0]-o.x(ts[-1])) <= tol
	assert np.fabs(y[0]-o.y(ts[-1])) <= tol
	assert np.fabs(z[0]-o.z(ts[-1])) <= tol
	assert np.fabs(vx[0]-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[0]-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[0]-o.vz(ts[-1])) <= tol

def test_tail_path_match(tol=0.1,ro=solar_ro,vo=solar_vo):

	tfinal=0.1
	nt=1000

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, 0.5 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), nt)
	o.integrate(ts,MWPotential2014)

	x,y,z=o.x(ts),o.y(ts),o.z(ts)
	vx,vy,vz=o.vx(ts),o.vy(ts),o.vz(ts)

	cluster=ctools.StarCluster(units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,sortstars=True,analyze=True)
	cluster.add_orbit(o.x(),o.y(),o.z(),o.vx(),o.vy(),o.vz(),ounits='kpckms')

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, -0.5 * tfinal / conversion.time_in_Gyr(ro=ro, vo=vo), nt)
	o.integrate(ts,MWPotential2014)

	x,y,z=o.x(ts),o.y(ts),o.z(ts)
	vx,vy,vz=o.vx(ts),o.vy(ts),o.vz(ts)
	cluster.add_stars(x,y,z,vx,vy,vz,sortstars=True,analyze=True)

	cluster.xc,cluster.yc,cluster.zc=0.,0.,0.
	cluster.vxc,cluster.vyc,cluster.vzc=0.,0.,0.

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=True,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=False,
    skypath=True,
    to_path=False,
    do_full=False,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=None,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
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

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=path,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=False,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=path,
    from_centre=False,
    skypath=False,
    to_path=True,
    do_full=False,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)

	t,dprog,dpath=cluster.tail_path_match(
    tfinal=tfinal,
    nt=10*nt,
	pot=MWPotential2014,
    path=path,
    from_centre=False,
    skypath=False,
    to_path=False,
    do_full=True,
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)