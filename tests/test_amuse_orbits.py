import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,IsochronePotential,KeplerPotential
try:
	from galpy.util import conversion
except:
	import galpy.util.bovy_conversion as conversion

try:
	import amuse.units.units as u
except:
	pass

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]

def test_initialize_orbit(tol=0.001):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_galaxy()
	cluster.to_amuse()

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	ocluster=cluster.initialize_orbit()

	assert np.fabs(o.x()-ocluster.x()) <= tol
	assert np.fabs(o.y()-ocluster.y()) <= tol
	assert np.fabs(o.z()-ocluster.z()) <= tol
	assert np.fabs(o.vx()-ocluster.vx()) <= tol
	assert np.fabs(o.vy()-ocluster.vy()) <= tol
	assert np.fabs(o.vz()-ocluster.vz()) <= tol

	ocluster=cluster.initialize_orbit(from_centre=True)

	assert np.fabs(o.x()-ocluster.x()+cluster.xc.value_in(u.kpc)) <= tol
	assert np.fabs(o.y()-ocluster.y()+cluster.yc.value_in(u.kpc)) <= tol
	assert np.fabs(o.z()-ocluster.z()+cluster.zc.value_in(u.kpc)) <= tol
	assert np.fabs(o.vx()-ocluster.vx()+cluster.vxc.value_in(u.kms)) <= tol
	assert np.fabs(o.vy()-ocluster.vy()+cluster.vyc.value_in(u.kms)) <= tol
	assert np.fabs(o.vz()-ocluster.vz()+cluster.vzc.value_in(u.kms)) <= tol

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
	cluster.to_amuse()


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

	np.testing.assert_allclose(x/1000.,ocluster.x()+cluster.xc.value_in(u.kpc),rtol=tol)
	np.testing.assert_allclose(y/1000.,ocluster.y()+cluster.yc.value_in(u.kpc),rtol=tol)
	np.testing.assert_allclose(z/1000.,ocluster.z()+cluster.zc.value_in(u.kpc),rtol=tol)

	np.testing.assert_allclose(vx,ocluster.vx()+cluster.vxc.value_in(u.kms),rtol=tol)
	np.testing.assert_allclose(vy,ocluster.vy()+cluster.vyc.value_in(u.kms),rtol=tol)
	np.testing.assert_allclose(vz,ocluster.vz()+cluster.vzc.value_in(u.kms),rtol=tol)

	cluster.to_cluster()
	ocluster=cluster.initialize_orbits()

	np.testing.assert_allclose(x/1000.,ocluster.x(),rtol=tol)
	np.testing.assert_allclose(y/1000.,ocluster.y(),rtol=tol)
	np.testing.assert_allclose(z/1000.,ocluster.z(),rtol=tol)

	np.testing.assert_allclose(vx,ocluster.vx(),rtol=tol)
	np.testing.assert_allclose(vy,ocluster.vy(),rtol=tol)
	np.testing.assert_allclose(vz,ocluster.vz(),rtol=tol)

def test_interpolate_orbit(tol=0.1,ro=solar_ro,vo=solar_vo):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy')
	cluster.to_amuse()

	x,y,z,vx,vy,vz=cluster.interpolate_orbit(pot=MWPotential2014,tfinal=1. | u.Gyr, nt=1000)

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo))
	o.integrate(ts,MWPotential2014)


	assert np.fabs(cluster.xgc.value_in(u.kpc)-o.x(ts[-1])) <= tol
	assert np.fabs(cluster.ygc.value_in(u.kpc)-o.y(ts[-1])) <= tol
	assert np.fabs(cluster.zgc.value_in(u.kpc)-o.z(ts[-1])) <= tol

	assert np.fabs(cluster.vxgc.value_in(u.kms)-o.vx(ts[-1])) <= tol
	assert np.fabs(cluster.vygc.value_in(u.kms)-o.vy(ts[-1])) <= tol
	assert np.fabs(cluster.vzgc.value_in(u.kms)-o.vz(ts[-1])) <= tol

def test_interpolate_orbits(tol=0.1,ro=solar_ro,vo=solar_vo):

	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy')
	cluster.to_amuse()

	xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.interpolate_orbit(pot=MWPotential2014,tfinal=1. | u.Gyr,nt=1000)

	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy',mbar=10)
	cluster.to_amuse()

	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)

	x,y,z,vx,vy,vz=cluster.interpolate_orbits(pot=MWPotential2014,tfinal=1. | u.Gyr,nt=1000)

	np.testing.assert_allclose(cluster.x.value_in(u.kpc),x.value_in(u.kpc),rtol=tol)
	np.testing.assert_allclose(cluster.y.value_in(u.kpc),y.value_in(u.kpc),rtol=tol)
	np.testing.assert_allclose(cluster.z.value_in(u.kpc),z.value_in(u.kpc),rtol=tol)
	np.testing.assert_allclose(cluster.vx.value_in(u.kms),vx.value_in(u.kms),rtol=tol)
	np.testing.assert_allclose(cluster.vy.value_in(u.kms),vy.value_in(u.kms),rtol=tol)
	np.testing.assert_allclose(cluster.vz.value_in(u.kms),vz.value_in(u.kms),rtol=tol)

	vxvv=np.column_stack([rad/1000.0/ro,vR/vo,vT/vo,zed/1000.0/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo))
	o.integrate(ts,MWPotential2014)

	np.testing.assert_allclose(cluster.x.value_in(u.kpc),o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.y.value_in(u.kpc),o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.z.value_in(u.kpc),o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vx.value_in(u.kms),o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vy.value_in(u.kms),o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vz.value_in(u.kms),o.vz(ts[-1]),rtol=tol)

	assert np.fabs(xgc.value_in(u.pc)-cluster.xgc.value_in(u.pc)) <= tol
	assert np.fabs(ygc.value_in(u.pc)-cluster.ygc.value_in(u.pc)) <= tol
	assert np.fabs(zgc.value_in(u.pc)-cluster.zgc.value_in(u.pc)) <= tol
	assert np.fabs(vxgc.value_in(u.kms)-cluster.vxgc.value_in(u.kms)) <= tol
	assert np.fabs(vygc.value_in(u.kms)-cluster.vygc.value_in(u.kms)) <= tol
	assert np.fabs(vzgc.value_in(u.kms)-cluster.vzgc.value_in(u.kms)) <= tol

	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='cluster',mbar=10)
	cluster.to_amuse()

	xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc
	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)

	cluster.to_galpy()
	pot=PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
	cluster.to_amuse()

	vxvv=np.column_stack([rad/ro/1000.0,vR/vo,vT/vo,zed/ro/1000.0,vzed/vo,phi])
	o=Orbit(vxvv,ro=solar_ro,vo=solar_vo)
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),1000)
	o.integrate(ts,pot)

	cluster.interpolate_orbits(pot=pot,tfinal=1 | u.Gyr,nt=1000)

	np.testing.assert_allclose(cluster.x.value_in(u.kpc),o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.y.value_in(u.kpc),o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.z.value_in(u.kpc),o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vx.value_in(u.kms),o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vy.value_in(u.kms),o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vz.value_in(u.kms),o.vz(ts[-1]),rtol=tol)

	assert np.fabs(xgc.value_in(u.pc)-cluster.xgc.value_in(u.pc)) <= tol
	assert np.fabs(ygc.value_in(u.pc)-cluster.ygc.value_in(u.pc)) <= tol
	assert np.fabs(zgc.value_in(u.pc)-cluster.zgc.value_in(u.pc)) <= tol
	assert np.fabs(vxgc.value_in(u.kms)-cluster.vxgc.value_in(u.kms)) <= tol
	assert np.fabs(vygc.value_in(u.kms)-cluster.vygc.value_in(u.kms)) <= tol
	assert np.fabs(vzgc.value_in(u.kms)-cluster.vzgc.value_in(u.kms)) <= tol


	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101',units='pckms',origin='centre',mbar=10)
	cluster.to_amuse()

	xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc

	cluster.to_kpckms()
	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)
	cluster.to_pckms()

	cluster.to_galpy()
	pot=PlummerPotential(cluster.mtot,b=cluster.rm/1.305,ro=ro,vo=vo)
	cluster.to_amuse()

	vxvv=np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=solar_ro,vo=solar_vo)
	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),1000)
	o.integrate(ts,pot)

	cluster.interpolate_orbits(pot=pot,tfinal=1000. | u.Myr,nt=1000)

	np.testing.assert_allclose(cluster.x.value_in(u.kpc),o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.y.value_in(u.kpc),o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.z.value_in(u.kpc),o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vx.value_in(u.kms),o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vy.value_in(u.kms),o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(cluster.vz.value_in(u.kms),o.vz(ts[-1]),rtol=tol)

	assert np.fabs(xgc.value_in(u.pc)-cluster.xgc.value_in(u.pc)) <= tol
	assert np.fabs(ygc.value_in(u.pc)-cluster.ygc.value_in(u.pc)) <= tol
	assert np.fabs(zgc.value_in(u.pc)-cluster.zgc.value_in(u.pc)) <= tol
	assert np.fabs(vxgc.value_in(u.kms)-cluster.vxgc.value_in(u.kms)) <= tol
	assert np.fabs(vygc.value_in(u.kms)-cluster.vygc.value_in(u.kms)) <= tol
	assert np.fabs(vzgc.value_in(u.kms)-cluster.vzgc.value_in(u.kms)) <= tol

def test_orbital_path(tol=0.1,ro=solar_ro,vo=solar_vo):
	cluster=ctools.load_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy',mbar=10)
	cluster.to_amuse()

	tfinal=0.1 | u.Gyr

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

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, (-1.0* tfinal.value_in(u.Gyr))/ conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[0].value_in(u.kpc)-o.x(ts[-1])) <= tol
	assert np.fabs(y[0].value_in(u.kpc)-o.y(ts[-1])) <= tol
	assert np.fabs(z[0].value_in(u.kpc)-o.z(ts[-1])) <= tol
	assert np.fabs(vx[0].value_in(u.kms)-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[0].value_in(u.kms)-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[0].value_in(u.kms)-o.vz(ts[-1])) <= tol

	assert np.fabs(cluster.xpath[0].value_in(u.kpc)-o.x(ts[-1])) <= tol
	assert np.fabs(cluster.ypath[0].value_in(u.kpc)-o.y(ts[-1])) <= tol
	assert np.fabs(cluster.zpath[0].value_in(u.kpc)-o.z(ts[-1])) <= tol
	assert np.fabs(cluster.vxpath[0].value_in(u.kms)-o.vx(ts[-1])) <= tol
	assert np.fabs(cluster.vypath[0].value_in(u.kms)-o.vy(ts[-1])) <= tol
	assert np.fabs(cluster.vzpath[0].value_in(u.kms)-o.vz(ts[-1])) <= tol

	assert np.fabs(cluster.orbit.x(ts[-1])-o.x(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.y(ts[-1])-o.y(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.z(ts[-1])-o.z(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.vx(ts[-1])-o.vx(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.vy(ts[-1])-o.vy(ts[-1])) <= tol
	assert np.fabs(cluster.orbit.vz(ts[-1])-o.vz(ts[-1])) <= tol

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, (1.0* tfinal.value_in(u.Gyr))/ conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[-1].value_in(u.kpc)-o.x(ts[-1])) <= tol
	assert np.fabs(y[-1].value_in(u.kpc)-o.y(ts[-1])) <= tol
	assert np.fabs(z[-1].value_in(u.kpc)-o.z(ts[-1])) <= tol
	assert np.fabs(vx[-1].value_in(u.kms)-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[-1].value_in(u.kms)-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[-1].value_in(u.kms)-o.vz(ts[-1])) <= tol

	#Check skypath
	t, x, y, z, vx, vy, vz=cluster.orbital_path(
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
	    ro=solar_ro,
	    vo=solar_vo,
	    solarmotion=solar_motion,
	    plot=False,
	)

	xgc,ygc,zgc=cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc
	vxgc,vygc,vzgc=cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc

	rad,phi,zed,vR,vT,vzed=ctools.cart_to_cyl(xgc.value_in(u.kpc),ygc.value_in(u.kpc),zgc.value_in(u.kpc),vxgc.value_in(u.kms),vygc.value_in(u.kms),vzgc.value_in(u.kms))
	vxvv=[rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi]

	o=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
	ts = np.linspace(0, (-1.0* tfinal.value_in(u.Gyr))/ conversion.time_in_Gyr(ro=ro, vo=vo), 1000)
	o.integrate(ts,MWPotential2014)

	assert np.fabs(x[0].value_in(u.kpc)-o.x(ts[-1])) <= tol
	assert np.fabs(y[0].value_in(u.kpc)-o.y(ts[-1])) <= tol
	assert np.fabs(z[0].value_in(u.kpc)-o.z(ts[-1])) <= tol
	assert np.fabs(vx[0].value_in(u.kms)-o.vx(ts[-1])) <= tol
	assert np.fabs(vy[0].value_in(u.kms)-o.vy(ts[-1])) <= tol
	assert np.fabs(vz[0].value_in(u.kms)-o.vz(ts[-1])) <= tol

def test_orbital_path_match(tol=0.1,ro=solar_ro,vo=solar_vo):

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

	cluster.to_amuse()

	t,dprog,dpath=cluster.orbital_path_match(
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

	t,dprog,dpath=cluster.orbital_path_match(
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

	t,dprog,dpath=cluster.orbital_path_match(
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

	t,dprog,dpath=cluster.orbital_path_match(
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
	cluster.to_amuse()

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
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
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
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
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
    ro=solar_ro,
    vo=solar_vo,
    solarmotion=solar_motion,
    plot=False,
    projected=False,
	)

	np.testing.assert_allclose(np.zeros(len(dpath)),dpath,rtol=tol,atol=1.)
	np.testing.assert_allclose(dprog,x,rtol=tol,atol=1.)
