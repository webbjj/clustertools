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

def test_integrate_orbit(tol=0.1,ro=solar_ro,vo=solar_vo):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',units='kpckms',origin='galaxy')
	ocluster=cluster.integrate_orbit(tfinal=1.,nt=1000)

	o=Orbit.from_name('NGC6101',ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)
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

def test_integrate_orbits(tol=0.0001,ro=solar_ro,vo=solar_vo):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101',mbar=10)
	cluster.to_galaxy()
	cluster.to_kpckms()

	ocluster=cluster.integrate_orbits(tfinal=1,nt=1000)

	rad,phi,zed,vR,vT,vzed=ctools.cyl_coords(cluster)
	vxvv=np.column_stack([rad/ro,vR/vo,vT/vo,zed/ro,vzed/vo,phi])
	o=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	ts=np.linspace(0,1./conversion.time_in_Gyr(ro=ro,vo=vo),1000)
	o.integrate(ts,MWPotential2014)

	np.testing.assert_allclose(ocluster.x(ts[-1]),o.x(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.y(ts[-1]),o.y(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.z(ts[-1]),o.z(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.vx(ts[-1]),o.vx(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.vy(ts[-1]),o.vy(ts[-1]),rtol=tol)
	np.testing.assert_allclose(ocluster.vz(ts[-1]),o.vz(ts[-1]),rtol=tol)
