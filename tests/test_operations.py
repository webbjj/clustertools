import clustertools as ctools
import numpy as np
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014,PlummerPotential,IsochronePotential,KeplerPotential
try:
	from galpy.util import conversion,coords
except:
	import galpy.util.bovy_conversion as conversion
	import galpy.util.bovy_coords as coords

from copy import copy

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]

def scale_test(init,final,tscale,mscale,xscale,vscale):
	#assume arrays are t,mass,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc

	for i in range(0,len(init)):
		if i==0:
			assert np.all(np.fabs(final[i]/init[i]-tscale) < 1e-10)
		elif i==1:
			assert np.all(np.fabs(final[i]/init[i]-mscale) < 1e-4)
		elif (i>=2 and i<=4) or (i>=8 and i<=10) or (i>=14 and i<=16):
			assert np.all(np.fabs(final[i]/init[i]-xscale) < 1e-10)
		elif (i>=5 and i<=7) or (i>=11 and i<=13) or (i>=17 and i<=19):
			assert np.all(np.fabs(final[i]/init[i]-vscale) < 1e-10)

def test_to_pckms(tol=0.0001,ro=solar_ro,vo=solar_vo):

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['kpckms','galpy','nbody','WDunits']
	tscales=[1000,(1000.*conversion.time_in_Gyr(ro=ro,vo=vo)),tbar,1000.0]
	xscales=[1000,ro*1000,rbar,1000.]
	vscales=[1.,vo,vbar,220.0/conversion.velocity_in_kpcGyr(220.0, 8.0)]
	mscales=[1.,conversion.mass_in_msol(ro=ro,vo=vo),zmbar,222288.4543021174]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_pckms()

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])


def test_to_kpckms(tol=0.0001,ro=solar_ro,vo=solar_vo):

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['pckms','galpy','nbody','WDunits']
	tscales=[1./1000,conversion.time_in_Gyr(ro=ro,vo=vo),tbar/1000,1.]
	xscales=[1./1000,ro,rbar/1000.,1.]
	vscales=[1.,vo,vbar,220.0/conversion.velocity_in_kpcGyr(220.0, 8.0)]
	mscales=[1.,conversion.mass_in_msol(ro=ro,vo=vo),zmbar,222288.4543021174]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_kpckms()

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

def test_to_galpy(tol=0.0001,ro=solar_ro,vo=solar_vo):

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['pckms','kpckms','nbody','WDunits']
	tscales=[1./1000/conversion.time_in_Gyr(ro=ro,vo=vo),1./conversion.time_in_Gyr(ro=ro,vo=vo),tbar/1000./conversion.time_in_Gyr(ro=ro,vo=vo),1./conversion.time_in_Gyr(ro=ro,vo=vo)]
	xscales=[1./1000/ro,1./ro,rbar/1000./ro,1./ro]
	vscales=[1./vo,1./vo,vbar/vo,220.0/conversion.velocity_in_kpcGyr(220.0, 8.0)/vo]
	mscales=[1./conversion.mass_in_msol(ro=ro,vo=vo),1./conversion.mass_in_msol(ro=ro,vo=vo),zmbar/conversion.mass_in_msol(ro=ro,vo=vo),222288.4543021174/conversion.mass_in_msol(ro=ro,vo=vo)]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_galpy()

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

def test_to_nbody(tol=0.0001,ro=solar_ro,vo=solar_vo):

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['pckms','kpckms','galpy','WDunits']
	tscales=[1./tbar,1000./tbar,1000.0*conversion.time_in_Gyr(ro=ro,vo=vo)/tbar,1000./tbar]
	xscales=[1./rbar,1000./rbar,ro*1000./rbar,1000./rbar]
	vscales=[1./vbar,1./vbar,vo/vbar,220.0/conversion.velocity_in_kpcGyr(220.0, 8.0)/vbar]
	mscales=[1./zmbar,1./zmbar,conversion.mass_in_msol(ro=ro,vo=vo)/zmbar,222288.4543021174/zmbar]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_nbody()

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

def test_to_wdunits(tol=0.01,ro=solar_ro,vo=solar_vo):
	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	vcon=220./conversion.velocity_in_kpcGyr(220.0, 8.0)

	init_units=['pckms','galpy','nbody','kpckms']
	tscales=[1./1000,conversion.time_in_Gyr(ro=ro,vo=vo),tbar/1000,1.]
	xscales=[1./1000,ro,rbar/1000.,1.]
	vscales=[1./vcon,vo/vcon,vbar/vcon,1./vcon]
	mscales=[1./222288.4543021174,conversion.mass_in_msol(ro=ro,vo=vo)/222288.4543021174,zmbar/222288.4543021174,1./222288.4543021174]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_WDunits()

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])
def test_to_radec(tol=0.0001,ro=solar_ro,vo=solar_vo):
	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x,y,z,vx,vy,vz)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	cluster.to_radec()

	assert np.all(np.fabs(ogc.ra()/cluster.ra_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.ra()/cluster.xgc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dec()/cluster.dec_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dec()/cluster.ygc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dist()/cluster.dist_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dist()/cluster.zgc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmra()/cluster.pmra_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmra()/cluster.vxgc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmdec()/cluster.pmdec_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmdec()/cluster.vygc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.vlos()/cluster.vlos_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.vlos()/cluster.vzgc-1.) < 1e-10)

	assert np.all(np.fabs(oc.ra()/cluster.ra_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.ra()/cluster.xc-1.) < 1e-10)
	assert np.all(np.fabs(oc.dec()/cluster.dec_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.dec()/cluster.yc-1.) < 1e-10)
	assert np.all(np.fabs(oc.dist()/cluster.dist_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.dist()/cluster.zc-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmra()/cluster.pmra_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmra()/cluster.vxc-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmdec()/cluster.pmdec_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmdec()/cluster.vyc-1.) < 1e-10)
	assert np.all(np.fabs(oc.vlos()/cluster.vlos_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.vlos()/cluster.vzc-1.) < 1e-10)

	assert np.all(np.fabs(os.ra()/cluster.ra-1.) < 1e-10)
	assert np.all(np.fabs(os.ra()/cluster.x-1.) < 1e-10)
	assert np.all(np.fabs(os.dec()/cluster.dec-1.) < 1e-10)
	assert np.all(np.fabs(os.dec()/cluster.y-1.) < 1e-10)
	assert np.all(np.fabs(os.dist()/cluster.dist-1.) < 1e-10)
	assert np.all(np.fabs(os.dist()/cluster.z-1.) < 1e-10)
	assert np.all(np.fabs(os.pmra()/cluster.pmra-1.) < 1e-10)
	assert np.all(np.fabs(os.pmra()/cluster.vx-1.) < 1e-10)
	assert np.all(np.fabs(os.pmdec()/cluster.pmdec-1.) < 1e-10)
	assert np.all(np.fabs(os.pmdec()/cluster.vy-1.) < 1e-10)
	assert np.all(np.fabs(os.vlos()/cluster.vlos-1.) < 1e-10)
	assert np.all(np.fabs(os.vlos()/cluster.vz-1.) < 1e-10)


	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x+xgc,y+ygc,z+zgc,vx+vxgc,vy+vygc,vz+vzgc)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	cluster.to_radec()

	newra=(os.ra()-ogc.ra())*np.cos(np.radians(ogc.dec()))
	newdec=os.dec()-ogc.dec()
	newdist=os.dist()-ogc.dist()
	newpmra=os.pmra()-ogc.pmra()
	newpmdec=os.pmdec()-ogc.pmdec()
	newvlos=os.vlos()-ogc.vlos()

	assert np.all(np.fabs(newra/cluster.x-1) < 1e-10)
	assert np.all(np.fabs(newdec/cluster.y-1) < 1e-10)
	assert np.all(np.fabs(newdist/cluster.z-1) < 1e-10)
	assert np.all(np.fabs(newpmra/cluster.vx-1) < 1e-10)
	assert np.all(np.fabs(newpmdec/cluster.vy-1) < 1e-10)
	assert np.all(np.fabs(newvlos/cluster.vz-1) < 1e-10)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x+xgc+xc,y+ygc+yc,z+zgc+zc,vx+vxgc+vxc,vy+vygc+vyc,vz+vzgc+vzc)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	cluster.to_radec()

	newra=(os.ra()-oc.ra())*np.cos(np.radians(oc.dec()))
	newdec=os.dec()-oc.dec()
	newdist=os.dist()-oc.dist()
	newpmra=os.pmra()-oc.pmra()
	newpmdec=os.pmdec()-oc.pmdec()
	newvlos=os.vlos()-oc.vlos()

	assert np.all(np.fabs(newra/cluster.x-1) < 1e-10)
	assert np.all(np.fabs(newdec/cluster.y-1) < 1e-10)
	assert np.all(np.fabs(newdist/cluster.z-1) < 1e-10)
	assert np.all(np.fabs(newpmra/cluster.vx-1) < 1e-10)
	assert np.all(np.fabs(newpmdec/cluster.vy-1) < 1e-10)
	assert np.all(np.fabs(newvlos/cluster.vz-1) < 1e-10)

def test_to_units(tol=0.01,ro=solar_ro,vo=solar_vo):

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['kpckms','galpy','nbody']
	tscales=[1000,(1000.*conversion.time_in_Gyr(ro=ro,vo=vo)),tbar]
	xscales=[1000,ro*1000,rbar]
	vscales=[1.,vo,vbar]
	mscales=[1.,conversion.mass_in_msol(ro=ro,vo=vo),zmbar]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_units('pckms')

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['pckms','galpy','nbody']
	tscales=[1./1000,conversion.time_in_Gyr(ro=ro,vo=vo),tbar/1000]
	xscales=[1./1000,ro,rbar/1000.]
	vscales=[1.,vo,vbar]
	mscales=[1.,conversion.mass_in_msol(ro=ro,vo=vo),zmbar]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_units('kpckms')

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['pckms','kpckms','nbody']
	tscales=[1./1000/conversion.time_in_Gyr(ro=ro,vo=vo),1./conversion.time_in_Gyr(ro=ro,vo=vo),tbar/1000./conversion.time_in_Gyr(ro=ro,vo=vo)]
	xscales=[1./1000/ro,1./ro,rbar/1000./ro]
	vscales=[1./vo,1./vo,vbar/vo]
	mscales=[1./conversion.mass_in_msol(ro=ro,vo=vo),1./conversion.mass_in_msol(ro=ro,vo=vo),zmbar/conversion.mass_in_msol(ro=ro,vo=vo)]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_units('galpy')

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

	tbar,rbar,vbar,zmbar=2.1,5.3,2.5,7553.7

	init_units=['pckms','kpckms','galpy']
	tscales=[1./tbar,1000./tbar,1000.0*conversion.time_in_Gyr(ro=ro,vo=vo)/tbar]
	xscales=[1./rbar,1000./rbar,ro*1000./rbar]
	vscales=[1./vbar,1./vbar,vo/vbar]
	mscales=[1./zmbar,1./zmbar,conversion.mass_in_msol(ro=ro,vo=vo)/zmbar]

	for i in range(0,len(init_units)):
		t=1.
		m=np.random.rand(100)
		x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)
		xc,yc,zc,vxc,vyc,vzc=np.random.rand(6)
		xgc,ygc,zgc,vxgc,vygc,vzgc=np.random.rand(6)

		cluster=ctools.StarCluster(tphys=t,units=init_units[i],origin='cluster')
		cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
		cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits=init_units[i])
		cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
		cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
		cluster.tbar,cluster.rbar,cluster.vbar,cluster.zmbar=tbar,rbar,vbar,zmbar

		init=[t,m,x,y,z,vx,vy,vz,xc,yc,zc,vxc,vyc,vzc,xgc,ygc,zgc,vxgc,vygc,vzgc]

		cluster.to_units('nbody')

		final=[cluster.tphys,cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc]
		scale_test(init,final,tscales[i],mscales[i],xscales[i],vscales[i])

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x,y,z,vx,vy,vz)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	cluster.to_units('radec')

	assert np.all(np.fabs(ogc.ra()/cluster.ra_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.ra()/cluster.xgc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dec()/cluster.dec_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dec()/cluster.ygc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dist()/cluster.dist_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.dist()/cluster.zgc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmra()/cluster.pmra_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmra()/cluster.vxgc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmdec()/cluster.pmdec_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.pmdec()/cluster.vygc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.vlos()/cluster.vlos_gc-1.) < 1e-10)
	assert np.all(np.fabs(ogc.vlos()/cluster.vzgc-1.) < 1e-10)

	assert np.all(np.fabs(oc.ra()/cluster.ra_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.ra()/cluster.xc-1.) < 1e-10)
	assert np.all(np.fabs(oc.dec()/cluster.dec_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.dec()/cluster.yc-1.) < 1e-10)
	assert np.all(np.fabs(oc.dist()/cluster.dist_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.dist()/cluster.zc-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmra()/cluster.pmra_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmra()/cluster.vxc-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmdec()/cluster.pmdec_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.pmdec()/cluster.vyc-1.) < 1e-10)
	assert np.all(np.fabs(oc.vlos()/cluster.vlos_c-1.) < 1e-10)
	assert np.all(np.fabs(oc.vlos()/cluster.vzc-1.) < 1e-10)

	assert np.all(np.fabs(os.ra()/cluster.ra-1.) < 1e-10)
	assert np.all(np.fabs(os.ra()/cluster.x-1.) < 1e-10)
	assert np.all(np.fabs(os.dec()/cluster.dec-1.) < 1e-10)
	assert np.all(np.fabs(os.dec()/cluster.y-1.) < 1e-10)
	assert np.all(np.fabs(os.dist()/cluster.dist-1.) < 1e-10)
	assert np.all(np.fabs(os.dist()/cluster.z-1.) < 1e-10)
	assert np.all(np.fabs(os.pmra()/cluster.pmra-1.) < 1e-10)
	assert np.all(np.fabs(os.pmra()/cluster.vx-1.) < 1e-10)
	assert np.all(np.fabs(os.pmdec()/cluster.pmdec-1.) < 1e-10)
	assert np.all(np.fabs(os.pmdec()/cluster.vy-1.) < 1e-10)
	assert np.all(np.fabs(os.vlos()/cluster.vlos-1.) < 1e-10)
	assert np.all(np.fabs(os.vlos()/cluster.vz-1.) < 1e-10)

def test_to_audays(tol=0.01):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.to_cluster(sortstars=True)
	cluster.reset_nbody_scale()
	cluster.to_nbody()
	cluster.bunits='nbody'

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

	cluster.to_audays()

	assert np.all(np.fabs(semi/cluster.semi-1./cluster.rbar_au) < 1e-10)
	assert np.all(np.fabs(pb/cluster.pb-1./cluster.tbar_days) < 1e-10)
	assert np.all(np.fabs(m1/cluster.m1-1./cluster.zmbar) < 1e-10)
	assert np.all(np.fabs(m2/cluster.m2-1./cluster.zmbar) < 1e-10)

	cluster.to_nbody()

	assert np.all(np.fabs(semi/cluster.semi-1) < 1e-10)
	assert np.all(np.fabs(pb/cluster.pb-1.) < 1e-10)
	assert np.all(np.fabs(m1/cluster.m1-1.) < 1e-10)
	assert np.all(np.fabs(m2/cluster.m2-1.) < 1e-10)

def test_to_sudays(tol=0.01):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')

	zmbar,rbar,vbar,tbar=ctools.reset_nbody_scale(cluster)
	cluster.add_nbody6(1.,1.,rbar,1.,1.,1.,1.,zmbar,vbar,tbar,1.,1.,1.,1.)
	cluster.to_cluster(sortstars=True)
	cluster.to_nbody()
	cluster.bunits='nbody'

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

	cluster.to_sudays()

	assert np.all(np.fabs(semi/cluster.semi-1./cluster.rbar_su) < 1e-10)
	assert np.all(np.fabs(pb/cluster.pb-1./cluster.tbar_days) < 1e-10)
	assert np.all(np.fabs(m1/cluster.m1-1./cluster.zmbar) < 1e-10)
	assert np.all(np.fabs(m2/cluster.m2-1./cluster.zmbar) < 1e-10)

	cluster.to_nbody()

	assert np.all(np.fabs(semi/cluster.semi-1) < 1e-10)
	assert np.all(np.fabs(pb/cluster.pb-1.) < 1e-10)
	assert np.all(np.fabs(m1/cluster.m1-1.) < 1e-10)
	assert np.all(np.fabs(m2/cluster.m2-1.) < 1e-10)


def test_to_cluster(tol=0.01):
	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_cluster(sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xgc)
	dec_gc = np.radians(cluster.ygc)

	x = np.cos(dec) * np.sin(ra - ra_gc)
	y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
	    dec_gc
	) * np.cos(ra - ra_gc)

	z = np.zeros(len(x))

	vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
	    ra - ra_gc
	)
	vy = pmra * np.sin(dec_gc) * np.sin(
	    ra - ra_gc
	) + pmdec_deg * (
	    np.cos(dec) * np.cos(dec_gc)
	    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	)

	z=cluster.z-cluster.zgc
	vz=cluster.vz-cluster.vzgc

	cluster.to_cluster(centre_method = "orthographic")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xgc)
	dec_gc = np.radians(cluster.ygc)

	x = (
	    (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
	)
	y = (
	    (10800.0 / np.pi)
	    * (
	        np.sin(dec) * np.cos(dec_gc)
	        - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	    )
	    / 60.0
	)

	z = cluster.z-cluster.zgc

	vx = cluster.vx-cluster.vxgc
	vy = cluster.vy-cluster.vygc
	vz = cluster.vz-cluster.vzgc

	cluster.to_cluster(centre_method = "VandeVen")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)


def test_to_centre(tol=0.01):
	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_centre(sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc+xc,np.random.rand(100)+ygc+yc,np.random.rand(100)+zgc+zc
	vx,vy,vz=np.random.rand(100)+vxgc+vxc,np.random.rand(100)+vygc+vyc,np.random.rand(100)+vzgc+vzc

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_centre(sortstars=True)

	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc+cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc+cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc+cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc+cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc+cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc+cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = np.cos(dec) * np.sin(ra - ra_gc)
	y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
	    dec_gc
	) * np.cos(ra - ra_gc)

	z = np.zeros(len(x))

	vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
	    ra - ra_gc
	)
	vy = pmra * np.sin(dec_gc) * np.sin(
	    ra - ra_gc
	) + pmdec_deg * (
	    np.cos(dec) * np.cos(dec_gc)
	    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	)

	z=cluster.z-cluster.zc
	vz=cluster.vz-cluster.vzc

	cluster.to_centre(centre_method = "orthographic")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = (
	    (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
	)
	y = (
	    (10800.0 / np.pi)
	    * (
	        np.sin(dec) * np.cos(dec_gc)
	        - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	    )
	    / 60.0
	)

	z = cluster.z-cluster.zc

	vx = cluster.vx-cluster.vxc
	vy = cluster.vy-cluster.vyc
	vz = cluster.vz-cluster.vzc

	cluster.to_centre(centre_method = "VandeVen")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

def test_to_center(tol=0.01):
	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_center(sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc+xc,np.random.rand(100)+ygc+yc,np.random.rand(100)+zgc+zc
	vx,vy,vz=np.random.rand(100)+vxgc+vxc,np.random.rand(100)+vygc+vyc,np.random.rand(100)+vzgc+vzc

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_center(sortstars=True)

	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc+cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc+cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc+cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc+cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc+cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc+cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = np.cos(dec) * np.sin(ra - ra_gc)
	y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
	    dec_gc
	) * np.cos(ra - ra_gc)

	z = np.zeros(len(x))

	vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
	    ra - ra_gc
	)
	vy = pmra * np.sin(dec_gc) * np.sin(
	    ra - ra_gc
	) + pmdec_deg * (
	    np.cos(dec) * np.cos(dec_gc)
	    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	)

	z=cluster.z-cluster.zc
	vz=cluster.vz-cluster.vzc

	cluster.to_center(centre_method = "orthographic")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = (
	    (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
	)
	y = (
	    (10800.0 / np.pi)
	    * (
	        np.sin(dec) * np.cos(dec_gc)
	        - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	    )
	    / 60.0
	)

	z = cluster.z-cluster.zc

	vx = cluster.vx-cluster.vxc
	vy = cluster.vy-cluster.vyc
	vz = cluster.vz-cluster.vzc

	cluster.to_center(centre_method = "VandeVen")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

def test_to_galaxy(tol=0.01, ro=solar_ro, vo=solar_vo):

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_galaxy(sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_galaxy(sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc+cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc+cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc+cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc+cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc+cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc+cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x,y,z,vx,vy,vz)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	cluster.to_radec()
	cluster.to_galaxy()

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

def test_to_sky(tol=0.01, ro=solar_ro, vo=solar_vo):

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x,y,z,vx,vy,vz)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
	cluster.to_sky()

	cluster2=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster2.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster2.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster2.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster2.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
	cluster2.to_radec()

	assert np.all(np.fabs(cluster2.x-cluster.x) < tol)
	assert np.all(np.fabs(cluster2.y-cluster.y) < tol)
	assert np.all(np.fabs(cluster2.z-cluster.z) < tol)
	assert np.all(np.fabs(cluster2.vx-cluster.vx) < tol)
	assert np.all(np.fabs(cluster2.vy-cluster.vy) < tol)
	assert np.all(np.fabs(cluster2.vz-cluster.vz) < tol)

def test_to_origin(tol=0.01, ro=solar_ro, vo=solar_vo):
	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('cluster',sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	cluster.to_origin('galaxy')
	cluster.to_origin('radec',sortstars=True)

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xgc)
	dec_gc = np.radians(cluster.ygc)

	x = np.cos(dec) * np.sin(ra - ra_gc)
	y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
	    dec_gc
	) * np.cos(ra - ra_gc)

	z = np.zeros(len(x))

	vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
	    ra - ra_gc
	)
	vy = pmra * np.sin(dec_gc) * np.sin(
	    ra - ra_gc
	) + pmdec_deg * (
	    np.cos(dec) * np.cos(dec_gc)
	    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	)

	z=cluster.z-cluster.zgc
	vz=cluster.vz-cluster.vzgc

	cluster.to_origin('cluster',centre_method = "orthographic")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	cluster.to_origin('galaxy')
	cluster.to_origin('radec')

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xgc)
	dec_gc = np.radians(cluster.ygc)

	x = (
	    (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
	)
	y = (
	    (10800.0 / np.pi)
	    * (
	        np.sin(dec) * np.cos(dec_gc)
	        - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	    )
	    / 60.0
	)

	z = cluster.z-cluster.zgc

	vx = cluster.vx-cluster.vxgc
	vy = cluster.vy-cluster.vygc
	vz = cluster.vz-cluster.vzgc

	cluster.to_origin('cluster',centre_method = "VandeVen")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('centre',sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc+xc,np.random.rand(100)+ygc+yc,np.random.rand(100)+zgc+zc
	vx,vy,vz=np.random.rand(100)+vxgc+vxc,np.random.rand(100)+vygc+vyc,np.random.rand(100)+vzgc+vzc

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('centre',sortstars=True)

	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc+cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc+cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc+cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc+cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc+cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc+cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = np.cos(dec) * np.sin(ra - ra_gc)
	y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
	    dec_gc
	) * np.cos(ra - ra_gc)

	z = np.zeros(len(x))

	vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
	    ra - ra_gc
	)
	vy = pmra * np.sin(dec_gc) * np.sin(
	    ra - ra_gc
	) + pmdec_deg * (
	    np.cos(dec) * np.cos(dec_gc)
	    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	)

	z=cluster.z-cluster.zc
	vz=cluster.vz-cluster.vzc

	cluster.to_origin('centre',centre_method = "orthographic")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = (
	    (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
	)
	y = (
	    (10800.0 / np.pi)
	    * (
	        np.sin(dec) * np.cos(dec_gc)
	        - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	    )
	    / 60.0
	)

	z = cluster.z-cluster.zc

	vx = cluster.vx-cluster.vxc
	vy = cluster.vy-cluster.vyc
	vz = cluster.vz-cluster.vzc

	cluster.to_origin('centre',centre_method = "VandeVen")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('center',sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc+xc,np.random.rand(100)+ygc+yc,np.random.rand(100)+zgc+zc
	vx,vy,vz=np.random.rand(100)+vxgc+vxc,np.random.rand(100)+vygc+vyc,np.random.rand(100)+vzgc+vzc

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('center',sortstars=True)

	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc+cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc+cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc+cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc+cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc+cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc+cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = np.cos(dec) * np.sin(ra - ra_gc)
	y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
	    dec_gc
	) * np.cos(ra - ra_gc)

	z = np.zeros(len(x))

	vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
	    ra - ra_gc
	)
	vy = pmra * np.sin(dec_gc) * np.sin(
	    ra - ra_gc
	) + pmdec_deg * (
	    np.cos(dec) * np.cos(dec_gc)
	    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	)

	z=cluster.z-cluster.zc
	vz=cluster.vz-cluster.vzc

	cluster.to_origin('center',centre_method = "orthographic")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	cluster.to_galaxy()
	cluster.to_radec()

	ra,dec,dist=copy(cluster.x),copy(cluster.y),copy(cluster.z)
	pmra,pmdec,vlos=copy(cluster.vx),copy(cluster.vy),copy(cluster.vz)

	pmdec_deg=copy(pmdec)

	ra = np.radians(ra)
	dec = np.radians(dec)
	pmra = np.radians(pmra / (1000.0 * 3600.0))
	pmdec = np.radians(pmdec / (1000.0 * 3600.0))
	ra_gc = np.radians(cluster.xc)
	dec_gc = np.radians(cluster.yc)

	x = (
	    (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
	)
	y = (
	    (10800.0 / np.pi)
	    * (
	        np.sin(dec) * np.cos(dec_gc)
	        - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
	    )
	    / 60.0
	)

	z = cluster.z-cluster.zc

	vx = cluster.vx-cluster.vxc
	vy = cluster.vy-cluster.vyc
	vz = cluster.vz-cluster.vzc

	cluster.to_origin('center',centre_method = "VandeVen")

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('galaxy',sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100),np.random.rand(100),np.random.rand(100)
	vx,vy,vz=np.random.rand(100),np.random.rand(100),np.random.rand(100)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	rorder=copy(cluster.rorder)

	cluster.to_origin('galaxy',sortstars=True)
	assert np.all(np.fabs(x-cluster.x)-np.fabs(cluster.xgc+cluster.xc) < tol)
	assert np.all(np.fabs(y-cluster.y)-np.fabs(cluster.ygc+cluster.yc) < tol)
	assert np.all(np.fabs(z-cluster.z)-np.fabs(cluster.zgc+cluster.zc) < tol)
	assert np.all(np.fabs(vx-cluster.vx)-np.fabs(cluster.vxgc+cluster.vxc) < tol)
	assert np.all(np.fabs(vy-cluster.vy)-np.fabs(cluster.vygc+cluster.vyc) < tol)
	assert np.all(np.fabs(vz-cluster.vz)-np.fabs(cluster.vzgc+cluster.vzc) < tol)
	assert np.all(np.not_equal(cluster.rorder,rorder)==1)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x,y,z,vx,vy,vz)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc

	cluster.to_radec()
	cluster.to_galaxy()

	assert np.all(np.fabs(x-cluster.x) < tol)
	assert np.all(np.fabs(y-cluster.y) < tol)
	assert np.all(np.fabs(z-cluster.z) < tol)
	assert np.all(np.fabs(vx-cluster.vx) < tol)
	assert np.all(np.fabs(vy-cluster.vy) < tol)
	assert np.all(np.fabs(vz-cluster.vz) < tol)

	t=1.
	m=np.random.rand(100)
	xc,yc,zc,vxc,vyc,vzc=0.01*np.random.rand(6)
	xgc,ygc,zgc,vxgc,vygc,vzgc=10.0*np.random.rand(6)

	x,y,z=np.random.rand(100)+xgc,np.random.rand(100)+ygc,np.random.rand(100)+zgc
	vx,vy,vz=np.random.rand(100)+vxgc,np.random.rand(100)+vygc,np.random.rand(100)+vzgc

	rads,phis,zeds,vrads,vphis,vzeds=ctools.cart_to_cyl(x,y,z,vx,vy,vz)
	vxvvs=np.column_stack([rads/ro,vrads/vo,vphis/vo,zeds/ro,vzeds/vo,phis])
	os=Orbit(vxvvs,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radgc,phigc,zedgc,vradgc,vphigc,vzedgc=ctools.cart_to_cyl(xgc,ygc,zgc,vxgc,vygc,vzgc)
	vxvv=[radgc/ro,vradgc/vo,vphigc/vo,zedgc/ro,vzedgc/vo,phigc]
	ogc=Orbit(vxvv,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	radc,phic,zedc,vradc,vphic,vzedc=ctools.cart_to_cyl(xgc+xc,ygc+yc,zgc+zc,vxgc+vxc,vygc+vyc,vzgc+vzc)
	vxvvc=[radc/ro,vradc/vo,vphic/vo,zedc/ro,vzedc/vo,phic]
	oc=Orbit(vxvvc,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion)

	cluster=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
	cluster.to_origin('sky')

	cluster2=ctools.StarCluster(tphys=t,units='kpckms',origin='galaxy')
	cluster2.add_stars(x,y,z,vx,vy,vz,m=m,analyze=True)
	cluster2.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits='kpckms')
	cluster2.xc,cluster.yc,cluster.zc=xc,yc,zc
	cluster2.vxc,cluster.vyc,cluster.vzc=vxc,vyc,vzc
	cluster2.to_radec()

	assert np.all(np.fabs(cluster2.x-cluster.x) < tol)
	assert np.all(np.fabs(cluster2.y-cluster.y) < tol)
	assert np.all(np.fabs(cluster2.z-cluster.z) < tol)
	assert np.all(np.fabs(cluster2.vx-cluster.vx) < tol)
	assert np.all(np.fabs(cluster2.vy-cluster.vy) < tol)
	assert np.all(np.fabs(cluster2.vz-cluster.vz) < tol)

def test_save_cluster():

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')

	units0,origin0,rorder0,rorder_origin0=ctools.save_cluster(cluster)
	cluster.save_cluster()

	assert units0==cluster.units0
	assert origin0==cluster.origin0
	assert np.all(rorder0==cluster.rorder0)
	assert rorder_origin0==cluster.rorder_origin0

def test_return_cluster():

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	units0,origin0,rorder0,rorder_origin0=ctools.save_cluster(cluster)

	cluster.to_galpy()
	cluster.to_galaxy()

	cluster.return_cluster(units0,origin0,rorder0,rorder_origin0)

	assert units0==cluster.units
	assert origin0==cluster.origin
	assert np.all(rorder0==cluster.rorder)
	assert rorder_origin0==cluster.rorder_origin

def test_reset_nbody_scale(tol=0.01):

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	cluster.reset_nbody_scale()
	cluster.to_nbody()
	cluster.virial_radius()

	assert np.fabs(cluster.mtot-1. <= tol)
	assert np.fabs(cluster.rv-1. <= tol)

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')

	mtot0=cluster.mtot

	cluster.reset_nbody_scale(mass=False)
	cluster.to_nbody()
	cluster.virial_radius()

	assert np.fabs(cluster.mtot-mtot0 <= tol)
	assert np.fabs(cluster.rv-1. <= tol)

	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')

	rv0=cluster.virial_radius()

	cluster.reset_nbody_scale(radius=False)
	cluster.to_nbody()
	cluster.virial_radius()

	assert np.fabs(cluster.mtot-1. <= tol)
	assert np.fabs(cluster.rv-rv0 <= tol)

def test_add_rotation(tol=0.1):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')
	vr, vtheta, vz = coords.rect_to_cyl_vec(
	    cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
	)

	negvtheta0=np.sum(vtheta<0)
	posvtheta0=np.sum(vtheta>=0)

	qrot=0.5

	cluster.add_rotation(qrot)

	vr, vtheta, vz = coords.rect_to_cyl_vec(
	    cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
	) 

	negvtheta=np.sum(vtheta<0)
	posvtheta=np.sum(vtheta>=0)

	assert np.fabs(negvtheta/negvtheta0-qrot) <= tol


def test_virialize(tol=0.01):
	cluster=ctools.setup_cluster(ctype='limepy',gcname='NGC6101')

	cluster.energies()

	qv=np.sqrt(-0.5/cluster.qvir)

	cluster.virialize()
	cluster.energies()

	assert cluster.qvir==-0.5
	assert cluster.qv==qv



