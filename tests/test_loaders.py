import clustertools as ctools
import numpy as np
from astropy.table import QTable
from amuse.lab import *
from amuse.units import nbody_system,units
from amuse.datamodel import Particles


from galpy.potential import KingPotential
from galpy.util import conversion
from galpy.df import kingdf
from limepy import limepy

solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]
solar_zo=0.0208 #Bennett, M. & Bovy, J. 2019, MNRAS, 483, 1417

mo=conversion.mass_in_msol(ro=solar_ro,vo=solar_vo)

def check_params(cluster,ctype,units,origin,projected,**kwargs):

	assert np.fabs(cluster.tphys-kwargs.get('tphys',0.)) < 0.001

	if 'nbody' in ctype:
		assert cluster.units=='nbody'
		assert cluster.origin=='cluster'
	else:
		assert cluster.units==units
		assert cluster.origin==origin

	assert cluster.bunits==cluster.units
	assert cluster.units_init==cluster.units
	assert cluster.origin_init==cluster.origin

	assert cluster.ctype==ctype

	assert cluster.projected==projected

	
	assert cluster.nsnap == int(kwargs.get("nsnap", 0))
	assert cluster.delimiter == kwargs.get("delimiter", None)
	assert cluster.wdir == kwargs.get("wdir", "./")
	assert cluster.nzfill == int(kwargs.get("nzfill", 5))
	assert cluster.snapbase == kwargs.get("snapbase", "")
	assert cluster.snapend == kwargs.get("snapend", ".dat")
	assert cluster.snapdir == kwargs.get("snapdir", "")
	assert cluster.skiprows == kwargs.get("skiprows", 0)
	if cluster.sfile is not None:
		if isinstance(cluster.sfile,str):
			assert cluster.sfile == kwargs.get("sfile")
		else:
			assert cluster.sfile.name == kwargs.get("sfile")
	if cluster.bfile is not None:
		assert cluster.bfile.name == kwargs.get("bfile")
	if cluster.ssefile !=None:
		assert cluster.ssefile.name == kwargs.get("ssefile")
	if cluster.bsefile != None:
		assert cluster.bsefile.name == kwargs.get("bsefile")
	if cluster.ofile is not None:
		assert cluster.ofile.name == kwargs.get("ofile", None)
	assert cluster.ofilename == kwargs.get("ofilename", None)
	assert cluster.orbit == kwargs.get("orbit", None)
	assert cluster.give==kwargs.get('give','mxv')
	assert cluster.centre_method == kwargs.get("centre_method", None)

	assert cluster._ro==kwargs.get('ro',solar_ro)
	assert cluster._vo==kwargs.get('vo',solar_vo)
	assert cluster._zo==kwargs.get('zo',solar_zo)
	assert cluster._solarmotion==kwargs.get('_solarmotion',solar_motion)


	# Total Number of Stars + Binaries in the cluster
	assert cluster.ntot == len(cluster.x)
	assert cluster.nb == len(cluster.id1)

	assert cluster.id.size != 0
	assert cluster.m.size != 0
	assert cluster.x.size != 0
	assert cluster.y.size != 0
	assert cluster.z.size != 0
	assert cluster.vx.size != 0
	assert cluster.vy.size != 0
	assert cluster.vz.size != 0
	assert cluster.m0.size != 0

	if 'nbody6pp' in ctype or 'nbody6++' in ctype:
		assert cluster.rhos.size != 0

	if 'nbody' in ctype:
		assert cluster.zmbar != 1.
		assert cluster.rbar != 1.
		assert cluster.vbar != 1.
		assert cluster.tbar != 1.
	else:
		assert cluster.zmbar == 1.
		assert cluster.rbar == 1.
		assert cluster.vbar == 1.
		assert cluster.tbar == 1.


	# variables for centre of cluster
	if origin == 'cluster':
		assert cluster.xc != 0.0
		assert cluster.yc != 0.0
		assert cluster.zc != 0.0
		if 'nbody' in ctype:
			assert cluster.vxc == 0.0
			assert cluster.vyc == 0.0
			assert cluster.vzc == 0.0
		else:
			assert cluster.vxc != 0.0
			assert cluster.vyc != 0.0
			assert cluster.vzc != 0.0
	else:
		assert np.fabs(cluster.xc) <= 1.e10
		assert np.fabs(cluster.yc) <= 1.e10
		assert np.fabs(cluster.zc) <= 1.e10
		assert np.fabs(cluster.vxc) <= 1.e10
		assert np.fabs(cluster.vyc ) <= 1.e10
		assert np.fabs(cluster.vzc) <= 1.e10

	# variable for galpy orbits
	assert cluster.orbits==None

	# variables for orbital position and kinematics
	if cluster.orbit is None and cluster.ofile is None:
		if cluster.origin=='galaxy':
			assert cluster.xgc != 0.0
			assert cluster.ygc != 0.0
			assert cluster.zgc != 0.0
			assert cluster.vxgc != 0.0
			assert cluster.vygc != 0.0
			assert cluster.vzgc != 0.0
		else:
			assert cluster.xgc == 0.0
			assert cluster.ygc == 0.0
			assert cluster.zgc == 0.0
			assert cluster.vxgc == 0.0
			assert cluster.vygc == 0.0
			assert cluster.vzgc == 0.0
	elif units=='radec' and cluster.orbit is not None:
	    assert cluster.xgc == cluster.orbit.ra()
	    assert cluster.ygc == cluster.orbit.dec()
	    assert cluster.zgc == cluster.orbit.dist()
	    assert cluster.vxgc == cluster.orbit.pmra()
	    assert cluster.vygc == cluster.orbit.pmdec()
	    assert cluster.vzgc == cluster.orbit.vlos()   
	elif cluster.orbit is not None:
	    assert cluster.xgc == cluster.orbit.x()
	    assert cluster.ygc == cluster.orbit.y()
	    assert cluster.zgc == cluster.orbit.z()
	    assert cluster.vxgc == cluster.orbit.vx()
	    assert cluster.vygc == cluster.orbit.vy()
	    assert cluster.vzgc == cluster.orbit.vz()
	elif origin=='galaxy' or cluster.ofile is not None:
		#Test clusters only have x and vy != zero
		assert cluster.xgc != 0.0
		#assert cluster.ygc != 0.0
		#assert cluster.zgc != 0.0
		#assert cluster.vxgc != 0.0
		assert cluster.vygc != 0.0
		#assert cluster.vzgc != 0.0

	# variable for cluster's on-sky coordinates


	assert cluster.ra.size == 0
	assert cluster.dec.size == 0
	assert cluster.dist.size == 0
	assert cluster.pmra.size == 0
	assert cluster.pmdec.size == 0
	assert cluster.vlos.size == 0


	if cluster.orbit is None:
	    assert cluster.ra_gc == 0.0
	    assert cluster.dec_gc == 0.0
	    assert cluster.dist_gc == 0.0
	    assert cluster.pmra_gc == 0.0
	    assert cluster.pmdec_gc == 0.0
	    assert cluster.vlos_gc == 0.0
	else:
	    assert cluster.ra_gc == cluster.orbit.ra()
	    assert cluster.dec_gc == cluster.orbit.dec()
	    assert cluster.dist_gc == cluster.orbit.dist()
	    assert cluster.pmra_gc == cluster.orbit.pmra()
	    assert cluster.pmdec_gc == cluster.orbit.pmdec()
	    assert cluster.vlos_gc == cluster.orbit.vlos()

	assert cluster.ra_c == 0.0
	assert cluster.dec_c == 0.0
	assert cluster.dist_c == 0.0
	assert cluster.pmra_c == 0.0
	assert cluster.pmdec_c == 0.0
	assert cluster.vlos_c == 0.0   

	if 'nbody' in ctype:
		# variables for add_nbody6
		# Number of stars in the core
		assert cluster.nc != 0
		# Core radius
		assert cluster.rc != 0
		# Distance scaling parameter
		assert cluster.rbar != 1.
		assert cluster.rbar_su!=1.
		assert cluster.rbar_au!=1.
		# Tidal limit from NBODY6 (not neccesarily a true tidal radius)
		assert cluster.rtide != 0.
		# Center of mass of cluster (x,yz)
		assert cluster.xc != 0.
		assert cluster.yc != 0.
		assert cluster.zc != 0.
		assert cluster.xcn != None
		assert cluster.ycn != None
		assert cluster.zcn != None
		# Mass scaling parameter
		assert cluster.zmbar != 1.
		# Velocity scaling parameter
		assert cluster.vbar != 1.
		# Time scaling parameter
		assert cluster.tbar != 1.
		assert cluster.tbar_days!=1.
		# Scale radius of cluster
		assert cluster.rscale != 1.
		# Number of single stars
		assert cluster.ns != 0
		# Number of binary stars
		assert cluster.nb == len(cluster.id1)
		# Number of particles (from NBODY6 when tidal tail is being integrated)
		assert cluster.n_p != 0
	else:
		# variables for add_nbody6
		# Number of stars in the core
		assert cluster.nc == 0
		# Core radius
		assert cluster.rc == 0
		# Distance scaling parameter
		assert cluster.rbar == 1.
		assert cluster.rbar_su==1.
		assert cluster.rbar_au==1.
		# Tidal limit from NBODY6 (not neccesarily a true tidal radius)
		assert cluster.rtide == 0.
		# Center of mass of cluster (x,yz)
		assert cluster.xcn == None
		assert cluster.ycn == None
		assert cluster.zcn == None
		# Mass scaling parameter
		assert cluster.zmbar == 1.
		# Velocity scaling parameter
		assert cluster.vbar == 1.
		# Time scaling parameter
		assert cluster.tbar == 1.
		assert cluster.tbar_days==1.
		# Scale radius of cluster
		assert cluster.rscale == 1.
		# Number of single stars
		assert cluster.ns == 0
		# Number of binary stars
		assert cluster.nb == 0
		# Number of particles (from NBODY6 when tidal tail is being integrated)
		assert cluster.n_p == 0

	# variables for add_sse (stellar evolution information)
	if 'nbody' in ctype and cluster.ssefile is not None:

		assert cluster.kw.size != 0.
		assert cluster.logl.size != 0.
		assert cluster.logr.size != 0.
		assert cluster.ep.size != 0.
		assert cluster.ospin.size != 0.
		assert cluster.lum.size != 0.
	else:
		assert cluster.kw.size == 0.
		assert cluster.logl.size == 0.
		assert cluster.logr.size == 0.
		assert cluster.ep.size == 0.
		assert cluster.ospin.size == 0.
		assert cluster.lum.size == 0.


	if 'nbody' in ctype and cluster.bsefile is not None:
		# variables for add_bse (binary star evolution information)
		assert cluster.id1.size != 0
		assert cluster.id2.size != 0
		assert cluster.kw1.size != 0
		assert cluster.kw2.size != 0
		assert cluster.kcm.size != 0
		assert cluster.ecc.size != 0
		assert cluster.pb.size != 0
		assert cluster.semi.size != 0
		assert cluster.m1.size != 0
		assert cluster.m2.size != 0
		assert cluster.logl1.size != 0
		assert cluster.logl2.size != 0
		assert cluster.logr1.size != 0
		assert cluster.logr2.size != 0
		assert cluster.ep1.size != 0
		assert cluster.ep2.size != 0
		assert cluster.ospin1.size != 0
		assert cluster.ospin2.size != 0
	else:
		assert cluster.id1.size == 0
		assert cluster.id2.size == 0
		assert cluster.kw1.size == 0
		assert cluster.kw2.size == 0
		assert cluster.kcm.size == 0
		assert cluster.ecc.size == 0
		assert cluster.pb.size == 0
		assert cluster.semi.size == 0
		assert cluster.m1.size == 0
		assert cluster.m2.size == 0
		assert cluster.logl1.size == 0
		assert cluster.logl2.size == 0
		assert cluster.logr1.size == 0
		assert cluster.logr2.size == 0
		assert cluster.ep1.size == 0
		assert cluster.ep2.size == 0
		assert cluster.ospin1.size == 0
		assert cluster.ospin2.size == 0

	# variables of energies
	if 'nbodypp' in ctype or 'nbody6++' in ctype:
		assert cluster.kin.size != 0
		assert cluster.pot.size != 0
		assert cluster.etot.size != 0
	else:
		assert cluster.kin.size == 0
		assert cluster.pot.size == 0
		assert cluster.etot.size == 0

	# Lagrange Radii,10% lagrage radius, half-mass radius, limiting radius, tidal radius, and virial radius
	assert cluster.rn == None
	assert cluster.r10 != None
	assert cluster.r10pro!=None
	assert cluster.rm != None
	assert cluster.rmpro != None

	if cluster.logl.size!=0:
		assert cluster.rh != None
		assert cluster.rhpro != None
		assert cluster.rl != None

	assert cluster.rt == None
	assert cluster.rv == None

	
	#3D and projected order of stars with respect to origin
	assert len(cluster.rorder) == len(cluster.x)
	assert cluster.rorder_origin == origin
	assert len(cluster.rproorder) == len(cluster.x)

	
	# Additional variables for operation and function calls
	assert cluster.trelax == None
	assert cluster.trh == None
	assert cluster.trc == None
	assert cluster.qv == None
	assert cluster.alpha == None
	assert cluster.eta == None
	assert cluster.rvmax == None
	assert cluster.vmax == None

	#For use with multiple populations
	assert cluster.npop.size == len(cluster.x)

	if 'gyrfalcon' in ctype or 'nemo' in ctype:
		#For use with extended nemo/gyrfalcon output
		if cluster.give == 'mxvpqael':
		    assert cluster.gyrpot.size != 0
		    assert cluster.gyrq.size != 0
		    assert cluster.gyracc.size != 0
		    assert cluster.eps.size != 0
		    assert cluster.gyrlev.size != 0
		elif cluster.give =='mxve':
		    assert cluster.eps.size != 0

	if 'nbody' in ctype and kwargs.get('hdf5',False):
		#For use with HDF5
		assert cluster.hdf5==kwargs.get('hdf5',False)
		assert cluster.ngroups==0
		assert cluster.ngroup==0

	return True

def test_nbody6():
	wdir='../docs/source/notebooks/nbody6_sim/'

	out3 = open("%sOUT3" % wdir, "rb")
	out33 = open("%sOUT33" % wdir, "rb")

	cluster=ctools.load_cluster('nbody6',wdir=wdir)

	assert check_params(cluster,'nbody6','nbody','cluster',False,sfile=out3.name,bfile=out33.name)

	cluster=ctools.advance_cluster(cluster)

	assert check_params(cluster,'nbody6','nbody','cluster',False,sfile=out3.name,bfile=out33.name,tphys=10)


def test_nbody6pp():
	wdir='../docs/source/notebooks/'

	conf3 = open("%sconf.3_0" % wdir, "rb")
	cluster=ctools.load_cluster('nbody6pp',wdir=wdir)
	assert check_params(cluster,'nbody6++','nbody','cluster',False,sfile=conf3.name,wdir=wdir)

	conf3 = open("%sconf.3_20" % wdir, "rb")
	cluster=ctools.advance_cluster(cluster,dtout=20)
	assert check_params(cluster,'nbody6++','nbody','cluster',False,sfile=conf3.name,wdir=wdir,tphys=20,nsnap=20)

def test_snapshot():
	wdir='../docs/source/notebooks/'
	cluster=ctools.load_cluster('snapshot',wdir=wdir,filename='00000.dat',units='pckms',origin='cluster')

	assert check_params(cluster,'snapshot','pckms','cluster',False,sfile='%s00000.dat' % wdir,wdir=wdir)

	cluster=ctools.load_cluster('snapshot',wdir=wdir,filename='00000.dat',units='pckms',origin='cluster',
                            col_names=["m", "x", "y", "z", "vx", "vy", "vz"],col_nums=[0, 1, 2, 3, 4, 5, 6],)

	assert check_params(cluster,'snapshot','pckms','cluster',False,sfile='%s00000.dat' % wdir,wdir=wdir)

	cluster=ctools.load_cluster('snapshot',wdir=wdir,filename='00000.dat',units='pckms',origin='cluster',
	                            col_names=["m", "x", "y", "z", "vx", "vy", "vz"],
	                            col_nums=[0, 1, 2, 3, 4, 5, 6], ofilename='orbit.dat')

	assert check_params(cluster,'snapshot','pckms','cluster',False,sfile='%s00000.dat' % wdir,wdir=wdir,ofile='%sorbit.dat' % wdir,ofilename='orbit.dat')

	cluster=ctools.advance_cluster(cluster,otime=False)

	assert check_params(cluster,'snapshot','pckms','cluster',False,sfile='%s00000.dat' % wdir,wdir=wdir,ofile='%sorbit.dat' % wdir,ofilename='orbit.dat',tphys=0,nsnap=1)

	cluster=ctools.load_cluster('snapshot',wdir=wdir,filename='00000.dat',units='pckms',origin='cluster',
	                            col_names=["m", "x", "y", "z", "vx", "vy", "vz"],
	                            col_nums=[0, 1, 2, 3, 4, 5, 6], ofilename='orbit.dat')
	cluster=ctools.advance_cluster(cluster,otime=True)

	assert check_params(cluster,'snapshot','pckms','cluster',False,sfile='%s00000.dat' % wdir,wdir=wdir,ofile='%sorbit.dat' % wdir,ofilename='orbit.dat',tphys=10,nsnap=1)

def test_gyrfalcon():
	wdir='../docs/source/notebooks/'
	cluster=ctools.load_cluster('gyrfalcon',filename='cluster.nemo.dat',wdir=wdir,units='WDunits',origin='centre')
	assert check_params(cluster,'nemo','WDunits','centre',False,sfile='%scluster.nemo.dat' % wdir,wdir=wdir,skiprows=13)

def test_new_gyrfalcon():
	wdir='../docs/source/notebooks/'
	cluster=ctools.load_cluster('new_gyrfalcon',filename='cluster.nemo.dat',wdir=wdir,units='WDunits',origin='centre')
	assert check_params(cluster,'new_nemo','WDunits','centre',False,sfile='%scluster.nemo.dat' % wdir,wdir=wdir,skiprows=13)

	cluster=ctools.advance_cluster(cluster,filename=cluster.sfile)


	assert check_params(cluster,'new_nemo','WDunits','centre',False,sfile='%scluster.nemo.dat' % wdir,wdir=wdir,skiprows=13,nsnap=1,tphys=0.484375)

	cluster=ctools.advance_cluster(cluster,filename=cluster.sfile)

	assert check_params(cluster,'new_nemo','WDunits','centre',False,sfile='%scluster.nemo.dat' % wdir,wdir=wdir,skiprows=13,nsnap=2,tphys=2.*0.484375)

def test_amuse():
	N=100
	Mcluster=100.0 | units.MSun
	Rcluster= 1.0 | units.parsec
	converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
	stars=new_plummer_sphere(N,converter)

	cluster=ctools.load_cluster('amuse',particles=stars,units='pckms',origin='centre')


	assert check_params(cluster,'amuse','pckms','centre',False)

def test_galpy(tol=0.01):
	kdf= kingdf(M=2.3,rt=1.4,W0=3.)
	sam= kdf.sample(n=1000)

	cluster=ctools.load_cluster('galpy',particles=sam,units='kpckms',origin='cluster')
	assert check_params(cluster,'galpy','kpckms','cluster',False)

	gcluster=ctools.load_cluster('galpy',particles=sam,units='galpy',origin='cluster',ro=1.,vo=1.,)
	assert check_params(gcluster,'galpy','galpy','cluster',False)

	assert np.all(np.fabs(cluster.r/gcluster.r-solar_ro) < tol)
	assert np.all(np.fabs(cluster.v/gcluster.v-solar_vo) < tol)

def test_astropy_table():
	wdir='../docs/source/notebooks/'
	data = QTable.read("%spal5_rp.dat" % wdir, format="ascii")
	cluster = ctools.load_cluster('astropy_table',particles=data, units='kpckms',origin='galaxy')
	assert check_params(cluster,'astropy_table','kpckms','galaxy',False)

def test_mcluster(tol=0.2):
	wdir='../docs/source/notebooks/'
	rvirp=1.20395
	rhp=1.0
	rplummer=0.76628

	m,x,y,z,vx,vy,vz=np.loadtxt(wdir+'N1k.dat.10',unpack=True)
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m)

	r_v=ctools.virial_radius(cluster)
	print(r_v,rvirp,cluster.rm,rhp)

	assert np.fabs(r_v-rvirp) <= tol
	assert np.fabs(cluster.rm-rhp) <= tol

	m,x,y,z,vx,vy,vz=np.loadtxt(wdir+'N1kb.dat.10',unpack=True)
	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m=m,nb=500)

	r_v=ctools.virial_radius(cluster)
	print(r_v,rvirp,cluster.rm,rhp)

	assert np.fabs(r_v-rvirp) <= tol
	assert np.fabs(cluster.rm-rhp) <= tol

