import clustertools as ctools
from clustertools.analysis.functions import tpl_func
from scipy.optimize import curve_fit

import numpy as np
from galpy.orbit import Orbit
from galpy.potential import NFWPotential,MWPotential2014,rtide,evaluateDensities
from galpy.df import isotropicNFWdf
from galpy.util import bovy_conversion

def test_rho_prof(tol=0.1):
	#Create a homogenous cluster so rho should be the same everywhere
	rho0=1.0

	rlower=np.arange(0,20,1)
	rupper=np.arange(1,21,1)
	rmid=(rupper+rlower)/2.

	vol=4.0*np.pi*(rupper**3.-rlower**3.)/3.
	m0=rho0*vol

	x=np.array([])
	for i,m in enumerate(m0):
		x=np.append(x,np.random.uniform(rlower[i],rupper[i],int(m)))

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vx,vy,vz=np.zeros(n),np.zeros(n),np.zeros(n)
	m=np.ones(n)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)

	#Test fixed bin
	rprof, pprof, nprof=ctools.rho_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	assert np.all(np.fabs(pprof/np.mean(pprof)-1) < tol)
	assert np.all(np.fabs(rprof/((rlower+rupper)/2.)-1) < tol)

	#test variable bins
	rprof, pprof, nprof=ctools.rho_prof(cluster)
	assert np.all(np.fabs(nprof/np.mean(nprof)-1) < tol)


	#Test Normalize
	rprof_norm, pprof_norm, nprof_norm=ctools.rho_prof(cluster,normalize=True)
	assert np.all(np.fabs(rprof_norm/(rprof/cluster.rm)-1) < tol)
	assert np.all(np.fabs(pprof/pprof_norm-1) < tol)

	#Test projected
	cluster.z=np.random.rand(len(cluster.x))
	rprof_pro, pprof_pro, nprof_pro=ctools.rho_prof(cluster,projected=True)

	assert np.all(np.fabs(rprof/rprof_pro-1) < tol)
	assert np.all(np.fabs(pprof_pro/pprof_pro-1) < tol)
	assert np.all(np.fabs(nprof_pro/nprof_pro-1) < tol)


	#test subcluster
	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=0.5,1.5
	vmin,vmax=-0.5,0.5

	mextra=np.append(np.random.uniform(0,mmin-0.1,500),np.random.uniform(mmax+0.1,2.0,500))

	xextra=np.random.uniform(np.amin(cluster.x),np.amax(cluster.x),len(mextra))
	yextra=np.zeros(len(mextra))
	zextra=np.zeros(len(mextra))


	vxextra=np.append(np.random.uniform(-1.,vmin-0.1,500),np.random.uniform(vmax+0.1,1.0,500))
	vyextra=np.zeros(len(mextra))
	vzextra=np.zeros(len(mextra))

	kwextra=np.ones(len(mextra))*12
	eextra=np.append(np.random.uniform(emin-1.,emin-0.1,500),np.random.uniform(emax+0.1,emax+1.0,500))

	cluster.add_stars(xextra,yextra,zextra,vxextra,vyextra,vzextra,mextra)
	cluster.kw=np.append(np.zeros(n),kwextra)
	cluster.etot=np.append(cluster.etot,eextra)
	cluster.analyze()

	print(cluster.ntot,len(cluster.etot),np.sum(cluster.etot<emin),np.sum(cluster.etot>emax))

	rprof_ex, pprof_ex, nprof_ex=ctools.rho_prof(cluster,mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprof/rprof_ex-1) < tol)
	assert np.all(np.fabs(pprof/pprof_ex-1) < tol)
	assert np.all(np.fabs(nprof/nprof_ex-1) < tol)

def test_m_prof(tol=0.1):
	#Create a homogenous cluster so rho should be the same everywhere
	rho0=1.0

	rlower=np.arange(0,20,1)
	rupper=np.arange(1,21,1)
	rmid=(rupper+rlower)/2.

	vol=4.0*np.pi*(rupper**3.-rlower**3.)/3.
	m0=rho0*vol

	x=np.array([])
	m=np.array([])
	for i,ms in enumerate(m0):
		x=np.append(x,np.random.uniform(rlower[i],rupper[i],int(ms)))
		mnew=np.ones(int(ms))*float(ms/int(ms))
		m=np.append(m,mnew)
		print(ms,np.sum(mnew))

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vx,vy,vz=np.zeros(n),np.zeros(n),np.zeros(n)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)

	#Test fixed bin
	rprof, mprof, nprof=ctools.m_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))


	assert np.all(np.fabs(mprof/m0-1) < tol)
	assert np.all(np.fabs(rprof/((rlower+rupper)/2.)-1) < tol)


	#test variable bins
	rprof, mprof, nprof=ctools.m_prof(cluster)
	assert np.all(np.fabs(nprof/np.mean(nprof)-1) < tol)


	#Test Normalize
	rprof_norm, mprof_norm, nprof_norm=ctools.m_prof(cluster,normalize=True)
	assert np.all(np.fabs(rprof_norm/(rprof/cluster.rm)-1) < tol)
	assert np.all(np.fabs(mprof/mprof_norm-1) < tol)


	#Test projected
	cluster.z=np.random.rand(len(cluster.x))
	rprof_pro, mprof_pro, nprof_pro=ctools.m_prof(cluster,projected=True)

	assert np.all(np.fabs(rprof/rprof_pro-1) < tol)
	assert np.all(np.fabs(mprof_pro/mprof_pro-1) < tol)
	assert np.all(np.fabs(nprof_pro/nprof_pro-1) < tol)

	

	#test subcluster
	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=-0.5,0.5

	mextra=np.append(np.random.uniform(0,mmin-0.1,500),np.random.uniform(mmax+0.1,2.0,500))

	xextra=np.random.uniform(np.amin(cluster.x),np.amax(cluster.x),len(mextra))
	yextra=np.zeros(len(mextra))
	zextra=np.zeros(len(mextra))


	vxextra=np.append(np.random.uniform(-1.,vmin-0.1,500),np.random.uniform(vmax+0.1,1.0,500))
	vyextra=np.zeros(len(mextra))
	vzextra=np.zeros(len(mextra))

	kwextra=np.ones(len(mextra))*12
	eextra=np.append(np.random.uniform(emin-1.,emin-0.1,500),np.random.uniform(emax+0.1,emax+1.0,500))

	cluster.add_stars(xextra,yextra,zextra,vxextra,vyextra,vzextra,mextra)
	cluster.kw=np.append(np.zeros(n),kwextra)
	cluster.etot=np.append(cluster.etot,eextra)
	cluster.analyze()

	print(cluster.ntot,len(cluster.etot),np.sum(cluster.etot<emin),np.sum(cluster.etot>emax))

	rprof_ex, mprof_ex, nprof_ex=ctools.m_prof(cluster,mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprof/rprof_ex-1) < tol)
	assert np.all(np.fabs(mprof/mprof_ex-1) < tol)
	assert np.all(np.fabs(nprof/nprof_ex-1) < tol)

	


