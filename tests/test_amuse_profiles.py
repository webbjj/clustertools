import clustertools as ctools
from clustertools.analysis.functions import tpl_func
from scipy.optimize import curve_fit

import numpy as np
from galpy.orbit import Orbit
from galpy.potential import NFWPotential,MWPotential2014,rtide,evaluateDensities
from galpy.df import isotropicNFWdf

import pytest

try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion

try:
	import amuse.units.units as u
	from amuse.units.quantities import ScalarQuantity,VectorQuantity
	noamuse=False
except:
	noamuse=True


solar_motion=[-11.1,12.24,7.25] #Schönrich, R., Binney, J., Dehnen, W., 2010, MNRAS, 403, 1829
solar_ro=8.275 #Gravity Collaboration, Abuter, R., Amorim, A., et al. 2020 ,A&A, 647, A59
solar_vo=solar_ro*30.39-solar_motion[1]

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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
	cluster.to_amuse()

	#Test fixed bin
	rprof, pprof, nprof=ctools.rho_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	assert np.all(np.fabs(pprof/np.mean(pprof)-1) < tol)
	assert np.all(np.fabs(rprof.value_in(u.parsec)/((rlower+rupper)/2.)-1) < tol)

	#test variable bins
	rprof, pprof, nprof=ctools.rho_prof(cluster)
	assert np.all(np.fabs(nprof/np.mean(nprof)-1) < tol)


	#Test Normalize
	rprof_norm, pprof_norm, nprof_norm=ctools.rho_prof(cluster,normalize=True)
	assert np.all(np.fabs(rprof_norm/(rprof/cluster.rm)-1) < tol)
	assert np.all(np.fabs(pprof/pprof_norm-1) < tol)

	#Test projected
	cluster.z=np.random.rand(len(cluster.x)) | u.parsec
	rprof_pro, pprof_pro, nprof_pro=ctools.rho_prof(cluster,projected=True)

	assert np.all(np.fabs(rprof/rprof_pro-1) < tol)
	assert np.all(np.fabs(pprof_pro/pprof_pro-1) < tol)
	assert np.all(np.fabs(nprof_pro/nprof_pro-1) < tol)


	#test subcluster
	cluster.to_pckms()

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

	cluster.to_amuse()
	cluster.analyze()

	print(cluster.etot,emin,emax)

	print(cluster.ntot,len(cluster.etot),np.sum(cluster.etot<emin),np.sum(cluster.etot>emax))

	rprof_ex, pprof_ex, nprof_ex=ctools.rho_prof(cluster,mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprof/rprof_ex-1) < tol)
	assert np.all(np.fabs(pprof/pprof_ex-1) < tol)
	assert np.all(np.fabs(nprof/nprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
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
	cluster.to_amuse()

	#Test fixed bin
	rprof, mprof, nprof=ctools.m_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))


	assert np.all(np.fabs(mprof.value_in(u.MSun)/m0-1) < tol)
	assert np.all(np.fabs(rprof.value_in(u.parsec)/((rlower+rupper)/2.)-1) < tol)


	#test variable bins
	rprof, mprof, nprof=ctools.m_prof(cluster)
	assert np.all(np.fabs(nprof/np.mean(nprof)-1) < tol)


	#Test Normalize
	rprof_norm, mprof_norm, nprof_norm=ctools.m_prof(cluster,normalize=True)
	assert np.all(np.fabs(rprof_norm/(rprof/cluster.rm)-1) < tol)
	assert np.all(np.fabs(mprof/mprof_norm-1) < tol)


	#Test projected
	cluster.z=np.random.rand(len(cluster.x)) | u.parsec
	rprof_pro, mprof_pro, nprof_pro=ctools.m_prof(cluster,projected=True)

	assert np.all(np.fabs(rprof/rprof_pro-1) < tol)
	assert np.all(np.fabs(mprof/mprof_pro-1) < tol)
	assert np.all(np.fabs(nprof/nprof_pro-1) < tol)

	

	#test subcluster
	cluster.to_pckms()
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

	cluster.to_amuse()

	cluster.analyze()

	print(cluster.ntot,len(cluster.etot),np.sum(cluster.etot<emin),np.sum(cluster.etot>emax))

	rprof_ex, mprof_ex, nprof_ex=ctools.m_prof(cluster,mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprof/rprof_ex-1) < tol)
	assert np.all(np.fabs(mprof/mprof_ex-1) < tol)
	assert np.all(np.fabs(nprof/nprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')	
def test_alpha_prof(tol=0.2):
	#Create a homogenous cluster so rho should be the same everywhere
	rho0=1.0

	mmin=0.1
	mmax=1.0
	mbins=np.arange(mmin,mmax+0.1,0.1)

	dm=0.1/2

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.

	alphas=np.arange(0,1,0.1)

	x=np.array([])
	m=np.array([])
	nbin=1000
	for i,a in enumerate(alphas):
		mfrac=mbins**a
		mfrac/=np.sum(mfrac)
		mtot=(mfrac*nbin).astype(int)
		mnew=np.array([])
		for j in range(0,len(mtot)):
			mnew=np.append(mnew,np.ones(mtot[j])*mbins[j])

		m=np.append(m,mnew+np.random.uniform(-dm,dm,len(mnew)))

		xnew=np.random.uniform(rlower[i],rupper[i],len(mnew))
		x=np.append(x,xnew)

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vx,vy,vz=np.zeros(n),np.zeros(n),np.zeros(n)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True,sortstars=True)

	da0=np.polyfit(np.log(rmid/cluster.rm),alphas,1)[0]

	cluster.to_amuse()
	cluster.analyze()

	#Test fixed mass and radius bins
	rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='fix')
	assert np.all(np.fabs(alphas-aprof) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.parsec)/((rlower+rupper)/2.)-1) < tol)
	assert np.fabs(dalpha/da0-1) < tol

	#Test number mass bin
	rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='num')
	assert np.all(np.fabs(alphas-aprof) < tol)

	#Test number radius bin
	rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = ctools.alpha_prof(cluster,bintype='num',nrad=10,nmass=10,mbintype='fix')
	assert np.all(np.fabs(alphas-aprof) < tol)

	#Test Normalize
	rprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='fix')
	rprofn_norm, aprof_norm, dalpha_norm, edalpha_norm, ydalpha_norm, eydalpha_norm=ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='fix',normalize=True)

	assert np.all(np.fabs(rprofn/(rprofn_norm*cluster.rm.value_in(u.pc))-1) < tol)
	assert np.all(np.fabs(aprof/aprof_norm-1) < tol)


	#test_mcorr
	mcorr=np.ones(cluster.ntot)
	rprofn_corr, aprof_corr, dalpha_corr, edalpha_corr, ydalpha_corr, eydalpha_corr=ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='fix',normalize=False,mcorr=mcorr)
	assert np.all(np.fabs(aprof/aprof_corr-1) < tol)

	#Test projected
	cluster.z=np.random.rand(len(cluster.x)) | u.parsec
	rprofn_pro, aprof_pro, dalpha_pro, edalpha_pro, ydalpha_pro, eydalpha_pro=ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='fix',projected=True)

	assert np.all(np.fabs(rprofn/rprofn_pro-1) < tol)
	assert np.all(np.fabs(aprof/aprof_pro-1) < tol)
	assert np.fabs(dalpha_pro/dalpha-1) < tol


	#test subcluster
	cluster.z=np.zeros(cluster.ntot) | u.parsec
	cluster.to_pckms()
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

	cluster.to_amuse()
	cluster.analyze()


	rprofn_ex, aprof_ex, dalpha_ex, edalpha_ex, ydalpha_ex, eydalpha_ex=ctools.alpha_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),nrad=10,nmass=10,mbintype='fix',mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(aprof/aprof_ex-1) < tol)
	assert np.fabs(dalpha_ex/dalpha-1) < tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_sigv_prof(tol=0.1):

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.
	stds=np.arange(1,11,1)

	print(len(rupper),len(stds))

	x=np.array([])
	vx=np.array([])

	n=100000
	num=int(n/len(rmid))

	for i in range(0,len(rmid)):

		x=np.append(x,np.random.uniform(rlower[i],rupper[i],num))
		vnew=np.random.normal(0.,stds[i],num)
		vx=np.append(vx,vnew)

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vy,vz=np.zeros(n),np.zeros(n)
	m=np.ones(n)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)
	cluster.to_amuse()

	#Test fixed bin
	rprofn, sigvprof=ctools.sigv_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx')

	assert np.all(np.fabs(sigvprof.value_in(u.kms)/stds-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test num bin
	rprofn, sigvprof=ctools.sigv_prof(cluster,bintype='num',nrad=10,coord='vx')

	assert np.all(np.fabs(sigvprof.value_in(u.kms)/stds-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test normalize
	rprofn_norm, sigvprof_norm=ctools.sigv_prof(cluster,bintype='num',nrad=10,coord='vx',normalize=True)

	assert np.all(np.fabs(sigvprof_norm.value_in(u.kms)/stds-1) < tol)
	assert np.all(np.fabs((cluster.rm.value_in(u.pc))*rprofn_norm.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#test projected
	cluster.z=np.random.rand(cluster.ntot) | u.parsec
	rprofn_pro, sigvprof_pro=ctools.sigv_prof(cluster,bintype='num',nrad=10,coord='vx',normalize=False,projected=True)

	assert np.all(np.fabs(rprofn/rprofn_pro-1) < tol)
	assert np.all(np.fabs(sigvprof/sigvprof_pro-1) < tol)

	#test subcluster
	cluster.z=np.zeros(cluster.ntot) | u.parsec
	rprofn, sigvprof=ctools.sigv_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx')

	cluster.to_pckms()
	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=np.amin(cluster.v),np.amax(cluster.v)

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

	cluster.to_amuse()
	cluster.analyze()

	rprofn_ex, sigvprof_ex=ctools.sigv_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx',mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	print(sigvprof)
	print(sigvprof_ex)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(sigvprof/sigvprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_beta_prof(tol=0.1):

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.
	stds=np.arange(1,11,1)

	print(len(rupper),len(stds))

	r=np.array([])
	vr=np.array([])

	n=100000
	num=int(n/len(rmid))

	for i in range(0,len(rmid)):

		r=np.append(r,np.random.uniform(rlower[i],rupper[i],num))
		vnew=np.random.normal(0.,stds[i],num)
		vr=np.append(vr,vnew)

	n=len(r)
	phi=2.0*np.pi*np.random.rand(n)
	theta=np.arccos(1.0-2.0*np.random.rand(n))
	vp,vt=vr,vr
	m=np.ones(n)

	x,y,z,vx,vy,vz=ctools.sphere_to_cart(r,phi,theta,vr,vp,vt)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)
	cluster.to_amuse()

	#Test fixed bin
	rprofn, bprof=ctools.beta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	assert np.all(np.fabs(bprof) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test num bin
	rprofn, bprof=ctools.beta_prof(cluster,bintype='num',nrad=10)

	assert np.all(np.fabs(bprof) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test normalize
	rprofn_norm, bprof_norm=ctools.beta_prof(cluster,bins=(rlower,rmid,rupper),normalize=True)
	assert np.all(np.fabs(cluster.rm.value_in(u.parsec)*rprofn_norm/((rlower+rupper)/2.)-1) < tol)

	#test subcluster
	rprofn, bprof=ctools.beta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	cluster.to_pckms()
	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=np.amin(cluster.v),np.amax(cluster.v)

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
	cluster.to_amuse()
	cluster.analyze()

	rprofn_ex, bprof_ex=ctools.beta_prof(cluster,bins=(rlower,rmid,rupper),mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(np.exp(bprof)/np.exp(bprof_ex)-1) < tol)
	
	#test projected
	n=len(r)
	phi=2.0*np.pi*np.random.rand(n)
	theta=np.arccos(1.0-2.0*np.random.rand(n))
	z=r*np.sin(theta)
	rad=r*np.cos(theta)
	vp,vz=vr,vz
	m=np.ones(n)

	x,y,z,vx,vy,vz=ctools.cyl_to_cart(rad,phi,z,vr,vp,vz)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)
	cluster.x=np.sqrt(cluster.x**2+cluster.z**2.)
	cluster.z=0
	cluster.to_amuse()
	cluster.analyze()


	rprofn, bprof=ctools.beta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),projected=True)

	assert np.all(np.fabs(bprof) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_v_prof(tol=0.1):

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.
	vs=np.arange(1,11,1)
	stds=np.ones(len(vs))*0.1

	print(len(rupper),len(stds))

	x=np.array([])
	vx=np.array([])

	n=100000
	num=int(n/len(rmid))

	for i in range(0,len(rmid)):

		x=np.append(x,np.random.uniform(rlower[i],rupper[i],num))
		vnew=np.random.normal(vs[i],stds[i],num)
		vx=np.append(vx,vnew)

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vy,vz=np.zeros(n),np.zeros(n)
	m=np.ones(n)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)
	cluster.to_amuse()

	#Test fixed bin
	rprofn, vprof=ctools.v_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx')

	assert np.all(np.fabs(vprof.value_in(u.kms)/vs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test num bin
	rprofn, vprof=ctools.v_prof(cluster,bintype='num',nrad=10,coord='vx')

	assert np.all(np.fabs(vprof.value_in(u.kms)/vs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test normalize
	rprofn_norm, vprof_norm=ctools.v_prof(cluster,bintype='num',nrad=10,coord='vx',normalize=True)

	assert np.all(np.fabs(vprof_norm.value_in(u.kms)/vs-1) < tol)
	assert np.all(np.fabs(cluster.rm.value_in(u.pc)*rprofn_norm/((rlower+rupper)/2.)-1) < tol)

	#test projected
	cluster.z=np.random.rand(cluster.ntot) | u.parsec
	rprofn_pro, vprof_pro=ctools.v_prof(cluster,bintype='num',nrad=10,coord='vx',normalize=False,projected=True)

	assert np.all(np.fabs(rprofn/rprofn_pro-1) < tol)
	assert np.all(np.fabs(vprof/vprof_pro-1) < tol)

	#test subcluster

	cluster.z=np.zeros(cluster.ntot) | u.parsec
	rprofn, vprof=ctools.v_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx')

	cluster.to_pckms()

	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=np.amin(cluster.v),np.amax(cluster.v)

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

	cluster.to_amuse()

	cluster.analyze()

	rprofn_ex, vprof_ex=ctools.v_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx',mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	print(rprofn)
	print(rprofn_ex)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(vprof/vprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_v2_prof(tol=0.1):

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.
	vs=np.arange(1,11,1)
	stds=np.ones(len(vs))*0.1

	print(len(rupper),len(stds))

	x=np.array([])
	vx=np.array([])

	n=100000
	num=int(n/len(rmid))

	for i in range(0,len(rmid)):

		x=np.append(x,np.random.uniform(rlower[i],rupper[i],num))
		vnew=np.random.normal(vs[i],stds[i],num)
		vx=np.append(vx,vnew)

	vs=vs**2.

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vy,vz=np.zeros(n),np.zeros(n)
	m=np.ones(n)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)

	cluster.to_amuse()

	#Test fixed bin
	rprofn, vprof=ctools.v2_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx')


	assert np.all(np.fabs(vprof.value_in(u.kms*u.kms)/vs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test num bin
	rprofn, vprof=ctools.v2_prof(cluster,bintype='num',nrad=10,coord='vx')

	assert np.all(np.fabs(vprof.value_in(u.kms*u.kms)/vs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/((rlower+rupper)/2.)-1) < tol)

	#Test normalize
	rprofn_norm, vprof_norm=ctools.v2_prof(cluster,bintype='num',nrad=10,coord='vx',normalize=True)

	assert np.all(np.fabs(vprof_norm.value_in(u.kms*u.kms)/vs-1) < tol)
	assert np.all(np.fabs(cluster.rm.value_in(u.parsec)*rprofn_norm/((rlower+rupper)/2.)-1) < tol)

	#test projected

	cluster.z=np.random.rand(cluster.ntot) | u.parsec
	rprofn_pro, vprof_pro=ctools.v2_prof(cluster,bintype='num',nrad=10,coord='vx',normalize=False,projected=True)

	assert np.all(np.fabs(rprofn/rprofn_pro-1) < tol)
	assert np.all(np.fabs(vprof/vprof_pro-1) < tol)

	#test subcluster
	cluster.z=np.zeros(cluster.ntot) | u.pc
	rprofn, vprof=ctools.v2_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx')

	cluster.to_pckms()
	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=np.amin(cluster.v),np.amax(cluster.v)

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

	cluster.to_amuse()
	cluster.analyze()

	rprofn_ex, vprof_ex=ctools.v2_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),coord='vx',mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	print(rprofn)
	print(vprof)
	print(rprofn_ex)
	print(vprof_ex)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(vprof/vprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_eta_prof(tol=0.1):

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.

	x=np.array([])
	vx=np.array([])
	m=np.array([])

	n=100000
	num=int(n/len(rmid))
	alphatest=-1

	#etatest=-0.25
	etatest=np.arange(-2,-0.,0.2)

	for i in range(0,len(rmid)):
		xnew=np.random.uniform(rlower[i],rupper[i],num)
		x=np.append(x,xnew)

		mnew=ctools.power_law_distribution_function(len(xnew),alphatest,0.1,1.)
		m=np.append(m,mnew)

		for j in range(0,len(xnew)):
			sig=mnew[j]**etatest[i]
			vx=np.append(vx,np.random.normal(0.,sig))

	y,z=np.zeros(n),np.zeros(n)
	vy,vz=np.zeros(n),np.zeros(n)

	cluster=ctools.StarCluster(units='pckms',origin='centre')
	cluster.add_stars(x,y,z,vx,vy,vz,m)
	cluster.analyze()

	deta0=np.polyfit(np.log(rmid/cluster.rm),etatest,1)[0]
	deta0_norm=np.polyfit(np.log(rmid),etatest,1)[0]

	cluster.to_amuse()
	cluster.analyze()

	rprofn, eprof, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	assert np.all(np.fabs(eprof/etatest-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.parsec)/((rlower+rupper)/2.)-1) < tol)
	assert np.fabs(deta0/deta-1) < tol

	#Test num bins

	rprofn, eprof, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,bintype='num',nrad=10)

	assert np.all(np.fabs(eprof/etatest-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.parsec)/((rlower+rupper)/2.)-1) < tol)
	assert np.fabs(deta0/deta-1) < tol

	#Test normalize
	rprofn, eprof, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))
	rprofn_norm, eprof_norm, deta_norm, edeta_norm, ydeta_norm, eydeta_norm=ctools.eta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),normalize=True)

	print(rprofn)
	print(rprofn_norm)
	print(cluster.rm)

	assert np.all(np.fabs(eprof/eprof_norm-1) < tol)
	assert np.all(np.fabs(cluster.rm.value_in(u.parsec)*rprofn_norm/((rlower+rupper)/2.)-1) < tol)
	assert np.fabs(deta0_norm/deta_norm-1) < tol

	#test projected
	cluster.z=np.random.rand(cluster.ntot) | u.parsec
	cluster.analyze()
	rprofn_pro, eprof_pro, deta_pro, edeta_pro, ydeta_pro, eydeta_pro=ctools.eta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),projected=True)

	assert np.all(np.fabs(rprofn/rprofn_pro-1) < tol)
	assert np.all(np.fabs(eprof/eprof_pro-1) < tol)

	#test subcluster
	cluster.z=np.zeros(cluster.ntot) | u.parsec
	rprofn, eprof, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	cluster.to_pckms()

	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=np.amin(cluster.v),np.amax(cluster.v)

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

	cluster.to_amuse()
	cluster.analyze()

	rprofn_ex, eprof_ex, deta_ex, edeta_ex, ydeta_ex, eydeta_ex=ctools.eta_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper),mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(eprof/eprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_vcirc_prof(tol=0.1):

	rlower=np.arange(0,20,2)
	rupper=rlower+2
	rmid=(rupper+rlower)/2.
	
	x=np.array([])

	n=100000
	num=int(n/len(rmid))

	for i in range(0,len(rmid)):

		x=np.append(x,np.random.uniform(rlower[i],rupper[i],num))

	n=len(x)
	y=np.zeros(n)
	z=np.zeros(n)
	vx,vy,vz=np.zeros(n),np.zeros(n),np.zeros(n)
	m=np.ones(n)

	grav=4.302e-3
	numsum=np.zeros(len(rmid))
	for i in range(0,len(numsum)):
		numsum[i]=np.sum(x<=rmid[i])


	vcs=np.sqrt(grav*numsum/rmid)

	cluster=ctools.StarCluster(units='pckms',origin='cluster')
	cluster.add_stars(x,y,z,vx,vy,vz,m,analyze=True)

	cluster.to_amuse()

	#Test fixed bin
	rprofn,vprof,rvmax,vmax =ctools.vcirc_prof(cluster,bintype='fix',bins=(rlower,rmid,rupper))

	assert np.all(np.fabs(vprof.value_in(u.kms)/vcs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/rmid-1) < tol)

	
	rprofn,vprof,rvmax,vmax =ctools.vcirc_prof(cluster,bintype='fix',nrad=10)
	assert np.all(np.fabs(vprof.value_in(u.kms)/vcs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.pc)/rmid-1) < tol)

	
	#Test num bin
	rprofn,vprof,rvmax,vmax =ctools.vcirc_prof(cluster,bintype='num',nrad=10)

	assert np.all(np.fabs(vprof.value_in(u.kms)/vcs-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.parsec)/((rlower+rupper)/2.)-1) < tol)

	#Test full
	rprofn,vprof,rvmax,vmax =ctools.vcirc_prof(cluster,nrad=None)

	assert np.all(np.fabs(vprof.value_in(u.kms)/np.sqrt(grav*np.cumsum(m)/np.sort(x))-1) < tol)
	assert np.all(np.fabs(rprofn.value_in(u.parsec)/(np.sort(x))-1) < tol)
	assert vmax.value_in(u.kms)==np.amax(np.sqrt(grav*np.cumsum(m)/np.sort(x)))
	assert rvmax.value_in(u.pc)==np.sort(x)[np.argmax(np.sqrt(grav*np.cumsum(m)/np.sort(x)))]

	#Test normalize
	rprofn,vprof,rvmax,vmax =ctools.vcirc_prof(cluster,bintype='fix',nrad=10)
	rprofn_norm,vprof_norm,rvmax_norm,vmax_norm =ctools.vcirc_prof(cluster,bintype='fix',nrad=10,normalize=True)

	assert np.all(np.fabs(vprof_norm/vprof-1) < tol)
	assert np.all(np.fabs(cluster.rm*rprofn_norm/rprofn-1) < tol)

	#test projected
	cluster.z=np.random.rand(cluster.ntot) | u.parsec
	cluster.analyze()
	rprofn_pro,vprof_pro,rvmax_pro,vmax_pro =ctools.vcirc_prof(cluster,bintype='fix',nrad=10,projected=True)

	assert np.all(np.fabs(rprofn/rprofn_pro-1) < tol)
	assert np.all(np.fabs(vprof/vprof_pro-1) < tol)

	#test subcluster
	cluster.z=np.zeros(cluster.ntot) | u.parsec
	cluster.analyze()
	rprofn,vprof,rvmax,vmax =ctools.vcirc_prof(cluster,bintype='fix',nrad=10)

	cluster.to_pckms()

	cluster.energies()
	emin=np.amin(cluster.etot)
	emax=np.amax(cluster.etot)

	mmin,mmax=np.amin(cluster.m),np.amax(cluster.m)
	vmin,vmax=np.amin(cluster.v),np.amax(cluster.v)

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

	cluster.to_amuse()
	cluster.analyze()

	rprofn_ex,vprof_ex,rvmax_ex,vmax_ex =ctools.vcirc_prof(cluster,bintype='fix',nrad=10, mmin=mmin,mmax=mmax,kwmax=1,emin=emin,emax=emax,vmin=vmin,vmax=vmax)

	assert np.all(np.fabs(rprofn/rprofn_ex-1) < tol)
	assert np.all(np.fabs(vprof/vprof_ex-1) < tol)

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_amuse_profile_units(tol=0.1):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_kpckms_galaxy_mbar1.dat',units='kpckms',origin='galaxy')
	cluster.add_orbit(-2.0052100994789871, -9.3481814843660959, -3.9454681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)

	cluster.m=np.random.rand(cluster.ntot)
	cluster.to_cluster()
	cluster.reset_nbody_scale(rvirial=False)

	init_units=['pckms','kpckms','pcmyr','kpcgyr','galpy','nbody','WDunits']
	mcons=[1.0,1.0,1.0,1.0,1./conversion.mass_in_msol(ro=solar_ro, vo=solar_vo),1./cluster.zmbar,1./222288.4543021174]
	pcons=[1.0,1.0,1.0,1.0,1./conversion.dens_in_msolpc3(ro=solar_ro, vo=solar_vo),(cluster.rbar**3.)/cluster.zmbar,1./222288.4543021174]
	rcons=[1.0,1.0,1.0,1.0,1./solar_ro,1.0/cluster.rbar,1.]
	vcons=[1.0,1.0,1.0,1.0,1./solar_vo,1.0/cluster.vbar,1.]

	#vcon=[1.,1.0,1.022712165045695,1.022712165045695,solar_vo,cluster.vbar,1.022712165045695]
	#rcon=[1.,1000.,1.,1000.,solar_ro,cluster.rbar,1000.0]

	projected=False

	for i in range(0,len(init_units)):
		print('DEBUG : ',i,init_units[i])

		if init_units[i]=='pckms':
			munit = u.MSun
			vunit = u.kms
			runit = u.pc
		elif init_units[i]=='kpckms':
			munit = u.MSun
			vunit = u.kms
			runit = u.kpc
		elif init_units[i]=='pcmyr':
			munit = u.MSun
			vunit = u.pc/u.Myr
			runit = u.pc
		elif init_units[i]=='kpcgyr':
			munit = u.MSun
			vunit = u.kpc/u.Gyr
			runit = u.kpc
		elif init_units[i]=='galpy':
			munit = u.MSun
			vunit = u.kms
			runit = u.kpc
		elif init_units[i]=='nbody':
			munit = u.MSun
			vunit = u.kms
			runit = u.pc
		elif init_units[i]=='WDunits':
			munit = u.MSun
			vunit = u.kpc/u.Gyr
			runit = u.kpc

		if init_units[i]=='galpy':
			punit=munit/((u.pc)**3.0)
		else:
			punit=munit/(runit**3.0)


		cluster.to_units(init_units[i])

		rprof0, prof0, nprof0=ctools.rho_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, nprof=ctools.rho_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(punit)*pcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0, nprof0=ctools.m_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, nprof=ctools.m_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(munit)*mcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0, dalpha, edalpha, ydalpha, eydalpha=ctools.alpha_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, dalpha, edalpha, ydalpha, eydalpha=ctools.alpha_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.sigv_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.sigv_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit)*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.beta_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.beta_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.v_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.v_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit)*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.v2_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.v2_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit*vunit)*vcons[i]*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])


		rprof0, prof0, rvmax0, vmax0=ctools.vcirc_prof(cluster,projected=projected)
		cluster.to_amuse(analyze=True)
		rprof, prof, rvmax, vmax=ctools.vcirc_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit)*vcons[i]
		rvmax=rvmax.value_in(runit)*rcons[i]
		vmax=vmax.value_in(vunit)*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		assert np.fabs(1.0-rvmax0/rvmax) < tol
		assert np.fabs(1.0-vmax0/vmax) < tol

@pytest.mark.skipif(noamuse, reason='amuse required to run this test')
def test_amuse_profile_units_projected(tol=0.1):
	cluster=ctools.load_cluster(ctype='snapshot',filename='ngc6101_kpckms_galaxy_mbar1.dat',units='kpckms',origin='galaxy')
	cluster.add_orbit(-2.0052100994789871, -9.3481814843660959, -3.9454681762489472, -296.18121334354328, 82.774301940161507, -190.84753679996979)
	cluster.m=np.random.rand(cluster.ntot)
	cluster.to_cluster()
	cluster.reset_nbody_scale(rvirial=False)

	init_units=['pckms','kpckms','pcmyr','kpcgyr','galpy','nbody','WDunits']
	mcons=[1.0,1.0,1.0,1.0,1./conversion.mass_in_msol(ro=solar_ro, vo=solar_vo),1./cluster.zmbar,1./222288.4543021174]
	pcons=[1.0,1.0,1.0,1.0,1./conversion.surfdens_in_msolpc2(ro=solar_ro, vo=solar_vo),(cluster.rbar**2.)/cluster.zmbar,1./222288.4543021174]
	rcons=[1.0,1.0,1.0,1.0,1./solar_ro,1.0/cluster.rbar,1.]
	vcons=[1.0,1.0,1.0,1.0,1./solar_vo,1.0/cluster.vbar,1.]

	#vcon=[1.,1.0,1.022712165045695,1.022712165045695,solar_vo,cluster.vbar,1.022712165045695]
	#rcon=[1.,1000.,1.,1000.,solar_ro,cluster.rbar,1000.0]

	projected=True

	for i in range(0,len(init_units)):
		print('DEBUG : ',i,init_units[i])

		if init_units[i]=='pckms':
			munit = u.MSun
			vunit = u.kms
			runit = u.pc
		elif init_units[i]=='kpckms':
			munit = u.MSun
			vunit = u.kms
			runit = u.kpc
		elif init_units[i]=='pcmyr':
			munit = u.MSun
			vunit = u.pc/u.Myr
			runit = u.pc
		elif init_units[i]=='kpcgyr':
			munit = u.MSun
			vunit = u.kpc/u.Gyr
			runit = u.kpc
		elif init_units[i]=='galpy':
			munit = u.MSun
			vunit = u.kms
			runit = u.kpc
		elif init_units[i]=='nbody':
			munit = u.MSun
			vunit = u.kms
			runit = u.pc
		elif init_units[i]=='WDunits':
			munit = u.MSun
			vunit = u.kpc/u.Gyr
			runit = u.kpc

		if init_units[i]=='galpy':
			punit=munit/((u.pc)**2.0)
		else:
			punit=munit/(runit**2.0)


		cluster.to_units(init_units[i])

		rprof0, prof0, nprof0=ctools.rho_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, nprof=ctools.rho_prof(cluster,projected=projected)

		print('DEBUG:',rprof)
		print(prof)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(punit)*pcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0, nprof0=ctools.m_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, nprof=ctools.m_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(munit)*mcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0, dalpha, edalpha, ydalpha, eydalpha=ctools.alpha_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, dalpha, edalpha, ydalpha, eydalpha=ctools.alpha_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.sigv_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.sigv_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit)*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.beta_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.beta_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.v_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.v_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit)*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0=ctools.v2_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof=ctools.v2_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit*vunit)*vcons[i]*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])

		rprof0, prof0, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,projected=projected)
		cluster.to_amuse()
		rprof, prof, deta, edeta, ydeta, eydeta=ctools.eta_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		cluster.to_units(init_units[i])


		rprof0, prof0, rvmax0, vmax0=ctools.vcirc_prof(cluster,projected=projected)
		cluster.to_amuse(analyze=True)
		rprof, prof, rvmax, vmax=ctools.vcirc_prof(cluster,projected=projected)

		rprof=rprof.value_in(runit)*rcons[i]
		prof=prof.value_in(vunit)*vcons[i]
		rvmax=rvmax.value_in(runit)*rcons[i]
		vmax=vmax.value_in(vunit)*vcons[i]

		assert np.all(np.fabs(rprof/rprof0-1) < tol)
		assert np.all(np.fabs(prof/prof0-1) < tol)
		assert np.fabs(1.0-rvmax0/rvmax) < tol
		assert np.fabs(1.0-vmax0/vmax) < tol

