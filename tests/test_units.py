import clustertools as ctools
from clustertools.util.units import _convert_length
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

def test_convert_length(tol=0.0001):

	init_units=['nbody','pckms','pcmyr','kpckms','kpcgyr','galpy','WDunits','radec']
	cluster_units=['nbody','pckms','pcmyr','kpckms','kpcgyr','galpy','WDunits','radec']

	rbar=10.0
	dist=10.0
	x0=1.0

	conversion=np.zeros((8,8))
	conversion[0]=np.array([1.,1./rbar,1./rbar,1000./rbar,1000./rbar,1000.0*ctools.solar_ro/rbar,1000./rbar,1.0/np.degrees(np.arctan2(rbar/1000.0,dist))])
	conversion[1]=np.array([rbar,1.,1.,1000.,1000.,1000.0*ctools.solar_ro,1000.,1.0/np.degrees(np.arctan2(1./1000.0,dist))])
	conversion[2]=np.array([rbar,1.,1.,1000.,1000.,1000.0*ctools.solar_ro,1000.,1.0/np.degrees(np.arctan2(1./1000.0,dist))])
	conversion[3]=np.array([rbar/1000.0,1./1000.0,1./1000.0,1.,1.,ctools.solar_ro,1.,1.0/np.degrees(np.arctan2(1.,dist))])
	conversion[4]=np.array([rbar/1000.0,1./1000.0,1./1000.0,1.,1.,ctools.solar_ro,1.,1.0/np.degrees(np.arctan2(1.,dist))])
	conversion[6]=np.array([rbar/1000.0,1./1000.0,1./1000.0,1.,1.,ctools.solar_ro,1.,1.0/np.degrees(np.arctan2(1.,dist))])
	conversion[5]=np.array([conversion[0,5]**-1.,conversion[1,5]**-1.,conversion[2,5]**-1.,conversion[3,5]**-1.,conversion[4,5]**-1.,1.0,conversion[6,5]**-1.,1.0/np.degrees(np.arctan2(ctools.solar_ro,dist))])
	conversion[7]=conversion[3]/(dist*np.tan(np.radians(x0)))
	conversion[7,-1]=1.

	for i,init_unit in enumerate(init_units):
		for j,cluster_unit in enumerate(cluster_units):
			cluster=ctools.StarCluster(units=cluster_unit,origin='cluster')
			cluster.add_orbit(0.,0.,dist,0.,0.,0.)
			cluster.rbar=rbar

			xnew=_convert_length(x0,init_unit,cluster)
			print(init_unit,cluster_unit,x0,xnew,x0/xnew,conversion[i,j])
			assert np.fabs(x0/xnew-conversion[i,j])<=tol


