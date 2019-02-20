#Change StarCluster units when necessary and recalulate key_params
#Real units count as Msun, km/s, and pc
from coordinates import *

#CAUTION: - Double Check that everything is in the same units when definining StarCluster and note that any operations performed on StarCluster will be in those units and any calculated variables will not have their units changed by the functions below.

import math

def nbody_to_realpc(cluster):
    if cluster.units=='nbody':
        cluster.m=cluster.m*cluster.zmbar
        cluster.x=cluster.x*cluster.rbar
        cluster.y=cluster.y*cluster.rbar
        cluster.z=cluster.z*cluster.rbar
        cluster.vx=cluster.vx*cluster.vstar
        cluster.vy=cluster.vy*cluster.vstar
        cluster.vz=cluster.vz*cluster.vstar
        
        cluster.xc=cluster.xc*cluster.rbar
        cluster.yc=cluster.yc*cluster.rbar
        cluster.zc=cluster.zc*cluster.rbar
        cluster.vxc=cluster.vxc*cluster.vstar
        cluster.vyc=cluster.vyc*cluster.vstar
        cluster.vzc=cluster.vzc*cluster.vstar
        
        cluster.xgc=cluster.xgc*cluster.rbar
        cluster.ygc=cluster.ygc*cluster.rbar
        cluster.zgc=cluster.zgc*cluster.rbar
        cluster.vxgc=cluster.vxgc*cluster.vstar
        cluster.vygc=cluster.vygc*cluster.vstar
        cluster.vzgc=cluster.vzgc*cluster.vstar
        

        cluster.units='realpc'

        if cluster.nb>0:
            yrs = (cluster.rbar*1296000./(2.0*np.pi))**1.5/np.sqrt(cluster.zmbar)
            days = 365.25*yrs
            pctoau=206265.

            cluster.pb*=days
            cluster.semi*=cluster.rbar*pctoau

        cluster.key_params()

    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units,' - NOTHING DONE')


def nbody_to_realkpc(cluster,subcluster=None):
    if cluster.units=='nbody':
        cluster.m=cluster.m*cluster.zmbar
        cluster.x=cluster.x*cluster.rbar/1000.0
        cluster.y=cluster.y*cluster.rbar/1000.0
        cluster.z=cluster.z*cluster.rbar/1000.0
        cluster.vx=cluster.vx*cluster.vstar
        cluster.vy=cluster.vy*cluster.vstar
        cluster.vz=cluster.vz*cluster.vstar
        
        cluster.xc=cluster.xc*cluster.rbar/1000.0
        cluster.yc=cluster.yc*cluster.rbar/1000.0
        cluster.zc=cluster.zc*cluster.rbar/1000.0
        cluster.vxc=cluster.vxc*cluster.vstar
        cluster.vyc=cluster.vyc*cluster.vstar
        cluster.vzc=cluster.vzc*cluster.vstar
        
        cluster.xgc=cluster.xgc*cluster.rbar/1000.0
        cluster.ygc=cluster.ygc*cluster.rbar/1000.0
        cluster.zgc=cluster.zgc*cluster.rbar/1000.0
        cluster.vxgc=cluster.vxgc*cluster.vstar
        cluster.vygc=cluster.vygc*cluster.vstar
        cluster.vzgc=cluster.vzgc*cluster.vstar

        cluster.units='realkpc'
        cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)


def realpc_to_nbody(cluster,subcluster=None):
    if cluster.units=='realpc':
        cluster.m=cluster.m/cluster.zmbar
        cluster.x=cluster.x/cluster.rbar
        cluster.y=cluster.y/cluster.rbar
        cluster.z=cluster.z/cluster.rbar
        cluster.vx=cluster.vx/cluster.vstar
        cluster.vy=cluster.vy/cluster.vstar
        cluster.vz=cluster.vz/cluster.vstar
        
        cluster.xc=cluster.xc/cluster.rbar
        cluster.yc=cluster.yc/cluster.rbar
        cluster.zc=cluster.zc/cluster.rbar
        cluster.vxc=cluster.vxc/cluster.vstar
        cluster.vyc=cluster.vyc/cluster.vstar
        cluster.vzc=cluster.vzc/cluster.vstar
        
        cluster.xgc=cluster.xgc/cluster.rbar
        cluster.ygc=cluster.ygc/cluster.rbar
        cluster.zgc=cluster.zgc/cluster.rbar
        cluster.vxgc=cluster.vxgc/cluster.vstar
        cluster.vygc=cluster.vygc/cluster.vstar
        cluster.vzgc=cluster.vzgc/cluster.vstar

        cluster.units='nbody'
        cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)


def nbody_to_galpy(cluster,r0,v0,subcluster=None):
    if cluster.units=='nbody':
        nbody_to_realkpc(cluster,subcluster)
        cluster.m=cluster.m
        cluster.x=cluster.x/r0
        cluster.y=cluster.y/r0
        cluster.z=cluster.z/r0
        cluster.vx=cluster.vx/v0
        cluster.vy=cluster.vy/v0
        cluster.vz=cluster.vz/v0

        cluster.xc=cluster.xc/r0
        cluster.yc=cluster.yc/r0
        cluster.zc=cluster.zc/r0
        cluster.vxc=cluster.vxc/v0
        cluster.vyc=cluster.vyc/v0
        cluster.vzc=cluster.vzc/v0
        
        cluster.xgc=cluster.xgc/r0
        cluster.ygc=cluster.ygc/r0
        cluster.zgc=cluster.zgc/r0
        cluster.vxgc=cluster.vxgc/v0
        cluster.vygc=cluster.vygc/v0
        cluster.vzgc=cluster.vzgc/v0
        

        cluster.units='galpy'
        cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)

def kpctopc(cluster):
    if cluster.units == 'realkpc':
        cluster.x*=1000.0
        cluster.y*=1000.0
        cluster.z*=1000.0
        
        cluster.xgc*=1000.0
        cluster.ygc*=1000.0
        cluster.zgc*=1000.0
        
        cluster.xc*=1000.0
        cluster.yc*=1000.0
        cluster.zc*=1000.0
      
        cluster.units='realpc'
        cluster.key_params()
    
    else:
        print('UNITS ALREADY IN PC - NOTHING DONE')

def pctokpc(cluster):
    if cluster.units=='realpc':
        cluster.x/=1000.0
        cluster.y/=1000.0
        cluster.z/=1000.0

        cluster.xgc/=1000.0
        cluster.ygc/=1000.0
        cluster.zgc/=1000.0

        cluster.xc/=1000.0
        cluster.yc/=1000.0
        cluster.zc/=1000.0

        cluster.units='realkpc'
        cluster.key_params()
    else:
        print('UNITS ALREADY IN KPC - NOTHING DONE')


