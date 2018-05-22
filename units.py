#Change StarCluster units when necessary and recalulate key_params
#Real units count as Msun, km/s, and pc
from coordinates import *

#CAUTION: - Double Check that everything is in the same units when definining StarCluster and note that any operations performed on StarCluster will be in those units and any calculated variables will not have their units changed by the functions below.

import math

def nbody_to_realpc(cluster):
    if cluster.units=='nbody':
        for i in range(0,len(cluster.x)):
            cluster.m[i]=cluster.m[i]*cluster.zmbar
            cluster.x[i]=cluster.x[i]*cluster.rbar
            cluster.y[i]=cluster.y[i]*cluster.rbar
            cluster.z[i]=cluster.z[i]*cluster.rbar
            cluster.vx[i]=cluster.vx[i]*cluster.vstar
            cluster.vy[i]=cluster.vy[i]*cluster.vstar
            cluster.vz[i]=cluster.vz[i]*cluster.vstar

        cluster.units='realpc'
        if cluster.keyparams:
            cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units,' - NOTHING DONE')


def nbody_to_realkpc(cluster,subcluster=None):
    if cluster.units=='nbody':
        for i in range(0,len(cluster.x)):
            cluster.m[i]=cluster.m[i]*cluster.zmbar
            cluster.x[i]=cluster.x[i]*cluster.rbar/1000.0
            cluster.y[i]=cluster.y[i]*cluster.rbar/1000.0
            cluster.z[i]=cluster.z[i]*cluster.rbar/1000.0
            cluster.vx[i]=cluster.vx[i]*cluster.vstar
            cluster.vy[i]=cluster.vy[i]*cluster.vstar
            cluster.vz[i]=cluster.vz[i]*cluster.vstar

        cluster.units='realkpc'
        if cluster.keyparams:
            cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)


def realpc_to_nbody(cluster,subcluster=None):
    if cluster.units=='realpc':
        for i in range(0,len(cluster.m)):
            cluster.m[i]=cluster.m[i]/cluster.zmbar
            cluster.x[i]=cluster.x[i]/cluster.rbar
            cluster.y[i]=cluster.y[i]/cluster.rbar
            cluster.z[i]=cluster.z[i]/cluster.rbar
            cluster.vx[i]=cluster.vx[i]/cluster.vstar
            cluster.vy[i]=cluster.vy[i]/cluster.vstar
            cluster.vz[i]=cluster.vz[i]/cluster.vstar

        cluster.units='nbody'
        if cluster.keyparams:
            cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)


def nbody_to_galpy(cluster,r0,v0,subcluster=None):
    if cluster.units=='nbody':
        nbody_to_realkpc(cluster,subcluster)
        for i in range(0,len(cluster.m)):
            cluster.m[i]=cluster.m[i]
            cluster.x[i]=cluster.x[i]/r0
            cluster.y[i]=cluster.y[i]/r0
            cluster.z[i]=cluster.z[i]/r0
            cluster.vx[i]=cluster.vx[i]/v0
            cluster.vy[i]=cluster.vy[i]/v0
            cluster.vz[i]=cluster.vz[i]/v0

        cluster.units='galpy'
        if cluster.keyparams:
            cluster.key_params()
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)

def kpctopc(cluster):
    if cluster.units == 'realkpc':
        if cluster.ntot==1:
            cluster.x*=1000.0
            cluster.y*=1000.0
            cluster.z*=1000.0
        else:
            for i in range(0,cluster.ntot):
                cluster.x[i]=cluster.x[i]*1000.0
                cluster.y[i]=cluster.y[i]*1000.0
                cluster.z[i]=cluster.z[i]*1000.0

        cluster.units='realpc'
        if cluster.keyparams:
            cluster.key_params()
    else:
        print('UNITS ALREADY IN PC - NOTHING DONE')

def pctokpc(cluster):
    if cluster.units=='realpc':
        if cluster.ntot==1:
            cluster.x/=1000.0
            cluster.y/=1000.0
            cluster.z/=1000.0
        else:
            for i in range(0,cluster.ntot):
                cluster.x[i]=cluster.x[i]/1000.0
                cluster.y[i]=cluster.y[i]/1000.0
                cluster.z[i]=cluster.z[i]/1000.0

        cluster.units='realkpc'
        if cluster.keyparams:
            cluster.key_params()
    else:
        print('UNITS ALREADY IN KPC - NOTHING DONE')


