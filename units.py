#Change StarCluster units when necessary
#Real units count as Msun, km/s, and pc

#CAUTION: - Double Check that everything is in the same units when definining StarCluster and note that any operations performed on StarCluster will be in those units and any calculated variables will not have their units changed by the functions below.

#Also added subroutine to shift positions and velocities and rename origin

#TO DO - CONVERT OTHER PARAMETERS AS WELL? RGC?

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
        
        cluster.key_params()
        cluster.units='realpc'
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
        
        cluster.key_params()
        cluster.units='realkpc'
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
        
        cluster.key_params()
        cluster.units='nbody'
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
        
        cluster.key_params()
        cluster.units='galpy'
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',cluster.units)

def xvscale(cluster,rscale=1.0,vscale=1.0,subcluster=None):
    for i in range(0,cluster.ntot):
        cluster.x[i]*=rscale
        cluster.y[i]*=rscale
        cluster.z[i]*=rscale
        cluster.vx[i]*=vscale
        cluster.vy[i]*=vscale
        cluster.vz[i]*=vscale
    cluster.key_params()


def xvshift(cluster,x,y,z,vx,vy,vz,origin=None,subcluster=None):
    if cluster.origin!=origin:
        for i in range(0,cluster.ntot):
            if origin=='cluster':
                cluster.x[i]-=x
                cluster.y[i]-=y
                cluster.z[i]-=z
                cluster.vx[i]-=vx
                cluster.vy[i]-=vy
                cluster.vz[i]-=vz
                cluster.origin=origin
            elif origin=='galaxy':
                cluster.x[i]+=x
                cluster.y[i]+=y
                cluster.z[i]+=z
                cluster.vx[i]+=vx
                cluster.vy[i]+=vy
                cluster.vz[i]+=vz
                cluster.origin=origin
            else:
                cluster.x[i]+=x
                cluster.y[i]+=y
                cluster.z[i]+=z
                cluster.vx[i]+=vx
                cluster.vy[i]+=vy
                cluster.vz[i]+=vz
                
        cluster.key_params()

    else:
        print(origin,' ALREADY EQUALS ',cluster.origin,' - NOTHING DONE')

def kpctopc(cluster):
    if cluster.units == 'realkpc':
        if cluster.ntot==1:
            cluster.x*=1000.0
            cluster.y*=1000.0
            cluster.z*=1000.0
            cluster.r*=1000.0
            cluster.rxy*=1000.0
        else:
            for i in range(0,cluster.ntot):
                cluster.x[i]=cluster.x[i]*1000.0
                cluster.y[i]=cluster.y[i]*1000.0
                cluster.z[i]=cluster.z[i]*1000.0
                cluster.rxy[i]=math.sqrt(cluster.x[i]**2.0+cluster.y[i]**2.0)
                cluster.r[i]=math.sqrt(cluster.x[i]**2.0+cluster.y[i]**2.0+cluster.z[i]**2.0)
        
        cluster.key_params()
        cluster.units='realpc'
    else:
        print('UNITS ALREADY IN PC - NOTHING DONE')

def pctokpc(cluster):
    if cluster.units=='realpc':
        if cluster.ntot==1:
            cluster.x/=1000.0
            cluster.y/=1000.0
            cluster.z/=1000.0
            cluster.r/=1000.0
            cluster.rxy/=1000.0
        else:
            for i in range(0,cluster.ntot):
                cluster.x[i]=cluster.x[i]/1000.0
                cluster.y[i]=cluster.y[i]/1000.0
                cluster.z[i]=cluster.z[i]/1000.0
                cluster.rxy[i]=math.sqrt(cluster.x[i]**2.0+cluster.y[i]**2.0)
                cluster.r[i]=math.sqrt(cluster.x[i]**2.0+cluster.y[i]**2.0+cluster.z[i]**2.0)

        cluster.key_params()
        cluster.units='realkpc'
    else:
        print('UNITS ALREADY IN KPC - NOTHING DONE')


