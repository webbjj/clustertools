#Change Star or StarCluster units when necessary
#Real units count as Msun, km/s, and pc

#CAUTION: - Double Check that everything is in the same units when definining Star or StarCluster and note that any operations performed on Star or StarCluster will be in those units and any calculated variables will not have their units changed by the functions below.
#Also added subroutine to shift positions and velocities and rename origin

#In what follows, object is typically a StarCluster while subobject will be a Star. Based on the current definitions of the two, Star sometimes needs information from StarCluster to change units. Not currently loving this setup, so this may change. Need to tweak how individual stars are defined.

import math

def nbody_to_realpc(object,subobject=None):
    if object.units=='nbody':
        if subobject!= None:
            subobject.m*=object.zmbar
            subobject.x*=object.rbar
            subobject.y*=object.rbar
            subobject.z*=object.rbar
            subobject.r*=object.rbar
            subobject.rxy*=object.rbar
            subobject.vx*=object.vstar
            subobject.vy*=object.vstar
            subobject.vz*=object.vstar
            subobject.v*=object.vstar
        else:
            for i in range(0,len(object.x)):
                object.m[i]=object.m[i]*object.zmbar
                object.x[i]=object.x[i]*object.rbar
                object.y[i]=object.y[i]*object.rbar
                object.z[i]=object.z[i]*object.rbar
                object.r[i]=object.r[i]*object.rbar
                object.rxy[i]=object.rxy[i]*object.rbar
                object.vx[i]=object.vx[i]*object.vstar
                object.vy[i]=object.vy[i]*object.vstar
                object.vz[i]=object.vz[i]*object.vstar
                object.v[i]=object.v[i]*object.vstar
            object.rtide=object.rtide*object.rbar
            object.rc=object.rc*object.rbar
            object.xc=object.xc*object.rbar
            object.yc=object.yc*object.rbar
            object.zc=object.zc*object.rbar
        object.units='realpc'
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',object.units,' - NOTHING DONE')


def nbody_to_realkpc(object,subobject=None):
    if object.units=='nbody':
        if subobject!= None:
            subobject.m*=object.zmbar
            subobject.x*=object.rbar/1000.0
            subobject.y*=object.rbar/1000.0
            subobject.z*=object.rbar/1000.0
            subobject.r*=object.rbar/1000.0
            subobject.rxy*=object.rbar/1000.0
            subobject.vx*=object.vstar
            subobject.vy*=object.vstar
            subobject.vz*=object.vstar
            subobject.v*=object.vstar
        else:
            for i in range(0,len(object.x)):
                object.m[i]=object.m[i]*object.zmbar
                object.x[i]=object.x[i]*object.rbar/1000.0
                object.y[i]=object.y[i]*object.rbar/1000.0
                object.z[i]=object.z[i]*object.rbar/1000.0
                object.r[i]=object.r[i]*object.rbar/1000.0
                object.rxy[i]=object.rxy[i]*object.rbar/1000.0
                object.vx[i]=object.vx[i]*object.vstar
                object.vy[i]=object.vy[i]*object.vstar
                object.vz[i]=object.vz[i]*object.vstar
                object.v[i]=object.v[i]*object.vstar
            object.rtide=object.rtide*object.rbar/1000.0
            object.rc=object.rc*object.rbar/1000.0
            object.xc=object.xc*object.rbar/1000.0
            object.yc=object.yc*object.rbar/1000.0
            object.zc=object.zc*object.rbar/1000.0
            object.units='realkpc'
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',object.units)


def realpc_to_nbody(object,subobject=None):
    if object.units=='realpc':
        if subobject!= None:
            subobject.m/=object.zmbar
            subobject.x/=object.rbar
            subobject.y/=object.rbar
            subobject.z/=object.rbar
            subobject.r/=object.rbar
            subobject.rxy/=object.rbar
            subobject.vx/=object.vstar
            subobject.vy/=object.vstar
            subobject.vz/=object.vstar
            subobject.v/=object.vstar
        else:
            for i in range(0,len(object.m)):
            
                object.m[i]=object.m[i]/object.zmbar
                object.x[i]=object.x[i]/object.rbar
                object.y[i]=object.y[i]/object.rbar
                object.z[i]=object.z[i]/object.rbar
                object.r[i]=object.r[i]/object.rbar
                object.rxy[i]=object.rxy[i]/object.rbar
                object.vx[i]=object.vx[i]/object.vstar
                object.vy[i]=object.vy[i]/object.vstar
                object.vz[i]=object.vz[i]/object.vstar
                object.v[i]=object.v[i]/object.vstar
            object.rtide=object.rtide/rbar/1000.0
            object.rc=object.rc/rbar
            object.xc=object.xc/rbar
            object.yc=object.yc/rbar
            object.zc=object.zc/rbar
        object.units='nbody'
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',object.units)


def nbody_to_galpy(object,r0,v0,subobject=None):
    if object.units=='nbody':
        if subobject!=None:
            nbody_to_realkpc(object)
            subobject.m/=subobject.m
            subobject.x/=r0
            subobject.y/=r0
            subobject.z/=r0
            subobject.r/=r0
            subobject.rxy/=r0
            subobject.vx/=v0
            subobject.vy/=v0
            subobject.vz/=v0
            subobject.v/=v0
        
        else:
            nbody_to_realkpc(object,subobject)
            for i in range(0,len(object.m)):
                object.m[i]=object.m[i]
                object.x[i]=object.x[i]/r0
                object.y[i]=object.y[i]/r0
                object.z[i]=object.z[i]/r0
                object.r[i]=object.r[i]/r0
                object.rxy[i]=object.rxy[i]/r0
                object.vx[i]=object.vx[i]/v0
                object.vy[i]=object.vy[i]/v0
                object.vz[i]=object.vz[i]/v0
                object.v[i]=object.v[i]/v0
            object.rtide=object.rtide/r0
            object.rc=object.rc/r0
            object.xc=object.xc/r0
            object.yc=object.yc/r0
            object.zc=object.zc/r0

        object.units='galpy'
    else:
        print('CANNOT RUN OPERATION WHEN UNITS = ',object.units)

def xvscale(object,rscale=1.0,vscale=1.0,subobject=None):
    if subobject !=None:
        subobject.x*=rscale
        subobject.y*=rscale
        subobject.z*=rscale
        subobject.vx*=vscale
        subobject.vy*=vscale
        subobject.vz*=vscale
        subobject.v=math.sqrt(subobject.vx**2.0+subobject.vy**2.0+subobject.vz**2.0)
        subobject.rxy=math.sqrt(subobject.x**2.0+subobject.y**2.0)
        subobject.r=math.sqrt(subobject.x**2.0+subobject.y**2.0+subobject.z**2.0)
    else:
        for i in range(0,object.ntot):
            object.x[i]*=rscale
            object.y[i]*=rscale
            object.z[i]*=rscale
            object.vx[i]*=vscale
            object.vy[i]*=vscale
            object.vz[i]*=vscale
            object.v[i]=math.sqrt(object.vx[i]**2.0+object.vy[i]**2.0+object.vz[i]**2.0)
            object.rxy[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0)
            object.r[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0+object.z[i]**2.0)

def xvshift(object,x,y,z,vx,vy,vz,origin=None,subobject=None):
    if object.origin!=origin:
        if subobject !=None:
            if origin=='cluster':
                subobject.x-=x
                subobject.y-=y
                subobject.z-=z
                subobject.vx-=vx
                subobject.vy-=vy
                subobject.vz-=vz
                subobject.v=math.sqrt(subobject.vx**2.0+subobject.vy**2.0+subobject.vz**2.0)
                subobject.rxy=math.sqrt(subobject.x**2.0+subobject.y**2.0)
                subobject.r=math.sqrt(subobject.x**2.0+subobject.y**2.0+subobject.z**2.0)
                subobject.origin=origin
            elif origin=='galaxy':
                subobject.x+=x
                subobject.y+=y
                subobject.z+=z
                subobject.vx+=vx
                subobject.vy+=vy
                subobject.vz+=vz
                subobfject.v=math.sqrt(subobject.vx**2.0+subobject.vy**2.0+subobject.vz**2.0)
                subobject.rxy=math.sqrt(subobject.x**2.0+subobject.y**2.0)
                subobject.r=math.sqrt(subobject.x**2.0+subobject.y**2.0+subobject.z**2.0)
                subobject.origin=origin
            else:
                subobject.x+=x
                subobject.y+=y
                subobject.z+=z
                subobject.vx+=vx
                subobject.vy+=vy
                subobject.vz+=vz
                subobfject.v=math.sqrt(subobject.vx**2.0+subobject.vy**2.0+subobject.vz**2.0)
                subobject.rxy=math.sqrt(subobject.x**2.0+subobject.y**2.0)
                subobject.r=math.sqrt(subobject.x**2.0+subobject.y**2.0+subobject.z**2.0)
        else:
            for i in range(0,object.ntot):
                if origin=='cluster':
                    object.x[i]-=x
                    object.y[i]-=y
                    object.z[i]-=z
                    object.vx[i]-=vx
                    object.vy[i]-=vy
                    object.vz[i]-=vz
                    object.v[i]=math.sqrt(object.vx[i]**2.0+object.vy[i]**2.0+object.vz[i]**2.0)
                    object.rxy[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0)
                    object.r[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0+object.z[i]**2.0)
                    object.origin=origin
                elif origin=='galaxy':
                    object.x[i]+=x
                    object.y[i]+=y
                    object.z[i]+=z
                    object.vx[i]+=vx
                    object.vy[i]+=vy
                    object.vz[i]+=vz
                    object.v[i]=math.sqrt(object.vx[i]**2.0+object.vy[i]**2.0+object.vz[i]**2.0)
                    object.rxy[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0)
                    object.r[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0+object.z[i]**2.0)
                    object.origin=origin
                else:
                    object.x[i]+=x
                    object.y[i]+=y
                    object.z[i]+=z
                    object.vx[i]+=vx
                    object.vy[i]+=vy
                    object.vz[i]+=vz
                    object.v[i]=math.sqrt(object.vx[i]**2.0+object.vy[i]**2.0+object.vz[i]**2.0)
                    object.rxy[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0)
                    object.r[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0+object.z[i]**2.0)
    else:
        print(origin,' ALREADY EQUALS ',object.origin,' - NOTHING DONE')

#TO DO - CONVERT OTHER PARAMETERS AS WELL? RGC?
def kpctopc(object):
    if object.units == 'realkpc':
        if object.ntot==1:
            object.x*=1000.0
            object.y*=1000.0
            object.z*=1000.0
            object.r*=1000.0
            object.rxy*=1000.0
        else:
            for i in range(0,object.ntot):
                object.x[i]=object.x[i]*1000.0
                object.y[i]=object.y[i]*1000.0
                object.z[i]=object.z[i]*1000.0
                object.rxy[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0)
                object.r[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0+object.z[i]**2.0)
        object.units='realpc'
    else:
        print('UNITS ALREADY IN PC - NOTHING DONE')

def pctokpc(object):
    if object.units=='realpc':
        if object.ntot==1:
            object.x/=1000.0
            object.y/=1000.0
            object.z/=1000.0
            object.r/=1000.0
            object.rxy/=1000.0
        else:
            for i in range(0,object.ntot):
                object.x[i]=object.x[i]/1000.0
                object.y[i]=object.y[i]/1000.0
                object.z[i]=object.z[i]/1000.0
                object.rxy[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0)
                object.r[i]=math.sqrt(object.x[i]**2.0+object.y[i]**2.0+object.z[i]**2.0)
        object.units='realkpc'
    else:
        print('UNITS ALREADY IN KPC - NOTHING DONE')


