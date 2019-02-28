from galpy.orbit import Orbit, Orbits
from galpy.util import bovy_coords,bovy_conversion
from galpy.potential import LogarithmicHaloPotential,MWPotential2014,rtide

import numpy as np
from cluster import StarCluster

def initialize_orbit(cluster,from_center=False,r0=8.,v0=220.):
 
    units0=cluster.units
    if units0!='galpy':
        cluster.to_galpy()

    if from_center:
       x,y,z=cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc
       vx,vy,vz=cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc
    else:
        x,y,z=cluster.xgc,cluster.ygc,cluster.zgc
        vx,vy,vz=cluster.vxgc,cluster.vygc,cluster.vzgc

    R,phi,z=bovy_coords.rect_to_cyl(x,y,z)
    vR,vT,vz=bovy_coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)
    o=Orbit([R,vR,vT,z,vz,phi],ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])

    if cluster.units!=units0:
        if units0=='realpc':
            cluster.to_realpc()
        elif units0=='realkpc':
            cluster.to_realkpc()
        elif units0=='nbody':
            cluster.to_nbody()

    return o

def initialize_orbits(cluster,r0=8.,v0=220.):
    units0=cluster.units
    origin0=cluster.origin
    center0=cluster.center
    
    if center0:
        cluster.from_center()
    
    if origin0!='galaxy':
        cluster.to_galaxy()
    
    if units0!='galpy':
        cluster.to_galpy()

    x,y,z=cluster.x,cluster.y,cluster.z
    vx,vy,vz=cluster.vx,cluster.vy,cluster.vz

    R,phi,z=bovy_coords.rect_to_cyl(x,y,z)
    vR,vT,vz=bovy_coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)

    vxvv=np.array([R,vR,vT,z,vz,phi])
    vxvv=np.rot90(vxvv)
    os=Orbits(vxvv,ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])

    if cluster.units!=units0:
        if units0=='realpc':
            cluster.to_realpc()
        elif units0=='realkpc':
            cluster.to_realkpc()
        elif units0=='nbody':
            cluster.to_nbody()

    if cluster.origin!=origin0:
        cluster.to_cluster()

    if center0:
        cluster.to_center()

    return os

def orbit_interpolate(cluster,dt,pot,from_center=False,do_tails=False,rmin=None,rmax=None,e_min=None,e_max=None,r0=8.,v0=220.):
    cluster.tphys+=dt
    
    if do_tails:
        
        cluster.to_cluster()
        if from_center:cluster.to_center()
        
        if rmin==None: rmin=np.min(cluster.r)
        if rmax==None: rmax=np.max(cluster.r)
        rindx=(cluster.r>=rmin) * (cluster.r<=rmax)
        
        if len(cluster.etot)==cluster.ntot:
            if e_min==None: e_min=np.min(cluster.etot)
            if e_max==None: e_max=np.max(cluster.etot)
            eindx=(cluster.etot>=e_min) * (cluster.etot<=e_max)
        else:
            eindx=cluster.id>-1
        
        indx=rindx * eindx
        tindx=np.invert(indx)
        
        if from_center:cluster.from_center()
        cluster.to_galaxy()
    
    else:
        indx=cluster.id>-1
    
    print('DO CLUSTER')
    
    cluster.orbit=initialize_orbit(cluster,from_center)
    ts=np.linspace(0,dt/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),10)
    print('INTEGRATE ORBIT')

    cluster.orbit.integrate(ts,pot)

    units0=cluster.units
    origin0=cluster.origin
    center0=cluster.center

    if units0!='realkpc':
        cluster.to_realkpc()
    if origin0!='galaxy':
        cluster.to_galaxy()

    if from_center:
        dx=cluster.orbit.x(ts[-1])-cluster.xc-cluster.xgc
        dy=cluster.orbit.y(ts[-1])-cluster.yc-cluster.ygc
        dz=cluster.orbit.z(ts[-1])-cluster.zc-cluster.zgc
        dvx=cluster.orbit.vx(ts[-1])-cluster.vxc-cluster.vxgc
        dvy=cluster.orbit.vy(ts[-1])-cluster.vyc-cluster.vygc
        dvz=cluster.orbit.vz(ts[-1])-cluster.vzc-cluster.vzgc
    else:
        dx=cluster.orbit.x(ts[-1])-cluster.xgc
        dy=cluster.orbit.y(ts[-1])-cluster.ygc
        dz=cluster.orbit.z(ts[-1])-cluster.zgc
        dvx=cluster.orbit.vx(ts[-1])-cluster.vxgc
        dvy=cluster.orbit.vy(ts[-1])-cluster.vygc
        dvz=cluster.orbit.vz(ts[-1])-cluster.vzgc
    
    print(dx,dy,dz,dvx,dvy,dvz)

    print('MOVING CLUSTER STARS')
    
    cluster.x[indx]+=dx
    cluster.y[indx]+=dy
    cluster.z[indx]+=dz
    cluster.vx[indx]+=dvx
    cluster.vy[indx]+=dvy
    cluster.vz[indx]+=dvz
    
    if from_center:
        cluster.xc,cluster.yc,cluster.zc=0.0,0.0,0.0
        cluster.vxc,cluster.vyc,cluster.vzc=0.0,0.0,0.0
    else:
        cluster.xc+=dx
        cluster.yc+=dy
        cluster.zc+=dz
        cluster.vxc+=dvx
        cluster.vyc+=dvy
        cluster.vzc+=dvz

    cluster.xgc,cluster.ygc,cluster.zgc=cluster.orbit.x(ts[-1]),cluster.orbit.y(ts[-1]),cluster.orbit.z(ts[-1])
    cluster.vxgc,cluster.vygc,cluster.vzgc=cluster.orbit.vx(ts[-1]),cluster.orbit.vy(ts[-1]),cluster.orbit.vz(ts[-1])

    if do_tails:
        print('DO TAILS')

        tail=StarCluster(np.sum(tindx),0.0,units=cluster.units,origin=cluster.origin,center=cluster.center)
        tail.add_stars(cluster.id[tindx],cluster.m[tindx],cluster.x[tindx],cluster.y[tindx],cluster.z[tindx],cluster.vx[tindx],cluster.vy[tindx],cluster.vz[tindx])

        otail=initialize_orbits(tail)
        ts=np.linspace(0,dt/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),10)
        print('INTEGRATE ORBITS')

        otail.integrate(ts,pot)
        
        print('MOVING TAIL STARS')

        cluster.x[tindx]=np.array(otail.x(ts[-1]))
        cluster.y[tindx]=np.array(otail.y(ts[-1]))
        cluster.z[tindx]=np.array(otail.z(ts[-1]))

        cluster.vx[tindx]=np.array(otail.vx(ts[-1]))
        cluster.vy[tindx]=np.array(otail.vy(ts[-1]))
        cluster.vz[tindx]=np.array(otail.vz(ts[-1]))

    if cluster.units!=units0:
        if units0=='realpc':
            cluster.to_realpc()
        elif units0=='realkpc':
            cluster.to_realkpc()
        elif units0=='nbody':
            cluster.to_nbody()

    if cluster.origin!=origin0:
        cluster.to_cluster()

def get_rtide(cluster,pot=MWPotential2014 ,rtiterate=0,rgc=None,r0=8.,v0=220.):
    
    units0=cluster.units
    origin0=cluster.origin
    center0=cluster.center

    if origin0!='cluster':
        cluster.to_cluster()
        cluster.to_center()
    elif not center0:
        cluster.to_center()

    if units0!='galpy':
        cluster.to_galpy()
    
    if rtide_rgc!=None:
        R=rgc/r0
        z=0.0
    elif cluster.units=='realkpc':
        R=numpy.sqrt(cluster.xgc**2.0+cluster.ygc**2.0)
        z=cluster.zgc

    #Calculate rtide
    rt=rtide(pot,R,z,M=cluster.mtot)

    for i in range(0,rtiterate):
        msum=0.0
        indx=(cluster.r < rt)
        msum=np.sum(cluster.m[indx])

        rtnew=rtide(pot,R/R0,z/R0,M=msum/bovy_conversion.mass_in_msol(ro=r0,vo=v0))*r0
        
        if cluster.units=='realpc':
            rtnew*=1000.0
        
        print(rt,rtnew,rtnew/rt,msum/cluster.mtot)

        if rtnew/rt>=0.9:
            break
        rt=rtnew

    if not center0:
        cluster.from_center()
    if units0!='cluster':
        cluster.to_galaxy()

    if units0=='realpc':
        rt*=1000.0*r0
        cluster.to_realpc()
    elif units0=='realkpc':
        rt*=r0
        cluster.to_realkpc()
    elif units=='nbody':
        rt*=(1000.0*r0/cluster.rbar)
        cluster.to_nbody()

    return rt
