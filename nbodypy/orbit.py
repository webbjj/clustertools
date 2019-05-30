#Functions and Operations that heavily use galpy and focus on the cluster's orbit

from galpy.orbit import Orbit, Orbits
from galpy.util import bovy_coords,bovy_conversion
from galpy import potential
from galpy.potential import LogarithmicHaloPotential,MWPotential2014,rtide
from galpy.actionAngle import actionAngleStaeckel
from galpy.actionAngle.actionAngleIsochroneApprox import actionAngleIsochroneApprox

import numpy as np

from ..util.recipes import rotate,interpolate,binmaker
from .operations import save_cluster,return_cluster
from .profiles import rho_prof
from ..util.plots import *

import astropy.coordinates as coord
import astropy.units as u

def initialize_orbit(cluster,from_centre=False,r0=8.,v0=220.):
 
    units0,origin0=save_cluster(cluster)
    cluster.to_galpy()

    if from_centre:
       x,y,z=cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc
       vx,vy,vz=cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc
    else:
        x,y,z=cluster.xgc,cluster.ygc,cluster.zgc
        vx,vy,vz=cluster.vxgc,cluster.vygc,cluster.vzgc

    R,phi,z=bovy_coords.rect_to_cyl(x,y,z)
    vR,vT,vz=bovy_coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)
    o=Orbit([R,vR,vT,z,vz,phi],ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])

    cluster.orbit=o

    return_cluster(cluster,units0,origin0)

    return o

def initialize_orbits(cluster,r0=8.,v0=220.):
    units0,origin0=save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_galpy()

    x,y,z=cluster.x,cluster.y,cluster.z
    vx,vy,vz=cluster.vx,cluster.vy,cluster.vz

    R,phi,z=bovy_coords.rect_to_cyl(x,y,z)
    vR,vT,vz=bovy_coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)

    vxvv=np.array([R,vR,vT,z,vz,phi])
    vxvv=np.rot90(vxvv)
    os=Orbits(vxvv,ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])

    return_cluster(cluster,units0,origin0)

    return os

def integrate_orbit(cluster,pot=MWPotential2014,tfinal=12.0,nt=1000,r0=8.,v0=220.,plot=False):
    cluster.orbit=initialize_orbit(cluster)
    ts=np.linspace(0,tfinal/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),nt)
    cluster.orbit.integrate(ts,pot)

    if plot:
        cluster.orbit.plot()

    return ts,cluster.orbit

def orbit_interpolate(cluster,dt,pot=MWPotential2014,from_centre=False,do_tails=False,rmin=None,rmax=None,emin=None,emax=None,r0=8.,v0=220.):
    cluster.tphys+=dt
    units0,origin0=save_cluster(cluster)
    cluster.to_galaxy()
    
    if do_tails:
        
        cluster.to_cluster()
        if from_centre:cluster.to_centre()
        
        if rmin==None: rmin=np.min(cluster.r)
        if rmax==None: rmax=np.max(cluster.r)
        rindx=(cluster.r>=rmin) * (cluster.r<=rmax)
        
        if len(cluster.etot)==cluster.ntot:
            if emin==None: emin=np.min(cluster.etot)
            if emax==None: emax=np.max(cluster.etot)
            eindx=(cluster.etot>=emin) * (cluster.etot<=emax)
        else:
            eindx=cluster.id>-1
        
        indx=rindx * eindx
        tindx=np.invert(indx)
        
        cluster.to_galaxy()
    
    else:
        indx=cluster.id>-1
    
    print('DO CLUSTER')
    
    cluster.orbit=initialize_orbit(cluster,from_centre)
    ts=np.linspace(0,dt/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),10)
    print('INTEGRATE ORBIT')

    cluster.orbit.integrate(ts,pot)


    cluster.to_realkpc()

    if from_centre:
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
    
    if from_centre:
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
        cluster.to_galaxy()
        cluster.to_galpy()

        x,y,z=cluster.x[tindx],cluster.y[tindx],cluster.z[tindx]
        vx,vy,vz=cluster.vx[tindx],cluster.vy[tindx],cluster.vz[tindx]

        R,phi,z=bovy_coords.rect_to_cyl(x,y,z)
        vR,vT,vz=bovy_coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)

        vxvv=np.array([R,vR,vT,z,vz,phi])
        vxvv=np.rot90(vxvv)
        otail=Orbits(vxvv,ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])

        cluster.to_realkpc()

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

    return_cluster(cluster,units0,origin0)


def orbital_path(cluster,dt=0.1,nt=100,pot=MWPotential2014,from_centre=False,skypath=False,r0=8.,v0=220.):
    o=initialize_orbit(cluster,from_centre=from_centre)

    ts=np.linspace(0,-1.*dt/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),nt)
    o.integrate(ts,pot)

    R,phi,z=bovy_coords.rect_to_cyl(o.x(ts[-1]),o.y(ts[-1]),o.z(ts[-1]))
    vR,vT,vz=bovy_coords.rect_to_cyl_vec(o.vx(ts[-1]),o.vy(ts[-1]),o.vz(ts[-1]),o.x(ts[-1]),o.y(ts[-1]),o.z(ts[-1]))
    o=Orbit([R/r0,vR/v0,vT/v0,z/r0,vz/v0,phi],ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])
    ts=np.linspace(0,2.*dt/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),2.*nt)
    o.integrate(ts,pot)

    if skypath:
        ra=np.array(o.ra(ts))
        dec=np.array(o.dec(ts))
        dist=np.array(o.dist(ts))
        pmra=np.array(o.pmra(ts))
        pmdec=np.array(o.pmdec(ts))
        vlos=np.array(o.vlos(ts))

        if cluster.units=='realpc':
            t=ts*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)
        elif cluster.units=='nbody':
            t=ts*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)/cluster.tstar
        elif cluster.units=='galpy':
            t=ts
        else:
            t=ts*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)

        return t,ra,dec,dist,pmra,pmdec,vlos,o
    else:
        x=np.array(o.x(ts))
        y=np.array(o.y(ts))
        z=np.array(o.z(ts))
        vx=np.array(o.vx(ts))
        vy=np.array(o.vy(ts))
        vz=np.array(o.vz(ts))

        if cluster.units=='realpc':
            x*=1000.
            y*=1000.
            z*=1000.
            t=ts*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)
        elif cluster.units=='nbody':
            x*=(1000./cluster.rbar)
            y*=(1000./cluster.rbar)
            z*=(1000./luster.rbar)
            vx/=cluster.vstar
            vy/=cluster.vstar
            vz/=cluster.vstar
            t=ts*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)/cluster.tstar

        elif cluster.units=='galpy':
            x/=r0
            y/=r0
            z/=r0
            vx/=v0 
            vy/=v0
            vz/=v0
            t=ts
        else:
            t=ts*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)

        return t,x,y,z,vx,vy,vz,o

def orbital_path_match(cluster,dt=0.1,nt=100,pot=MWPotential2014,from_centre=False,r0=8.,v0=220.):

    units0,origin0=save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_realkpc()

    t,x,y,z,vx,vy,vz,o=orbital_path(cluster,dt=dt,nt=nt,pot=pot,from_centre=False,r0=8.,v0=220.)

    ts=np.linspace(t[0],t[-1],10*nt)/bovy_conversion.time_in_Gyr(ro=r0,vo=v0)
    x=o.x(ts)
    y=o.y(ts)
    z=o.z(ts)
    vx=o.vx(ts)
    vy=o.vy(ts)
    vz=o.vz(ts)

    pindx=np.argmin(np.fabs(ts-dt))
    print('DEBUG: ',pindx,ts[pindx])

    dx=np.tile(np.array(o.x(ts)),cluster.ntot).reshape(cluster.ntot,len(ts))-np.repeat(cluster.x,len(ts)).reshape(cluster.ntot,len(ts))
    dy=np.tile(np.array(o.y(ts)),cluster.ntot).reshape(cluster.ntot,len(ts))-np.repeat(cluster.y,len(ts)).reshape(cluster.ntot,len(ts))
    dz=np.tile(np.array(o.z(ts)),cluster.ntot).reshape(cluster.ntot,len(ts))-np.repeat(cluster.z,len(ts)).reshape(cluster.ntot,len(ts))
    dr=np.sqrt(dx**2.+dy**2.+dz**2.)
    
    indx=np.argmin(dr,axis=1)
    dpath=np.min(dr,axis=1)
    tstar=ts[indx]*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)


    dprog=[]
    for i in range(0,cluster.ntot):
        if indx[i]==pindx:
            dprog.append(0.)
        elif indx[i] > pindx:

            xdiff=x-x[indx[i]]
            ydiff=y-y[indx[i]]
            zdiff=z-z[indx[i]]

            dx=xdiff[pindx:indx[i]]-xdiff[pindx+1:indx[i]+1]
            dy=ydiff[pindx:indx[i]]-ydiff[pindx+1:indx[i]+1]
            dz=zdiff[pindx:indx[i]]-zdiff[pindx+1:indx[i]+1]
            dprog.append(np.sqrt(np.cumsum(np.fabs(dx))[-1]**2.+np.cumsum(np.fabs(dy))[-1]**2.+np.cumsum(np.fabs(dz))[-1]**2.))            
        else:   

            xdiff=x-x[pindx]
            ydiff=y-y[pindx]
            zdiff=z-z[pindx]

            dx=xdiff[indx[i]:pindx]-xdiff[indx[i]+1:pindx+1]
            dy=ydiff[indx[i]:pindx]-ydiff[indx[i]+1:pindx+1]
            dz=zdiff[indx[i]:pindx]-zdiff[indx[i]+1:pindx+1]
            dprog.append(-1.*np.sqrt(np.cumsum(np.fabs(dx))[-1]**2.+np.cumsum(np.fabs(dy))[-1]**2.+np.cumsum(np.fabs(dz))[-1]**2.))

    #Assign negative to stars with position vectors in opposite direction as local angular momentum vector
    rgc=np.column_stack([o.x(ts[indx]),o.y(ts[indx]),o.z(ts[indx])])
    vgc=np.column_stack([o.vx(ts[indx]),o.vy(ts[indx]),o.vz(ts[indx])])
    lz=np.cross(rgc,vgc)

    rstar=np.column_stack([cluster.x-o.x(ts[indx]),cluster.y-o.y(ts[indx]),cluster.z-o.z(ts[indx])])

    ldot=np.sum(rstar*lz,axis=1)
    dpath[ldot<0]*=-1 

    return_cluster(cluster,units0,origin0)

    return np.array(tstar),np.array(dprog),np.array(dpath),o

def stream_path(cluster,dt=0.1,nt=100,pot=MWPotential2014,from_centre=False,r0=8.,v0=220.):

    units0,origin0=save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_realkpc()

    tstar,dprog,dpath,o=orbital_path_match(cluster=cluster,dt=dt,nt=nt,pot=pot,from_centre=from_centre,r0=r0,v0=v0)

    t_lower,t_mid,t_upper,t_hist=binmaker(tstar,nbin=50)
    xstream=np.array([])
    ystream=np.array([])
    zstream=np.array([])
    vxstream=np.array([])
    vystream=np.array([])
    vzstream=np.array([])

    for i in range(0,len(t_mid)):
        indx=(tstar>=t_lower[i]) * (tstar<=t_upper[i])
        xstream=np.append(xstream,np.mean(cluster.x[indx]))
        ystream=np.append(ystream,np.mean(cluster.y[indx]))
        zstream=np.append(zstream,np.mean(cluster.z[indx]))
        vxstream=np.append(vxstream,np.mean(cluster.vx[indx]))
        vystream=np.append(vystream,np.mean(cluster.vy[indx]))
        vzstream=np.append(vzstream,np.mean(cluster.vz[indx]))

    return_cluster(cluster,units0,origin0)

    return t_mid,xstream,ystream,zstream,vxstream,vystream,vzstream

def stream_path_match(cluster,dt=0.1,nt=100,pot=MWPotential2014,from_centre=False,r0=8.,v0=220.):

    units0,origin0=save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_realkpc()

    ts,x,y,z,vx,vy,vz=stream_path(cluster,dt=dt,nt=nt,pot=pot,from_centre=False,r0=8.,v0=220.)

    pindx=np.argmin(np.fabs(ts-dt))
    print('DEBUG: ',pindx,ts[pindx])

    dx=np.tile(x,cluster.ntot).reshape(cluster.ntot,len(ts))-np.repeat(cluster.x,len(ts)).reshape(cluster.ntot,len(ts))
    dy=np.tile(y,cluster.ntot).reshape(cluster.ntot,len(ts))-np.repeat(cluster.y,len(ts)).reshape(cluster.ntot,len(ts))
    dz=np.tile(z,cluster.ntot).reshape(cluster.ntot,len(ts))-np.repeat(cluster.z,len(ts)).reshape(cluster.ntot,len(ts))
    dr=np.sqrt(dx**2.+dy**2.+dz**2.)
    
    indx=np.argmin(dr,axis=1)
    dstream=np.min(dr,axis=1)
    tstar=ts[indx]*bovy_conversion.time_in_Gyr(ro=r0,vo=v0)


    dprog=[]
    for i in range(0,cluster.ntot):
        if indx[i]==pindx:
            dprog.append(0.)
        elif indx[i] > pindx:

            xdiff=x-x[indx[i]]
            ydiff=y-y[indx[i]]
            zdiff=z-z[indx[i]]

            dx=xdiff[pindx:indx[i]]-xdiff[pindx+1:indx[i]+1]
            dy=ydiff[pindx:indx[i]]-ydiff[pindx+1:indx[i]+1]
            dz=zdiff[pindx:indx[i]]-zdiff[pindx+1:indx[i]+1]
            dprog.append(np.sqrt(np.cumsum(np.fabs(dx))[-1]**2.+np.cumsum(np.fabs(dy))[-1]**2.+np.cumsum(np.fabs(dz))[-1]**2.))            
        else:   

            xdiff=x-x[pindx]
            ydiff=y-y[pindx]
            zdiff=z-z[pindx]

            dx=xdiff[indx[i]:pindx]-xdiff[indx[i]+1:pindx+1]
            dy=ydiff[indx[i]:pindx]-ydiff[indx[i]+1:pindx+1]
            dz=zdiff[indx[i]:pindx]-zdiff[indx[i]+1:pindx+1]
            dprog.append(-1.*np.sqrt(np.cumsum(np.fabs(dx))[-1]**2.+np.cumsum(np.fabs(dy))[-1]**2.+np.cumsum(np.fabs(dz))[-1]**2.))

    #Assign negative to stars with position vectors in opposite direction as local angular momentum vector
    rgc=np.column_stack([x[indx],y[indx],z[indx]])
    vgc=np.column_stack([vx[indx],vy[indx],vz[indx]])
    lz=np.cross(rgc,vgc)

    rstar=np.column_stack([cluster.x-x[indx],cluster.y-y[indx],cluster.z-z[indx]])

    ldot=np.sum(rstar*lz,axis=1)
    dstream[ldot<0]*=-1 

    return_cluster(cluster,units0,origin0)

    return np.array(tstar),np.array(dprog),np.array(dstream)

def rtidal(cluster,pot=MWPotential2014 ,rtiterate=0,rgc=None,r0=8.,v0=220.):
    
    units0,origin0=save_cluster(cluster)

    cluster.to_centre()
    cluster.to_galpy()

    if rgc!=None:
        R=rgc/r0
        z=0.0
    else:
        R=np.sqrt(cluster.xgc**2.0+cluster.ygc**2.0)
        z=cluster.zgc

    #Calculate rtide
    rt=rtide(pot,R,z,M=cluster.mtot)
    nit=0
    for i in range(0,rtiterate):
        msum=0.0
        
        indx=(cluster.r < rt)
        msum=np.sum(cluster.m[indx])

        rtnew=rtide(pot,R,z,M=msum)
        
        print(rt,rtnew,rtnew/rt,msum/cluster.mtot)

        if rtnew/rt>=0.9:
            break
        rt=rtnew
        nit+=1

    print('FINAL RT: ',rt*r0*1000.0, 'pc after',nit,' of ',rtiterate,' iterations')

    if units0=='realpc':
        rt*=1000.0*r0
    elif units0=='realkpc':
        rt*=r0
    elif units0=='nbody':
        rt*=(1000.0*r0/cluster.rbar)

    cluster.rt=rt

    return_cluster(cluster,units0,origin0)

    return rt

def rlimiting(cluster,pot=MWPotential2014 ,rgc=None,r0=8.,v0=220.,nrad=20,projected=False,obs_cut=False,plot=False,**kwargs):

    units0,origin0=save_cluster(cluster)

    cluster.to_centre()
    cluster.to_galpy()

    if rgc!=None:
        R=rgc/r0
        z=0.0
    else:
        R=np.sqrt(cluster.xgc**2.0+cluster.ygc**2.0)
        z=cluster.zgc

    #Calculate local density:
    rho_local=potential.evaluateDensities(pot,R,z,ro=r0,vo=v0)/bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)

    rprof,pprof,nprof=rho_prof(cluster,nrad=nrad,projected=projected,obs_cut=obs_cut)

    if pprof[-1] > rho_local:
        rl=rprof[-1]
    elif pprof[0] < rho_local:
        rl=0.0
    else:
        indx=np.argwhere(pprof < rho_local)[0][0]
        r1=(rprof[indx-1],pprof[indx-1])
        r2=(rprof[indx],pprof[indx])

        rl=interpolate(r1,r2,y=rho_local)

    print('FINAL RL: ',rl*r0*1000.0, 'pc')

    if units0=='realpc':
        rl*=(1000.0*r0)
    elif units0=='realkpc':
        rl*=r0
    elif units=='nbody':
        rl*=(1000.0*r0/cluster.rbar)

    cluster.rl=rl

    return_cluster(cluster,units0,origin0)

    if plot:
        print('LOCAL DENSITY = ',rho_local)    

        filename=kwargs.pop('filename',None)   
        overplot=kwargs.pop('overplot',False)        
     
        if cluster.units=='nbody':
            rprof*=(r0*1000.0/cluster.rbar)
            pprof*=(bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)*(cluster.rbar**3.)/cluster.zmbar)
            rho_local*=(bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)*(cluster.rbar**3.)/cluster.zmbar)
            xunits=' (NBODY)'
            yunits=' (NBODY)'
        elif cluster.units=='realpc':
            rprof*=(r0*1000.0)
            pprof*=bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)
            rho_local*=bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)
            xunits=' (pc)'
            if projected:
                yunits=' Msun/pc^2'
            else:
                yunits=' Msun/pc^3'
        elif cluster.units=='realkpc':
            rprof*=r0
            pprof*=bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)*(1000.0**3.)
            rho_local*=bovy_conversion.dens_in_msolpc3(ro=r0,vo=v0)*(1000.0**3.)

            xunits=' (kpc)'
            if projected:
                yunits=' Msun/kpc^2'
            else:
                yunits=' Msun/kpc^3'
        elif cluster.units=='galpy':
            xunits=' (GALPY)'
            yunits=' (GALPY)'


        else:
            xunits=''
            yunits=''

        x,y,n=rprof,pprof,nprof
        nlplot(x,y,xlabel='R'+xunits,ylabel='rho'+yunits,title='Time = %f' % cluster.tphys,log=True,overplot=overplot,filename=filename)
        nlplot(x,np.ones(len(x))*rho_local,'--',overplot=True)
        nlplot(np.ones(len(y))*rl,y,'--',overplot=True)

        if filename!=None:
            plt.savefig(filename)

    return rl

def get_cluster_orbit(gcname='list',names=False,r0=8.,v0=220.):

    data=np.loadtxt('/Users/webbjj/Codes/nbodypy/tables/orbits.dat',str)
    i_d=data[:,0].astype('int')

    name=data[:,1]
    name2=data[:,2]
    ra=data[:,3].astype('float64')
    dec=data[:,4].astype('float64')
    dist=data[:,5].astype('float64')
    vlos=data[:,6].astype('float64')
    evlos=data[:,7].astype('float64')
    pmra=data[:,8].astype('float64')
    pmdec=data[:,9].astype('float64')
    epmra=data[:,10].astype('float64')
    epmdec=data[:,11].astype('float64')
    corr=data[:,12].astype('float64')
    rscale=data[:,13].astype('float64')
    nstar=data[:,14].astype('float64')
    simbad_name=data[:,15]

    if 'list' in gcname:
        print('SELECT CLUSTER:')
        for i in range(0,len(name)):
            print('%s %s' % (name[i],name2[i]))

        indx=np.ones(len(name),bool)
        if names:
            return name[indx]
        else:
            return -1
    elif 'all' in gcname or 'ALL' in gcname:
        indx=(dist>0.) * (i_d<=151)
    elif isinstance(gcname,str):
        gcname=gcname.upper()
        indx=np.logical_or(name==gcname,name2==gcname)
    else:
        for i in range(0,len(gcname)):
            gcname[i]=gcname[i].upper()
        indx=np.logical_or(np.in1d(name,gcname),np.in1d(name2,gcname))

    if np.sum(indx)==0:
        return -1
    elif np.sum(indx)==1:
        print('GETTING ORBIT: ',name[indx])
        vxvv=[float(ra[indx]),float(dec[indx]),float(dist[indx]),float(pmra[indx]),float(pmdec[indx]),float(vlos[indx])]
        o=Orbit(vxvv,ro=r0,vo=v0,radec=c,solarmotion=[-11.1,24.,7.25])
    else:
        print('GETTING ORBIT: ',name[indx])
        vxvv = coord.SkyCoord(ra=ra[indx]*u.degree, dec=dec[indx]*u.degree, distance=dist[indx]*u.kpc, pm_ra_cosdec=pmra[indx]*u.mas/u.yr, pm_dec=pmdec[indx]*u.mas/u.yr, radial_velocity=vlos[indx]*u.km/u.s)
        o=Orbits(vxvv,ro=r0,vo=v0,solarmotion=[-11.1,24.,7.25])

    if names:
        return o,name[indx]
    else:
        return o

def calc_actions(cluster,pot=MWPotential2014,r0=8.,v0=220.,**kwargs):
    """
    NAME:

       calc_actions

    PURPOSE:

       Calculate action angle values for each star 

    INPUT:

       cluster - StarCluster instance

       pot - GALPY potential used to calculate actions

       r0,v0 - GALPY scaling parameters

    KWARGS:

       type - method for calculating actions (default: staeckel)

       delta - focus for staeckel method (default: 0.45 - optimal for MWPotential2014)

       c - if True, always use C for calculations

       Additional KWARGS can be included for other action angle calculation methods in galpy

    OUTPUT:

        JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz

    HISTORY:

       2019 - Written - Webb (UofT)

    """  


    os=initialize_orbits(cluster,r0,v0)
    atype=kwargs.pop('type','staeckel')
    delta=kwargs.pop('delta',0.45)
    c=kwargs.pop('c',True)

    JR=os.jr(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    Jphi=os.jp(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    Jz=os.jz(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    OR=os.Or(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    Ophi=os.Op(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    Oz=os.Oz(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    TR=os.Tr(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    Tphi=os.Tp(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)
    Tz=os.Tz(pot=pot,type=atype,delta=delta,c=c,ro=r0,vo=v0,**kwargs)

    return JR,Jphi,Jz,OR,Ophi,Oz,TR,Tphi,Tz

