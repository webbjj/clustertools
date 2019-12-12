#MAIN VERSION 02/28/19
import sys
import os, os.path,glob
import csv
import numpy as np
import math
import nbodypy as npy
from galpy.util import bovy_conversion
from galpy.potential import LogarithmicHaloPotential,MWPotential2014,NFWPotential,rtide


#Maximum timestep to read to
tphys_max=12100.0
#Calculate energies
energies=False
#Track center of cluster via stars within tidal or limiting racius
track_centre=False
#Use 3D or projected values?
projected=False
#Calculate GALPY orbit
gorbit=False
r0=8.
v0=220.

#Extrct key parameters
extrct=True
rtextrct=False
dvprof=False
vprof=False
pprof=False
dalpha=False
aprof=False

# Coordinate for profiles:
coord='r'

#Mass Function:
mmin=0.3
mmax=0.8
rmin=None
rmax=None
kwmin=0
kwmax=1
obscut=None #WIP

#Rtide calculation
calc_rtide=False
calc_rlim=False
if calc_rlim:
    rprefix='rl'
else:
    rprefix='rt'

pot=MWPotential2014
#pot=NFWPotential(a=2.,normalize=0.35, ro=r0, vo=v0)
#pot=LogarithmicHaloPotential(normalize=1.,q=1.0)
rtiterate=0
#Specificy rgc or rtmax if want to keep rt fixed
rgc=None
rtmax=None
#Specify number of radial bins for calculating rlim
nrad=50

#Extrct snaps shops
snaps=False
snapdir= 'snaps/'
rtsnaps=False
rtsnapdir= 'rtsnaps/'
#Set origin (centre,cluster,galaxy)
set_origin='galaxy'
#Set units (nbody,realpc,realkpc,galpy)
set_units='realkpc'

def tsep(ctype,**kwargs):
 
    nstep=0

    print('CTYPE = ',ctype)
    print(kwargs)    
    wdir=kwargs.get('wdir','./')
    units=kwargs.pop('units','realpc')
    origin=kwargs.pop('origin','cluster')

    #Open extrct files for writing
    if extrct: exfile=open(wdir+'extrct.npy','w')
    if dvprof:dvfile=open(wdir+'dvprof.npy','w')
    if dalpha:dalpha_file=open(wdir+'dalpha_prof.npy','w')
    if aprof:aprof_file=open(wdir+'alpha_prof.npy','w')
    if pprof: pfile=open(wdir+'pprof.npy','w')
    if vprof: vfile=open(wdir+'vprof.npy','w')
    
    if rtextrct:
        rtexfile=open(wdir+rprefix+'extrct.npy','w')
        if dvprof:rtdvfile=open(wdir+rprefix+'dvprof.npy','w')
        if dalpha:rtdalpha_file=open(wdir+rprefix+'dalpha_prof.npy','w')
        if aprof:rtaprof_file=open(wdir+rprefix+'alpha_prof.npy','w')
        if pprof:rtpfile=open(wdir+rprefix+'pprof.npy','w')
        if vprof:rtdvfile=open(wdir+rprefix+'vprof.npy','w')


    list_file=kwargs.get('list',None)
    if list_file!=None:
        filenames=np.loadtxt(list_file,str)
        do_filenames=True
    else:
        do_filenames=False

    if snaps and not os.path.exists(wdir+snapdir):
        os.makedirs(wdir+snapdir)

    if rtsnaps and not os.path.exists(wdir+rtsnapdir):
        os.makedirs(wdir+rtsnapdir)

    #Output file for orbital and COM information
    norbit=open(wdir+'norbit.npy','w')

    #See if orbit file is provided
    if 'ofilename' in kwargs:
        ofile=open(kwargs.pop('ofilename'),'r')
    elif os.path.isfile('gc_orbit.dat'):
        ofile=open('gc_orbit.dat','r')
    elif os.path.isfile('orbit.dat'):
        ofile=open('orbit.dat','r')
    elif os.path.isfile('gorbit.dat'):
        ofile=open('gorbit.dat','r')
    else:  
        ofile=None

    #Initialize cluster
    if do_filenames:
        print(filenames[nstep])
        cluster=npy.load_cluster(ctype,units=units,origin=origin,ofile=ofile,filename=filenames[nstep],**kwargs)
    else:
        cluster=npy.load_cluster(ctype,units=units,origin=origin,ofile=ofile,**kwargs)

    units0,origin0=npy.save_cluster(cluster)

    if gorbit:
        go=npy.initialize_orbit(cluster)
        ts=np.linspace(0.,(tphys_max/1000.0)/bovy_conversion.time_in_Gyr(ro=r0,vo=v0),1200)
        go.integrate(ts,pot)

    bound_id=cluster.id

    #Begin reading in timesteps
    while cluster.ntot>0 and cluster.tphys<=tphys_max:
                
        print('READING: ',nstep, ctype, cluster.tphys, cluster.ntot, cluster.nsnap)

        if gorbit:
            to=(cluster.tphys/1000.)/bovy_conversion.time_in_Gyr(ro=r0,vo=v0)
            xgc,ygc,zgc=go.x(to),go.y(to),go.z(to)
            vxgc,vygc,vzgc=go.vx(to),go.vy(to),go.vz(to)
            cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,'realkpc')

        #Shift to clustercentric units for analysis
        cluster.to_cluster()
        if track_centre:
            bindx=np.in1d(cluster.id,bound_id)
            print('DEBUG BINDX: ',np.sum(bindx))
            cluster.find_centre(indx=bindx)
        else:
            bindx=np.ones(cluster.ntot,bool)

        print('CLUSTER ORBIT: ',cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc)
        print('MOVE TO CENTER: ',cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc)

        cluster.to_realpc()
        cluster.to_centre()

        if extrct:npy.extrct_out(cluster,exfile,projected=projected)
        if dvprof:npy.sigv_out(cluster,dvfile,projected=projected)
        if pprof:npy.p_prof_out(cluster,pfile,projected=projected)
        if vprof:npy.v_out(cluster,vfile,coord=coord,projected=projected)
        if dalpha: npy.dalpha_out(cluster,dalpha_file,projected=projected)
        if aprof: npy.alpha_prof_out(cluster,aprof_file,mmin=mmin,mmax=mmax,rmin=rmin,rmax=rmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obscut)

        if energies:
            ektot,pottot=npy.energies(cluster,specific=True)
            nbound=np.sum((cluster.etot<0.))

        if calc_rtide or calc_rlim:
            if energies:
                rtcluster=npy.sub_cluster(cluster,emax=0.0,indx=bindx)
            if calc_rtide:
                rtide=npy.rtidal(cluster,pot,rtiterate,rgc=rgc)
                rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rtide,indx=bindx)
            elif calc_rlim:
                rlim=npy.rlimiting(cluster,pot,rgc,nrad=nrad,projected=projected,obs_cut=obscut)
                rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rlim,indx=bindx)
            else:
                rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rtmax,indx=bindx)
           
            bound_id=rtcluster.id

        if rtextrct and rtcluster.ntot>0:

            npy.extrct_out(rtcluster,rtexfile,projected=projected)

            if dvprof: npy.sigv_out(rtcluster,rtdvfile,projected=projected)
            if pprof: npy.p_prof_out(rtcluster,rtpfile,projected=projected)
            if vprof: npy.v_out(rtcluster,rtvfile,coord=coord,projected=projected)
            if dalpha: npy.dalpha_out(rtcluster,rtdalpha_file,projected=projected)
            if aprof: npy.alpha_prof_out(rtcluster,rtaprof_file,mmin=mmin,mmax=mmax,rmin=rmin,rmax=rmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obscut)

        if snaps:
            if cluster.units!=set_units:
                if set_units=='realkpc':
                    cluster.to_realkpc()
                elif set_units=='realpc':
                    cluster.to_realpc()
                elif set_units=='nbody':
                    cluster.to_nbody()

            #Split snapshots and save to separate files
            if set_origin=='cluster':
                cluster.to_cluster()
            elif set_origin=='galaxy':
                cluster.to_galaxy()

            npy.snapout(cluster,os.path.join(snapdir,'%s.dat' % str(nstep).zfill(5)))

        if rtsnaps and rtcluster.ntot>0:
            if rtcluster.units!=set_units:
                if set_units=='realkpc':
                    rtcluster.to_realkpc()
                elif set_units=='realpc':
                    rtcluster.to_realpc()
                elif set_units=='nbody':
                    rtcluster.to_nbody()

            #Split snapshots and save to separate files
            if set_origin=='cluster':
                rtcluster.to_cluster()
            elif set_origin=='galaxy':
                rtcluster.to_galaxy()

            npy.snapout(rtcluster,os.path.join(rtsnapdir,'%s.dat' % str(nstep).zfill(5)))
    
        norbit.write("%f %f %f %f %f %f %f %f %f %f %f %f %f\n" % (cluster.tphys,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc))
        
        npy.return_cluster(cluster,units0,origin0)

        if do_filenames and nstep+1<len(filenames):
            print(filenames[nstep+1])
            cluster=npy.advance_cluster(cluster,ofile,filename=filenames[nstep+1],**kwargs)
        elif do_filenames:
            print('END OF LIST')
            break
        else:
            cluster=npy.advance_cluster(cluster,ofile,**kwargs)

        nstep+=1
    return None

if __name__ == '__main__':
    
    ctype=sys.argv[1]
    
    if len(sys.argv) < 2:
        raise SyntaxError("Insufficient arguments")
    elif len(sys.argv)==2:
        tsep(ctype)
    else:
        # If there are keyword arguments
        tsep(ctype,**dict(arg.split('=') for arg in sys.argv[2:]))
