import sys
import os, os.path
import csv
import numpy
import math
import nbodypy as npy
from galpy.util import bovy_conversion
from galpy.potential import MWPotential2014,rtide

debug=False

tphys_max=12000.0

#Correct for Center of Mass
cmcorr=False
#Use galpy orbit integrator
use_galpy_orbit=False
#Extrct key parameters
extrct=True
rtextrct=True
calc_rtide=False
rtiterate=1.
rtide_rgc=None
rtmax=201.188696061
#Extrct snaps shops
snaps=True
#Calculate boundedness
boundedness=False
#Set origin
set_origin='galaxy'
#Set units
set_units='realkpc'
#Swap y and z to compare with gyrfalcon
yzswap=False
#If there is mass correction (e.g. gyrfalcon):
mcon=222288.4543021174

def get_rtide(cluster,rtiterate,rtide_rgc=None):
    R0=8.
    V0=220.
    #pot=MWPotential2014
    pot=potential.LogarithmicHaloPotential(normalize=1.,q=1.0)
    
    if rtide_rgc!=None:
        R=rtide_rgc
        z=0.0
    else:
        R=numpy.sqrt(cluster.xgc**2.0+cluster.ygc**2.0)
        z=cluster.zgc
    
    rt=rtide(pot,R/R0,z/R0,M=cluster.mtot/bovy_conversion.mass_in_msol(ro=R0,vo=V0))*R0

    if cluster.units=='realpc':
        rt*=1000.0
    
    for i in range(0,rtiterate):
        msum=0.0
        for j in range(0,cluster.ntot):
            if cluster.r[j]<=rt:
                msum+=cluster.m[j]

        rtnew=rtide(pot,R/R0,z/R0,M=msum/bovy_conversion.mass_in_msol(ro=R0,vo=V0))*R0
        if cluster.units=='realpc':
            rtnew*=1000.0
        print(rt,rtnew,rtnew/rt,msum/cluster.mtot)
        if rtnew/rt>=0.9:
            break
        rt=rtnew

    return rt

def get_cluster(type=None,f82=None,f83=None):
    
    if type=='OUT34':
        cluster=npy.get_nbody6_out34(f83,do_keyparams=False)
    elif type=='NBODY6':
        cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)
    elif type=='GYRFALCON':
        vcon=220.0/bovy_conversion.velocity_in_kpcGyr(220.0,8.0)
        cluster=npy.get_gyrfalcon(f83,r0=8.0,v0=220.0,vcon=vcon,mcon=mcon,do_keyparams=False,find_center=True)
    elif type=='SNAPAUTO':
        cluster=npy.get_nbody6_snapauto(f83,do_keyparams=False)
    return cluster

def get_cluster_orbit(orbit):

    data=orbit.readline().split()
    xgc=float(data[1])
    ygc=float(data[2])
    zgc=float(data[3])
    vxgc=float(data[4])
    vygc=float(data[5])
    vzgc=float(data[6])
    rgc=math.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
    return xgc,ygc,zgc,vxgc,vygc,vzgc,rgc

def tsep(type,snapfile82=None,snapfile83=None):
    #Directories
    snapdir= 'snaps/'
    if type=='SNAPAUTO':
        nsnap=1
        snapfile83=('./snapshot/snap.dat.%s' % str(nsnap).zfill(4))
    
    #Open extrct files for writing
    if extrct:
        exfile=open('extrct.npy','w')
    if rtextrct:
        rtexfile=open('rtextrct.npy','w')
    
    if type=='OUT34':
        basefilename='fort'
    else:
        basefilename= snapfile83.split('.')[0]

    #OPEN FILES AND READ FIRST TIMESTEPS
    f83=open(snapfile83,'r')
    if snapfile82!=None:
        f82=open(snapfile82,'r')
        cluster=get_cluster(type=type,f82=f82,f83=f83)
    else:
        cluster=get_cluster(type=type,f82=None,f83=f83)

    #READ IN OR WRITE ORBIT DATA
    if use_galpy_orbit:
        gorbit=open('gorbit.dat','r')
    elif type=='OUT34':
        orbit=open('Norbit.dat','w')
    elif type=='GYRFALCON':
        oname=os.path.isfile('orbit.dat')
        if oname:
            orbit=open('orbit.dat','r')
        else:
            orbit=open('gorbit.dat','w')
    elif type=='NBODY6' or type=='SNAPAUTO':
        oname=os.path.isfile('orbit.dat')
        if oname:
            orbit=open('orbit.dat','r')
        else:
            orbit=open('orbit.dat','w')

    #Begin reading in timesteps
    nstep=0
    while cluster.ntot>0 and cluster.tphys<=tphys_max:
        
        if not debug:
            print('READING: ',nstep)

        if type=='OUT34':
            orbit.write("%f %f %f %f %f %f %f\n" % (cluster.tphys,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc))
        elif type=='GYRFALCON':
            if oname:
                odata=orbit.readline().split()
                cluster.add_orbit(float(odata[1])*8.,float(odata[2])*8.,float(odata[3])*8.,float(odata[4])*220.,float(odata[5])*220.,float(odata[6])*220.)
            else:
                orbit.write("%f %f %f %f %f %f %f\n" % (cluster.tphys,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc))

        elif type=='NBODY6' or type=='SNAPAUTO':
            if oname:
                xgc,ygc,zgc,vxgc,vygc,vzgc,rgc=get_cluster_orbit(orbit)
            else:
                xgc,ygc,zgc,vxgc,vygc,vzgc=numpy.zeros(6)
                cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

        if type!='GYRFALCON':

            if cmcorr:
                npy.xvshift(cluster,cluster.xc,cluster.yc,cluster.zc,0.0,0.0,0.0,cluster.origin)
            
            npy.nbody_to_realpc(cluster)

        if extrct:
            npy.extrct_out(cluster,exfile)

        if rtextrct:

            if cluster.origin=='galaxy':
                npy.xvshift(cluster,-cluster.xgc,-cluster.ygc,-cluster.zgc,-cluster.vxgc,-cluster.vygc,-cluster.vzgc,'cluster')
                npy.kpctopc(cluster)

                if calc_rtide:
                    rtide=get_rtide(cluster,rtiterate,rtide_rgc)
                    rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rtide)
                else:
                    rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rtmax)

                npy.pctokpc(cluster)
                npy.xvshift(cluster,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,'galaxy')
            else:
                if calc_rtide:
                    rtide=get_rtide(cluster,rtiterate,rtide_rgc)
                    rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rtide)
                else:
                    rtcluster=npy.sub_cluster(cluster,rmin=0.0,rmax=rtmax)
            npy.extrct_out(rtcluster,rtexfile)


        if debug:
            print('DEBUG ',numpy.min(cluster.x),cluster.tphys,cluster.ntot,rtmax)

        if boundedness:
            energy,bound=cluster.bound('nbody')
        else:
            energy=[0.0]*cluster.ntot
            bound=[True]*cluster.ntot
        
        if use_galpy_orbit:
            gdata=gorbit.readline().split()
            cluster.xgc=float(gdata[1])*8.
            cluster.ygc=float(gdata[2])*8.
            cluster.zgc=float(gdata[3])*8.
            cluster.vxgc=float(gdata[4])*220.
            cluster.vygc=float(gdata[5])*220.
            cluster.vzgc=float(gdata[6])*220.

        if snaps:
            csvfile= open(os.path.join(snapdir,basefilename+'_%s.dat' % str(nstep).zfill(5)),'w')
            writer= csv.writer(csvfile,
                           delimiter=',')

            if cluster.units!=set_units and set_units=='realkpc':
                if cluster.units=='nbody':
                    npy.nbody_to_realkpc(cluster)
                elif cluster.units=='realpc':
                    npy.pctokpc(cluster)

            if cluster.units!=set_units and set_units=='realpc':
                if cluster.units=='nbody':
                    npy.nbody_to_realpc(cluster)
                elif cluster.units=='realkpc':
                    npy.kpctopc(cluster)

            if cluster.units!=set_units and set_units=='nbody':
                if cluster.units=='realkpc':
                    npy.realkpc_to_nbody(cluster)
                elif cluster.units=='realpc':
                    npy.realpc_to_nbody(cluster)

            #Split snapshots and save to separate files
            if cluster.origin!=set_origin and set_origin=='galaxy':
                npy.xvshift(cluster,cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc,'galaxy')
            elif cluster.origin!=set_origin and set_origin=='cluster':
                npy.xvshift(cluster,-cluster.xgc,-cluster.ygc,-cluster.zgc,-cluster.vxgc,-cluster.vygc,-cluster.vzgc,'cluster')

            for j in range(cluster.ntot):
                if type=='GYRFALCON' and yzswap:
                    writer.writerow([cluster.m[j],cluster.x[j],cluster.y[j],cluster.z[j],
                                 (cluster.vx[j]),(cluster.vy[j]),(cluster.vz[j]),cluster.id[j]])
                elif yzswap:
                    writer.writerow([cluster.m[j],cluster.x[j],cluster.z[j],cluster.y[j],
                             (cluster.vx[j]),(cluster.vz[j]),(cluster.vy[j]),cluster.id[j]])
                else:
                    writer.writerow([cluster.m[j],cluster.x[j],cluster.y[j],cluster.z[j],
                             (cluster.vx[j]),(cluster.vy[j]),(cluster.vz[j]),cluster.id[j]])
            csvfile.close()

        if snapfile82!=None:
            cluster=get_cluster(type=type,f82=f82,f83=f83)
        elif type=='SNAPAUTO':
            nsnap+=1
            f83.close()
            snapfile83=('./snapshot/snap.dat.%s' % str(nsnap).zfill(4))
            f83=open(snapfile83,'r')
            cluster=get_cluster(type=type,f82=None,f83=f83)
        else:
            cluster=get_cluster(type=type,f82=None,f83=f83)

        if cluster.ntot==0:
            f83.close()
            snapfile83='./cont/'+snapfile83
            if os.path.isfile(snapfile83):
                f83=open(snapfile83,'r')
                if snapfile82!=None:
                    snapfile82='./cont/'+snapfile82
                    f82=open(snapfile83,'r')
                    cluster=get_cluster(type=type,f82=f82,f83=f83)
                else:
                    cluster=get_cluster(type=type,f82=None,f83=f83)

        nstep+=1
    return None

if __name__ == '__main__':
    
    print(len(sys.argv),sys.argv[1])
    
    if 'OUT34' in sys.argv[1]:
        type='OUT34'
    elif 'fort' in sys.argv[1] or 'bnd' in sys.argv[1] or 'esc' in sys.argv[1]:
        type='NBODY6'
    elif 'snapauto' in sys.argv[1]:
        type='SNAPAUTO'
    else:
        type='GYRFALCON'


    #If not in current directory, check for a /start directory
    if not os.path.isfile(sys.argv[1]):
        sys.argv[1]='./start/'+sys.argv[1]
    
    if len(sys.argv) ==2:
        tsep(type,snapfile82=None,snapfile83=sys.argv[1])
    elif len(sys.argv) ==3:
        if 'start' in sys.argv[1]:
            sys.argv[2]='./start/'+sys.argv[2]
        tsep(type,snapfile82=sys.argv[1],snapfile83=sys.argv[2])

