import sys
import os, os.path
import csv
import numpy
import math
import nbodypy as npy
from galpy.util import bovy_conversion
from galpy.potential import LogarithmicHaloPotential,MWPotential2014,rtide

debug=False

#Correct for Center of Mass
cmcorr=True
#Calculate boundedness
boundedness=True

def virialize(cluster):

    pot=abs(numpy.sum(cluster.pot))
    kin=abs(numpy.sum(cluster.kin))

    q = 0.5
    qv = numpy.sqrt(q*pot/kin)

    print(pot,kin,q,qv)

    cluster.vx*=qv
    cluster.vy*=qv
    cluster.vz*=qv


def get_cluster(type=None,f82=None,f83=None):

    cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)

    return cluster

def get_cluster_orbit(orbit):

    data=orbit.readline().split()
    
    xgc=float(data[8])
    ygc=float(data[9])
    zgc=float(data[10])
    vxgc=float(data[11])
    vygc=float(data[12])
    vzgc=float(data[13])

    rgc=math.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
    return xgc,ygc,zgc,vxgc,vygc,vzgc,rgc

def tsep(type,snapfile82=None,snapfile83=None,tphys_max=0.0):

    #OPEN FILES AND READ FIRST TIMESTEPS
    f83=open(snapfile83,'r')
    f82=open(snapfile82,'r')
    cluster=get_cluster(type=type,f82=f82,f83=f83)

    oname=os.path.isfile('gc_orbit.dat')
    if oname:
        orbit=open('gc_orbit.dat','r')
    else:
        oname=os.path.isfile('orbit.dat')
        if oname:
            orbit=open('orbit.dat','r')

    #Begin reading in timesteps
    nstep=0
    snapsearch=True
    while cluster.ntot>0 and snapsearch:
        
        if not debug:
            print('READING: ',nstep, oname)

        if oname:
            xgc,ygc,zgc,vxgc,vygc,vzgc,rgc=get_cluster_orbit(orbit)
        else:
            xgc,ygc,zgc,vxgc,vygc,vzgc=numpy.zeros(6)
        cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

        if cmcorr:
            npy.xvshift(cluster,cluster.xc,cluster.yc,cluster.zc,0.0,0.0,0.0,cluster.origin)

        xshift,yshift,zshift,vxshift,vyshift,vzshift=cluster.find_center(0.0,0.0,0.0)
        print('FIND CENTER: ',xshift,yshift,zshift)
        npy.xvshift(cluster,-xshift,-yshift,-zshift,-vxshift,-vyshift,-vzshift)

        if cluster.tphys==0.0:
            m0=numpy.zeros(cluster.ntot)
            for i in range(0,cluster.ntot):
                mindx=int(cluster.id[i])-1
                m0[mindx]=cluster.m[i]*cluster.zmbar

        if cluster.tphys>=tphys_max:
            cluster.key_params() 
            if boundedness:
                pot,ektot=npy.get_energies(cluster,True)
            
            snapsearch=False
            csvfile= open('fort.10.realpc','w')
            writer= csv.writer(csvfile,
                           delimiter=' ')

            npy.nbody_to_realpc(cluster)

            if boundedness:
                indx=cluster.etot < 0.
            else:
                indx=cluster.id > 0.

            print('MEAN MASS = ',(numpy.sum(cluster.m[indx]))/len(cluster.m[indx]))
            print('MTOT = ',numpy.sum(cluster.m[indx]))
            print('ORBIT: ',xgc+xshift/1000.,ygc+yshift/1000.,zgc+zshift/1000.,vxgc+vxshift,vygc+vyshift,vzgc+vzshift)

            new_zmbar=numpy.sum(cluster.m[indx])
            new_rbar=cluster.rbar
            new_vstar=0.06557*numpy.sqrt(new_zmbar/new_rbar)

            #cluster.m/=numpy.sum(cluster.m[indx]) 
            #cluster.x/=new_rbar
            #cluster.y/=new_rbar
            #cluster.z/=new_rbar
            #cluster.vx/=new_vstar
            #cluster.vy/=new_vstar
            #cluster.vz/=new_vstar

            for j in range(cluster.ntot):
                if boundedness and cluster.etot[j]<0:
                    writer.writerow([cluster.m[j],cluster.x[j],cluster.y[j],cluster.z[j],
                             (cluster.vx[j]),(cluster.vy[j]),(cluster.vz[j])])
                elif not boundedness:
                    writer.writerow([cluster.m[j],cluster.x[j],cluster.y[j],cluster.z[j],
                                     (cluster.vx[j]),(cluster.vy[j]),(cluster.vz[j])])

            csvfile.close()

            csvfile= open('fort.12','w')
            writer= csv.writer(csvfile,
                           delimiter=' ')

            for j in range(cluster.ntot):
                if boundedness and cluster.etot[j]<0:
                    mindx=int(cluster.id[j])-1
                    writer.writerow([cluster.m[j],cluster.kw[j],m0[mindx],cluster.ep[j],cluster.ospin[j]])
                elif not boundedness:
                    mindx=int(cluster.id[j])-1
                    writer.writerow([cluster.m[j],cluster.kw[j],m0[mindx],cluster.ep[j],cluster.ospin[j]])

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
    
    tphys_max=float(sys.argv[1])
    tsep(type,snapfile82=sys.argv[2],snapfile83=sys.argv[3],tphys_max=tphys_max)

