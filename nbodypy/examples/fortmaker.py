#MAIN VERSION 02/19/19
import sys
import os, os.path
import csv
import numpy
import math
import nbodypy as npy
from galpy.util import bovy_conversion
from galpy.potential import LogarithmicHaloPotential,MWPotential2014,rtide

debug=False

#Calculate boundedness
boundedness=True

def fmaker(ctype,tphys_max,**kwargs):

    #See if orbit file is provided
    if 'oname' in kwargs:
        ofile=open(kwargs.pop('oname'),'r')
    elif os.path.isfile('Norbit.dat'):
        ofile=open('Norbit.dat','r')
    elif os.path.isfile('norbit.npy'):
        ofile=open('norbit.npy','r')
    elif os.path.isfile('gc_orbit.dat'):
        ofile=open('gc_orbit.dat','r')
    elif os.path.isfile('orbit.dat'):
        ofile=open('orbit.dat','r')
    elif os.path.isfile('gorbit.dat'):
        ofile=open('gorbit.dat','r')
    else:
        ofile=None

    cluster=npy.load_cluster(ctype,ofile=ofile,**kwargs)
    
    while cluster.ntot > 0.:

        print(cluster.tphys,cluster.ctype)

        if cluster.tphys==0.0:
            m0=numpy.zeros(cluster.ntot)
            for i in range(0,cluster.ntot):
                mindx=int(cluster.id[i])-1
                m0[mindx]=cluster.m[i]*cluster.zmbar
        
        if cluster.tphys>=tphys_max:

            cluster.to_realpc()

            print('CLUSTER ORBIT: ',cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc)
            print('CLUSTER CENTER: ',cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc)
            cluster.to_centre()

            if 'rmax' in kwargs:
                rmax=float(kwargs.get('rmax'))
                print('GET RADII: ',rmax)
                new_cluster=npy.sub_cluster(cluster,rmax=rmax,new_centre=True)
            elif boundedness:
                print('GET ENERGIES: ',boundedness)
                npy.energies(cluster,specific=True)
                new_cluster=npy.sub_cluster(cluster,e_max=0.,new_centre=True)
            else:
                new_cluster=cluster

            print('NEW ORBIT: ',new_cluster.xgc,new_cluster.ygc,new_cluster.zgc,new_cluster.vxgc,new_cluster.vygc,new_cluster.vzgc)
            print('NEW CENTER: ',new_cluster.xc,new_cluster.yc,new_cluster.zc)
            print('MEAN MASS = ',(numpy.sum(new_cluster.m))/len(new_cluster.m))
            print('MTOT = ',numpy.sum(new_cluster.m))

            new_zmbar=numpy.sum(new_cluster.m)
            new_rbar=cluster.rbar
            new_vstar=0.06557*numpy.sqrt(new_zmbar/new_rbar)

            npy.snapout(new_cluster,'fort.10.realpc.t%i' % int(tphys_max))

            csvfile= open('fort.12.hrplot.t%i' % int(tphys_max),'w')
            writer= csv.writer(csvfile,delimiter=' ')
                
            for j in range(0,new_cluster.ntot):
                mindx=int(new_cluster.id[j])-1
                writer.writerow([new_cluster.m[j],new_cluster.kw[j],m0[mindx],new_cluster.ep[j],new_cluster.ospin[j]])
           
            csvfile.close()

            break

        cluster=npy.advance_cluster(cluster,ofile,**kwargs)

    return None

if __name__ == '__main__':
    
    ctype=sys.argv[1]
    tphys_max=float(sys.argv[2])
    
    if len(sys.argv) < 3:
        raise SyntaxError("Insufficient arguments.")
    elif len(sys.argv)==3:
        fmaker(ctype,tphys_max)
    else:
        # If there are keyword arguments
        fmaker(ctype,tphys_max,**dict(arg.split('=') for arg in sys.argv[3:]))

