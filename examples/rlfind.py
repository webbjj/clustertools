import numpy as np
import matplotlib.pyplot as plt
import nbodypy as npy
from galpy import potential
from galpy.potential import MWPotential2014
import os
r0=8.
v0=220.

#for i in range(0,16):
for i in range(11,12):

    wdir='./nfw%i/' % i
    pot=potential.NFWPotential(a=2.,normalize=0.35,ro=r0,vo=v0)


    if os.path.isfile(wdir+'norbit.npy'):


        wfile=open(wdir+'rldata.dat','w')
        print(i,wdir)

        for j in range(0,13):
            nsnap=j
            cluster=npy.load_cluster('nbodypy',units0='realkpc',origin0='galaxy',nsnap=nsnap,wdir=wdir,oname='norbit.npy')
            cluster.rlimiting(pot=pot,nrad=50)
            bcluster=npy.sub_cluster(cluster,rmax=cluster.rl)
            if bcluster.ntot==0:
                break
            else:
                print(j,bcluster.mtot,bcluster.rm)

                wfile.write('%i %f %f %f\n' % (j,bcluster.mtot,bcluster.rm,cluster.rl))

        wfile.close()

#for i in range(0,16):
for i in range(11,12):
    wdir='./mwpot%i/' % i
    pot=MWPotential2014

    if os.path.isfile(wdir+'norbit.npy'):


        wfile=open(wdir+'rldata.dat','w')
        print(i,wdir)

        for j in range(0,13):
            print(j)
            nsnap=j
            cluster=npy.load_cluster('nbodypy',units0='realkpc',origin0='galaxy',nsnap=nsnap,wdir=wdir,oname='norbit.npy')
            cluster.rlimiting(pot=pot,nrad=50)
            bcluster=npy.sub_cluster(cluster,rmax=cluster.rl)
            if bcluster.ntot==0:
                break
            else:
                print(j,bcluster.mtot,bcluster.rm)

                wfile.write('%i %f %f %f\n' % (j,bcluster.mtot,bcluster.rm,cluster.rl))

        wfile.close()
