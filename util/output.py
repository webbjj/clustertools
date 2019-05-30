#Output files containing key cluster properties 
#Only functions and profiles should be called here

import numpy as np

from .coordinates import sky_coords
from ..nbodypy.cluster import sub_cluster
from ..nbodypy.functions import *
from ..nbodypy.profiles import *
from ..nbodypy.operations import *

def extrct_out(cluster,fileout,projected=False):
    #Write key properties (extrct.npy)
    
    units0=cluster.units
    origin0=cluster.origin

    if cluster.ntot==0:
        nb=0
        cluster.mtot=0.0
        trh=0.0
        rn=np.zeros(10)
    else:
        if origin0=='galaxy':
            cluster.to_cluster()
        if units0!='realpc':
            cluster.to_realpc()

        if cluster.nb>0:
            nb=len(cluster.m2)
        else:
            nb=0

        trh=relaxation_time(cluster,local=False,multimass=True,projected=projected)

        if cluster.ntot > 10:    
            if cluster.rn==None:
                rn=rlagrange(cluster,nlagrange=10,projected=projected)
        else:
            rn=np.zeros(10)
                
    fileout.write("%i %i %f %f %f " % (cluster.ntot,nb,cluster.tphys,trh,cluster.mtot))

    for i in range(0,len(rn)):
        fileout.write("%f " % rn[i])

    fileout.write("%f " % cluster.rmpro)

    if len(cluster.logl)>0:
        fileout.write("%f " % cluster.rhpro)

    fileout.write("\n")

    if units0=='realkpc':
        cluster.to_realkpc()
    elif units0=='nbody':
        cluster.to_nbody()
    elif units0=='galpy':
        cluster.to_galpy()
    if origin0=='galaxy':
        cluster.to_galaxy()

def snapout(cluster,filename,energies=False,observations=False):
    #Output a snapshot in nbodypy format

    if observations:

        ra,dec,d0,pmra,pmdec,vr0=sky_coords(cluster)

        if energies:
            np.savetxt(filename,np.olumn_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,cluster.kin,cluster.pot,cluster.etot,ra,dec,d0,pmra,pmdec,vr0]))
        else:
            np.savetxt(filename,np.column_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,ra,dec,d0,pmra,pmdec,vr0]))

    else:
        if energies:
            np.savetxt(filename,np.olumn_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw,cluster.kin,cluster.pot,cluster.etot]))
        else:
            np.savetxt(filename,np.column_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz,cluster.id,cluster.kw]))

    return 0

def fortout(cluster,filename='fort.10',reset_nbody_scale=False,reset_nbody_mass=True,reset_nbody_radii=True):
    #Output a snapshot in fort.10 format so it can be used to star a new Nbody6 simulation
    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    if reset_nbody_scale:
        reset_nbody_scale(cluster,mass=reset_nbody_mass,radii=reset_nbody_radii)

    cluster.to_nbody()

    np.savetxt(filename,np.column_stack([cluster.m,cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz]))

    return_cluster(cluster,units0,origin0)


    return 0
