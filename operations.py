#Perform an operation on a cluster and return a new cluster

import numpy as np
from recipes import rotate
from galpy.util import bovy_conversion

def save_cluster(cluster):
    return cluster.units,cluster.origin,cluster.center

def return_cluster(cluster,units0,origin0,center0):

    if not (center0 and cluster.center):
        if center0:
            cluster.to_center()
        else:
            cluster.from_center()

    if cluster.origin!=origin0:
        if origin0=='cluster':
            cluster.to_cluster()
        else:
            cluster.to_galaxy()

    if cluster.units!=units0:
        if units0=='realpc':
            cluster.to_realpc()
        elif units0=='realkpc':
            cluster.to_realkpc()
        elif units0=='nbody':
            cluster.to_nbody()
        elif units0=='galpy':
            cluster.to_galpy()

def rotate_to_stream(cluster):
    if cluster.origin!='cluster':
        cluster.to_cluster()
    v=np.array([cluster.vxgc,cluster.vygc,cluster.vzgc])
    thetax=np.arccos(np.dot([0.,0.,1.],v)/np.linalg.norm(v))
    thetay=np.arccos(np.dot([0.,1.,0.],v)/np.linalg.norm(v))
    thetaz=np.arccos(np.dot([1.,0.,0.],v)/np.linalg.norm(v))

    x,y,z=rotate(cluster.x,cluster.y,cluster.z,thetax,thetay,thetaz)

    cluster.x=x
    cluster.y=y
    cluster.z=z


