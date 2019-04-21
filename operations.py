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

def add_rotation(cluster,qrot):
    r,theta,z=bovy_coords.rect_to_cyl(cluster.x,cluster.y,cluster.z)
    vr,vtheta,vz=bovy_coords.rect_to_cyl_vec(cluster.vx,cluster.vy,cluster.vz,cluster.x,cluster.y,cluster.z)
    
    indx=(vtheta < 0.)
    rindx=(np.random.rand(cluster.ntot) < qrot)

    vtheta[indx*rindx]=np.sqrt(vtheta[indx*rindx]*vtheta[indx*rindx])
    cluster.x,cluster.y,cluster.z=bovy_coords.cyl_to_rect(r,theta,z)
    cluster.vx,cluster.vy,cluster.vz=bovy_coords.cyl_to_rect_vec(vr,vtheta,vz,theta)

def nbody_scale(cluster,mass=True,radii=True):
    units0,origin0,center0=save_cluster(cluster)
    cluster.to_center()
    cluster.to_realpc()

    if mass:
        cluster.zmbar=cluster.mtot
    if radii:
        cluster.rbar=4.0*cluster.rm/3.0

    cluster.vstar=0.06557*np.sqrt(cluster.zmbar/cluster.rbar)
    cluster.tstar=cluster.rbar/cluster.vstar

    return_cluster(cluster,units0,origin0,center0)
    




