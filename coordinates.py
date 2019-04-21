"""
COORDINATES

Designed to accept StarCluster instance as in put to
calculate stellar positions and velocities in different coordinates

"""

import numpy as np
from galpy.util import bovy_coords

def rect_to_sphere(cluster):
    """
    NAME:

       rect_to_sphere

    PURPOSE:

       Get the spherical coordinates of every stars in the cluster

    INPUT:

       cluster - a StarCluster-class object

    OUTPUT:

       r,theta,phi,vr,vtheta,vphi

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    
    r=np.sqrt(cluster.x**2.+cluster.y**2.+cluster.z**2.)
    theta=np.arccos(cluster.z/r)
    phi=np.arctan2(cluster.y,cluster.x)

    vr=cluster.vx*np.sin(theta)*np.cos(phi)+cluster.vy*np.sin(theta)*np.sin(phi)+cluster.vz*np.cos(theta)
    vtheta=cluster.vx*np.cos(theta)*np.cos(phi)+cluster.vy*np.cos(theta)*np.sin(phi)-cluster.vz*np.sin(theta)
    vphi=cluster.vx*-np.sin(phi)+cluster.vy*np.cos(phi)

    return r,theta,phi,vr,vtheta,vphi

def sky_coords(cluster):
    """
    NAME:

       sky_coords

    PURPOSE:

       Get the sky coordinates of every stars in the cluster

    INPUT:

       cluster - a StarCluster-class object

    OUTPUT:

       ra,dec,d0,pmra,pmdec,vr0

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    origin0=cluster.origin
    center0=cluster.center

    if origin0!='galaxy':
 	    cluster.to_galaxy()

    x0,y0,z0=bovy_coords.galcenrect_to_XYZ(cluster.x,cluster.y,cluster.z,Xsun=8.,Zsun=0.025).T
    vx0,vy0,vz0=bovy_coords.galcenrect_to_vxvyvz(cluster.vx,cluster.vy,cluster.vz,Xsun=8.,Zsun=0.025,vsun=[-11.1,244.,7.25]).T

    l0,b0,d0=bovy_coords.XYZ_to_lbd(x0,y0,z0,degree=True).T
    ra,dec=bovy_coords.lb_to_radec(l0,b0,degree=True).T

    vr0,pmll0,pmbb0=bovy_coords.vxvyvz_to_vrpmllpmbb(vx0,vy0,vz0,l0,b0,d0,degree=True).T
    pmra,pmdec=bovy_coords.pmllpmbb_to_pmrapmdec(pmll0,pmbb0,l0,b0,degree=True).T

    if center0:
        cluster.to_center()
    elif origin0=='cluster':
    	 cluster.to_cluster()

    return ra,dec,d0,pmra,pmdec,vr0

def cyl_coords(cluster):
    r,theta,z=bovy_coords.rect_to_cyl(cluster.x,cluster.y,cluster.z)
    vr,vtheta,vz=bovy_coords.rect_to_cyl_vec(cluster.vx,cluster.vy,cluster.vz,cluster.x,cluster.y,cluster.z)
    return r,theta,z,vr,vtheta,z
