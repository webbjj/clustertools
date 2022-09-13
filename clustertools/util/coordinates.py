"""For changing coordinate systems

Designed to accept StarCluster instance as in put to
calculate stellar positions and velocities in different coordinates

"""
__author__ = "Jeremy J Webb"


__all__ = [
    "sphere_coords",
    "cart_to_sphere",
    "sphere_to_cart",
    "cyl_coords",
    "cart_to_cyl",
    "cyl_to_cart",
    "sky_coords",
    "cart_to_sky"
]

import numpy as np
try:
    from galpy.util import coords
except:
    import galpy.util.bovy_coords as coords
from galpy.orbit import Orbit

import copy

from ..util.constants import *

def sphere_coords(cluster):
    """Get the spherical coordinates of every star in the cluster

    Parameters
    ----------
    cluster : class
      StarCluster

    Returns
    -------
    r,theta,phi,vr,vtheta,vphi : float
      stellar positions and velocities in spherical coordinates

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if origin0=='sky':
        cluster.to_galaxy(ro=ro,vo=vo,solarmotion=solar_motion)

    if cluster.units=='amuse':
        cluster.to_pckms()

    x,y,z,vx,vy,vz=cart_to_sphere(cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)


    return x,y,z,vx,vy,vz


def cart_to_sphere(x,y,z,vx,vy,vz):
    """Convert cartesian coordinates to spherical coordinates

    Parameters
    ----------
    x,y,z,vx,vy,vz : float
      positions and velocities in cartesian coordinates

    Returns
    -------
    r,theta,phi,vr,vtheta,vphi : float
      positions and velocities in spherical coordinates

    History
    -------
    2018 - Written - Webb (UofT)
    """

    r=np.sqrt(x**2.+y**2.+z**2.)
    phi=np.arctan2(y,x)
    theta=np.arccos(z/r)
    
    rhatx=np.cos(phi)*np.sin(theta)
    rhaty=np.sin(phi)*np.sin(theta)
    rhatz=np.cos(theta)
    
    phatx=-1.*np.sin(phi)
    phaty=np.cos(phi)
    phatz=0.

    thatx=np.cos(phi)*np.cos(theta)
    thaty=np.sin(phi)*np.cos(theta)
    thatz=-1.*np.sin(theta)
    
    vr=vx*rhatx+vy*rhaty+vz*rhatz
    vphi=vx*phatx+vy*phaty+vz*phatz
    vtheta=vx*thatx+vy*thaty+vz*thatz
    
    return r,phi,theta,vr,vphi,vtheta

def sphere_to_cart(r,phi,theta,vr,vphi,vtheta):
    """Convert cartesian coordinates to spherical coordinates

    Parameters
    ----------
    x,y,z,vx,vy,vz : float
      positions and velocities in cartesian coordinates

    Returns
    -------
    r,theta,phi,vr,vtheta,vphi : float
      positions and velocities in spherical coordinates

    History
    -------
    2018 - Written - Webb (UofT)
    """

    x=r*np.cos(phi)*np.sin(theta)
    y=r*np.sin(phi)*np.sin(theta)
    z=r*np.cos(theta)

    xhatr=np.cos(phi)*np.sin(theta)
    xhatp=-1.*np.sin(phi)
    xhatt=np.cos(phi)*np.cos(theta)

    yhatr=np.sin(phi)*np.sin(theta)
    yhatp=np.cos(phi)
    yhatt=np.sin(phi)*np.cos(theta)

    zhatr=np.cos(theta)
    zhatt=-1.*np.sin(theta)
    zhatp=0.

    vx=vr*xhatr+vtheta*xhatt+vphi*xhatp
    vy=vr*yhatr+vtheta*yhatt+vphi*yhatp
    vz=vr*zhatr+vtheta*zhatt+vphi*zhatp

    return x,y,z,vx,vy,vz

def cyl_coords(cluster):
    """Get the cylindrical coordinates of every star in the cluster

    Parameters
    ----------
    cluster : class
      StarCluster

    Returns
    -------
    r, theta, z, vr, vtheta, vz : float
      stellar positions and velocities in cylindrical coordinates

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if origin0=='sky':
        cluster.to_galaxy(ro=ro,vo=vo,solarmotion=solar_motion)

    if cluster.units=='amuse':
        cluster.to_pckms()
        
    R,phi,z,vR,vT,vz=cart_to_cyl(cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz)

    cluster.return_cluster(units0, origin0, rorder0, rorder_origin0)

    return R,phi,z,vR,vT,vz

def cart_to_cyl(x,y,z,vx,vy,vz):
    """Convert cartesian coordinates to cylindrical coordinates

    Parameters
    ----------
    x,y,z,vx,vy,vz : float
      positions and velocities in cartesian coordinates

    Returns
    -------
    r, theta, z, vr, vtheta, vz : float
      stellar positions and velocities in cylindrical coordinates

    History
    -------
    2018 - Written - Webb (UofT)
    """
    r, theta, zed = coords.rect_to_cyl(x, y, z)
    vr, vtheta, vzed = coords.rect_to_cyl_vec(
        vx, vy, vz, x, y, z
    )
    return r, theta, copy.copy(zed), vr, vtheta, copy.copy(vzed)

def cyl_to_cart(r,theta,zed,vr,vtheta,vz):
    """Convert cylindrical coordinates to cartesian coordinates 

    Parameters
    ----------
    r, theta, z, vr, vtheta, vz : float
      stellar positions and velocities in cylindrical coordinates

    Returns
    -------
    x,y,z,vx,vy,vz : float
      positions and velocities in cartesian coordinates

    History
    -------
    2022 - Written - Webb (UofT)
    """

    x,y,z=coords.cyl_to_rect(r,theta,zed)
    vx,vy,vz=coords.cyl_to_rect_vec(vr,vtheta,vz,theta)


    return x, y, copy.copy(z), vx, vy, copy.copy(vz)

def sky_coords(cluster,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion,zo=solar_zo):
    """Get the sky coordinates of every star in the cluster

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    ra,dec,d0,pmra,pmdec,vr0 : float
      on-sky positions and velocities of cluster stars
      
    History
    -------
    2018 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    cluster.to_radec(ro=ro,vo=vo,solarmotion=solar_motion)

    ra,dec,d0,pmra,pmdec,vr0=cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    return ra, dec, d0, pmra, pmdec, vr0

def cart_to_sky(x,y,z,vx,vy,vz,ro=solar_ro,vo=solar_vo,solarmotion=solar_motion):
    """Convert cartesian coordinates to sky coordinates

    Parameters
    ----------
    x,y,z,vx,vy,vz : float
      positions and velocities in cartesian coordinates

    Returns
    -------
    ra,dec,d0,pmra,pmdec,vr0 : float
      on-sky positions and velocities 
      
    History
    -------
    2021 - Written - Webb (UofT)
    """

    R, phi, z = coords.rect_to_cyl(x,y,z)
    vR, vT, vz = coords.rect_to_cyl_vec(vx,vy,vz,x,y,z)

    if isinstance(R,float):
        o = Orbit(
            [R/ro, vR/vo, vT/vo, z/ro, vz/vo, phi], ro=ro, vo=vo, solarmotion=solarmotion
        )
    else:

        vxvv=np.column_stack([R/ro, vR/vo, vT/vo, z/ro, vz/vo, phi])

        o = Orbit(
            vxvv, ro=ro, vo=vo, solarmotion=solarmotion
        )

    return o.ra(),o.dec(),o.dist(),o.pmra(),o.pmdec(),o.vlos()
