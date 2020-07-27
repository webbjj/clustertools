"""For changing coordinate systems

Designed to accept StarCluster instance as in put to
calculate stellar positions and velocities in different coordinates

"""
__author__ = "Jeremy J Webb"


__all__ = [
    "sphere_coords",
    "cart_to_sphere",
    "cyl_coords",
    "cart_to_cyl",
    "sky_coords"
]

import numpy as np
from galpy.util import bovy_coords


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
    return cart_to_sphere(cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz)


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
    theta=np.arccos(z/r)
    phi=np.arctan2(y,x)
    
    rhatx=x/r
    rhaty=y/r
    rhatz=z/r
    
    thatx=np.cos(theta)*np.cos(phi)
    thaty=np.cos(theta)*np.sin(phi)
    thatz=-1.*np.sin(theta)
    
    phatx=-1.*np.sin(phi)
    phaty=np.cos(phi)
    phatz=0.
    
    vr=vx*rhatx+vy*rhaty+vz*rhatz
    vtheta=vx*thatx+vy*thaty+vz*thatz
    vphi=vx*phatx+vy*phaty+vz*phatz
    
    return r,phi,theta,vr,vphi,vtheta

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
    return cart_to_cyl(cluster.x,cluster.y,cluster.z,cluster.vx,cluster.vy,cluster.vz)

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
    r, theta, z = bovy_coords.rect_to_cyl(x, y, z)
    vr, vtheta, vz = bovy_coords.rect_to_cyl_vec(
        vx, vy, vz, x, y, z
    )
    return r, theta, z, vr, vtheta, vz

def sky_coords(cluster):
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

    if origin0 != "galaxy":
        cluster.to_galaxy(starsort=False)

    x0, y0, z0 = bovy_coords.galcenrect_to_XYZ(
        cluster.x, cluster.y, cluster.z, Xsun=8.0, Zsun=0.025
    ).T
    vx0, vy0, vz0 = bovy_coords.galcenrect_to_vxvyvz(
        cluster.vx,
        cluster.vy,
        cluster.vz,
        Xsun=8.0,
        Zsun=0.025,
        vsun=[-11.1, 244.0, 7.25],
    ).T

    l0, b0, d0 = bovy_coords.XYZ_to_lbd(x0, y0, z0, degree=True).T
    ra, dec = bovy_coords.lb_to_radec(l0, b0, degree=True).T

    vr0, pmll0, pmbb0 = bovy_coords.vxvyvz_to_vrpmllpmbb(
        vx0, vy0, vz0, l0, b0, d0, degree=True
    ).T
    pmra, pmdec = bovy_coords.pmllpmbb_to_pmrapmdec(pmll0, pmbb0, l0, b0, degree=True).T

    cluster.return_cluster()

    return ra, dec, d0, pmra, pmdec, vr0
