# Perform an operation on a cluster and return a new cluster

import numpy as np
from ..util.recipes import rotate
from galpy.util import bovy_conversion


def save_cluster(cluster):
    """
    NAME:

       save_cluster

    PURPOSE:

       Save cluster's units and origin

    INPUT:

       cluster - StarCluster instance

    OUTPUT:

       units, origin

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    return cluster.units, cluster.origin


def return_cluster(cluster, units0, origin0, do_order=False, do_key_params=False):
    """
    NAME:

       return_cluster

    PURPOSE:

       return cluster to a specific combination of units and origin

    INPUT:

       cluster - StarCluster instance

       units0 - units that StarCluster will be changed to

       origin0 - origin that StarCluster will be changed to


    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)
    """
    if cluster.units != units0:
        cluster.to_units(units0)
    if cluster.origin != origin0:
        cluster.to_origin(origin0, do_order=do_order, do_key_params=do_key_params)


def rotate_to_stream(cluster):
    """
    NAME:

       rotate_to_stream

    PURPOSE:

       Rotate the coordinate system to be aligned with cluster's orbit (Work in progress)

    INPUT:

       cluster - StarCluster instance

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    if cluster.origin != "cluster":
        cluster.to_cluster()
    v = np.array([cluster.vxgc, cluster.vygc, cluster.vzgc])
    thetax = np.arccos(np.dot([0.0, 0.0, 1.0], v) / np.linalg.norm(v))
    thetay = np.arccos(np.dot([0.0, 1.0, 0.0], v) / np.linalg.norm(v))
    thetaz = np.arccos(np.dot([1.0, 0.0, 0.0], v) / np.linalg.norm(v))

    x, y, z = rotate(cluster.x, cluster.y, cluster.z, thetax, thetay, thetaz)

    cluster.x = x
    cluster.y = y
    cluster.z = z


def add_rotation(cluster, qrot):
    """
    NAME:

       add_rotation

    PURPOSE:

       Add a degree of rotation to an already generated StarCluster

    INPUT:

       cluster - StarCluster instance

       qrot - the fraction of stars with v_phi < 0 that are switched to vphi > 0

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    r, theta, z = bovy_coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
    vr, vtheta, vz = bovy_coords.rect_to_cyl_vec(
        cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
    )

    indx = vtheta < 0.0
    rindx = np.random.rand(cluster.ntot) < qrot

    vtheta[indx * rindx] = np.sqrt(vtheta[indx * rindx] * vtheta[indx * rindx])
    cluster.x, cluster.y, cluster.z = bovy_coords.cyl_to_rect(r, theta, z)
    cluster.vx, cluster.vy, cluster.vz = bovy_coords.cyl_to_rect_vec(
        vr, vtheta, vz, theta
    )


def reset_nbody_scale(cluster, mass=True, radii=True, rvirial=False, **kwargs):
    """
    NAME:

       reset_nbody_scale

    PURPOSE:

       Assign new conversions for real mass, size, and velocity to Nbody units

    INPUT:

       cluster - StarCluster instance

       mass - find new mass conversion (Default: True)

       radii - find new radius conversion (Default: True)

       rvirial - use virial radius to set conversion rate for radii as opposed to the approximation that rbar=4/3 rm

    KWARGS:

       same as rvirial

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    units0, origin0 = save_cluster(cluster)
    cluster.to_centre()
    cluster.to_realpc()

    if mass:
        cluster.zmbar = cluster.mtot

    if radii:
        if rvirial:

            H = kwargs.get("H", 70.0)
            Om = kwargs.get("Om", 0.3)
            overdens = kwargs.get("overdens", 200.0)
            nrad = kwargs.get("nrad", 20.0)
            projected = kwargs.get("projected", False)

            cluster.rvirial(
                H=H, Om=Om, overdens=overdens, nrad=nrad, projected=projected
            )
        else:
            cluster.rbar = 4.0 * cluster.rm / 3.0

    cluster.vstar = 0.06557 * np.sqrt(cluster.zmbar / cluster.rbar)
    cluster.tstar = cluster.rbar / cluster.vstar

    return_cluster(cluster, units0, origin0)
