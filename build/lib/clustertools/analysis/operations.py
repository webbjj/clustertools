""" Perform an operation on a cluster and return a new cluster

"""

__author__ = "Jeremy J Webb"


__all__ = [
    "find_centre",
    "find_centre_of_density",
    "find_centre_of_mass",
    "to_pckms",
    "to_kpckms",
    "to_nbody",
    "to_radec",
    "to_galpy",
    "to_units",
    "to_centre",
    "to_cluster",
    "to_galaxy",
    "to_sky",
    "to_tail",
    "to_origin",
    "save_cluster",
    "return_cluster",
    "rotate_to_tail",
    "reset_nbody_scale",
    "convert_binary_units",
    "add_rotation",
]

import numpy as np
from ..util.recipes import rotate
from galpy.util import bovy_conversion


def find_centre(
    self,
    xstart=0.0,
    ystart=0.0,
    zstart=0.0,
    vxstart=0.0,
    vystart=0.0,
    vzstart=0.0,
    indx=None,
    nsigma=1.0,
    nsphere=100,
    density=True,
    rmin=0.1,
    nmax=100,
    ro=8.0,
    vo=220.0,
):
    """
    NAME:

       find_centre

    PURPOSE:

       Find the centre of mass of the cluster
       Notes:
        - The routine first works to identify a sphere of nsphere stars around the centre in which
        to perform the C.O.M calculation. This step prevents long tidal tails from affecting the 
        calculation

    INPUT:
       xstart,ystart,zstart - starting position for centre
       vxstart,vystart,vzstart - starting velocity for centre
       indx - subset of stars to use when finding center
       nsigma - number of standard deviations to within which to keep stars
       nsphere - number of stars in centre sphere (default:100)
       density - use Yohai Meiron's centre of density calculator instead (Default: True)
       if density:
           - rmin - minimum radius to start looking for stars
           - nmax - maximum number of iterations to find centre
       ro,vo - For converting to and from galpy units (Default: 8., 220.)

    OUTPUT:

       xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    if indx is None:
        indx = np.ones(cluster.ntot, bool)
    elif np.sum(indx) == 0.0:
        print("NO SUBSET OF STARS GIVEN")
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    if density:
        xc,yc,zc,vxc,vyc,vzc=cluster.find_centre_of_density(
            xstart=xstart,
            ystart=ystart,
            zstart=zstart,
            vxstart=vxstart,
            vystart=vystart,
            vzstart=vzstart,
            indx=indx,
            rmin=rmin,
            nmax=nmax,
            ro=ro,
            vo=vo,
        )
    else:

        x = cluster.x[indx] - xstart
        y = cluster.y[indx] - ystart
        z = cluster.z[indx] - zstart
        r = np.sqrt(x ** 2.0 + y ** 2.0 + z ** 2.0)
        i_d = cluster.id[indx]

        while len(r) > nsphere:
            sigma = nsigma * np.std(r)
            indx = r < sigma

            if len(r[indx]) > nsphere:
                i_d = i_d[indx]
                x = x[indx] - np.mean(x[indx])
                y = y[indx] - np.mean(y[indx])
                z = z[indx] - np.mean(z[indx])
                r = np.sqrt(x * x + y * y + z * z)
            else:
                break

        # Find centre of mass and velocity of inner stars:
        indx = np.in1d(cluster.id, i_d)

        xc = np.sum(cluster.m[indx] * cluster.x[indx]) / np.sum(cluster.m[indx])
        yc = np.sum(cluster.m[indx] * cluster.y[indx]) / np.sum(cluster.m[indx])
        zc = np.sum(cluster.m[indx] * cluster.z[indx]) / np.sum(cluster.m[indx])

        vxc = np.sum(cluster.m[indx] * cluster.vx[indx]) / np.sum(cluster.m[indx])
        vyc = np.sum(cluster.m[indx] * cluster.vy[indx]) / np.sum(cluster.m[indx])
        vzc = np.sum(cluster.m[indx] * cluster.vz[indx]) / np.sum(cluster.m[indx])

        return xc, yc, zc, vxc, vyc, vzc

def find_centre_of_density(
    self,
    xstart=0.0,
    ystart=0.0,
    zstart=0.0,
    vxstart=0.0,
    vystart=0.0,
    vzstart=0.0,
    indx=None,
    rmin=0.1,
    nmax=100,
    ro=8.0,
    vo=220.0,
):
    """
    NAME:

       find_centre_of_density

    PURPOSE:

       Find the centre of density of the cluster
       Notes - the general framework for this code was taken from Yohai Meiron,
       who made a python adaptation of the centre of density finder in PhiGrape

    INPUT:
       xstart,ystart,zstart - starting position for centre
       vxstart,vystart,vzstart - starting velocity for centre
       rmin - minimum radius of sphere around which to estimate density centre (default: 0.1 pc)
       nmax - maximum number of iterations (default:100)

    OUTPUT:

       xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

    HISTORY:

       2019 - Written - Webb (UofT)
       - with kudos to Yohai Meiron (UofT)

    """

    # Need to change rmin for pc to cluster.units:
    if cluster.units != "nbody":
        rmin /= cluster.rbar
    elif cluster.units == "kpckms":
        rmin /= 1000.0
    elif cluster.units == "galpy":
        rmin /= 1000.0 * ro

    if indx is None:
        indx = np.ones(cluster.ntot, bool)

    m = cluster.m[indx]
    x = cluster.x[indx] - xstart
    y = cluster.y[indx] - ystart
    z = cluster.z[indx] - zstart
    vx = cluster.vx[indx] - vxstart
    vy = cluster.vy[indx] - vystart
    vz = cluster.vz[indx] - vzstart

    r = np.sqrt(x ** 2.0 + y ** 2.0 + z ** 2.0)
    rlim = np.amax(r)

    xdc, ydc, zdc = xstart, ystart, zstart
    vxdc, vydc, vzdc = vxstart, vystart, vzstart

    n = 0

    while (rlim > rmin) and (n < nmax):
        r2 = x ** 2.0 + y ** 2.0 + z ** 2.0
        indx = r2 < rlim ** 2
        nc = np.sum(indx)
        mc = np.sum(m[indx])

        if mc == 0:
            xc, yc, zc = 0.0, 0.0, 0.0
            vxc, vyc, vzc = 0.0, 0.0, 0.0
        else:

            xc = np.sum(m[indx] * x[indx]) / mc
            yc = np.sum(m[indx] * y[indx]) / mc
            zc = np.sum(m[indx] * z[indx]) / mc

            vxc = np.sum(m[indx] * vx[indx]) / mc
            vyc = np.sum(m[indx] * vy[indx]) / mc
            vzc = np.sum(m[indx] * vz[indx]) / mc

        if (mc > 0) and (nc > 100):
            x -= xc
            y -= yc
            z -= zc
            xdc += xc
            ydc += yc
            zdc += zc

            vx -= vxc
            vy -= vyc
            vz -= vzc
            vxdc += vxc
            vydc += vyc
            vzdc += vzc

        else:
            break
        rlim *= 0.8
        n += 1

    return xdc, ydc, zdc,vxdc, vydc, vzdc


def find_centre_of_mass(self):
    """
    NAME:

       find_centre_of_mass

    PURPOSE:

       Find the centre of mass of the cluster

    INPUT:
       None

    OUTPUT:

       xc,yc,zc,vxc,vyc,vzc - coordinates of centre of mass

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    xc = np.sum(cluster.m * cluster.x) / cluster.mtot
    yc = np.sum(cluster.m * cluster.y) / cluster.mtot
    zc = np.sum(cluster.m * cluster.z) / cluster.mtot

    vxc = np.sum(cluster.m * cluster.vx) / cluster.mtot
    vyc = np.sum(cluster.m * cluster.vy) / cluster.mtot
    vzc = np.sum(cluster.m * cluster.vz) / cluster.mtot

    return xdc, ydc, zdc,vxdc, vydc, vzdc

def to_pckms(cluster, do_key_params=False):
    """
    NAME:

       to_pckms

    PURPOSE:

       Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and km/s

    INPUT:

       None

    OUTPUT:

        None

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    if cluster.units == "galpy" or cluster.units == "radec":
        cluster.to_kpckms()

    if cluster.units == "nbody":
        cluster.m *= cluster.zmbar
        cluster.x *= cluster.rbar
        cluster.y *= cluster.rbar
        cluster.z *= cluster.rbar
        cluster.vx *= cluster.vstar
        cluster.vy *= cluster.vstar
        cluster.vz *= cluster.vstar

        cluster.xc *= cluster.rbar
        cluster.yc *= cluster.rbar
        cluster.zc *= cluster.rbar
        cluster.vxc *= cluster.vstar
        cluster.vyc *= cluster.vstar
        cluster.vzc *= cluster.vstar

        cluster.xgc *= cluster.rbar
        cluster.ygc *= cluster.rbar
        cluster.zgc *= cluster.rbar
        cluster.vxgc *= cluster.vstar
        cluster.vygc *= cluster.vstar
        cluster.vzgc *= cluster.vstar

        cluster.units = "pckms"

        #if cluster.nb > 0:
        #    yrs = (cluster.rbar * 1296000.0 / (2.0 * np.pi)) ** 1.5 / np.sqrt(
        #        cluster.zmbar
        #    )
        #    days = 365.25 * yrs
        #    pctoau = 206265.0

        #    cluster.pb *= days
        #    cluster.semi *= cluster.rbar * pctoau

    elif cluster.units == "kpckms":
        cluster.x *= 1000.0
        cluster.y *= 1000.0
        cluster.z *= 1000.0

        cluster.xgc *= 1000.0
        cluster.ygc *= 1000.0
        cluster.zgc *= 1000.0

        cluster.xc *= 1000.0
        cluster.yc *= 1000.0
        cluster.zc *= 1000.0

        cluster.units = "pckms"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_kpckms(self, do_key_params=False, ro=8.0, vo=220.0):
    """
    NAME:

       to_kpckms

    PURPOSE:

       Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and km/s

    INPUT:

       None

    OUTPUT:

        None

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    if cluster.units == "radec":
        cluster.from_radec()

    if cluster.units == "galpy":
        cluster.m *= bovy_conversion.mass_in_msol(ro=ro, vo=vo)
        cluster.x *= ro
        cluster.y *= ro
        cluster.z *= ro
        cluster.vx *= vo
        cluster.vy *= vo
        cluster.vz *= vo

        cluster.xc *= ro
        cluster.yc *= ro
        cluster.zc *= ro
        cluster.vxc *= vo
        cluster.vyc *= vo
        cluster.vzc *= vo

        cluster.xgc *= ro
        cluster.ygc *= ro
        cluster.zgc *= ro
        cluster.vxgc *= vo
        cluster.vygc *= vo
        cluster.vzgc *= vo

        cluster.units = "kpckms"

    elif cluster.units == "nbody":
        cluster.m *= cluster.zmbar
        cluster.x *= cluster.rbar / 1000.0
        cluster.y *= cluster.rbar / 1000.0
        cluster.z *= cluster.rbar / 1000.0
        cluster.vx *= cluster.vstar
        cluster.vy *= cluster.vstar
        cluster.vz *= cluster.vstar

        cluster.xc *= cluster.rbar / 1000.0
        cluster.yc *= cluster.rbar / 1000.0
        cluster.zc *= cluster.rbar / 1000.0
        cluster.vxc *= cluster.vstar
        cluster.vyc *= cluster.vstar
        cluster.vzc *= cluster.vstar

        cluster.xgc *= cluster.rbar / 1000.0
        cluster.ygc *= cluster.rbar / 1000.0
        cluster.zgc *= cluster.rbar / 1000.0
        cluster.vxgc *= cluster.vstar
        cluster.vygc *= cluster.vstar
        cluster.vzgc *= cluster.vstar

        cluster.units = "kpckms"

    elif cluster.units == "pckms":
        cluster.x /= 1000.0
        cluster.y /= 1000.0
        cluster.z /= 1000.0

        cluster.xgc /= 1000.0
        cluster.ygc /= 1000.0
        cluster.zgc /= 1000.0

        cluster.xc /= 1000.0
        cluster.yc /= 1000.0
        cluster.zc /= 1000.0

        cluster.units = "kpckms"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_nbody(self, do_key_params=False, ro=8.0, vo=220.0):
    """
    NAME:

       to_nbody

    PURPOSE:

       Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Nbody units
       Notes:
        - requires that cluster.zmbar, cluster.rbar, cluster.vstar are set (defaults are 1 in add_nbody6)

    INPUT:

       None

    OUTPUT:

        None

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    if cluster.units != "pckms":
        cluster.to_pckms(do_key_params=False)

    if cluster.units == "pckms":
        cluster.m /= cluster.zmbar
        cluster.x /= cluster.rbar
        cluster.y /= cluster.rbar
        cluster.z /= cluster.rbar
        cluster.vx /= cluster.vstar
        cluster.vy /= cluster.vstar
        cluster.vz /= cluster.vstar

        cluster.xc /= cluster.rbar
        cluster.yc /= cluster.rbar
        cluster.zc /= cluster.rbar
        cluster.vxc /= cluster.vstar
        cluster.vyc /= cluster.vstar
        cluster.vzc /= cluster.vstar

        cluster.xgc /= cluster.rbar
        cluster.ygc /= cluster.rbar
        cluster.zgc /= cluster.rbar
        cluster.vxgc /= cluster.vstar
        cluster.vygc /= cluster.vstar
        cluster.vzgc /= cluster.vstar

        cluster.units = "nbody"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_radec(self, do_key_params=False, ro=8.0, vo=220.0):
    """
    NAME:

       to_radec

    PURPOSE:

       Convert to on-sky position, proper motion, and radial velocity of cluster

    INPUT:

       None

    OUTPUT:

        ra,dec,dist,pmra,pmdec,vlos

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    try:
        cluster.x = copy(cluster.ra)
        cluster.y = copy(cluster.dec)
        cluster.z = copy(cluster.dist)
        cluster.vx = copy(cluster.pmra)
        cluster.vy = copy(cluster.pmdec)
        cluster.vz = copy(cluster.vlos)

        cluster.units = "radec"
        cluster.origin = "sky"

    except:

        units0, origin0 = cluster.units, cluster.origin

        cluster.to_galaxy()
        cluster.to_kpckms()

        x0, y0, z0 = bovy_coords.galcenrect_to_XYZ(
            cluster.x, cluster.y, cluster.z, Xsun=8.0, Zsun=0.025
        ).T

        cluster.dist = np.sqrt(x0 ** 2.0 + y0 ** 2.0 + z0 ** 2.0)

        vx0, vy0, vz0 = bovy_coords.galcenrect_to_vxvyvz(
            cluster.vx,
            cluster.vy,
            cluster.vz,
            Xsun=8.0,
            Zsun=0.025,
            vsun=[-11.1, 244.0, 7.25],
        ).T

        cluster.vlos = (vx0 * x0 + vy0 * y0 + vz0 * z0) / np.sqrt(
            x0 ** 2.0 + y0 ** 2.0 + z0 ** 2.0
        )

        l0, b0, cluster.dist = bovy_coords.XYZ_to_lbd(x0, y0, z0, degree=True).T
        cluster.ra, cluster.dec = bovy_coords.lb_to_radec(l0, b0, degree=True).T

        vr0, pmll0, pmbb0 = bovy_coords.vxvyvz_to_vrpmllpmbb(
            vx0, vy0, vz0, l0, b0, cluster.dist, degree=True
        ).T
        cluster.pmra, cluster.pmdec = bovy_coords.pmllpmbb_to_pmrapmdec(
            pmll0, pmbb0, l0, b0, degree=True
        ).T

        x0, y0, z0 = bovy_coords.galcenrect_to_XYZ(
            cluster.xgc, cluster.ygc, cluster.zgc, Xsun=8.0, Zsun=0.025
        )
        vx0, vy0, vz0 = bovy_coords.galcenrect_to_vxvyvz(
            cluster.vxgc,
            cluster.vygc,
            cluster.vzgc,
            Xsun=8.0,
            Zsun=0.025,
            vsun=[-11.1, 244.0, 7.25],
        )

        cluster.vlos_gc = (vx0 * x0 + vy0 * y0 + vz0 * z0) / np.sqrt(
            x0 ** 2.0 + y0 ** 2.0 + z0 ** 2.0
        )

        l0, b0, cluster.dist_gc = bovy_coords.XYZ_to_lbd(x0, y0, z0, degree=True)
        cluster.ra_gc, cluster.dec_gc = bovy_coords.lb_to_radec(l0, b0, degree=True)

        vr0, pmll0, pmbb0 = bovy_coords.vxvyvz_to_vrpmllpmbb(
            vx0, vy0, vz0, l0, b0, cluster.dist_gc, degree=True
        )
        cluster.pmra_gc, cluster.pmdec_gc = bovy_coords.pmllpmbb_to_pmrapmdec(
            pmll0, pmbb0, l0, b0, degree=True
        )

        cluster.x = copy(cluster.ra)
        cluster.y = copy(cluster.dec)
        cluster.z = copy(cluster.dist)
        cluster.vx = copy(cluster.pmra)
        cluster.vy = copy(cluster.pmdec)
        cluster.vz = copy(cluster.vlos)

        cluster.xgc = copy(cluster.ra_gc)
        cluster.ygc = copy(cluster.dec_gc)
        cluster.zgc = copy(cluster.dist_gc)
        cluster.vxgc = copy(cluster.pmra_gc)
        cluster.vygc = copy(cluster.pmdec_gc)
        cluster.vzgc = copy(cluster.vlos_gc)

        cluster.units = "radec"
        cluster.origin = "sky"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params(do_order=do_order)

def from_radec(self, do_order=False, do_key_params=False):
    """
    NAME:

       from_radec

    PURPOSE:

       Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster

    INPUT:

       None

    OUTPUT:

        x,y,z,vx,vy,vz

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    if cluster.units == "radec" and cluster.origin == "sky":

        origin0 = cluster.origin

        l, b = bovy_coords.radec_to_lb(cluster.ra, cluster.dec, degree=True).T
        x0, y0, z0 = bovy_coords.lbd_to_XYZ(l, b, cluster.dist, degree=True).T
        cluster.x, cluster.y, cluster.z = bovy_coords.XYZ_to_galcenrect(
            x0, y0, z0, Xsun=8.0, Zsun=0.025
        ).T

        pml, pmb = bovy_coords.pmrapmdec_to_pmllpmbb(
            cluster.pmra, cluster.pmdec, cluster.ra, cluster.dec, degree=True
        ).T
        vx0, vy0, vz0 = bovy_coords.vrpmllpmbb_to_vxvyvz(
            cluster.vlos, pml, pmb, l, b, cluster.dist, degree=True
        ).T
        cluster.vx, cluster.vy, cluster.vz = bovy_coords.vxvyvz_to_galcenrect(
            vx0,
            vy0,
            vz0,
            vsun=[0.0, 220.0, 0.0],
            Xsun=8.0,
            Zsun=0.025,
            _extra_rot=True,
        ).T

        l_gc, b_gc = bovy_coords.radec_to_lb(cluster.ra_gc, cluster.dec_gc, degree=True)
        x0_gc, y0_gc, z0_gc = bovy_coords.lbd_to_XYZ(
            l_gc, b_gc, cluster.dist_gc, degree=True
        )
        cluster.xgc, cluster.ygc, cluster.zgc = bovy_coords.XYZ_to_galcenrect(
            x0_gc, y0_gc, z0_gc, Xsun=8.0, Zsun=0.025
        )

        pml_gc, pmb_gc = bovy_coords.pmrapmdec_to_pmllpmbb(
            cluster.pmra_gc, cluster.pmdec_gc, cluster.ra_gc, cluster.dec_gc, degree=True
        )
        vx0_gc, vy0_gc, vz0_gc = bovy_coords.vrpmllpmbb_to_vxvyvz(
            cluster.vlos_gc, pml_gc, pmb_gc, l_gc, b_gc, cluster.dist_gc, degree=True
        )
        cluster.vx_gc, cluster.vy_gc, cluster.vz_gc = bovy_coords.vxvyvz_to_galcenrect(
            vx0_gc,
            vy0_gc,
            vz0_gc,
            vsun=[0.0, 220.0, 0.0],
            Xsun=8.0,
            Zsun=0.025,
            _extra_rot=True,
        )

        cluster.origin = "galaxy"
        cluster.units = "kpckms"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params(do_order=do_order)

def to_galpy(self, do_key_params=False, ro=8.0, vo=220.0):
    """
    NAME:

       to_galpy

    PURPOSE:

       Convert stellar positions/velocities, centre of mass, and orbital position and velocity to galpy units

    INPUT:

       None

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    if cluster.units != "kpckms" and cluster.units != "galpy":
        cluster.to_kpckms(do_key_params=False)

    if cluster.units == "kpckms":
        cluster.m = cluster.m / bovy_conversion.mass_in_msol(ro=ro, vo=vo)
        cluster.x /= ro
        cluster.y /= ro
        cluster.z /= ro
        cluster.vx /= vo
        cluster.vy /= vo
        cluster.vz /= vo

        cluster.xc /= ro
        cluster.yc /= ro
        cluster.zc /= ro
        cluster.vxc /= vo
        cluster.vyc /= vo
        cluster.vzc /= vo

        cluster.xgc /= ro
        cluster.ygc /= ro
        cluster.zgc /= ro
        cluster.vxgc /= vo
        cluster.vygc /= vo
        cluster.vzgc /= vo

        cluster.units = "galpy"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params()

def to_units(self, units, do_order=False, do_key_params=False, ro=8.0, vo=220.0):
    """
    NAME:

       to_units

    PURPOSE:

       Convert stellar positions/velocities, centre of mass, and orbital position and velocity to user defined units

    INPUT:

       units - 'nbody','pckms','kpckms','galpy'

    OUTPUT:

        None

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    if units == "nbody":
        cluster.to_nbody(do_key_params=do_key_params)
    elif units == "galpy":
        cluster.to_galpy(do_key_params=do_key_params, ro=ro, vo=vo)
    elif units == "pckms":
        cluster.to_pckms(do_key_params=do_key_params)
    elif units == "kpckms":
        cluster.to_kpckms(do_key_params=do_key_params)
    elif units == "radec":
        origin0 = cluster.origin
        cluster.to_radec(do_key_params=do_key_params)
        cluster.to_origin(origin0, do_order=do_order, do_key_params=do_key_params)

def rotate_to_tail(cluster):
    """
    NAME:

       rotate_to_tail

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

def to_centre(self, do_order=False, do_key_params=False, centre_method=None):
    """
    NAME:

       to_centre

    PURPOSE:

       Shift coordinates such that origin is the centre of mass

    INPUT:

       None

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    if cluster.origin != "centre":

        if cluster.origin != "cluster":
            cluster.to_cluster(do_key_params=False, centre_method=centre_method)

        cluster.x -= cluster.xc
        cluster.y -= cluster.yc
        cluster.z -= cluster.zc
        cluster.vx -= cluster.vxc
        cluster.vy -= cluster.vyc
        cluster.vz -= cluster.vzc

        cluster.origin = "centre"

        cluster.rv3d()

    if do_key_params:
        cluster.key_params(do_order=do_order)

def to_cluster(self, do_order=False, do_key_params=False, centre_method=None):
    """
    NAME:

       to_cluster

    PURPOSE:

       Shift coordinates to clustercentric reference frame

    INPUT:

       do_order - re-sort cluster radii (default: False)

       do_key_params - call key_params after shift (default: False)

       centre_method - method for shifting to clustercentric coordinates

    OUTPUT:

        None

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    if centre_method is not None:
        cluster.centre_method = centre_method

    if cluster.origin != "cluster":
        if cluster.units == "radec" and cluster.origin == "sky":
            ra = np.radians(cluster.x)
            dec = np.radians(cluster.y)
            pmra = np.radians(cluster.vx / (1000.0 * 3600.0))
            pmdec = np.radians(cluster.vy / (1000.0 * 3600.0))
            ra_gc = np.radians(cluster.xgc)
            dec_gc = np.radians(cluster.ygc)

            if cluster.centre_method == "orthographic":
                cluster.x = np.cos(dec) * np.sin(ra - ra_gc)
                cluster.y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
                    dec_gc
                ) * np.cos(ra - ra_gc)
                cluster.z = np.zeros(len(cluster.x))

                cluster.vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
                    ra - ra_gc
                )
                cluster.vy = pmra * np.sin(dec_gc) * np.sin(
                    ra - ra_gc
                ) + cluster.pmdec * (
                    np.cos(dec) * np.cos(dec_gc)
                    + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
                )
                cluster.vz = np.zeros(len(cluster.x))

            else:
                if cluster.centre_method == "VandeVen":
                    cluster.x = (
                        (10800.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc) / 60.0
                    )
                    cluster.y = (
                        (10800.0 / np.pi)
                        * (
                            np.sin(dec) * np.cos(dec_gc)
                            - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
                        )
                        / 60.0
                    )
                else:
                    cluster.x = (cluster.x - cluster.xgc) * np.cos(np.radians(cluster.ygc))
                    cluster.y = cluster.y - cluster.ygc

                cluster.z = np.zeros(len(cluster.x))

                cluster.vx -= -cluster.vxgc
                cluster.vy -= cluster.vygc
                cluster.vz -= cluster.vzgc

        elif cluster.origin == "centre":
            cluster.x += cluster.xc
            cluster.y += cluster.yc
            cluster.z += cluster.zc
            cluster.vx += cluster.vxc
            cluster.vy += cluster.vyc
            cluster.vz += cluster.vzc
        elif cluster.origin == "galaxy":
            cluster.x -= cluster.xgc
            cluster.y -= cluster.ygc
            cluster.z -= cluster.zgc
            cluster.vx -= cluster.vxgc
            cluster.vy -= cluster.vygc
            cluster.vz -= cluster.vzgc

        cluster.rv3d()

        cluster.origin = "cluster"
        if do_key_params:
            cluster.key_params(do_order=do_order)

def to_galaxy(self, do_order=False, do_key_params=False):
    """
    NAME:

       to_galaxy

    PURPOSE:

       Shift coordinates to galactocentric reference frame

    INPUT:

       None

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)
    """
    if cluster.units == "radec" and cluster.origin == "sky":
        cluster.from_radec(do_key_params=False)

    elif cluster.origin != "galaxy":
        if cluster.origin == "centre":
            cluster.to_cluster(do_key_params=False)

        cluster.x += cluster.xgc
        cluster.y += cluster.ygc
        cluster.z += cluster.zgc
        cluster.vx += cluster.vxgc
        cluster.vy += cluster.vygc
        cluster.vz += cluster.vzgc

        cluster.origin = "galaxy"

    cluster.rv3d()

    if do_key_params:
        cluster.key_params(do_order=do_order)

def to_sky(self, do_order=False, do_key_params=False):
    """
    NAME:

       to_sky

    PURPOSE:

       Calculate on-sky position, proper motion, and radial velocity of cluster
        --> Also changes units to radec
    INPUT:

       None

    OUTPUT:

        None

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    cluster.to_radec(do_key_params=do_key_params)

def from_sky(self, do_order=False, do_key_params=False):
    """
    NAME:

       from_sky

    PURPOSE:

       Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster
        -> Also changes units to kpckms
    INPUT:

       None

    OUTPUT:

        None

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    cluster.from_radec(do_order=do_order, do_key_params=do_key_params)

def to_tail(self, plot=False):
    """
    NAME:

       to_tail

    PURPOSE:

       Calculate positions and velocities of stars when rotated such that clusters velocity vector
       points along x-axis

    INPUT:

       None

    OUTPUT:

        x_rot,y_rot,z_rot,vx_rot,vy_rot,vz_rot

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    units0, origin0 = cluster.units, cluster.origin

    cluster.to_centre()

    v_vec = np.array([cluster.vxgc, cluster.vygc, cluster.vzgc])
    new_v_vec = np.array([1.0, 0.0, 0.0])

    rot = _rotate_to_arbitrary_vector(
        np.atleast_2d(v_vec), new_v_vec, inv=False, _dontcutsmall=False
    )

    x_tail = (
        cluster.x * rot[:, 0, 0] + cluster.y * rot[:, 1, 0] + cluster.z * rot[:, 2, 0]
    )
    y_tail = (
        cluster.x * rot[:, 0, 1] + cluster.y * rot[:, 1, 1] + cluster.z * rot[:, 2, 1]
    )
    z_tail = (
        cluster.x * rot[:, 0, 2] + cluster.y * rot[:, 1, 2] + cluster.z * rot[:, 2, 2]
    )
    vx_tail = (
        cluster.vx * rot[:, 0, 0] + cluster.vy * rot[:, 1, 0] + cluster.vz * rot[:, 2, 0]
    )
    vy_tail = (
        cluster.vx * rot[:, 0, 1] + cluster.vy * rot[:, 1, 1] + cluster.vz * rot[:, 2, 1]
    )
    vz_tail = (
        cluster.vx * rot[:, 0, 2] + cluster.vy * rot[:, 1, 2] + cluster.vz * rot[:, 2, 2]
    )

    cluster.to_origin(origin0)

    return x_tail,y_tail,z_tail,vx_tail,vy_tail,vz_tail

def to_origin(self, origin, do_order=False, do_key_params=False):
    """
    NAME:

       to_origin

    PURPOSE:

       Shift cluster to origin as defined by keyword

    INPUT:

       origin - accepts 'cluster','centre','galaxy'

    OUTPUT:

        None

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    if origin == "centre":
        cluster.to_centre(do_order=do_order, do_key_params=do_key_params)
    elif origin == "cluster":
        cluster.to_cluster(do_order=do_order, do_key_params=do_key_params)
    elif origin == "galaxy":
        cluster.to_galaxy(do_order=do_order, do_key_params=do_key_params)
    elif origin == "sky":
        cluster.to_sky(do_order=do_order, do_key_params=do_key_params)


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
    cluster.to_pckms()

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

def convert_binary_units(self,param,from_units,to_units):
    yrs = (cluster.rbar * 1296000.0 / (2.0 * np.pi)) ** 1.5 / np.sqrt(cluster.zmbar)
    days = 365.25 * yrs
    au = 1.49597870700e13
    pc = 1296000.0e0/(2.0*np.pi)*au
    rsun=6.960e10
    su=pc/rsun*cluster.rbar


    param=np.array(param)
    from_units=np.array(from_units)
    tp_units=np.array(to_units)

    for i in range(0,len(param)):
        p=param[i]

        if p=='pb':
            #Convert to nbody first
            if from_units[i]=='days':
                cluster.pb/=days
            elif from_units[i]=='years':
                cluster.pb/=yrs
            elif from_units[i]=='nbody':
                pass
            else:
                print('UNIT %s NOT FOUND' % from_units[i])

            if to_units[i]=='days':
                cluster.pb*=days
            elif to_units[i]=='years':
                cluster.pb*=yrs
            elif to_units[i]=='nbody':
                pass
            else:
                print('UNIT %s NOT FOUND' % from_units[i])

        elif p=='semi':
            #Convert to nbody first
            if from_units[i]=='pc':
                cluster.semi/=cluster.rbar
            elif from_units[i]=='su':
                cluster.semi/=su
            elif from_units[i]=='au':
                cluster.semi/=(pc/au)*cluster.rbar
            elif from_units[i]=='nbody':
                pass
            else:
                print('UNIT %s NOT FOUND' % from_units[i])

            if to_units[i]=='pc':
                cluster.semi*=cluster.rbar
            elif to_units[i]=='su':
                cluster.semi*=su
            elif to_units[i]=='au':
                cluster.semi*=(pc/au)*cluster.rbar
            elif to_units[i]=='nbody':
                pass
            else:
                print('UNIT %s NOT FOUND' % to_units[i])

        elif p=='mass':
            if from_units[i]=='Msun' or from_units[i]=='msun' :
                cluster.m1/=cluster.zmbar
                cluster.m2/=cluster.zmbar
            elif from_units[i]=='nbody':
                pass

            if to_units=='Msun' or to_units[i]=='msun' :
                cluster.m1*=cluster.zmbar
                cluster.m2*=cluster.zmbar
            elif to_units[i]=='nbody':
                pass

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
