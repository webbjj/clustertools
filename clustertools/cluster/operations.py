""" Perform an operation on a cluster and return a new cluster

"""

__author__ = "Jeremy J Webb"

__all__ = [
    "to_pckms",
    "to_kpckms",
    "to_nbody",
    "to_pcmyr",
    "to_kpcgyr",
    "to_radec",
    "to_galpy",
    "to_WDunits",
    "to_amuse",
    "from_amuse",
    "to_units",
    "to_sudays",
    "to_audays",
    "to_centre",
    "to_cluster",
    "to_galaxy",
    "to_sky",
    "to_origin",
    "save_cluster",
    "return_cluster",
    "reset_nbody_scale",
    "add_rotation",
    "virialize",
]

from ..util.constants import *
import numpy as np
from galpy.orbit import Orbit
from galpy.util import _rotate_to_arbitrary_vector
try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion
    
from copy import copy
try:
    import numba
    full_default=True
except:
    full_default=False
    pass

try:    
    import amuse.units.units as u
except:
    pass

def to_pckms(cluster,analyze=True):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and km/s

    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)
    Returns
    -------
    None

    History
    -------
    2018 - Written - Webb (UofT)

    """

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units == "galpy" or cluster.units == "radec" or cluster.units =="WDunits" or cluster.units=='kpcgyr':
            cluster.to_kpckms(analyze=False)
        elif cluster.units == 'amuse':
            from_amuse(cluster)

        if cluster.units == "nbody":
            cluster.tphys *= cluster.tbar
            cluster.m *= cluster.zmbar
            cluster.x *= cluster.rbar
            cluster.y *= cluster.rbar
            cluster.z *= cluster.rbar
            cluster.vx *= cluster.vbar
            cluster.vy *= cluster.vbar
            cluster.vz *= cluster.vbar

            cluster.xc *= cluster.rbar
            cluster.yc *= cluster.rbar
            cluster.zc *= cluster.rbar
            cluster.vxc *= cluster.vbar
            cluster.vyc *= cluster.vbar
            cluster.vzc *= cluster.vbar

            cluster.xgc *= cluster.rbar
            cluster.ygc *= cluster.rbar
            cluster.zgc *= cluster.rbar
            cluster.vxgc *= cluster.vbar
            cluster.vygc *= cluster.vbar
            cluster.vzgc *= cluster.vbar

            cluster.units = "pckms"


        elif cluster.units == "kpckms":

            cluster.tphys *= 1000.0

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

        elif cluster.units == 'pcmyr':
            cluster.vx/=1.022712165045695
            cluster.vy/=1.022712165045695
            cluster.vz/=1.022712165045695

            cluster.vxgc/=1.022712165045695
            cluster.vygc/=1.022712165045695
            cluster.vzgc/=1.022712165045695

            cluster.vxc/=1.022712165045695
            cluster.vyc/=1.022712165045695
            cluster.vzc/=1.022712165045695

            cluster.units = "pckms"

        if analyze: cluster.analyze()

def to_kpckms(cluster,analyze=True):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and km/s

    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units == "radec":
            from_radec(cluster)

        elif cluster.units == 'pcmyr':
            cluster.to_pckms(analyze=False)
        elif cluster.units == 'amuse':
            from_amuse(cluster)

        if cluster.units == "galpy":
            cluster.tphys *= conversion.time_in_Gyr(ro=ro,vo=vo)
            cluster.m *= conversion.mass_in_msol(ro=ro, vo=vo)
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
            cluster.tphys *= cluster.tbar / 1000.0
            cluster.m *= cluster.zmbar
            cluster.x *= cluster.rbar / 1000.0
            cluster.y *= cluster.rbar / 1000.0
            cluster.z *= cluster.rbar / 1000.0
            cluster.vx *= cluster.vbar
            cluster.vy *= cluster.vbar
            cluster.vz *= cluster.vbar

            cluster.xc *= cluster.rbar / 1000.0
            cluster.yc *= cluster.rbar / 1000.0
            cluster.zc *= cluster.rbar / 1000.0
            cluster.vxc *= cluster.vbar
            cluster.vyc *= cluster.vbar
            cluster.vzc *= cluster.vbar

            cluster.xgc *= cluster.rbar / 1000.0
            cluster.ygc *= cluster.rbar / 1000.0
            cluster.zgc *= cluster.rbar / 1000.0
            cluster.vxgc *= cluster.vbar
            cluster.vygc *= cluster.vbar
            cluster.vzgc *= cluster.vbar

            cluster.units = "kpckms"

        elif cluster.units == "pckms":

            cluster.tphys /= 1000.0

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

        elif cluster.units == 'kpcgyr':
            cluster.vx/=1.022712165045695
            cluster.vy/=1.022712165045695
            cluster.vz/=1.022712165045695

            cluster.vxgc/=1.022712165045695
            cluster.vygc/=1.022712165045695
            cluster.vzgc/=1.022712165045695

            cluster.vxc/=1.022712165045695
            cluster.vyc/=1.022712165045695
            cluster.vzc/=1.022712165045695

            cluster.units = "kpckms"

        elif cluster.units == 'WDunits':
            vcon = vo / conversion.velocity_in_kpcGyr(vo=vo, ro=ro)
            mcon = 222288.4543021174

            cluster.m*=mcon
            cluster.vx*=vcon
            cluster.vy*=vcon
            cluster.vz*=vcon
            cluster.vxc*=vcon
            cluster.vyc*=vcon
            cluster.vzc*=vcon
            cluster.vxgc*=vcon
            cluster.vygc*=vcon
            cluster.vzgc*=vcon
            cluster.units = "kpckms"


        if analyze: cluster.analyze()

def to_pcmyr(cluster,analyze=True):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to pc and pc/Myr
       
    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2022 - Written - Webb (UofT)

    """
    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units!='pckms' and cluster.units!='pcmyr':
            cluster.to_pckms(analyze=False)


        if cluster.units!='pcmyr':
            cluster.vx*=1.022712165045695
            cluster.vy*=1.022712165045695
            cluster.vz*=1.022712165045695

            cluster.vxgc*=1.022712165045695
            cluster.vygc*=1.022712165045695
            cluster.vzgc*=1.022712165045695

            cluster.vxc*=1.022712165045695
            cluster.vyc*=1.022712165045695
            cluster.vzc*=1.022712165045695

            cluster.units='pcmyr'
            if analyze: cluster.analyze()

def to_kpcgyr(cluster,analyze=True):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to kpc and kpc/Gyr
       
    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2022 - Written - Webb (UofT)

    """
    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units!='kpckms' and cluster.units!='kpcgyr':
            cluster.to_kpckms(analyze=False)

        if cluster.units!='kpcgyr':
            cluster.vx*=1.022712165045695
            cluster.vy*=1.022712165045695
            cluster.vz*=1.022712165045695

            cluster.vxgc*=1.022712165045695
            cluster.vygc*=1.022712165045695
            cluster.vzgc*=1.022712165045695

            cluster.vxc*=1.022712165045695
            cluster.vyc*=1.022712165045695
            cluster.vzc*=1.022712165045695

            cluster.units='kpcgyr'
            if analyze: cluster.analyze()


def to_nbody(cluster,analyze=True):
    """Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Nbody units
       
    - requires that cluster.zmbar, cluster.rbar, cluster.vbar are set (defaults are 1)

    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units != "pckms" and cluster.units != "nbody":
            cluster.to_pckms(analyze=False)

        if cluster.units != "nbody":
            cluster.tphys /= cluster.tbar
            cluster.m /= cluster.zmbar
            cluster.x /= cluster.rbar
            cluster.y /= cluster.rbar
            cluster.z /= cluster.rbar
            cluster.vx /= cluster.vbar
            cluster.vy /= cluster.vbar
            cluster.vz /= cluster.vbar

            cluster.xc /= cluster.rbar
            cluster.yc /= cluster.rbar
            cluster.zc /= cluster.rbar
            cluster.vxc /= cluster.vbar
            cluster.vyc /= cluster.vbar
            cluster.vzc /= cluster.vbar

            cluster.xgc /= cluster.rbar
            cluster.ygc /= cluster.rbar
            cluster.zgc /= cluster.rbar
            cluster.vxgc /= cluster.vbar
            cluster.vygc /= cluster.vbar
            cluster.vzgc /= cluster.vbar

            cluster.units = "nbody"

        if cluster.bunits == 'audays' and cluster.nb>0:
            cluster.semi/=cluster.rbar_au
            cluster.pb/=cluster.tbar_days
            cluster.m1/=cluster.zmbar
            cluster.m2/=cluster.zmbar
            cluster.bunits='nbody'
        elif cluster.bunits == 'sudays' and cluster.nb>0:
            cluster.semi/=cluster.rbar_su
            cluster.pb/=cluster.tbar_days
            cluster.m1/=cluster.zmbar
            cluster.m2/=cluster.zmbar
            cluster.bunits='nbody'

        cluster.analyze()

def to_radec(cluster, sortstars=True,centre_method=None,analyze=True):
    """Convert to on-sky position, proper motion, and radial velocity of cluster
    
    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)
    """

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:

        origin0=cluster.origin

        if cluster.units != "radec":

            cluster.to_kpckms(analyze=False)
            cluster.to_galaxy()

            if cluster.orbit is None: cluster.initialize_orbit(ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
            if cluster.orbits is None: cluster.initialize_orbits(ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)

            R, phi, z = coords.rect_to_cyl(cluster.xgc+cluster.xc, cluster.ygc+cluster.yc, cluster.zgc+cluster.zc)
            vR, vT, vz = coords.rect_to_cyl_vec(cluster.vxgc+cluster.vxc, cluster.vygc+cluster.vyc, cluster.vzgc+cluster.vzc,cluster.xgc+cluster.xc, cluster.ygc+cluster.yc, cluster.zgc+cluster.zc)

            oc = Orbit(
                [R/ro, vR/vo, vT/vo, z/ro, vz/vo, phi], ro=ro, vo=vo, zo=zo, solarmotion=solarmotion
            )

            cluster.ra=cluster.orbits.ra()
            cluster.dec=cluster.orbits.dec()
            cluster.dist=cluster.orbits.dist()
            cluster.pmra=cluster.orbits.pmra()
            cluster.pmdec=cluster.orbits.pmdec()
            cluster.vlos=cluster.orbits.vlos()

            cluster.ra_gc=cluster.orbit.ra()
            cluster.dec_gc=cluster.orbit.dec()
            cluster.dist_gc=cluster.orbit.dist()
            cluster.pmra_gc=cluster.orbit.pmra()
            cluster.pmdec_gc=cluster.orbit.pmdec()
            cluster.vlos_gc=cluster.orbit.vlos()

            cluster.ra_c=oc.ra()
            cluster.dec_c=oc.dec()
            cluster.dist_c=oc.dist()
            cluster.pmra_c=oc.pmra()
            cluster.pmdec_c=oc.pmdec()
            cluster.vlos_c=oc.vlos()

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

            cluster.xc = copy(cluster.ra_c)
            cluster.yc = copy(cluster.dec_c)
            cluster.zc = copy(cluster.dist_c)
            cluster.vxc = copy(cluster.pmra_c)
            cluster.vyc = copy(cluster.pmdec_c)
            cluster.vzc = copy(cluster.vlos_c)

            cluster.units = "radec"
            cluster.origin = "sky"

            if origin0=='cluster':
                cluster.to_cluster(centre_method=centre_method)
            elif origin0=='centre':
                cluster.to_centre(centre_method=centre_method)

            if analyze: cluster.analyze(sortstars=sortstars)

def from_radec(cluster, sortstars=True):
    """Calculate galactocentric coordinates from on-sky position, proper motion, and radial velocity of cluster

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units=='radec' and cluster.origin!='sky':
            cluster.to_sky()

        if cluster.units == "radec" and cluster.origin == "sky":

            if cluster.orbit is None: cluster.initialize_orbit(ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)
            if cluster.orbits is None: cluster.initialize_orbits(ro=ro,vo=vo,zo=zo,solarmotion=solarmotion)


            if cluster.xcn is not None:
                cluster.xc,cluster.yc,cluster.zc=cluster.xcn,cluster.ycn,cluster.zcn
                cluster.xc*=cluster.rbar/1000.
                cluster.yc*=cluster.rbar/1000.
                cluster.zc*=cluster.rbar/1000.
                cluster.vxc=0.
                cluster.vyc=0.
                cluster.vzc=0.

            else:
                oc = Orbit(
                    [cluster.ra_c,cluster.dec_c,cluster.dist_c,cluster.pmra_c,cluster.pmdec_c,cluster.vlos_c], radec=True,ro=ro, vo=vo, zo=zo, solarmotion=solarmotion
                )

                cluster.xc=oc.x()
                cluster.yc=oc.y()
                cluster.zc=oc.z()
                cluster.vxc=oc.vx()
                cluster.vyc=oc.vy()
                cluster.vzc=oc.vz()

            cluster.x = cluster.orbits.x()
            cluster.y = cluster.orbits.y()
            cluster.z = cluster.orbits.z()
            cluster.vx = cluster.orbits.vx()
            cluster.vy = cluster.orbits.vy()
            cluster.vz = cluster.orbits.vz()

            cluster.xgc = cluster.orbit.x()
            cluster.ygc = cluster.orbit.y()
            cluster.zgc = cluster.orbit.z()
            cluster.vxgc = cluster.orbit.vx()
            cluster.vygc = cluster.orbit.vy()
            cluster.vzgc = cluster.orbit.vz()

            cluster.origin = "galaxy"
            cluster.units = "kpckms"

def to_galpy(cluster,analyze=True):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to galpy units
    
    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if cluster.units != "kpckms" and cluster.units != "galpy":
            cluster.to_kpckms(analyze=False)

        if cluster.units == "kpckms":
            cluster.tphys /= conversion.time_in_Gyr(ro=ro, vo=vo)
            cluster.m = cluster.m / conversion.mass_in_msol(ro=ro, vo=vo)
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

        cluster.analyze()

def to_WDunits(cluster,analyze=True):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to Walter Dehnen Units

    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2022 - Written - Webb (UofT)

    """

    ro,vo,zo,solarmotion=cluster._ro,cluster._vo,cluster._zo,cluster._solarmotion

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        vcon = vo / conversion.velocity_in_kpcGyr(vo=vo,ro=ro)
        mcon = 222288.4543021174

        if cluster.units!='kpckms' and cluster.units!="WDunits":
            cluster.to_kpckms(analyze=False)

        if cluster.units!="WDunits":

            cluster.m/=mcon
            cluster.vx/=vcon
            cluster.vy/=vcon
            cluster.vz/=vcon
            cluster.vxc/=vcon
            cluster.vyc/=vcon
            cluster.vzc/=vcon
            cluster.vxgc/=vcon
            cluster.vygc/=vcon
            cluster.vzgc/=vcon

            cluster.units = "WDunits"

            cluster.analyze()

def to_amuse(cluster,analyze=True):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to AMUSE Vector and Scalar Quantities

    Parameters
    ----------
    cluster : class
        StarCluster
    analyze : bool
        run analysis function (default: True)

    Returns
    -------
    None

    History:
    -------
    2022 - Written - Webb (UofT)

    """

    if cluster.units != 'amuse':

        cluster.to_pckms(analyze=False)

        cluster.tphys = cluster.tphys | u.Myr

        cluster.x=cluster.x | u.parsec
        cluster.y=cluster.y | u.parsec
        cluster.z=cluster.z | u.parsec
        cluster.vx=cluster.vx | u.kms
        cluster.vy=cluster.vy | u.kms
        cluster.vz=cluster.vz | u.kms
        cluster.m=cluster.m | u.MSun

        cluster.xc=cluster.xc | u.parsec
        cluster.yc=cluster.yc | u.parsec
        cluster.zc=cluster.zc | u.parsec
        cluster.vxc=cluster.vxc | u.kms
        cluster.vyc=cluster.vyc | u.kms
        cluster.vzc=cluster.vzc | u.kms

        cluster.xgc=cluster.xgc | u.parsec
        cluster.ygc=cluster.ygc | u.parsec
        cluster.zgc=cluster.zgc | u.parsec
        cluster.vxgc=cluster.vxgc | u.kms
        cluster.vygc=cluster.vygc | u.kms
        cluster.vzgc=cluster.vzgc | u.kms


        cluster.units = "amuse"

        if analyze: cluster.analyze()

def from_amuse(cluster):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity fronm AMUSE Vector and Scalar Quantities to pckms

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    None

    History:
    -------
    2022 - Written - Webb (UofT)

    """

    if cluster.units == 'amuse':

        cluster.tphys = cluster.tphys.value_in(u.Myr)

        cluster.x=cluster.x.value_in(u.parsec)
        cluster.y=cluster.y.value_in(u.parsec)
        cluster.z=cluster.z.value_in(u.parsec)
        cluster.vx=cluster.vx.value_in(u.kms)
        cluster.vy=cluster.vy.value_in(u.kms)
        cluster.vz=cluster.vz.value_in(u.kms)
        cluster.m=cluster.m.value_in(u.MSun)

        cluster.xc=cluster.xc.value_in(u.parsec)
        cluster.yc=cluster.yc.value_in(u.parsec)
        cluster.zc=cluster.zc.value_in(u.parsec)
        cluster.vxc=cluster.vxc.value_in(u.kms)
        cluster.vyc=cluster.vyc.value_in(u.kms)
        cluster.vzc=cluster.vzc.value_in(u.kms)

        cluster.xgc=cluster.xgc.value_in(u.parsec)
        cluster.ygc=cluster.ygc.value_in(u.parsec)
        cluster.zgc=cluster.zgc.value_in(u.parsec)
        cluster.vxgc=cluster.vxgc.value_in(u.kms)
        cluster.vygc=cluster.vygc.value_in(u.kms)
        cluster.vzgc=cluster.vzgc.value_in(u.kms)


        cluster.units = "pckms"


def to_units(cluster, units):
    """ Convert stellar positions/velocities, centre of mass, and orbital position and velocity to user defined units

    Parameters
    ----------
    cluster : class
        StarCluster
    units : str
        units to be converted to (nbody,galpy,pckms,kpckms,radec)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """

    if cluster.units is None:
        print('NO UNITS SPECIFIED')
    else:
        if units == "nbody":
            cluster.to_nbody()
        elif units == "galpy":
            cluster.to_galpy()
        elif units == "pckms":
            cluster.to_pckms()
        elif units == "kpckms":
            cluster.to_kpckms()
        elif units == "radec":
            cluster.to_radec(sortstars=False,)
        elif units == "WDunits":
            cluster.to_WDunits()
        elif units == "pcmyr":
            cluster.to_pcmyr()
        elif units == "kpcgyr":
            cluster.to_kpcgyr()
        elif units == 'amuse':
            cluster.to_amuse()

def to_sudays(cluster):
    """ Convert binary star semi-major axis and periods to solar radii and days
        Note: Masses will be converted to solar masses

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    None

    History
    -------
    2022 - Written - Webb (UofT)

    """
    if cluster.bunits is not None:
        if cluster.bunits == 'audays':
            cluster.semi/=cluster.rbar_au
            cluster.pb/=cluster.tbar_days
            cluster.m1/=cluster.zmbar
            cluster.m2/=cluster.zmbar
            cluster.bunits='nbody'

        if cluster.bunits == 'nbody':
            cluster.semi*=cluster.rbar_su
            cluster.pb*=cluster.tbar_days
            cluster.m1*=cluster.zmbar
            cluster.m2*=cluster.zmbar
        cluster.bunits='sudays'
    else:
        print("NO UNITS SPECIFIED")

def to_audays(cluster):
    """ Convert binary star semi-major axis and periods to solar radii and days

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    None

    History
    -------
    2022 - Written - Webb (UofT)

    """
    if cluster.bunits is not None:
        if cluster.bunits == 'sudays':
            cluster.semi/=cluster.rbar_su
            cluster.pb/=cluster.tbar_days
            cluster.m1/=cluster.zmbar
            cluster.m2/=cluster.zmbar
            cluster.bunits='nbody'

        if cluster.bunits == 'nbody':
            cluster.semi*=cluster.rbar_au
            cluster.pb*=cluster.tbar_days
            cluster.m1*=cluster.zmbar
            cluster.m2*=cluster.zmbar
            cluster.bunits='audays'
    else:
        print("NO UNITS SPECIFIED")

def to_centre(cluster, sortstars=True, centre_method=None):
    """Shift coordinates such that origin is the centre of mass

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:

        if centre_method is not None:
            cluster.centre_method = centre_method

        if cluster.origin != 'centre':

            if cluster.units=='radec':
                if cluster.origin!='sky':
                    cluster.to_sky()

                if cluster.centre_method == "orthographic":
                    _get_orthographic_centre(cluster,from_centre=True)
                else:
                    if cluster.centre_method == "VandeVen":
                        _get_vandeven_centre(cluster,from_centre=True)
                    else:
                        cluster.x = (cluster.x - cluster.xc) * np.cos(np.radians(cluster.yc))
                        cluster.y = cluster.y - cluster.yc

                        cluster.z -= cluster.zc

                        cluster.vx -= cluster.vxc
                        cluster.vy -= cluster.vyc
                        cluster.vz -= cluster.vzc


            else:
                if cluster.origin != "cluster":
                    cluster.to_cluster(centre_method=centre_method)

                cluster.x -= cluster.xc
                cluster.y -= cluster.yc
                cluster.z -= cluster.zc
                cluster.vx -= cluster.vxc
                cluster.vy -= cluster.vyc
                cluster.vz -= cluster.vzc

            cluster.origin = "centre"

            cluster.analyze(sortstars=sortstars)

def to_center(cluster, sortstars=True, centre_method=None):
    """Shift coordinates such that origin is the centre of mass

    - this is the same to to_centre, but allows for alternate spelling

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    to_centre(cluster, sortstars=sortstars, centre_method=centre_method)


def to_cluster(cluster, sortstars=True, centre_method=None):
    """Shift coordinates to clustercentric reference frame

    - When units='radec' and origin='sky', the cluster will be shifted to clustercentric coordinates using either:
    --centre_method=None: angular distance between each star's RA/DEC and the RA/DEC of the cluster's centre with proper motions directly subtracted off
    --centre_method='orthographic' , positions and velocities changed to orthnormal coordinates (Helmi et al. 2018)
    -- centre_method='VandeVen' , positions and velocities changed to clustercentric coordinates using method outlined by Van de Ven et al. 2005

    - Note the the line of site positions and velocities will just have the galactocentric coordinates of the cluster
    subtracted off. Be sure to set projected=True when making any calculations to use only x and y coordinates

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates. (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:

        if centre_method is not None:
            cluster.centre_method = centre_method

        if cluster.origin != "cluster":
            if cluster.units == "radec":

                if cluster.origin!='sky':
                    cluster.to_sky()

                if cluster.centre_method == "orthographic":
                    _get_orthographic_centre(cluster)
                else:
                    if cluster.centre_method == "VandeVen":
                        _get_vandeven_centre(cluster)
                    else:
                        cluster.x = (cluster.x - cluster.xgc) * np.cos(np.radians(cluster.ygc))
                        cluster.y = cluster.y - cluster.ygc

                        cluster.z -= cluster.zgc

                        cluster.vx -= cluster.vxgc
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


            cluster.origin = "cluster"

            cluster.analyze(sortstars=sortstars)

def _get_orthographic_centre(cluster,from_centre=False):

    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:


        if cluster.units!='radec' and cluster.origin!='sky':
            cluster.to_radec()
        elif cluster.units=='radec' and cluster.origin!='sky':
            cluster.to_sky()

        ra = copy(np.radians(cluster.x))
        dec = copy(np.radians(cluster.y))
        pmra = copy(cluster.vx) # / (1000.0 * 3600.0)))
        pmdec = copy(cluster.vy) # / (1000.0 * 3600.0)))

        if from_centre:
            ra_gc = copy(np.radians(cluster.xc))
            dec_gc = copy(np.radians(cluster.yc))
        else:
            ra_gc = copy(np.radians(cluster.xgc))
            dec_gc = copy(np.radians(cluster.ygc))

        cluster.x = np.cos(dec) * np.sin(ra - ra_gc)
        cluster.y = np.sin(dec) * np.cos(dec_gc) - np.cos(dec) * np.sin(
            dec_gc
        ) * np.cos(ra - ra_gc)


        cluster.x*=(180.0/np.pi)
        cluster.y*=(180.0/np.pi)

        cluster.vx = pmra * np.cos(ra - ra_gc) - pmdec * np.sin(dec) * np.sin(
            ra - ra_gc
        )
        cluster.vy = pmra * np.sin(dec_gc) * np.sin(
            ra - ra_gc
        ) + cluster.pmdec * (
            np.cos(dec) * np.cos(dec_gc)
            + np.sin(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
        )

        cluster.vx-=cluster.vxgc
        cluster.vy-=cluster.vygc

        if from_centre:
            cluster.z -= cluster.zc
            cluster.vz -= cluster.vzc
        else:
            cluster.z -= cluster.zgc
            cluster.vz -= cluster.vzgc

def _get_vandeven_centre(cluster,from_centre=False):

    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:

        if cluster.units!='radec' and cluster.origin!='sky':
            cluster.to_radec()
        elif cluster.units=='radec' and cluster.origin!='sky':
            cluster.to_sky()

        ra = np.radians(cluster.x)
        dec = np.radians(cluster.y)

        if from_centre:
            ra_gc = np.radians(cluster.xc)
            dec_gc = np.radians(cluster.yc)
        else:
            ra_gc = np.radians(cluster.xgc)
            dec_gc = np.radians(cluster.ygc)

        cluster.x = (
            (180.0 / np.pi) * np.cos(dec) * np.sin(ra - ra_gc)
        )
        cluster.y = (
            (180.0 / np.pi)
            * (
                np.sin(dec) * np.cos(dec_gc)
                - np.cos(dec) * np.sin(dec_gc) * np.cos(ra - ra_gc)
            )
        )

        if from_centre:
            cluster.z -= cluster.zc

            cluster.vx -= cluster.vxc
            cluster.vy -= cluster.vyc
            cluster.vz -= cluster.vzc
        else:
            cluster.z -= cluster.zgc

            cluster.vx -= cluster.vxgc
            cluster.vy -= cluster.vygc
            cluster.vz -= cluster.vzgc

def to_galaxy(cluster, sortstars=True):
    """Shift coordinates to galactocentric reference frame

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """
    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:

        if cluster.units == "radec":
            if cluster.origin!='sky':
                cluster.to_sky()

            from_radec(cluster)

        elif cluster.origin != "galaxy":
            if cluster.origin == "centre":
                cluster.to_cluster(sortstars=False)

            cluster.x += cluster.xgc
            cluster.y += cluster.ygc
            cluster.z += cluster.zgc
            cluster.vx += cluster.vxgc
            cluster.vy += cluster.vygc
            cluster.vz += cluster.vzgc

            cluster.origin = "galaxy"

        cluster.analyze(sortstars=sortstars)


def to_sky(cluster, sortstars=True,centre_method=None):
    """Calculate on-sky position, proper motion, and radial velocity of cluster
        
    - Also changes units to radec

    Parameters
    ----------
    cluster : class
        StarCluster
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates. (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)
    """

    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:

        if cluster.units=='radec' and cluster.origin!='sky':
            cluster.x=copy(cluster.ra)
            cluster.y=copy(cluster.dec)
            cluster.z=copy(cluster.dist)
            cluster.vx=copy(cluster.pmra)
            cluster.vy=copy(cluster.pmdec)
            cluster.vz=copy(cluster.vlos)

            cluster.origin='sky'

        elif cluster.units!='radec':
            cluster.to_galaxy()
            cluster.to_radec(sortstars=sortstars,centre_method=centre_method)

        cluster.analyze(sortstars=sortstars)

def to_origin(cluster, origin, sortstars=True, centre_method=None):
    """Shift cluster to origin as defined by keyword

    Parameters
    ----------
    cluster : class
        StarCluster
    origin : str
        origin of coordinate system to be shifted to (centre, cluster, galaxy, sky)
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)

    """

    if cluster.origin is None:
        print('NO ORIGIN SPECIFIED')
    else:
        if origin == "centre" or origin=="center":
            cluster.to_centre(sortstars=sortstars,centre_method=centre_method)
        elif origin == "cluster":
            cluster.to_cluster(sortstars=sortstars,centre_method=centre_method)
        elif origin == "galaxy":
            cluster.to_galaxy(sortstars=sortstars)
        elif origin == "sky":
            cluster.to_sky(sortstars=sortstars,centre_method=centre_method)
        elif origin == "radec":
            cluster.to_radec(sortstars=sortstars,centre_method=centre_method)


        cluster.analyze(sortstars=sortstars)


def save_cluster(cluster):
    """Save cluster's units and origin

    Parameters
    ----------
    cluster : class
        StarCluster

    Returns
    -------
    cluster.units, cluster.origin : str
        units and origin of StarCluster

    History:
    -------
    2018 - Written - Webb (UofT)
    """

    return cluster.units, cluster.origin, cluster.rorder, cluster.rorder_origin


def return_cluster(cluster, units0, origin0, rorder0, rorder_origin0 ,sortstars=False, centre_method=None):
    """ return cluster to a specific combination of units and origin

    Parameters
    ----------

    cluster : class
        StarCluster
    units0 : str
        units that StarCluster will be changed to
    origin0 : str
        origin that StarCluster will be changed to
    sortstars : bool
        sort star by radius after coordinate change (default: False)
    centre_method : str
        method for shifting coordinates to clustercentric coordinates (see to_cluster). (default: None)

    Returns
    -------
    None

    History:
    -------
    2018 - Written - Webb (UofT)
    """
    if cluster.units != units0:
        cluster.to_units(units0)
    if cluster.origin != origin0:
        cluster.to_origin(origin0, sortstars=sortstars, centre_method=centre_method)

    cluster.rorder=rorder0
    cluster.rorder_origin=rorder_origin0
    cluster.analyze()

def reset_nbody_scale(cluster, mass=True, radii=True, rvirial=True, projected=False,**kwargs):
    """ Assign new conversions for real mass, size, and velocity to Nbody units
    
    - kwargs are passed to the virial_radius function. See the virial_radius documenation in functions.py

    Parameters
    ----------

    cluster : class
        StarCluster instance
    mass : bool
        find new mass conversion (default: True)
    radii : bool
        find new radius conversion (default: True)
    rvirial : bool (default: True)
        use virial radius to set conversion rate for radii as opposed to the approximation that rbar=4/3 rm
    projected : bool
        use projected values to calculate radius conversion (default: False)

    Returns
    -------

    zmbar : float
        mass conversion
    rbar : float
        radius conversion
    vbar : float
        velocity conversion
    tbar : float
        time conversion

    History:

    2018 - Written - Webb (UofT)
    """
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)

    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre(sortstars=True)

    cluster.to_pckms()

    if mass:
        zmbar = np.sum(cluster.m)
    else:
        zmbar = cluster.zmbar

    if radii:
        if rvirial:

            method=kwargs.get("method", 'inverse_distance')
            full=kwargs.get("full", True)
            H = kwargs.get("H", 70.0)
            Om = kwargs.get("Om", 0.3)
            overdens = kwargs.get("overdens", 200.0)
            nrad = kwargs.get("nrad", 20.0)
            plot = kwargs.get("plot", False)

            cluster.virial_radius(method=method,full=full,
                H=H,Om=Om,overdens=overdens,nrad=nrad,projected=projected,
                plot=plot,**kwargs)

            rbar=cluster.rv
        else:
            rbar = 4.0 * cluster.rm / 3.0
    else:
        rbar = cluster.rbar

    vbar = 0.06557 * np.sqrt(zmbar / rbar)

    tbar=rbar/(1.023*vbar)

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return zmbar,rbar,vbar,tbar

def virialize(cluster, qvir=0.5, specific=True, full=full_default, projected=False, softening=0.0):
    """ Adjust stellar velocities so cluster is in virial equilibrium

    Parameters
    ----------
    cluster : class
        StarCluster
    qvir : float
        value you wish to virial parameter to be (default: 0.5)
    specific : bool
        find specific energies (default: True)
    full: bool
        do full array of stars at once with numba (default: full_default)
    projected : bool
        use projected values when calculating energies (default: False)
    softening : float
      Plummer softening length in cluster.units (default: 0.0)
    Returns
    -------
    qv : float
        scaling factor used to adjust velocities

    History
    -------
    2019 - Written - Webb (UofT)
    """
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre()

    try:
        qv = np.sqrt(np.abs(qvir / cluster.qvir))
    except:
        print("NEED TO CALCULATE ENERGIES FIRST")
        cluster.energies(specific=specific, full=full, projected=projected,softening=softening)
        qv = np.sqrt(np.abs(qvir / cluster.qvir))

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return qv

def add_rotation(cluster, qrot):
    """Add a degree of rotation to an already generated StarCluster

    Parameters
    ----------
    cluster : class
        StarCluster
    qrot : float
        fraction of stars with v_phi < 0 that are switched to vphi > 0

    Returns
    -------
    x,y,z,vx,vy,vz : float
        stellar positions and velocities (now with rotation)

    History
    -------
    2019 - Written - Webb (UofT)
    """
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)

    if origin0 != 'cluster' and origin0 != 'centre':
        cluster.to_centre()

    r, theta, z = coords.rect_to_cyl(cluster.x, cluster.y, cluster.z)
    vr, vtheta, vz = coords.rect_to_cyl_vec(
        cluster.vx, cluster.vy, cluster.vz, cluster.x, cluster.y, cluster.z
    )

    indx = vtheta < 0.0
    rindx = np.random.rand(cluster.ntot) < qrot

    vtheta[indx * rindx] = np.sqrt(vtheta[indx * rindx] * vtheta[indx * rindx])
    x,y,z = coords.cyl_to_rect(r, theta, z)
    vx,vy,vz = coords.cyl_to_rect_vec(
        vr, vtheta, vz, theta
    )

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return x,y,z,vx,vy,vz
