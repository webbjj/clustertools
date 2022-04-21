import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

import os, struct
from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from ..util.constants import *
from .orbit import _get_cluster_orbit

def _get_galpy_orbits(
    orbits, units=None, origin=None, ofile=None,ro=solar_ro,vo=solar_vo,**kwargs
):
    """Convert galpy orbits to a StarCluster instance

    Parameters
    ----------
    orbits : orbits
        GALPY orbits instance
    masses : float
        mass of stars (default: None)
    units : str
        units of input data (default: None)
    origin : str
        origin of input data (default: None)
    ofile : file
        opened file containing orbital information
    Returns
    -------
    cluster : class
        StarCluster

    Other Parameters
    ----------------

    m : float
        masses of individual stars (in galpy units)

    Same as load_cluster

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster = StarCluster(
        tphys=0.0,
        units=units,
        origin=origin,
        ctype="galpy",
        **kwargs
    )

    if units=='galpy':
        mo=conversion.mass_in_msol(ro=ro, vo=vo)
    else:
        mo=1.

    i_d = np.linspace(1, len(orbits), len(orbits), dtype="int")
    masses=kwargs.get('m',None)

    if masses is None:
        masses = np.ones(len(orbits))
    elif isinstance(masses,float):
        masses = np.ones(len(orbits))*masses

    if units=='galpy' and (ro != 1. or vo != 1.):
            print('WARNING - YOU HAVE SELECTED GALPY UNITS BUT RO AND/OR VO ARE NOT 1')

    cluster.add_stars(orbits.x(ro=ro,vo=vo), orbits.y(ro=ro,vo=vo), orbits.z(ro=ro,vo=vo), orbits.vx(ro=ro,vo=vo), orbits.vy(ro=ro,vo=vo), orbits.vz(ro=ro,vo=vo), masses/mo, i_d, sortstars=False)

    if origin == "galaxy":
        if ofile == None:
            cluster.find_centre()
        else:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            cluster.to_cluster(sortstars=False)
            cluster.find_centre()
            cluster.to_centre(sortstars=sortstars)
            cluster.to_galaxy()

    elif origin == "cluster" or origin=='centre':
        cluster.find_centre()
        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            cluster.analyze(sortstars=sortstars)

        if ofile != None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    return cluster

