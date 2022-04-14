import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

import os, struct
from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from ..planet.clusterwplanets import StarClusterwPlanets
from .orbit import _get_cluster_orbit

def _get_galpy_orbits(
    orbits, units="kpckms", origin="galaxy", ofile=None,ro=8.,vo=220.,**kwargs
):
    """Convert galpy orbits to a StarCluster instance

    Parameters
    ----------
    orbits : orbits
        GALPY orbits instance
    masses : float
        mass of stars (default: None)
    units : str
        units of input data (default: kpckms)
    origin : str
        origin of input data (default: cluster)
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

    mo=conversion.mass_in_msol(ro=ro, vo=vo)

    i_d = np.linspace(1, len(orbits), len(orbits), dtype="int")

    masses=kwargs.get('m',None)

    if masses is None:
        masses = np.ones(len(orbits))/mo
    elif isinstance(masses,float):
        masses = np.ones(len(orbits))*masses

    if units=='kpckms':
        cluster.add_stars(orbits.x(), orbits.y(), orbits.z(), orbits.vx(), orbits.vy(), orbits.vz(), masses*mo, i_d, sortstars=False)
    elif units=='pckms':
        cluster.add_stars(orbits.x()*1000.0, orbits.y()*1000.0, orbits.z()*1000.0, orbits.vx(), orbits.vy(), orbits.vz(), masses*mo, i_d, sortstars=False)
    elif units=='radec' and origin=='sky':
        cluster.add_stars(orbits.ra(), orbits.dec(), orbits.dist(), orbits.pmra(), orbits.pmdec(), orbits.vlos(), masses*mo, i_d, sortstars=False)
    elif units=='galpy':
        cluster.add_stars(orbits.x(use_physical=False), orbits.y(use_physical=False), orbits.z(use_physical=False), orbits.vx(use_physical=False), orbits.vy(use_physical=False), orbits.vz(use_physical=False), masses, i_d, sortstars=False)

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
        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            cluster.analyze(sortstars=sortstars)

        if ofile != None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    return cluster

