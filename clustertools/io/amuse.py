import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from .orbit import _get_cluster_orbit
from ..util.units import _convert_amuse

# Try Importing AMUSE. Only necessary for _get_amuse_particles
try:
    import amuse.units.units as u
    from amuse.io import read_set_from_file
except:
    pass

#Try importing hdf5. Only necessary with Nbody6++ and hdf5 output
try: 
    import h5py
except:
    pass


def _get_amuse_particles(
    particles, units=None, origin=None, ofile=None, **kwargs
):
    """Convert AMUSE particle dataset to a StarCluster instance

    Parameters
    ----------
    particles : particles
        AMUSE particle dataset
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
    Same as load_cluster

    History
    -------
    2018 - Written - Webb (UofT)
    """

    cluster = StarCluster(
        tphys=0.0,
        units=units,
        origin=origin,
        ctype="amuse",
        **kwargs
    )

    if units is None or units=='amuse':
        cluster.add_stars(particles.x, particles.y, particles.z, particles.vx, particles.vy, particles.vz, particles.mass, particles.key, sortstars=False)
    else:
        x,y,z,vx,vy,vz,m,id=_convert_amuse(particles,cluster)
        cluster.add_stars(x,y,z,vx,vy,vz,m,id)

    if origin == "galaxy":
        if ofile is not None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)
        else:
            cluster.find_centre()

        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            if cluster.origin is not None:
                cluster.to_cluster(sortstars=False)
                cluster.find_centre()
                cluster.to_centre(sortstars=sortstars)
                cluster.to_origin(origin)
            else:
                cluster.find_centre()
                cluster.analyze(sortstars=sortstars)

    elif origin == "cluster" or origin=='centre':
        cluster.find_centre()
        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            cluster.analyze(sortstars=sortstars)

        if ofile != None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    return cluster
