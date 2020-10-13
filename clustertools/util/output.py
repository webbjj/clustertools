""" For quickly writing key cluster properties to previously opened output files

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "snapout",
    "fortout",
    "gyrout",
]

import numpy as np
from galpy.util import bovy_conversion


from .coordinates import sky_coords
from ..analysis.functions import *
from ..analysis.profiles import *
from ..analysis.operations import *


def snapout(cluster, filename, energies=False, radec=False):
    """Output a snapshot in clustertools format
        - clustertools column format is mass,x,y,z,vx,vy,vz,id,kwtype where
        positions and velocities are in cartiesian coordinates, id is the star id, and
        kwtype is the stars stellar evolution type
        - units and coordinate system will be cluster.units and cluster.origin.
        - if energies==True, columns for kinetic energy, potential energy, and total energy will be added
        - if radec==True, columns for RA, Dec, distance, Proper Motion in RA, Proper Motion in Dec, and radial velocity will be added

    Parameters
    ----------
    cluster : class
        StarCluster
    filenamw : str
        name of file to be written to
    energies : bool 
        include energies in output (default: False)
    radec : bool
        include sky coordinates of stars (default: False) 

    Returns
    -------
    None

    History
    -------
    2018 - Written - Webb (UofT)
    """
    if radec:

        ra, dec, d0, pmra, pmdec, vr0 = sky_coords(cluster)

        if energies:
            np.savetxt(
                filename,
                np.column_stack(
                    [
                        cluster.m,
                        cluster.x,
                        cluster.y,
                        cluster.z,
                        cluster.vx,
                        cluster.vy,
                        cluster.vz,
                        cluster.id,
                        cluster.kw,
                        cluster.kin,
                        cluster.pot,
                        cluster.etot,
                        ra,
                        dec,
                        d0,
                        pmra,
                        pmdec,
                        vr0,
                    ]
                ),
            )
        else:
            np.savetxt(
                filename,
                np.column_stack(
                    [
                        cluster.m,
                        cluster.x,
                        cluster.y,
                        cluster.z,
                        cluster.vx,
                        cluster.vy,
                        cluster.vz,
                        cluster.id,
                        cluster.kw,
                        ra,
                        dec,
                        d0,
                        pmra,
                        pmdec,
                        vr0,
                    ]
                ),
            )

    else:
        if energies:
            np.savetxt(
                filename,
                np.column_stack(
                    [
                        cluster.m,
                        cluster.x,
                        cluster.y,
                        cluster.z,
                        cluster.vx,
                        cluster.vy,
                        cluster.vz,
                        cluster.id,
                        cluster.kw,
                        cluster.kin,
                        cluster.pot,
                        cluster.etot,
                    ]
                ),
            )
        else:
            np.savetxt(
                filename,
                np.column_stack(
                    [
                        cluster.m,
                        cluster.x,
                        cluster.y,
                        cluster.z,
                        cluster.vx,
                        cluster.vy,
                        cluster.vz,
                        cluster.id,
                        cluster.kw,
                    ]
                ),
            )

    return 0


def fortout(
    cluster,
    filename="fort.10",
    reset_nbody_scale=False,
    reset_nbody_mass=True,
    reset_nbody_radii=True,
):
    """Output a snapshot in NBODY6 fort.10 format

    Parameters
    ----------
    cluster : class
        StarCluster
    filename : str
        name of file to be written to (default: 'fort.10')
    reset_nbody_scale : bool
        reset nbody scaling parameters (default: False)
    reset_nbody_mass : bool
        reset nbody mass scaling parameter (default: False)
    reset_nbody_radii : bool
        reset nbody radii scaling parameter (default: False)

    Returns
    -------
    None

    History
    -------
    2019 - Written - Webb (UofT)

    """
    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_centre(sortstars=False)

    if reset_nbody_scale:
        reset_nbody_scale(cluster, mass=reset_nbody_mass, radii=reset_nbody_radii)

    cluster.to_nbody()

    np.savetxt(
        filename,
        np.column_stack(
            [
                cluster.m,
                cluster.x,
                cluster.y,
                cluster.z,
                cluster.vx,
                cluster.vy,
                cluster.vz,
            ]
        ),
    )

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return 0


def gyrout(cluster, filename="init.nemo.dat",eps=None,epsunits=None,ro=8.):
    """Output a snapshot in gyrfalcon/NEMO format

    - Note that units are converted to Walter Dehnen units (WDunits in GYRFALCON)

    Parameters
    ----------
    cluster : class
        StarCluster
    filename : str
        name of file to be written to (default: 'init.nemo.dat')
    eps : array
        values for softening length of individual particles (default : None)
    epsunits : str
        units for softening lengths (default: cluster.units)
    ro : float
        galpy spaital scaling parameter, in case epsunits=='gal[y' (default: 8.)

    Returns
    -------
    None

    History
    -------
    2019 - Written - Webb (UofT)

    """
    vcon = 220.0 / bovy_conversion.velocity_in_kpcGyr(220.0, 8.0)
    mcon = 222288.4543021174

    units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)
    cluster.to_galaxy(sortstars=False)
    cluster.to_kpckms()

    if eps is not None:

        if epsunits is None:
            epsunits=units0

        if epsunits=='pckms' or epsunits=='pc':
            eps/=1000.0
        elif epsunits=='nbody':
            eps*=(cluster.rbar/1000.0)
        elif epsunits=='galpy':
            eps*=ro

        np.savetxt(
                filename,
                np.column_stack(
                    [
                        cluster.m / mcon,
                        cluster.x,
                        cluster.y,
                        cluster.z,
                        cluster.vx / vcon,
                        cluster.vy / vcon,
                        cluster.vz / vcon,
                        eps,
                    ]
                ),
            )


    else:


        np.savetxt(
            filename,
            np.column_stack(
                [
                    cluster.m / mcon,
                    cluster.x,
                    cluster.y,
                    cluster.z,
                    cluster.vx / vcon,
                    cluster.vy / vcon,
                    cluster.vz / vcon,
                ]
            ),
        )

    return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)

    return 0
