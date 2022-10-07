""" For quickly writing key cluster properties to previously opened output files

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "snapout",
    "sseout",

    "fortout",
    "gyrout",
]

import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

from .coordinates import sky_coords
from .constants import *
from ..util.units import _convert_length,_convert_time,_convert_velocity



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

            if len(cluster.kw)==cluster.ntot:
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
                        ]
                    ),
                )

    return 0

def sseout(cluster, filename):
    """Output simple stellar evolution information in clustertools format
        - clustertools column format is initial mass, kw type, current mass, epoch and ospin
        masses are in solar masses, epoch and spin are in Nbody6 units
       
    Parameters
    ----------
    cluster : class
        StarCluster
    filenamw : str
        name of file to be written to

    Returns
    -------
    None

    History
    -------
    2021 - Written - Webb (UofT)
    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    cluster.to_pckms()

    np.savetxt(
        filename,
        np.column_stack(
            [
                cluster.m0,
                cluster.kw,
                cluster.m,
                cluster.ep,
                cluster.ospin,
            ]
        ),fmt=('%f %i %f %f %f')
    )

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

def fortout(
    cluster,
    filename="fort.10",
    sse=False,
    sse_filename='fort.12',
    reset_nbody=False,
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
    sse : bool
        include file containing stellar evolution information
    sse_filename : str
        name of stellar evolution file to be written to (default: 'fort.12')
    reset_nbody : bool
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
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    cluster.to_centre(sortstars=False)

    if reset_nbody:
        cluster.reset_nbody_scale(mass=reset_nbody_mass, radii=reset_nbody_radii)

    cluster.to_nbody()

    np.savetxt(
        cluster.wdir+filename,
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

    cluster.to_pckms()

    if sse:
        np.savetxt(
            cluster.wdir+sse_filename,
            np.column_stack(
                [
                    cluster.m0,
                    cluster.kw,
                    cluster.m,
                    cluster.ep,
                    cluster.ospin
                ]
            ),fmt=('%f %i %f %f %f'),
        )       

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    return 0


def gyrout(cluster, filename="init.nemo.dat",eps=None,epsunits=None,ro=solar_ro,vo=solar_vo):
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
        galpy spaital scaling parameter, in case epsunits=='galpy' (default: solar_ro)
    vo : float
        galpy velocity scaling parameter, in case epsunits=='galpy' (default: solar_vo)

    Returns
    -------
    None

    History
    -------
    2019 - Written - Webb (UofT)

    """
    vcon = 220.0 / conversion.velocity_in_kpcGyr(vo=solar_vo, ro=solar_ro)
    mcon = 222288.4543021174

    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    cluster.to_galaxy(sortstars=False)
    cluster.to_kpckms()

    if eps is not None:

        if epsunits is None:
            epsunits=units0

        eps=_convert_length(eps,epsunits,cluster)

        if isinstance(eps,float):
            eps=np.ones(cluster.ntot)*eps

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

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    return 0
