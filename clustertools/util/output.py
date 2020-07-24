""" For quickly writing key cluster properties to previously opened output files

"""
__author__ = "Jeremy J Webb"
__all__ = [
    "extrct_out",
    "rho_prof_out",
    "alpha_prof_out",
    "sigv_prof_out",
    "eta_prof_out",
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


def extrct_out(cluster, fileout, projected=False):
    """Extrct key cluster properties and write to file
        - Number of stars, number of binary stars, half-mass relaxation time, total mass, 10-100% Lagrange radii, 
        projected 50% radius, projected half-light radius, virial radius, limiting radius, tidal radius, 
        radius of maximal circular velocity, maximal circular velocity
        - Note that limiting radius and tidal radius are only calculated if cluster.rl and cluster.rt are != None

    Parameters
    ----------
    cluster : class
        StarCluster
    fileout : file
        opened file to write data to
    projected : bool
        calculate projected values when possible (default: False)

    Returns
    -------
    None

    History
    -------
    2018 - Written - Webb (UofT)
    """
    units0, origin0 = save_cluster(cluster)

    if cluster.ntot == 0:
        nb = 0
        cluster.mtot = 0.0
        trh = 0.0
        rn = np.zeros(10)
    else:
        cluster.to_pckms()
        cluster.to_centre( do_key_params=True)

        if cluster.nb > 0:
            nb = len(cluster.m2)
        else:
            nb = 0

        trh = half_mass_relaxation_time(cluster, multimass=True, projected=projected)

        if cluster.ntot > 10:
            if cluster.rn == None or (
                origin0 != cluster.origin or units0 != cluster.units
            ):
                rn = rlagrange(cluster, nlagrange=10, projected=projected)
        else:
            rn = np.zeros(10)

    fileout.write(
        "%i %i %f %f %f " % (cluster.ntot, nb, cluster.tphys, trh, cluster.mtot)
    )

    for i in range(0, len(rn)):
        fileout.write("%f " % rn[i])

    fileout.write("%f " % cluster.rmpro)

    if len(cluster.logl) > 0:
        fileout.write("%f " % cluster.rhpro)
    else:
        fileout.write("-1. ")

    # Write additional parameters if they have been calculated:
    if cluster.rv is not None:
        fileout.write("%f " % cluster.rv)
    if cluster.rl is not None:
        fileout.write("%f " % cluster.rl)
    if cluster.rt is not None:
        fileout.write("%f " % cluster.rt)

    try:
        fileout.write("%f %f " % (cluster.rvmax, cluster.vmax))
    except:
        pass

    fileout.write("\n")

    return_cluster(cluster, units0, origin0)

def rho_prof_out(cluster, fileout, **kwargs):
    """Write density profile to file

    Parameters
    ----------
    cluster : class
        StarCluster
    fileout : file
        opened file to write data to

    Returns
    -------
    None

    Other Parameters
    ----------------
    kwargs : str
        key word arguments for rho_prof

    History
    -------
    2018 - Written - Webb (UofT)
    """

    fileout.write("%f " % (cluster.tphys))
    
    rprof, pprof, nprof=rho_prof(cluster,**kwargs)

    for r in rprof:
        fileout.write("%f " % (r))

    for p in pprof:
        fileout.write("%f " % (p))

    for n in nprof:
        fileout.write("%f " % (n))

    fileout.write("\n")


def alpha_prof_out(
    cluster,
    fileout,
    **kwargs,
):
    """Write alpha profile and delta_alpha to file

    Parameters
    ----------
    cluster : class
        StarCluster
    fileout : file
        opened file to write data to

    Returns
    -------
    None

    Other Parameters
    ----------------
    kwargs : str
        key word arguments for alpha_prof

    History
    -------
    2018 - Written - Webb (UofT)
    """

    m_mean, m_hist, dm, alpha, ealpha, yalpha, eyalpha = mass_function(
        cluster,
        **kwargs,
    )
    lrprofn, aprof, dalpha, edalpha, ydalpha, eydalpha = alpha_prof(
        cluster,
        **kwargs,
    )

    fileout.write("%f %f %f %f %f " % (cluster.tphys, alpha, ealpha, yalpha, eyalpha))
    for i in range(0, len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0, len(dm)):
        fileout.write("%f " % dm[i])
    for i in range(0, len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0, len(aprof)):
        fileout.write("%f " % aprof[i])

    fileout.write("%f %f %f %f\n" % (dalpha, edalpha, ydalpha, eydalpha))


def sigv_prof_out(cluster, fileout,**kwargs,
):
    """Write velocity dispersion profile to file

    Parameters
    ----------
    cluster : class
        StarCluster
    fileout : file
        opened file to write data to

    Returns
    -------
    None

    Other Parameters
    ----------------
    kwargs : str
        key word arguments for sigv_prof

    History
    -------
    2018 - Written - Webb (UofT)
    """
    fileout.write("%f %f " % (cluster.tphys, cluster.mtot))

    lrprofn, sigvprof=sigv_prof(cluster, **kwargs)

    for lr in lrprofn:
        fileout.write("%f " % lr)
    for sig in sigvprof:
        fileout.write("%f " % sig)

    fileout.write("%f\n" % (cluster.rm))


def eta_prof_out(
    cluster,
    fileout,
    **kwargs,
):
    """Write eta profile to file

    Parameters
    ----------
    cluster : class
        StarCluster
    fileout : file
        opened file to write data to

    Returns
    -------
    None

    Other Parameters
    ----------------
    kwargs : str
        key word arguments for eta_function and eta_prof

    History
    -------
    2018 - Written - Webb (UofT)
    """

    m_mean, sigvm, eta, eeta, yeta, eyeta = eta_function(
        cluster,
        **kwargs,
    )
    lrprofn, eprof, deta, edeta, ydeta, eydeta = eta_prof(
        cluster,
        **kwargs,
    )

    fileout.write("%f %f %f %f %f " % (cluster.tphys, eta, eeta, yeta, eyeta))
    for i in range(0, len(m_mean)):
        fileout.write("%f " % m_mean[i])
    for i in range(0, len(sigvm)):
        fileout.write("%f " % sigvm[i])
    for i in range(0, len(lrprofn)):
        fileout.write("%f " % lrprofn[i])
    for i in range(0, len(eprof)):
        fileout.write("%f " % eprof[i])

    fileout.write("%f %f %f %f\n" % (deta, edeta, ydeta, eydeta))

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
                np.olumn_stack(
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
    units0, origin0 = save_cluster(cluster)
    cluster.to_centre()

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

    return_cluster(cluster, units0, origin0)

    return 0


def gyrout(cluster, filename="init.nemo.dat"):
    """Output a snapshot in gyrfalcon/NEMO format

    - Note that units are converted to Walter Dehnen units (WDunits in GYRFALCON)

    Parameters
    ----------
    cluster : class
        StarCluster
    filename : str
        name of file to be written to (default: 'init.nemo.dat')
   
    Returns
    -------
    None

    History
    -------
    2019 - Written - Webb (UofT)

    """
    vcon = 220.0 / bovy_conversion.velocity_in_kpcGyr(220.0, 8.0)
    mcon = 222288.4543021174

    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_kpckms()

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

    return_cluster(cluster, units0, origin0)

    return 0
