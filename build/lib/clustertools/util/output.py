# Output files containing key cluster properties
# Only functions and profiles should be called here

import numpy as np
from galpy.util import bovy_conversion


from .coordinates import sky_coords
from ..main.cluster import *
from ..main.functions import *
from ..main.profiles import *
from ..main.operations import *


def extrct_out(cluster, fileout, projected=False):
    """
    NAME:

       extrct_out

    PURPOSE:

       Extrct key cluster properties and write to file
       --> N, Nbinary, trh, mtot, rlagrange_1-10, rmpro, rhpro, rv, rl, rt, rvmax, vmax

    INPUT:

       cluster - a StarCluster-class object

       fileout - opened file to write data to

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    units0, origin0 = save_cluster(cluster)

    if cluster.ntot == 0:
        nb = 0
        cluster.mtot = 0.0
        trh = 0.0
        rn = np.zeros(10)
    else:
        cluster.to_realpc()
        cluster.to_centre(do_order=True, do_key_params=True)

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


def snapout(cluster, filename, energies=False, observations=False):
    """
    NAME:

       snapout

    PURPOSE:

       Output a snapshot in clustertools format

    INPUT:

       cluster - a StarCluster-class object

       filename - name of file to be written to

       energies - include energies in output (Default: False)

       observations - include sky values of stars (Default: False) 

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    if observations:

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
    """
    NAME:

       fortout

    PURPOSE:

       Output a snapshot in NBODY6 fort.10 format

    INPUT:

       cluster - a StarCluster-class object

       filename - name of file to be written to (Default: 'fort.10')

       reset_nbody_scale - reset nbody scaling parameters (Default: False)

       reset_nbody_mass - reset nbody mass scaling parameter (Default: False)

       reset_nbody_radii - reset nbody radii scaling parameter (Default: False)

    OUTPUT:

       None

    HISTORY:

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
    """
    NAME:

       gyrout

    PURPOSE:

       Output a snapshot in gyrfalcon/NEMO format

    INPUT:

       cluster - a StarCluster-class object

       filename - name of file to be written to (Default: 'init.nemo.dat')

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)

    """
    vcon = 220.0 / bovy_conversion.velocity_in_kpcGyr(220.0, 8.0)
    mcon = 222288.4543021174

    units0, origin0 = save_cluster(cluster)
    cluster.to_galaxy()
    cluster.to_realkpc()

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
