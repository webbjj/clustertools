import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from .orbit import _get_cluster_orbit

import os, struct

def _get_snapshot(
    filename=None,
    tphys=0.0,
    ctype="snapshot",
    col_names=["m", "x", "y", "z", "vx", "vy", "vz"],
    col_nums=[0, 1, 2, 3, 4, 5, 6],
    units=None,
    origin=None,
    ofile=None,
    advance=False,
    **kwargs
):
    """Load a generic snapshot where column names and numbers can be manually assigned

    Parameters
    ----------
    ctype : str
        code used to generate data (nbody6/nbody6se/nemo/gyrfalcon/...)
    filename : str
        name of file
    col_names : str
        names corresponding to mass, position, and velocity
    col_nums : int
        column numbers corresponding to each column name
    units : str
        units of input data (default: kpckms)
    origin : str
        origin of input data (default: cluster)
    ofile : file
        opened file containing orbital information
    advance : bool
        is this a snapshot that has been advanced to from initial  load_cluster? (default: False)

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

    col_names = np.array(col_names)
    col_nums = np.array(col_nums)

    nsnap = int(kwargs.get("nsnap", "0"))
    nzfill = int(kwargs.get("nzfill", 5))
    delimiter = kwargs.get("delimiter", None)
    wdir = kwargs.get("wdir", "./")
    snapdir = kwargs.get("snapdir", "snaps/")
    snapbase = kwargs.get("snapbase", "")
    snapend = kwargs.get("snapend", ".dat")
    skiprows = kwargs.get("skiprows", 0)

    if filename != None:
        if os.path.isfile("%s%s%s" % (wdir, snapdir, filename)):
            data = np.loadtxt(
                "%s%s%s" % (wdir, snapdir, filename),
                delimiter=delimiter,
                skiprows=skiprows,
            )
        elif os.path.isfile("%s%s" % (wdir, filename)):
            data = np.loadtxt(
                "%s%s" % (wdir, filename), delimiter=delimiter, skiprows=skiprows
            )
        else:
            print("NO FILE FOUND: %s, %s, %s" % (wdir, snapdir, filename))
            cluster = StarCluster( 0., ctype=ctype, **kwargs)
            return cluster
    elif os.path.isfile(
        "%s%s%s%s%s" % (wdir, snapdir, snapbase, str(nsnap).zfill(nzfill), snapend)
    ):
        filename = "%s%s%s%s%s" % (
            wdir,
            snapdir,
            snapbase,
            str(nsnap).zfill(nzfill),
            snapend,
        )
        data = np.loadtxt(
            (
                "%s%s%s%s%s"
                % (wdir, snapdir, snapbase, str(nsnap).zfill(nzfill), snapend)
            ),
            delimiter=delimiter,
            skiprows=skiprows,
        )
    elif os.path.isfile(
        "%s%s%s%s" % (wdir, snapbase, str(nsnap).zfill(nzfill), snapend)
    ):
        filename = "%s%s%s%s" % (wdir, snapbase, str(nsnap).zfill(nzfill), snapend)
        data = np.loadtxt(
            ("%s%s%s%s" % (wdir, snapbase, str(nsnap).zfill(nzfill), snapend)),
            delimiter=delimiter,
            skiprows=skiprows,
        )
    else:
        print(
            "NO FILE FOUND - %s%s%s%s%s"
            % (wdir, snapdir, snapbase, str(nsnap).zfill(nzfill), snapend)
        )
        filename = "%s%s%s%s%s" % (
            wdir,
            snapdir,
            snapbase,
            str(nsnap).zfill(nzfill),
            snapend,
        )
        cluster = StarCluster( 0., ctype=ctype, **kwargs)
        return cluster

    if "m" in col_names:
        mindx = np.argwhere(col_names == "m")[0][0]
        m = data[:, col_nums[mindx]]
    else:
        m=1.

    if "x" in col_names:
        xindx = np.argwhere(col_names == "x")[0][0]
        x = data[:, col_nums[xindx]]
    else:
        x=0.

    if "y" in col_names:
        yindx = np.argwhere(col_names == "y")[0][0]
        y = data[:, col_nums[yindx]]
    else:
        y=0.

    if "z" in col_names:
        zindx = np.argwhere(col_names == "z")[0][0]
        z = data[:, col_nums[zindx]]
    else:
        z=0.

    if "vx" in col_names:
        vxindx = np.argwhere(col_names == "vx")[0][0]
        vx = data[:, col_nums[vxindx]]
    else:
        vx=0.

    if "vy" in col_names:
        vyindx = np.argwhere(col_names == "vy")[0][0]
        vy = data[:, col_nums[vyindx]]
    else:
        vy=0.

    if "vz" in col_names:
        vzindx = np.argwhere(col_names == "vz")[0][0]
        vz = data[:, col_nums[vzindx]]
    else:
        vz=0.

    if "id" in col_names:
        idindx = np.argwhere(col_names == "id")[0][0]
        i_d = data[:, col_nums[idindx]]
        idindx=True
    else:
        idindx=False

    if "kw" in col_names:
        kwindx = np.argwhere(col_names == "kw")[0][0]
        kw = data[:, col_nums[kwindx]]
        kwindx=True
    else:
        kwindx=False

    cluster = StarCluster(
        tphys, units=units, origin=origin, ctype=ctype, **kwargs
    )

    if idindx:
        cluster.add_stars(x, y, z, vx, vy, vz, m, i_d,sortstars=False)
    else:
        cluster.add_stars(x, y, z, vx, vy, vz, m,sortstars=False)

    if kwindx: 
        cluster.kw = kw

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
            cluster.find_centre()
            sortstars=kwargs.get("sortstars", True)
            cluster.analyze(sortstars=sortstars)

        if ofile != None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    return cluster

