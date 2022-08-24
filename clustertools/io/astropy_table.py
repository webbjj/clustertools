import numpy as np

from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from .orbit import _get_cluster_orbit

def _get_astropy_table(
    table,
    column_mapper = None,
    units  = None,
    origin = None,
    ofile = None,
    advance = False,
    verbose = False,
    **kwargs,
):
    """Convert astropy table to a StarCluster instance

    - :class:`~astropy.table.Table`.

    Parameters
    ----------
    table : `~astropy.table.Table` instance
    column_mapper: dict, optional
        Map the ``StarCluster.add_stars`` input to column names in `table`
        If not None,
        the mandatory keys are: "x", "y", "z", "vx", "vy", "vz", and
        the optional keys are "m", "id".

        If None, then performs a basic search of the table column names,
        lowercasing all names for fuzzy matching. The search options
        are different if `units` is "radec" and `origin` is "sky".
        In this case, the search parameters are, ordered by preference and
        lowercased:

            - x : "right ascension" or "ra"
            - y : "declination" or "dec"
            - z : "distance" or "dist" or "dist"
            - vx : "pm_ra_cosdec" or "pm_ra" or "pmra"
            - vy : "pm_dec" or "pmdec"
            - vz : "radial_velocity" or "rvel" or "v_los" or "vlos"
    units : str
        units of input data (default: None)
    origin : str
        origin of input data (default: None)
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

    Raises
    ------
    ValueError
        If `table` missing mandatory argument and
        cannot be found with `column_mapper`
    KeyError
        If missing a mandatory key in `column_mapper`.
    History
    -------
    2020 - Written - Starkman (UofT)
    """

    cm = column_mapper or {}  # None -> {}

    if column_mapper is None:
        # lower-case colum names
        colnames = [n.lower() for n in table.colnames]

        def _helper(*vs: str, to: str, optional=False):
            for v in vs:  # equivalent to if/elif series
                if v in colnames:
                    cm[to] = table.colnames[colnames.index(v)]
                    return table[cm[to]]

            if not optional:
                raise ValueError(
                    (
                        f"Table missing input {to} "
                        f"with searched column names {vs}."
                    )
                )

        # /def

        ID = _helper("id", to="id", optional=True)
        m = _helper("m", "mass", to="m", optional=True)

        if units == "radec" and origin == "sky":  # TOOD lowercase
            # positions
            x = _helper("right ascension", "ra", to="x", optional=False)
            y = _helper("declination", "dec", to="y", optional=False)
            z = _helper("distance", "dist", "d", to="z", optional=False)
            # velocities
            vx = _helper(
                "pm_ra_cosdec", "pm_ra", "pmra", to="vx", optional=False
            )
            vy = _helper("pm_dec", "pmdec", to="vy", optional=False)
            vz = _helper(
                "radial_velocity",
                "rvel",
                "v_los",
                "vlos",
                to="vz",
                optional=False,
            )

        else:
            # positions
            x = _helper("x", to="x", optional=False)
            y = _helper("y", to="y", optional=False)
            z = _helper("z", to="z", optional=False)
            # velocities
            vx = _helper("v_x", "vx", to="vx", optional=False)
            vy = _helper("v_x", "vy", to="vy", optional=False)
            vz = _helper("v_x", "vz", to="vz", optional=False)

    else:  # column_mapper not None
        x = table[cm.pop("x")]
        y = table[cm.pop("y")]
        z = table[cm.pop("z")]
        vx = table[cm.pop("vx")]
        vy = table[cm.pop("xy")]
        vz = table[cm.pop("vz")]

        m = table[cm.pop("m")] if "m" in cm else None
        ID = table[cm.pop("id")] if "id" in cm else None

    cluster = StarCluster(
        0., units=units, origin=origin, ctype='astropy_table', **kwargs
    )

    cluster.add_stars(np.array(x), np.array(y), np.array(z), np.array(vx), np.array(vy), np.array(vz), np.array(m), np.array(ID), sortstars=False)

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

    if verbose:
        print(cm)

    return cluster
