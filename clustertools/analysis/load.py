""" Read in cluster from Nbody simulations or generate an Nbody cluster

"""

__author__ = "Jeremy J Webb"


__all__ = [
    "load_cluster",
    "advance_cluster"
]

import numpy as np
from galpy.util import bovy_conversion
import os
from .cluster import StarCluster
from .operations import *
from .orbit import initialize_orbit

# Try Importing AMUSE. Only necessary for get_amuse_particles
try:
    import amuse.units.units as u
except:
    pass


def load_cluster(
    ctype="snapshot",
    units="pckms",
    origin="cluster",
    ofile=None,
    orbit=None,
    filename=None,
    **kwargs
):
    """
    NAME:

       load_cluster

    PURPOSE:

       Load a StarCluster snapshot from a generic code output
       --> Note user can easily create their own def my_code() for new code output types

    INPUT:

       ctype - Type of file being loaded (Currently supports nbody6, nbody6se, gyrfalcon, snaptrim,snapauto, clustertools, snapshot)

       units - units of input data (default: kpckms)

       origin - origin of input data (default: cluster)

       ofile - an already opened file containing orbit information (default: None)

       orbit - a galpy orbit to be used for the StarCluster's orbital information (default: None)

       filename - name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)

    KWARGS:

        ofilename - orbit filename if ofile is not given

        ounits - units of orbital information (else assumed equal to StarCluster.units)

        nsnap - if a specific snapshot is to be read in instead of starting from zero
    
        nzfill - value for zfill when reading and writing snapshots (Default: 5)
    
        delimiter - choice of delimiter when reading ascii/csv files (Default: ',')
    
        wdir - working directory of snapshots if not current directory

        intialize - initialize a galpy orbit after reading in orbital information (default: False)

        kwfile - open file containing stellar evolution type (kw) of individual stars (used if snapshot file does not contain kw and kw is needed)

        projected - calculate projected values as well as 3D values (Default: True)

        do_key_params - calculate key parameters

        do_rorder - sort stars in order from closes to the origin to the farthest


    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    wdir = kwargs.get("wdir", "./")
    initialize = kwargs.get("initialize", False)

    if "ofilename" in kwargs and ofile == None:
        ofile = open(wdir + kwargs.get("ofilename"), "r")

    if ctype == "nbody6se":
        # When stellar evolution is turned on, read in fort.82 and fort.83 and if possible gc_orbit.dat
        fort82 = open("%sfort.82" % wdir, "r")
        fort83 = open("%sfort.83" % wdir, "r")
        cluster = get_nbody6_jarrod(
            fort82, fort83, ofile=ofile, advance=False, **kwargs
        )
    elif ctype == "nbody6":
        # With stellar evolution turned off, read in OUT9 and OUT34. Orbit data already in OUT34
        if os.path.isfile("%sOUT9" % wdir):
            out9 = open("%sOUT9" % wdir, "r")
        else:
            out9 = None
        out34 = open("%sOUT34" % wdir, "r")
        cluster = get_nbody6_out(out9, out34, advance=False, **kwargs)
    elif ctype == "snapauto":
        # Read in snapshot produced from snapauto.f which reads binary files from either NBODY6 or NBODY6++
        cluster = get_nbody6_snapauto(
            filename=filename,
            units=units,
            origin=origin,
            ofile=ofile,
            advance=False,
            **kwargs
        )
    elif ctype == "gyrfalcon":
        # Read in snapshot from gyrfalcon.
        filein = open(wdir + filename, "r")
        cluster = get_gyrfalcon(filein, "WDunits", "galaxy", advance=False, **kwargs)
    elif ctype == "snaptrim":
        # Read in snaptrim snapshot from gyrfalcon.
        cluster = get_snaptrim(
            filename=filename, units="WDunits", origin="galaxy", advance=False, **kwargs
        )
    elif ctype == "clustertools":
        # Read in standard clustertools snapshot
        cluster = get_clustertools_snapshot(
            filename=filename,
            units=units,
            origin=origin,
            ofile=ofile,
            advance=False,
            **kwargs
        )
    elif ctype == "snapshot":
        # Read in standard generic snapshot
        col_names = kwargs.pop("col_names", ["m", "x", "y", "z", "vx", "vy", "vz"])
        col_nums = kwargs.pop("col_nums", [0, 1, 2, 3, 4, 5, 6])
        cluster = get_snapshot(
            filename=filename,
            col_names=col_names,
            col_nums=col_nums,
            units=units,
            origin=origin,
            ofile=ofile,
            advance=False,
            **kwargs
        )
    elif ctype == "mycode":
        # Read in new cluster type
        cluster = get_mycode()
    else:
        print("NO CTYPE GIVEN")
        return 0

    kwfile = kwargs.get("kwfile", None)
    if kwfile != None:
        cluster.kw[cluster.id - 1] = get_kwtype(cluster, kwfile)

    # Add galpy orbit if given
    if orbit != None:
        cluster.orbit = orbit
        if cluster.units == "pckms":
            t = (cluster.tphys / 1000.0) / bovy_conversion.time_in_Gyr(ro=8.0, vo=220.0)
            cluster.add_orbit(
                orbit.x(t) * 1000.0,
                orbit.y(t) * 1000.0,
                orbit.z(t) * 1000.0,
                orbit.vx(t),
                orbit.vy(t),
                orbit.vz(t),
            )
        if cluster.units == "kpckms":
            cluster.add_orbit(
                orbit.x(t),
                orbit.y(t),
                orbit.z(t),
                orbit.vx(t),
                orbit.vy(t),
                orbit.vz(t),
            )
        elif cluster.units == "nbody":
            t = (cluster.tphys * cluster.tstar / 1000.0) / bovy_conversion.time_in_Gyr(
                ro=8.0, vo=220.0
            )
            cluster.add_orbit(
                orbit.x(t) * 1000.0 / cluster.rbar,
                orbit.y(t) * 1000.0 / cluster.rbar,
                orbit.z(t) * 1000.0 / cluster.rbar,
                orbit.vx(t) / cluster.vstar,
                orbit.vy(t) / cluster.vstar,
                orbit.vz(t) / cluster.vstar,
            )
        elif cluster.units == "galpy":
            t = cluster.tphys
            cluster.add_orbit(
                orbit.x(t) / 8.0,
                orbit.y(t) / 8.0,
                orbit.z(t) / 8.0,
                orbit.vx(t) / 220.0,
                orbit.vy(t) / 220.0,
                orbit.vz(t) / 220.0,
            )

        units0, origin0 = save_cluster(cluster)

        cluster.to_cluster()
        cluster.find_centre()

        return_cluster(cluster, units0, origin0)
    elif initialize:
        initialize_orbit(cluster)

    # cluster.key_params()

    return cluster


def advance_cluster(cluster, ofile=None, orbit=None, filename=None, **kwargs):
    """
    NAME:

       advance_cluster

    PURPOSE:

       Advance a loaded StarCluster snapshot to the next timestep
       --> ofile or orbit need to be provded again, same as load_cluster.
       --> Be sure that advance is set to True so next line of orbit file is read in

    INPUT:

       filename - name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)

    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)

    """
    advance_kwargs = get_advanced_kwargs(cluster, **kwargs)

    if "kwfile" in kwargs:
        kw0 = cluster.kw

    # Continue reading in cluster opened in get_cluster()
    if cluster.ctype == "nbody6se":
        cluster = get_nbody6_jarrod(
            cluster.bfile, cluster.sfile, ofile=ofile, advance=True, **advance_kwargs
        )
    elif cluster.ctype == "nbody6":
        cluster = get_nbody6_out(
            cluster.bfile, cluster.sfile, advance=True, **advance_kwargs
        )
    elif cluster.ctype == "snapauto":
        cluster = get_nbody6_snapauto(
            filename=filename,
            units=cluster.units,
            origin=cluster.origin,
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )
    elif cluster.ctype == "gyrfalcon":

        cluster = get_gyrfalcon(
            cluster.sfile,
            units="WDunits",
            origin="galaxy",
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )
    elif cluster.ctype == "snaptrim":
        cluster = get_snaptrim(
            filename=filename,
            units="WDunits",
            origin="galaxy",
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )
    elif cluster.ctype == "clustertools":
        cluster = get_clustertools_snapshot(
            filename=filename,
            units=cluster.units,
            origin=cluster.origin,
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )
    elif cluster.ctype == "snapshot":
        col_names = kwargs.pop("col_names", ["m", "x", "y", "z", "vx", "vy", "vz"])
        col_nums = kwargs.pop("col_nums", [0, 1, 2, 3, 4, 5, 6])

        cluster = get_snapshot(
            filename=filename,
            col_names=col_names,
            col_nums=col_nums,
            units=cluster.units,
            origin=cluster.origin,
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )
    elif cluster.ctype == "mycode":
        cluster = get_mycode()
    else:
        cluster = StarCuster(ctype=cluster.ctype)

    # Check for restart
    if cluster.ntot == 0.0:
        print('NTOT = 0',cluster.wdir,advance_kwargs.get('wdir','./'))
        try:
            wdir = cluster.wdir + "cont/"
        except:
            print("WDIR NOT SET")
            wdir = "./cont/"

        try:
            ofilename = ofile.name
        except:
            print("OFILE NOT SET")
            ofile = None

        if os.path.exists(wdir):
            old_wdir=advance_kwargs.pop('wdir')
            cluster = load_cluster(
                ctype=cluster.ctype, ofile=ofile, wdir=wdir, **advance_kwargs
            )

    if cluster.ntot != 0.0:
        if "kwfile" in kwargs:
            cluster.kw[cluster.id - 1] = kw0

        # Add galpy orbit if given
        if orbit != None:
            cluster.orbit - orbit
            if cluster.units == "pckms" or cluster.units == "kpckms":
                t = (cluster.tphys / 1000.0) / bovy_conversion.time_in_Gyr(
                    ro=8.0, vo=220.0
                )
            elif cluster.units == "nbody":
                t = (
                    cluster.tphys * cluster.tstar / 1000.0
                ) / bovy_conversion.time_in_Gyr(ro=8.0, vo=220.0)
            elif cluster.units == "galpy":
                t = cluster.tphys

            cluster.add_orbit(
                orbit.x(t),
                orbit.y(t),
                orbit.z(t),
                orbit.vx(t),
                orbit.vy(t),
                orbit.vz(t),
            )

        cluster.key_params()

    return cluster


def get_advanced_kwargs(cluster, **kwargs):
    """
    NAME:

       get_advanced_kwargs

    PURPOSE:

       get **kwargs from current cluster before advancing

    INPUT:

       cluster - StarCluster instance

    KWARGS:

        same a load_cluster

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    kwfile = kwargs.get("kwfile", None)

    nsnap = np.maximum(int(kwargs.get("nsnap", 0)), cluster.nsnap) + 1
    delimiter = kwargs.get("delimiter", cluster.delimiter)
    wdir = kwargs.get("wdir", cluster.wdir)
    nzfill = int(kwargs.get("nzfill", cluster.nzfill))
    snapbase = kwargs.get("snapbase", cluster.snapbase)
    snapend = kwargs.get("snapend", cluster.snapend)
    snapdir = kwargs.get("snapdir", cluster.snapdir)
    skiprows = kwargs.get("skiprows", cluster.skiprows)

    projected = kwargs.get("projected", cluster.projected)

    return {
        "kwfile": kwfile,
        "nsnap": nsnap,
        "delimiter": delimiter,
        "wdir": wdir,
        "nzfill": nzfill,
        "snapbase": snapbase,
        "snapend": snapend,
        "snapdir": snapdir,
        "skiprows": skiprows,
        "projected": projected,
    }  # ,"sfile":sfile,"bfile":bfile}


def get_cluster_orbit(cluster, ofile, advance=False, **kwargs):
    """
    NAME:

       get_cluster_orbit

    PURPOSE:

       Read in cluster oribit from an ascii file and apply it to StarCluster instance
        -->Columns assumed to be time,x,y,z,vx,vy,vz.
        -->As I have hardcoded for gc_orbit.dat, if a certain filename is commonly used
        with specific columns, an if statement can easily be added.
            --> Note that gc_orbit.dat comes from running 'grep 'CLUSTER ORBIT' logfile > gc_orbit.dat' on 
            the standard Nbody6 logfile

    INPUT:

       cluster - StarCluster instance

       ofile - opened orbit file

    KWARGS:

        advance - Is this a continuation from a previous timestep, in which case read next line (default: False)
        nsnap - if nsnap is provided, read line # nsnap from the orbit file
        ounits - if you units are not the same as StarCluster units, provide them and they will be converted

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    nsnap = int(kwargs.get("nsnap", cluster.nsnap))
    ounits = kwargs.get("ounits", None)

    # Read in orbital information from orbit
    if nsnap != 0 and not advance:
        for i in range(0, int(nsnap) + 1):
            data = ofile.readline().split()
    else:
        data = ofile.readline().split()

    if "gc_orbit.dat" in ofile.name:
        # Saved orbit from doing a grep of NBODY6 or NBODY6++ logfile

        if len(data) == 18:
            xgc = float(data[9])
            ygc = float(data[10])
            zgc = float(data[11])
            vxgc = float(data[12])
            vygc = float(data[13])
            vzgc = float(data[14])
        else:
            xgc = float(data[8])
            ygc = float(data[9])
            zgc = float(data[10])
            vxgc = float(data[11])
            vygc = float(data[12])
            vzgc = float(data[13])
    else:
        tphys = float(data[0])
        xgc = float(data[1])
        ygc = float(data[2])
        zgc = float(data[3])
        vxgc = float(data[4])
        vygc = float(data[5])
        vzgc = float(data[6])

        if cluster.tphys == 0.0:
            cluster.tphys = tphys

    if ounits == None and "gc_orbit.dat" in ofile.name:
        ounits = "kpckms"

    cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc, ounits)

    return


def get_kwtype(cluster, kwfile):
    """
    NAME:

       get_kwtype

    PURPOSE:

       Get stellar evolution type of the star (kw in Nbody6 syntax) if information is in a separate file
       --> Columns assumed to be ID, KW
       --> see def kwtypes() in cluster.py for what different kw types mean

    INPUT:

       cluster - StarCluster instance

       kwfile - opened file containing kw type

    KWARGS:

        advance - Is this a continuation from a previous timestep, in which case read next line (default: False)
        nsnap - if nsnap is provided, read line # nsnap from the orbit file
        ounits - if you units are not the same as StarCluster units, provide them and they will be converted

    OUTPUT:

       None

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    kw = np.loadtxt(kwfile)[:, 1]
    indx = cluster.id - 1
    kw0 = kw[indx]
    return kw0


# Get StarCluster from Gyrfalcon output
def get_gyrfalcon(
    filein, units="WDunits", origin="galaxy", ofile=None, advance=False, **kwargs
):
    """
    NAME:

       get_gyrfalcon

    PURPOSE:

       Extract a single snapshot from an ascii file output from a gyrfalcon simulation

    INPUT:

       filein - opened gyrfalcon file

       units - units of data (default:'WDunits')

       ofile - opened file containing orbital information

       advance - is this a snapshot that has been advanced to from initial load_cluster?


    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)

    """

    if units == "WDunits":
        vcon = 220.0 / bovy_conversion.velocity_in_kpcGyr(220.0, 8.0)
        mcon = 222288.4543021174
        units = "kpckms"
    else:
        vcon = 1.0
        mcon = 1.0

    # Default **kwargs
    skiprows = kwargs.pop("skiprows", 13)

    i_d = []
    m = []
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []

    over_head = False
    ntot = 0
    tphys = 0.0

    for j in range(0, skiprows):
        data = filein.readline().split()
        if "#" not in data:
            over_head = True
            print("OVER HEAD")
            break
        if len(data) == 0:
            print("END OF FILE")
            return StarCluster(0.0,ctype="gyrfalcon",**kwargs)
        if any("Ntot" in dat for dat in data):
            sntot = data[2]
            ntot = int(sntot[:-1])
        if any("time" in dat for dat in data):
            tphys = float(data[2]) * 1000.0

    cluster = StarCluster(
        ntot,
        tphys,
        units=units,
        origin=origin,
        ctype="gyrfalcon",
        sfile=filein,
        bfile=None,
        skiprows=skiprows,
        **kwargs
    )

    for j in range(ntot):
        if over_head:
            over_head = False
        else:
            data = filein.readline().split()
        if "#" in data:
            break

        i_d.append(j + 1)
        m.append(float(data[0]) * mcon)
        x.append(float(data[1]))
        y.append(float(data[2]))
        z.append(float(data[3]))
        vx.append(float(data[4]) * vcon)
        vy.append(float(data[5]) * vcon)
        vz.append(float(data[6]) * vcon)

    if ntot > 0:

        cluster.add_stars(x, y, z, vx, vy, vz, m, i_d)

        if ofile == None:
            cluster.find_centre()
        else:
            get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

        if kwargs.get("do_key_params", True):
            do_order=kwargs.get("do_key_params", True)
            cluster.to_cluster()
            cluster.find_centre()
            cluster.to_centre(do_key_params=True, do_order=do_order)
            cluster.to_galaxy()

    return cluster


def get_snaptrim(
    filename=None, units="WDunits", origin="galaxy", ofile=None, advance=False, **kwargs
):
    """
    NAME:

       get_snaptrim

    PURPOSE:

       Load a gyrfalcon snapshot as produced by snaptrim

    INPUT:

       filename = name of file

       units - units of input data (default: WDunits)

       origin - origin of input data (default: galaxy)

       advance - is this a snapshot that has been advanced to from initial load_cluster?

       
    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    # Default **kwargs

    nzfill = int(kwargs.pop("nzfill", 1))
    skiprows = kwargs.pop("skiprows", 13)
    delimiter = kwargs.pop("delimiter", None)
    nsnap = int(kwargs.get("nsnap", "0"))
    wdir = kwargs.get("wdir", "./")
    snapdir = kwargs.get("snapdir", "snaps/")
    snapbase = kwargs.get("snapbase", "")
    snapend = kwargs.get("snapend", ".dat")

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
            cluster = StarCluster( 0.0, ctype='snaptrim', **kwargs)
            print(cluster.ntot)
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
    elif os.path.isfile(
        "%s%s%s%s" % (wdir, snapbase, str(nsnap).zfill(nzfill), snapend)
    ):
        filename = "%s%s%s%s" % (wdir, snapbase, str(nsnap).zfill(nzfill), snapend)
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
        cluster = StarCluster( 0.,ctype='snaptrim', sfile=filename, **kwargs)
        print(cluster.ntot)
        return cluster

    ntot = 0
    tphys = 0.0

    filein = open(filename, "r")

    for j in range(0, skiprows):
        data = filein.readline().split()
        if "#" not in data:
            print("OVER HEAD")
            break
        if len(data) == 0:
            print("END OF FILE")
            return StarCluster( 0.0, ctype="snaptrim",**kwargs)
        if any("Ntot" in dat for dat in data):
            sntot = data[2]
            ntot = int(sntot[:-1])
        if any("time" in dat for dat in data):
            tphys = float(data[2]) * 1000.0

    filein.close()

    cluster = get_snapshot(
        filename=filename,
        tphys=tphys,
        col_names=["m", "x", "y", "z", "vx", "vy", "vz"],
        col_nums=[0, 1, 2, 3, 4, 5, 6],
        units=units,
        origin=origin,
        ofile=ofile,
        advance=advance,
        ctype="snaptrim",
        nzfill=nzfill,
        skiprows=skiprows,
        delimiter=delimiter,
        **kwargs
    )

    return cluster


def get_nbody6_jarrod(fort82, fort83, ofile=None, advance=False, **kwargs):
    """
    NAME:

       get_nbody6_jarrod

    PURPOSE:

       Extract a single snapshot from the fort.82 and fort.83 files output by Jarrod Hurleys version of Nbody6
       --> Difference between Jarrod's output and the public version is found in /Node/hrplot.f of Nbody6
       --> Called for Nbody6 simulations with stellar evolution

    INPUT:

       fort82 - opened fort.82 file containing binary star information

       fort83 - opened fort.83 file containing single star information

       ofile - opened file containing orbital information

       advance - is this a snapshot that has been advanced to from initial load_cluster?

    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    line1 = fort83.readline().split()
    if len(line1) == 0:
        print("END OF FILE")
        return StarCluster( 0.0,ctype='nbody6se',**kwargs)

    line2 = fort83.readline().split()
    line3 = fort83.readline().split()
    line1b = fort82.readline().split()

    ns = int(line1[0])
    tphys = float(line1[1])
    nc = int(line2[0])
    rc = max(float(line2[1]), 0.01)
    rbar = float(line2[2])
    rtide = float(line2[3])
    xc = float(line2[4])
    yc = float(line2[5])
    zc = float(line2[6])
    zmbar = float(line3[0])
    vstar = 0.06557 * np.sqrt(zmbar / rbar)
    rscale = float(line3[2])
    nb = int(line1b[0])
    ntot = ns + nb

    nsbnd = 0
    nbbnd = 0
    nbnd = 0

    i_d = []
    id1 = []
    id2 = []
    kw = []
    kw1 = []
    kw2 = []
    kcm = []
    ecc = []
    pb = []
    semi = []
    m1 = []
    m2 = []
    m = []
    logl1 = []
    logl2 = []
    logl = []
    logr1 = []
    logr2 = []
    logr = []
    x = []
    y = []
    z = []
    rxy = []
    r = []
    vx = []
    vy = []
    vz = []
    v = []
    ep = []
    ep1 = []
    ep2 = []
    ospin = []
    ospin1 = []
    ospin2 = []
    kin = []
    pot = []
    etot = []

    data = fort82.readline().split()

    while int(data[0]) > 0 and len(data) > 0:
        id1.append(int(data[0]))
        id2.append(int(data[1]))
        i_d.append(id1[-1])
        kw1.append(int(data[2]))
        kw2.append(int(data[3]))
        kw.append(max(kw1[-1], kw2[-1]))
        kcm.append(float(data[4]))
        ecc.append(float(data[5]))
        pb.append(10.0**float(data[6]))
        semi.append(10.0**float(data[7]))
        m1.append(float(data[8]) / zmbar)
        m2.append(float(data[9]) / zmbar)
        m.append(m1[-1] + m2[-1])
        logl1.append(float(data[10]))
        logl2.append(float(data[11]))
        logl.append(max(logl1[-1], logl2[-1]))
        logr1.append(float(data[12]))
        logr2.append(float(data[13]))
        logr.append(max(logr1, logr2))
        x.append(float(data[14]))
        y.append(float(data[15]))
        z.append(float(data[16]))
        vx.append(float(data[17]))
        vy.append(float(data[18]))
        vz.append(float(data[19]))

        if "bnd" in fort82.name or "esc" in fort82.name:
            kin.append(float(data[20]))
            pot.append(float(data[21]))
            etot.append(float(data[23]))
        else:
            kin.append(0.0)
            pot.append(0.0)
            etot.append(0.0)
            ep1.append(float(data[20]))
            ep2.append(float(data[21]))
            ospin1.append(float(data[22]))
            ospin2.append(float(data[23]))

            ep.append(ep1[-1])
            ospin.append(ospin1[-1])

        nbbnd += 1
        data = fort82.readline().split()

    data = fort83.readline().split()
    while int(data[0]) > 0 and len(data) > 0:
        i_d.append(int(data[0]))
        kw.append(int(data[1]))
        m.append(float(data[2]) / zmbar)
        logl.append(float(data[3]))
        logr.append(float(data[4]))
        x.append(float(data[5]))
        y.append(float(data[6]))
        z.append(float(data[7]))
        vx.append(float(data[8]))
        vy.append(float(data[9]))
        vz.append(float(data[10]))

        if "bnd" in fort83.name or "esc" in fort83.name:
            kin.append(float(data[11]))
            pot.append(float(data[12]))
            etot.append(float(data[14]))
        else:
            kin.append(0.0)
            pot.append(0.0)
            etot.append(0.0)
            ep.append(float(data[11]))
            ospin.append(float(data[12]))

        nsbnd += 1
        data = fort83.readline().split()

    nbnd = nsbnd + nbbnd

    if nbnd > 0:
        cluster = StarCluster(
            nbnd,
            tphys,
            units="nbody",
            origin="cluster",
            ctype="nbody6se",
            sfile=fort83,
            bfile=fort82,
            **kwargs
        )
        cluster.add_nbody6(
            nc, rc, rbar, rtide, xc, yc, zc, zmbar, vstar, rscale, nsbnd, nbbnd
        )
        cluster.add_stars(x, y, z, vx, vy, vz, m, i_d)
        cluster.add_sse(kw, logl, logr, ep, ospin)
        cluster.add_bse(
            id1,
            id2,
            kw1,
            kw2,
            kcm,
            ecc,
            pb,
            semi,
            m1,
            m2,
            logl1,
            logl2,
            logr1,
            logr2,
            ep1,
            ep2,
            ospin1,
            ospin2,
        )
        cluster.add_energies(kin, pot, etot)

        if ofile != None:
            get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

        if kwargs.get("do_key_params", True):
            do_order=kwargs.get("do_order", True)
            # Estimate centre
            cluster.find_centre()
            cluster.to_centre(do_key_params=True, do_order=do_order)
            cluster.to_cluster()

    else:
        cluster = StarCluster( tphys, ctype="nbody6se", **kwargs)

    return cluster


def get_nbody6_out(out9, out34, advance=False, **kwargs):
    """
    NAME:

       get_nbody6_out34

    PURPOSE:

       Extract a single snapshot from the OUT9 and OUT34 output of Jeremy Webb's version of Nbody6
       --> Difference between Jeremy's output and the public version is found in /Node/output.f of Nbody6
       --> Called for Nbody6 simulations without stellar evolution

    INPUT:

       out9 - opened OUT34 file containing binary star information

       out34 - opened OUT34 file containing single star information

       ofile - opened file containing orbital information

       advance - is this a snapshot that has been advanced to from initial load_cluster?

    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    line1 = out34.readline().split()
    if len(line1) == 0:
        print("END OF FILE")
        return StarCluster( 0.0, ctype='nbody6',**kwargs)

    line2 = out34.readline().split()
    line3 = out34.readline().split()

    ns = int(line1[0])
    tphys = float(line1[1])
    n_p = int(line1[4])
    if len(line1) > 11:
        nb = int(float(line1[11]))
    else:
        nb = 0

    if out9 != None:
        line1b = out9.readline().split()
        line2b = out9.readline().split()
        line3b = out9.readline().split()

        if nb != int(line1b[0]):
            print("ERROR: NUMBER OF BINARIES DO NOT MATCH - ",nb,int(line1b[0]))

    nc = int(line2[0])
    rc = max(float(line2[1]), 0.01)
    rbar = float(line2[2])
    rtide = float(line2[3])
    xc = float(line2[4])
    yc = float(line2[5])
    zc = float(line2[6])
    zmbar = float(line3[0])
    vstar = 0.06557 * np.sqrt(zmbar / rbar)
    rscale = float(line3[2])
    ntot = ns + nb

    # Orbital Properties
    xgc = float(line1[5])
    ygc = float(line1[6])
    zgc = float(line1[7])
    vxgc = float(line1[8])
    vygc = float(line1[9])
    vzgc = float(line1[10])

    nsbnd = 0
    nbbnd = 0
    nbnd = 0

    i_d = []
    kw = []
    m = []
    logl = []
    logr = []
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    kin = []
    pot = []
    etot = []

    if out9 != None:

        yrs = (rbar * 1296000.0 / (2.0 * np.pi)) ** 1.5 / np.sqrt(zmbar)
        days = 365.25 * yrs

        id1 = []
        kw1 = []
        m1 = []
        logl1 = []
        logr1 = []
        id2 = []
        kw2 = []
        m2 = []
        logl2 = []
        logr2 = []
        pb = []
        kcm = []
        ecc = []
        semi = []

        for i in range(0, nb):
            data = out9.readline().split()

            #Ignore massless ghost particles ouput by NBODY6
            if (float(data[4])+float(data[5])) > 0:

                nbbnd += 1

                ecc.append(float(data[1]))
                m1.append(float(data[4]) / zmbar)
                m2.append(float(data[5]) / zmbar)
                pb.append(float(data[6]) / days)
                id1.append(int(float(data[7])))
                id2.append(int(float(data[8])))
                kw1.append(int(data[9]))
                kw2.append(int(data[10]))
                kcm.append(int(data[11]))

                logl1.append(1.0)
                logl2.append(1.0)
                logr1.append(1.0)
                logr2.append(1.0)

                x1 = float(data[12])
                y1 = float(data[13])
                z1 = float(data[14])
                vx1 = float(data[15])
                vy1 = float(data[16])
                vz1 = float(data[17])
                x2 = float(data[18])
                y2 = float(data[19])
                z2 = float(data[20])
                vx2 = float(data[21])
                vy2 = float(data[22])
                vz2 = float(data[23])

                ''' It seems binary COM information is included in OUT34
                x.append((x1 * m1[-1] + x2 * m2[-1]) / (m1[-1] + m2[-1]) + xc)
                y.append((y1 * m1[-1] + y2 * m2[-1]) / (m1[-1] + m2[-1]) + yc)
                z.append((z1 * m1[-1] + z2 * m2[-1]) / (m1[-1] + m2[-1]) + zc)
                vx.append((vx1 * m1[-1] + vx2 * m2[-1]) / (m1[-1] + m2[-1]))
                vy.append((vy1 * m1[-1] + vy2 * m2[-1]) / (m1[-1] + m2[-1]))
                vz.append((vz1 * m1[-1] + vz2 * m2[-1]) / (m1[-1] + m2[-1]))
                m.append(m1[-1] + m2[-1])
                i_d.append(id1[-1])
                kw.append(max(kw1[-1], kw2[-1]))
                logl.append(1.0)
                logr.append(1.0)
                

                r1 = np.sqrt((x1 - x[-1]) ** 2.0 + (y1 - y[-1]) ** 2.0 + (z1 - z[-1]) ** 2.0)
                r2 = np.sqrt((x2 - x[-1]) ** 2.0 + (y2 - y[-1]) ** 2.0 + (z2 - z[-1]) ** 2.0)
                '''
                mb=(m1[-1] + m2[-1])

                semi.append((pb[-1]**2.*mb)**(1./3.))

    data = out34.readline().split()

    while int(float(data[0])) >= -999:
        # IGNORE GHOST PARTICLES
        if float(data[2]) == 0.0:
            ns -= 1
            ntot -= 1
        else:
            i_d.append(int(float(data[0])))
            kw.append(int(data[1]))
            m.append(float(data[2]))
            logl.append(float(data[3]))
            logr.append(float(data[4]))
            x.append(float(data[5]) + xc)
            y.append(float(data[6]) + yc)
            z.append(float(data[7]) + zc)
            vx.append(float(data[8]))
            vy.append(float(data[9]))
            vz.append(float(data[10]))

            if len(data) > 14:
                kin.append(float(data[13]))
                pot.append(float(data[14]))
                etot.append(float(data[15]))
            else:
                kin.append(0.0)
                pot.append(0.0)
                etot.append(0.0)

            nsbnd += 1
        data = out34.readline().split()

        if len(data)==0:
            break

    nbnd = nsbnd + nbbnd

    cluster = StarCluster(
        nsbnd,
        tphys,
        units="nbody",
        origin="cluster",
        ctype="nbody6",
        sfile=out34,
        bfile=out9,
        **kwargs
    )
    cluster.add_nbody6(
        nc, rc, rbar, rtide, xc, yc, zc, zmbar, vstar, rscale, nsbnd, nbbnd, n_p
    )
    # Add back on the centre of mass which has been substracted off by NBODY6
    cluster.add_stars(x, y, z, vx, vy, vz, m, i_d)
    cluster.add_sse(kw, logl, logr, np.zeros(nbnd), np.zeros(nbnd))
    cluster.add_energies(kin, pot, etot)
    if out9 != None:
        cluster.add_bse(
            id1, id2, kw1, kw2, kcm, ecc, pb, semi, m1, m2, logl1, logl2, logr1, logr2
        )

    if kwargs.get("do_key_params", True):
        do_order=kwargs.get("do_key_params", True)
        # Estimate centre of distribution
        cluster.find_centre()
        cluster.to_centre(do_key_params=True, do_order=do_order)
        cluster.to_cluster()

    cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc)

    return cluster


def get_nbody6_out34(out34, advance=False, **kwargs):
    """
    NAME:

       get_nbody6_out34

    PURPOSE:

       Extract a single snapshot from the OUT34 output of Jeremy Webb's version of Nbody6
       --> Difference between Jeremy's output and the public version is found in /Node/output.f of Nbody6
       --> Called for Nbody6 simulations without stellar evolution

    INPUT:

       out34 - opened OUT34 file containing single star information

       advance - is this a snapshot that has been advanced to from initial load_cluster?

    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    line1 = out34.readline().split()
    if len(line1) == 0:
        print("END OF FILE")
        return StarCluster( 0.0, ctype='nbody6',**kwargs)

    line2 = out34.readline().split()
    line3 = out34.readline().split()

    ns = int(line1[0])
    tphys = float(line1[1])
    n_p = int(line1[4])

    nc = int(line2[0])
    rc = max(float(line2[1]), 0.01)
    rbar = float(line2[2])
    rtide = float(line2[3])
    xc = float(line2[4])
    yc = float(line2[5])
    zc = float(line2[6])
    zmbar = float(line3[0])
    vstar = 0.06557 * np.sqrt(zmbar / rbar)
    rscale = float(line3[2])
    nb = 0
    ntot = ns + nb

    # Orbit Properties
    xgc = float(line1[5])
    ygc = float(line1[6])
    zgc = float(line1[7])
    vxgc = float(line1[8])
    vygc = float(line1[9])
    vzgc = float(line1[10])

    nsbnd = 0
    nbbnd = 0
    nbnd = 0

    i_d = []
    kw = []
    m = []
    logl = []
    logr = []
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    kin = []
    pot = []
    etot = []

    data = out34.readline().split()
    while int(float(data[0])) >= -999 and len(data) > 0:
        # IGNORE GHOST PARTICLES
        if float(data[2]) == 0.0:
            ns -= 1
            ntot -= 1
        else:
            i_d.append(int(float(data[0])))
            kw.append(int(data[1]))
            m.append(float(data[2]))
            logl.append(float(data[3]))
            logr.append(float(data[4]))
            x.append(float(data[5]))
            y.append(float(data[6]))
            z.append(float(data[7]))
            vx.append(float(data[8]))
            vy.append(float(data[9]))
            vz.append(float(data[10]))

            if len(data) > 14:
                kin.append(float(data[13]))
                pot.append(float(data[14]))
                etot.append(float(data[15]))
            else:
                kin.append(0.0)
                pot.append(0.0)
                etot.append(0.0)

            nsbnd += 1
        data = out34.readline().split()

    nbnd = nsbnd + nbbnd

    cluster = StarCluster(
        nbnd,
        tphys,
        units="nbody",
        origin="cluster",
        ctype="nbody6",
        sfile=out34,
        bfile=None,
        **kwargs
    )
    cluster.add_nbody6(
        nc, rc, rbar, rtide, xc, yc, zc, zmbar, vstar, rscale, nsbnd, nbbnd, n_p
    )
    # Add back on the centre of mass which has been substracted off by NBODY6
    cluster.add_stars(x + xc, y + yc, z + zc, vx, vy, vz, m, i_d)
    cluster.add_sse(kw, logl, logr, np.zeros(nbnd), np.zeros(nbnd))
    cluster.add_energies(kin, pot, etot)

    if kwargs.get("do_key_params", True):
        do_order=kwargs.get("do_key_params", True)

        # Estimate centre of distribution using median function
        cluster.find_centre()
        cluster.to_centre(do_key_params=True, do_order=do_order)
        cluster.to_cluster()

    cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc)

    return cluster


def get_snapshot(
    filename=None,
    tphys=0.0,
    ctype="snapshot",
    col_names=["m", "x", "y", "z", "vx", "vy", "vz"],
    col_nums=[0, 1, 2, 3, 4, 5, 6],
    units="pckms",
    origin="cluster",
    ofile=None,
    advance=False,
    **kwargs
):
    """
    NAME:

       get_snapshot

    PURPOSE:

       Load a generic snapshot where column names and numbers can be manually assigned

    INPUT:

       ctype - code used to generate data (nbody6/nbody6se/gyrfalcon/...)

       filename - name of file

       col_names - names corresponding to mass, position, and velocity

       col_nums - column numbers corresponding to each column name

       units - units of input data (default: kpckms)

       origin - origin of input data (default: cluster)

       ofile - opened file containing orbital information

       advance - is this a snapshot that has been advanced to from initial load_cluster?

       
    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

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

    if units == "WDunits":
        vcon = 220.0 / bovy_conversion.velocity_in_kpcGyr(220.0, 8.0)
        mcon = 222288.4543021174
        units = "kpckms"
    else:
        vcon = 1.0
        mcon = 1.0

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
            print(cluster.ntot)
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
        print(cluster.ntot)
        return cluster

    mindx = np.argwhere(col_names == "m")[0][0]
    m = data[:, col_nums[mindx]] * mcon

    xindx = np.argwhere(col_names == "x")[0][0]
    x = data[:, col_nums[xindx]]
    yindx = np.argwhere(col_names == "y")[0][0]
    y = data[:, col_nums[yindx]]
    zindx = np.argwhere(col_names == "z")[0][0]
    z = data[:, col_nums[zindx]]
    vxindx = np.argwhere(col_names == "vx")[0][0]
    vx = data[:, col_nums[vxindx]] * vcon
    vyindx = np.argwhere(col_names == "vy")[0][0]
    vy = data[:, col_nums[vyindx]] * vcon
    vzindx = np.argwhere(col_names == "vz")[0][0]
    vz = data[:, col_nums[vzindx]] * vcon

    if "id" in col_names:
        idindx = np.argwhere(col_names == "id")[0][0]
        i_d = data[:, col_nums[idindx]]
    else:
        i_d = np.linspace(0, len(x), len(x)) + 1

    if "kw" in col_names:
        kwindx = np.argwhere(col_names == "kw")[0][0]
        kw = data[:, col_nums[kwindx]]
    else:
        kw = np.zeros(len(x))

    nbnd = len(m)

    cluster = StarCluster(
        nbnd, tphys, units=units, origin=origin, ctype=ctype, **kwargs
    )
    cluster.add_stars(x, y, z, vx, vy, vz, m, i_d)
    cluster.kw = kw

    if origin == "galaxy":
        if ofile == None:
            cluster.find_centre()
        else:
            get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

        if kwargs.get("do_key_params", True):
            do_order=kwargs.get("do_key_params", True)

            cluster.to_cluster()
            cluster.find_centre()
            cluster.to_centre(do_key_params=True, do_order=do_order)
            cluster.to_galaxy()

    elif origin == "cluster":
        if kwargs.get("do_key_params", True):
            do_order=kwargs.get("do_key_params", True)
            # Estimate centre of distribution
            cluster.find_centre()
            cluster.to_centre(do_key_params=True, do_order=do_order)
            cluster.to_cluster()

        if ofile != None:
            get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    return cluster


def get_nbody6_snapauto(
    filename=None, units="pckms", origin="cluster", ofile=None, advance=False, **kwargs
):
    """
    NAME:

       get_nbody6_snapauto

    PURPOSE:

       Load a single snapshot as produced by Jongsuk Hong's snpauto.f routine for extracting snapshots from
       Nbody6 and Nbody6++ binary files 

    INPUT:

       filename = name of file

       units - units of input data (default: kpckms)

       origin - origin of input data (default: cluster)

       ofile - opened file containing orbital information

       advance - is this a snapshot that has been advanced to from initial load_cluster?

    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    cluster = get_snapshot(
        ctype="snapauto",
        filename=filename,
        col_names=["m", "x", "y", "z", "vx", "vy", "vz", "id"],
        col_nums=[0, 1, 2, 3, 4, 5, 6, 7],
        units=units,
        origin=origin,
        ofile=ofile,
        advance=advance,
        **kwargs
    )

    return cluster


def get_clustertools_snapshot(
    filename=None, units="pckms", origin="cluster", ofile=None, advance=False, **kwargs
):
    """
    NAME:

       get_clustertools_snapshot

    PURPOSE:

       Load a single snapshot as produced by clustertools snapout routine

    INPUT:

       filename = name of file

       units - units of input data (default: kpckms)

       origin - origin of input data (default: cluster)

       ofile - opened file containing orbital information

       advance - is this a snapshot that has been advanced to from initial load_cluster?
       
    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    cluster = get_snapshot(
        ctype="clustertools",
        filename=filename,
        col_names=["m", "x", "y", "z", "vx", "vy", "vz", "id", "kw"],
        col_nums=[0, 1, 2, 3, 4, 5, 6, 7, 8],
        units=units,
        origin=origin,
        ofile=ofile,
        advance=advance,
        **kwargs
    )

    return cluster


def get_amuse_particles(
    particles, units="kpckms", origin="galaxy", ofile=None, **kwargs
):
    """
    NAME:

       get_amuse_particles

    PURPOSE:

       Convert AMUSE particle dataset to a StarCluster instance

    INPUT:

       particles - AMUSE particle dataset

       units - units of input data (default: kpckms)

       origin - origin of input data (default: cluster)
       
    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    cluster = StarCluster(
        len(particles),
        tphys=0.0,
        units=units,
        origin=origin,
        ctype="amuse",
        **kwargs
    )
    i_d = np.linspace(1, len(particles), len(particles), dtype="int")

    m = particles.mass.value_in(u.MSun)

    if units == "pckms":
        x = particles.x.value_in(u.parsec)
        y = particles.y.value_in(u.parsec)
        z = particles.z.value_in(u.parsec)
        vx = particles.vx.value_in(u.kms)
        vy = particles.vy.value_in(u.kms)
        vz = particles.vz.value_in(u.kms)

    elif units == "kpckms":
        x = particles.x.value_in(u.kpc)
        y = particles.y.value_in(u.kpc)
        z = particles.z.value_in(u.kpc)
        vx = particles.vx.value_in(u.kms)
        vy = particles.vy.value_in(u.kms)
        vz = particles.vz.value_in(u.kms)

    else:
        print("PLEASE SPECIFY UNITS")
        return 0

    cluster.add_stars(x, y, z, vx, vy, vz, m, i_d, do_key_params=True)

    if origin == "galaxy":
        if ofile == None:
            cluster.find_centre()
        else:
            get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

        if kwargs.get("do_key_params", True):
            do_order=kwargs.get("do_key_params", True)
            cluster.to_cluster()
            cluster.find_centre()
            cluster.to_centre(do_key_params=True, do_order=do_order)
            cluster.to_galaxy()

    elif origin == "cluster":
        if kwargs.get("do_key_params", True):
            do_order=kwargs.get("do_key_params", True)
            # Estimate centre of distribution
            cluster.find_centre()
            cluster.to_centre(do_key_params=True, do_order=do_order)
            cluster.to_cluster()

        if ofile != None:
            get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    return cluster


def get_mycode():
    """
    NAME:

       get_mycode

    PURPOSE:

       Empty shell for a new user to write their own routine for reading in data and generating a StarCluster instance
       --> Easiest to build around get_snapshot if files are separated
    INPUT:

       TBD
       
    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       XXXX - Written - XXXX
    """

    cluster = StarCluster(0.0)
    return cluster
