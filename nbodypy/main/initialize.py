import numpy as np
from galpy.util import bovy_conversion
import os

# import limepy
from ..community.limepy import *

from .cluster import StarCluster
from .profiles import m_prof
from .orbit import get_cluster_orbit


def setup_cluster(ctype, units="realpc", origin="cluster", orbit=None, **kwargs):
    """
    NAME:

       setup_cluster

    PURPOSE:

       Setup an N-body realization of a StarCluster with specific parameters
       --> Relies heavily on LIMEPY/SPES models (REFERENCE)

    INPUT:

       ctype - Type of model used to generate cluster ('SPES',LIMEPY','WOOLEY','KING','WILSON',Name/List of Galactic Clusters)
             --> To add: PLUMMER

       units - units of generated model (Default: 'realpc')

       origin - origin of generated model (Default: 'cluster')

       orbit - Galpy orbit of cluster to be generate


    KWARGS:

       g - model type for LIMEPY

       M/rh - Mass and half-mass radius of cluster (for Plummer)

       source - Source for extracting Galactic GC parameters (Harris 1996 (2010 Edition), de Boer et al. 2019). The default checks 
                de Boer et al. 2019 first and then pulls from Harris 1996 (2010 Edition) if no cluster found

       mbar - mean mass of stars in the cluster (only single mass models available at the moment)

       names - return names of Galactic clusters generate (Default: False)

       -- Additional KWARGS passed to get_limepy, and get_spes

    OUTPUT:

        StarCluster instance

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    if isinstance(ctype, str):

        if ctype == "plummer":
            M = float(kwargs.get("M"))
            rh = float(kwargs.get("rh"))
            cluster = get_plummer(M, rh)
        elif ctype == "spes":
            cluster = get_spes(**kwargs)
        elif ctype == "limepy":
            g = kwargs.pop("g")
            cluster = get_limepy(g=g, **kwargs)
        elif ctype == "woolley":
            g = kwargs.pop("g", 0)
            cluster = get_limepy(g=g, **kwargs)
        elif ctype == "king":
            g = kwargs.pop("g", 1)
            cluster = get_limepy(g=g, **kwargs)
        elif ctype == "wilson":
            g = kwargs.pop("g", 2)
            cluster = get_limepy(g=g, **kwargs)
        else:
            source = kwargs.pop("source", "default")
            mbar = kwargs.pop("mbar", 0.4)
            names = kwargs.pop("names", False)
            cluster = get_cluster(ctype, source, mbar, names)
    else:
        source = kwargs.pop("source", "default")
        mbar = kwargs.pop("mbar", 0.4)
        names = kwargs.pop("names", False)
        cluster = get_cluster(ctype, source, mbar, names)

    # Add galpy orbit if given
    if orbit != None:
        cluster.orbit = orbit
        t = (cluster.tphys / 1000.0) / bovy_conversion.time_in_Gyr(ro=8.0, vo=220.0)
        cluster.add_orbit(
            orbit.x(t),
            orbit.y(t),
            orbit.z(t),
            orbit.vx(t),
            orbit.vy(t),
            orbit.vz(t),
            "realkpc",
        )

    cluster.key_params(do_order=True)

    return cluster


def get_limepy(g=1, **kwargs):
    """
    NAME:

       get_limepy

    PURPOSE:

       Get an Nbody realization of a LIMEPY model cluster (REFERENCE)

    INPUT:

       g - model type for LIMEPY

    KWARGS:

       phi0 - Central potential

       project - return projected values

       M - cluster mass

       r0/rh/rv/rt - radius used for scaling population from model units (scale,half-mass,virial,limitig radius)

       N - number of stars to generate in model cluster

    OUTPUT:

        StarCluster instance

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    phi0 = float(kwargs.get("phi0"))
    project = bool(kwargs.get("project", False))

    if "M" in kwargs:
        units = "realpc"
        M = float(kwargs.get("M"))
        if "rt" in kwargs:
            rt = float(kwargs.get("rt"))
            lmodel = limepy(phi0, g, M=M, rt=rt, project=project)
        elif "rv" in kwargs:
            rv = float(kwargs.get("rv"))
            lmodel = limepy(phi0, g, M=M, rv=rv, project=project)
        elif "rh" in kwargs:
            rh = float(kwargs.get("rh"))
            lmodel = limepy(phi0, g, M=M, rh=rh, project=project)
        elif "r0" in kwargs:
            r0 = float(kwargs.get("r0"))
            lmodel = limepy(phi0, g, M=M, r0=r0, project=project)
        else:
            lmodel = limepy(phi0, g, M=M, r0=1.0, project=project)
    else:
        units = "nbody"
        lmodel = limepy(phi0, g, G=1, M=1, rv=1, project=project)

    N = int(kwargs.get("N", 1000))

    ldata = sample(lmodel, N=N)

    cluster = StarCluster(N, units=units, origin="cluster")
    cluster.ctype = "limepy"
    cluster.add_stars(
        np.linspace(1, N, N, dtype=int),
        ldata.m,
        ldata.x,
        ldata.y,
        ldata.z,
        ldata.vx,
        ldata.vy,
        ldata.vz,
    )
    cluster.find_centre()
    cluster.key_params()

    return cluster


def get_spes(**kwargs):
    """
    NAME:

       get_spes

    PURPOSE:

       Get an Nbody realization of a SPES model cluster (REFERENCE)

    INPUT:

       None

    KWARGS:

       phi0 - Central potential

       project - return projected values

       B - B parameter

       eta - eta parameter

       fpe - fpe parameter

       M - cluster mass

       r0/rh/rv/rt - radius used for scaling population from model units (scale,half-mass,virial,limitig radius)

       N - number of stars to generate in model cluster

    OUTPUT:

        StarCluster instance

    HISTORY:

       2019 - Written - Webb (UofT)

    """
    phi0 = float(kwargs.get("phi0"))
    B = float(kwargs.get("B"))
    eta = float(kwargs.get("eta"))
    fpe = float(kwargs.get("fpe"))

    project = bool(kwargs.get("project", False))

    if "M" in kwargs:
        units = "realpc"
        M = float(kwargs.get("M"))
        if "rt" in kwargs:
            rt = float(kwargs.get("rt"))
            smodel = spes(phi0, B=B, eta=eta, fpe=fpe, M=M, rt=rt, project=project)
        elif "rv" in kwargs:
            rv = float(kwargs.get("rv"))
            smodel = spes(phi0, B=B, eta=eta, fpe=fpe, M=M, rv=rv, project=project)
        elif "rh" in kwargs:
            rh = float(kwargs.get("rh"))
            smodel = spes(phi0, B=B, eta=eta, fpe=fpe, M=M, rh=rh, project=project)
        elif "r0" in kwargs:
            r0 = float(kwargs.get("r0"))
            smodel = spes(phi0, B=B, eta=eta, fpe=fpe, M=M, r0=r0, project=project)
        else:
            smodel = spes(phi0, B=B, eta=eta, fpe=fpe, M=M, r0=1.0, project=project)
    else:
        units = "nbody"
        smodel = spes(phi0, B=B, eta=eta, fpe=fpe, G=1, M=1, rv=1, project=project)

    N = int(kwargs.get("N", 1000))

    sdata = sample(smodel, N=N)

    cluster = StarCluster(N, units=units, origin="cluster")
    cluster.ctype = "spes"
    cluster.add_stars(
        np.linspace(1, N, N, dtype=int),
        sdata.m,
        sdata.x,
        sdata.y,
        sdata.z,
        sdata.vx,
        sdata.vy,
        sdata.vz,
    )
    cluster.find_centre()
    cluster.key_params()

    return cluster


def get_plummer(M, rm):
    # WIP
    cluster = StarCluster()

    return cluster


def c_to_w0(c, invert=False):
    """
    NAME:

       c_to_w0

    PURPOSE:

       Convert King central concentration (c) values to central potential (W0) values 

    INPUT:

       c - central concentration

       invert - convert from W0 to c instead, in which case the input c is the central potential (Default: False)

    OUTPUT:

        W0 (or c if invert==True)

    HISTORY:

       2019 - Written - Webb (UofT)

    """
    # From gridfit (Dean McLaughlin)
    w0 = np.array(
        [
            0.300000,
            0.400000,
            0.500000,
            0.600000,
            0.700000,
            0.800000,
            0.900000,
            1.000000,
            1.100000,
            1.200000,
            1.300000,
            1.400000,
            1.500000,
            1.600000,
            1.700000,
            1.800000,
            1.900000,
            2.000000,
            2.100000,
            2.200000,
            2.300000,
            2.400000,
            2.500000,
            2.600000,
            2.700000,
            2.800000,
            2.900000,
            3.000000,
            3.100000,
            3.200000,
            3.300000,
            3.400000,
            3.500000,
            3.600000,
            3.700000,
            3.800000,
            3.900000,
            4.000000,
            4.100000,
            4.200000,
            4.300000,
            4.400000,
            4.500000,
            4.600000,
            4.700000,
            4.800000,
            4.900000,
            5.000000,
            5.100000,
            5.200000,
            5.300000,
            5.400000,
            5.500000,
            5.600000,
            5.700000,
            5.800000,
            5.900000,
            6.000000,
            6.100000,
            6.200000,
            6.300000,
            6.400000,
            6.500000,
            6.600000,
            6.700000,
            6.800000,
            6.900000,
            7.000000,
            7.100000,
            7.200000,
            7.300000,
            7.400000,
            7.500000,
            7.600000,
            7.700000,
            7.800000,
            7.900000,
            8.000000,
            8.100000,
            8.200000,
            8.300000,
            8.400000,
            8.500000,
            8.600000,
            8.700000,
            8.800000,
            8.900000,
            9.000000,
            9.100000,
            9.200000,
            9.300000,
            9.400000,
            9.500000,
            9.600000,
            9.700000,
            9.800000,
            9.900000,
            10.000000,
            10.100000,
            10.200000,
            10.300000,
            10.400000,
            10.500000,
            10.600000,
            10.700000,
            10.800000,
            10.900000,
            11.000000,
            11.100000,
            11.200000,
            11.300000,
            11.400000,
            11.500000,
            11.600000,
            11.700000,
            11.800000,
            11.900000,
            12.000000,
            12.100000,
            12.200000,
            12.300000,
            12.400000,
            12.500000,
            12.600000,
            12.700000,
            12.800000,
            12.900000,
            13.000000,
            13.100000,
            13.200000,
            13.300000,
            13.400000,
            13.500000,
            13.600000,
            13.700000,
            13.800000,
            13.900000,
            14.000000,
            14.100000,
            14.200000,
            14.300000,
            14.400000,
            14.500000,
            14.600000,
            14.700000,
            14.800000,
            14.900000,
            15.000000,
            15.100000,
            15.200000,
            15.300000,
            15.400000,
            15.500000,
            15.600000,
            15.700000,
            15.800000,
            15.900000,
            16.000000,
            16.100000,
            16.200000,
            16.300000,
            16.400000,
            16.500000,
            16.600000,
            16.700000,
            16.800000,
            16.900000,
            17.000000,
            17.100000,
            17.200000,
            17.300000,
            17.400000,
            17.500000,
            17.600000,
            17.700000,
            17.800000,
            17.900000,
            18.000000,
            18.100000,
            18.200000,
            18.300000,
            18.400000,
            18.500000,
            18.600000,
            18.700000,
            18.800000,
            18.900000,
            19.000000,
            19.100000,
            19.200000,
            19.300000,
            19.400000,
            19.500000,
        ]
    )

    conc = np.array(
        [
            1.004710,
            1.171350,
            1.322660,
            1.463770,
            1.597760,
            1.726690,
            1.851980,
            1.974730,
            2.095780,
            2.215830,
            2.335470,
            2.455190,
            2.575460,
            2.696690,
            2.819260,
            2.943540,
            3.069880,
            3.198640,
            3.330170,
            3.464800,
            3.602910,
            3.744850,
            3.891000,
            4.041760,
            4.197530,
            4.358750,
            4.525890,
            4.699410,
            4.879840,
            5.067720,
            5.263660,
            5.468270,
            5.682240,
            5.906290,
            6.141230,
            6.387900,
            6.647220,
            6.920200,
            7.207920,
            7.511580,
            7.832460,
            8.171960,
            8.531630,
            8.913140,
            9.318330,
            9.749200,
            10.208000,
            10.697100,
            11.219100,
            11.777000,
            12.374000,
            13.013600,
            13.699700,
            14.436500,
            15.228700,
            16.081400,
            17.000300,
            17.991600,
            19.062000,
            20.218900,
            21.470400,
            22.825300,
            24.293000,
            25.883700,
            27.608600,
            29.479300,
            31.508300,
            33.708600,
            36.093800,
            38.678000,
            41.475300,
            44.499900,
            47.765800,
            51.286400,
            55.074300,
            59.140700,
            63.495500,
            68.146800,
            73.100800,
            78.361400,
            83.930800,
            89.808800,
            95.993700,
            102.482000,
            109.270000,
            116.352000,
            123.724000,
            131.381000,
            139.319000,
            147.537000,
            156.034000,
            164.810000,
            173.871000,
            183.222000,
            192.871000,
            202.828000,
            213.106000,
            223.721000,
            234.690000,
            246.032000,
            257.769000,
            269.924000,
            282.522000,
            295.592000,
            309.162000,
            323.263000,
            337.930000,
            353.196000,
            369.099000,
            385.679000,
            402.977000,
            421.038000,
            439.907000,
            459.634000,
            480.270000,
            501.871000,
            524.493000,
            548.199000,
            573.053000,
            599.122000,
            626.479000,
            655.201000,
            685.368000,
            717.064000,
            750.382000,
            785.415000,
            822.265000,
            861.038000,
            901.847000,
            944.812000,
            990.059000,
            1037.720000,
            1087.940000,
            1140.860000,
            1196.650000,
            1255.460000,
            1317.470000,
            1382.870000,
            1451.860000,
            1524.620000,
            1601.400000,
            1682.400000,
            1767.870000,
            1858.070000,
            1953.250000,
            2053.690000,
            2159.690000,
            2271.550000,
            2389.600000,
            2514.160000,
            2645.600000,
            2784.280000,
            2930.590000,
            3084.940000,
            3247.740000,
            3419.440000,
            3600.510000,
            3791.430000,
            3992.700000,
            4204.840000,
            4428.420000,
            4664.010000,
            4912.200000,
            5173.620000,
            5448.940000,
            5738.830000,
            6044.000000,
            6365.210000,
            6703.230000,
            7058.880000,
            7433.020000,
            7826.530000,
            8240.350000,
            8675.460000,
            9132.880000,
            9613.680000,
            10119.000000,
            10650.000000,
            11207.900000,
            11794.100000,
            12409.800000,
            13056.600000,
            13735.900000,
            14449.300000,
            15198.400000,
            15985.100000,
            16811.200000,
            17678.500000,
            18589.200000,
            19545.300000,
            20549.200000,
            21603.100000,
            22709.700000,
        ]
    )

    if invert:

        w = c
        indx = np.argmin(abs(w0 - w))

        if w0[indx] < w:
            m = (conc[indx + 1] - conc[indx]) / (w0[indx + 1] - w0[indx])
        else:
            m = (conc[indx] - conc[indx - 1]) / (w0[indx] - w0[indx - 1])

        b = conc[indx] - m * w0[indx]

        return np.log10(m * w + b)

    else:

        c = 10.0 ** c

        indx = np.argmin(abs(conc - c))

        if conc[indx] < c:
            m = (w0[indx + 1] - w0[indx]) / (conc[indx + 1] - conc[indx])
        else:
            m = (w0[indx] - w0[indx - 1]) / (conc[indx] - conc[indx - 1])

        b = w0[indx] - m * conc[indx]

        return m * c + b


def w0_to_c(w0):
    """
    NAME:

       w0_to_c

    PURPOSE:

       Convert central potential (W0) values to King central concentration (c)

    INPUT:

       W0 - central potential


    OUTPUT:

       c

    HISTORY:

       2019 - Written - Webb (UofT)

    """
    return c_to_w0(w0, invert=True)


def get_cluster(gcname="list", source="default", mbar=0.4, names=False, params=False):
    """
    NAME:

       get_cluster

    PURPOSE:

       Generate a StarCluster instance based on a Galactic Globular Cluster

    INPUT:

       gcname - name of cluster (or clusters) to be modelled. Also accepts 'list' to simply print out
                the list of clusters available and 'all' to generate each Galactic Globular Cluster

       source - source of model parameters to generate cluster. Default looks to LIMEPY models from
                de Boer et al. 2019 first before checking Harris 1996 (2010) edition second.
              - Can also specifiy de Boer et al. 2019 or Harris 1996 (2010) only

        mbar - mean mass of stars in model

        names - return names of clusters genereated (Default: False)

        params -  return mass and size of clusters generate (Default: False)

    OUTPUT:

       StarCluster (if gcname != 'list')

       Name (if names==True)

       Mtot and rm (if params ==True)

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    ddata = np.loadtxt(
        "/Users/webbjj/Codes/nbodypy/tables/deBoer2019.dat", str, skiprows=1
    )
    dname = ddata[:, 0]
    dmass = ddata[:, 7].astype(float)
    drad = ddata[:, 5].astype(float)

    hdata = np.loadtxt(
        "/Users/webbjj/Codes/nbodypy/tables/harris2010.dat", str, skiprows=2
    )
    hname = hdata[:, 0]
    hname2 = hdata[:, 1]
    hmass = hdata[:, 2].astype(float)
    hrad = hdata[:, 4].astype(float)

    name_list = []
    mass_list = []
    rm_list = []

    if isinstance(gcname, str):

        if (gcname == "list" or gcname == "all") and (
            "deboer" in source or "deBoer" in source
        ):
            for i in range(0, len(dname)):
                print(dname[i])
                name_list.append(dname[i])
                mass_list.append(dmass[i])
                rm_list.append(drad[i])
        elif (gcname == "list" or gcname == "all") and (
            "harris" in source or "Harris" in source
        ):
            for i in range(0, len(hname)):
                print(hname[i], hname2[i])
                name_list.append(hname[i])
                mass_list.append(hmass[i])
                rm_list.append(hrad[i])

        elif (gcname == "list" or gcname == "all") and source == "default":
            for i in range(0, len(dname)):
                print("DEBOER: ", dname[i])
                name_list.append(dname[i])
                mass_list.append(dmass[i])
                rm_list.append(drad[i])

            indx = np.in1d(hname, dname, invert=True) * np.in1d(
                hname2, dname, invert=True
            )
            for i in range(0, np.sum(indx)):
                print("HARRIS: ", hname[indx][i])
                name_list.append(hname[indx][i])
                mass_list.append(hmass[indx][i])
                rm_list.append(hrad[indx][i])

        else:
            gcname = gcname.upper()
            if (
                source == "default" or "deboer" in source or "deBoer" in source
            ) and gcname in dname:
                cluster = get_deBoer_cluster(ddata, gcname, mbar, names)
            elif (source == "default" or "harris" in source or "Harris" in source) and (
                gcname in hname or gcname in hname2
            ):
                cluster = get_harris_cluster(hdata, gcname, mbar, names)

            if names and params:
                return cluster, gcname, cluster.mtot, cluster.rm
            elif params:
                return cluster, cluster.mtot, cluster.rm
            elif names:
                return cluster, gcname
            else:
                return cluster
    else:
        name_list = gcname

    if len(name_list) > 0 and "list" not in gcname:
        cluster = []
        cluster_name = []
        cluster_mass = []
        cluster_rm = []
        for i in range(0, len(name_list)):
            name_list[i] = name_list[i].upper()
            if (
                source == "default" or "deboer" in source or "deBoer" in source
            ) and name_list[i] in dname:
                cluster.append(get_deBoer_cluster(ddata, name_list[i], mbar, names))
                cluster[-1].ctype = name_list[i]
                cluster_name.append(name_list[i])
                cluster_mass.append(cluster[-1].mtot)
                cluster_rm.append(cluster[-1].rm)
            elif (source == "default" or "harris" in source or "Harris" in source) and (
                name_list[i] in hname or name_list[i] in hname2
            ):
                cluster.append(get_harris_cluster(hdata, name_list[i], mbar, names))

                cluster[-1].ctype = name_list[i]
                cluster_name.append(name_list[i])
                cluster_mass.append(cluster[-1].mtot)
                cluster_rm.append(cluster[-1].rm)

            else:
                print("COULD NOT FIND CLUSTER %s" % name_list[i])

        if names and params:
            return cluster, cluster_name, cluster_mass, cluster_rm
        elif params:
            return cluster, cluster_mass, cluster_rm
        elif names:
            return cluster, cluster_name
        else:
            return cluster
    elif "list" in gcname:
        if names and params:
            return name_list, mass_list, rm_list
        elif params:
            return mass_list, rm_list
        elif names:
            return name_list
    else:
        return


def get_deBoer_cluster(data, gcname, mbar=0.4, names=False):
    """
    NAME:

       get_deBoer_cluster

    PURPOSE:

       Generate a StarCluster instance based on a measurements of Galactic Globular Clusters by de Boer et al. 2019
       --> Cluster is also assigned an orbit based on Vasiliev 2019

    INPUT:

       data - table of parameters from de Boer et al. 2019 (see ..tables)

       gcname - name of cluster (or clusters) to be modelled. Also accepts 'list' to simply print out
                the list of clusters available and 'all' to generate each Galactic Globular Cluster

       mbar - mean mass of stars in model

       names - return names of clusters genereated (Default: False)

    OUTPUT:

       StarCluster

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    spes = False  # Not yet implemented into LIMEPY
    limepy = True

    name = data[:, 0]

    gcname = gcname.upper()
    indx = name == gcname
    i_d = np.argwhere(indx == True)[0]

    if spes:
        # W_pe e_W_pe eta_pe e_eta_pe log1minB_pe e_log1minB_pe rt_pe e_rt_pe M_pe e_M_pe log_fpe e_log_fpe
        W_pe = data[i_d, 9].astype(float)
        eta_pe = data[i_d, 11].astype(float)
        B_pe = 10.0 ** data[i_d, 13].astype(float)
        rt_pe = data[i_d, 15].astype(float)
        M_pe = data[i_d, 17].astype(float)
        fpe = 10.0 ** data[i_d, 19].astype(float)
        N = M_pe / mbar
        cluster = get_spes(
            phi0=W_pe, B=B_pe, eta=eta_pe, fpe=fpe, M=M_pe, rt=rt_pe, N=N
        )

    else:
        W_lime = data[i_d, 1].astype(float)
        g_lime = data[i_d, 3].astype(float)
        rt_lime = data[i_d, 5].astype(float)
        M_lime = data[i_d, 7].astype(float)
        N = M_lime / mbar

        cluster = get_limepy(g=g_lime, phi0=W_lime, M=M_lime, rt=rt_lime, N=N)

    cluster.orbit = get_cluster_orbit(name[i_d])

    if cluster.orbit != -1:
        cluster.add_orbit(
            cluster.orbit.x(),
            cluster.orbit.y(),
            cluster.orbit.z(),
            cluster.orbit.vx(),
            cluster.orbit.vy(),
            cluster.orbit.vz(),
            ounits="realkpc",
        )
    else:
        cluster.orbit = None

    return cluster


def get_harris_cluster(data, gcname, mbar=0.4, names=False):
    """
    NAME:

       get_harris_cluster

    PURPOSE:

       Generate a StarCluster instance based on the Harris 1996 (2010 Edition) catalogue of Galactic Globular Clusters 
       --> Cluster is also assigned an orbit based on Vasiliev 2019

    INPUT:

       data - table of parameters from Harris 1996 (2010 Edition) (see ..tables)

       gcname - name of cluster (or clusters) to be modelled. Also accepts 'list' to simply print out
                the list of clusters available and 'all' to generate each Galactic Globular Cluster

       mbar - mean mass of stars in model

       names - return names of clusters genereated (Default: False)

    OUTPUT:

       StarCluster

    HISTORY:

       2019 - Written - Webb (UofT)

    """

    name = data[:, 0]
    name2 = data[:, 1]

    indx = np.logical_or(np.in1d(name, gcname), np.in1d(name2, gcname))
    i_d = np.argwhere(indx == True)[0]
    mgc = data[i_d, 2].astype(float)
    rc = data[i_d, 3].astype(float)
    rh = data[i_d, 4].astype(float)
    rl = data[i_d, 5].astype(float)
    c = np.log10(rl / rc)
    w0 = c_to_w0(c)
    N = mgc / mbar
    cluster = get_limepy(g=1.0, phi0=w0, M=mgc, rt=rl, N=N)
    cluster.orbit = get_cluster_orbit(name[i_d])

    if cluster.orbit != -1:
        cluster.add_orbit(
            cluster.orbit.x(),
            cluster.orbit.y(),
            cluster.orbit.z(),
            cluster.orbit.vx(),
            cluster.orbit.vy(),
            cluster.orbit.vz(),
            ounits="realkpc",
        )
    else:
        cluster.orbit = None

    return cluster
