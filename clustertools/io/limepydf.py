""" Initialize a StarCluster with limepy

"""

__author__ = "Jeremy J Webb"
__all__ = [
    "setup_cluster",
    "c_to_w0",
    "w0_to_c",
]

import numpy as np
from ..util.constants import *

try:
    from galpy.util import coords,conversion
except:
    import galpy.util.bovy_coords as coords
    import galpy.util.bovy_conversion as conversion

from galpy import potential
from galpy.orbit import Orbit
import os

try:
    import limepy
    from limepy import limepy,spes,sample
except:
    pass

from ..cluster.cluster import StarCluster

def _get_limepy(units="pckms", origin="cluster", orbit=None, ofile=None, advance=False, **kwargs):
    """ Setup an N-body realization of a StarCluster with specific parameters
    
    -Relies heavily on LIMEPY/SPES models (Woolley 1954, King 1966, Wilson, 1975, Gieles & Zocchi 2015, Claydon et al. 2019)
    - When setting up a specific Galactic cluster, makes use of de Boer et al. 2019 and Harris 1996 (2010 Edition). Cluster is also assigned an orbit based on Vasiliev 2019

    de Boer, T. J. L., Gieles, M., Balbinot, E., Hénault-Brunet, V., Sollima, A., Watkins, L. L., Claydon, I. 2019, MNRAS, 485, 4906
    Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576
    Harris, W.E. 1996 (2010 Edition), AJ, 112, 1487
    Vasiliev E., 2019, MNRAS, 484,2832  

    Parameters
    ----------
    units : str
        units of generated model (default: 'pckms')
    origin : str
        origin of generated model (default: 'cluster')
    orbit : class
        Galpy orbit of cluster to be generated
    ofile : file
        opened file containing orbital information
    advance : bool
        is this a snapshot that has been advanced to from initial  load_cluster? (default: False)

    Returns
    -------
    cluster: class
        StarCluster

    Other Parameters
    ----------------
    N : int
        number of stars in the cluster (default: 1000)
    model : str/object
        limepy model object
    gcname : str
        name of globular cluster to generate model for
    source : str
        Source for extracting Galactic GC parameters (Harris 1996 (2010 Edition) or de Boer et al. 2019). The default checks 
            de Boer et al. 2019 first and then pulls from Harris 1996 (2010 Edition) if no cluster found
    mbar : float
        mean mass of stars in the cluster (only single mass models available at the moment)
    citation : boolean
        print the proper citaion for setting up a mock Galactic GC (Default: True)
    kwargs : str
        Additional key word arguments needed by limepy and spes models can be passed. See https://readthedocs.org/projects/limepy/

    History
    -------
    2019 - Written - Webb (UofT)
    """

    gcname=kwargs.pop("gcname",None)
    model=kwargs.pop("model",None)
    citation=kwargs.pop("citation",True)

    if gcname is not None:
        source = kwargs.pop("source", "default")
        mbar = kwargs.pop("mbar", 1.)
        cluster = _get_cluster(gcname, source, mbar, **kwargs)

        if cluster.units!=units:
            if units=='nbody': cluster.reset_nbody_scale(rvirial=True)
            cluster.to_units(units)

        if cluster.origin!=origin:
            cluster.to_origin(origin)

    else:
        cluster = _sample_limepy(model=model, **kwargs)
        cluster.units=units

        if cluster.origin!=origin:
            cluster.to_origin(origin)

    cluster.ctype='limepy'

    # Add galpy orbit if given
    if orbit != None:
        cluster.orbit = orbit
        t = (cluster.tphys / 1000.0) / conversion.time_in_Gyr(ro=solar_ro, vo=solar_vo)
        cluster.add_orbit(
            orbit.x(t),
            orbit.y(t),
            orbit.z(t),
            orbit.vx(t),
            orbit.vy(t),
            orbit.vz(t),
            "kpckms",
        )

    elif ofile != None:
        _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)

    if origin=='galaxy':
        cluster.to_galaxy()

    cluster.analyze(sortstars=True)

    if gcname is not None and citation: 
        print('LOAD_CLUSTER MADE USE OF:')
        print("Gieles, M. & Zocchi, A. 2015, MNRAS, 454, 576")
        print("Vasiliev E., 2019, MNRAS, 484,2832 ")
        if source=="default":
            print("de Boer, T. J. L., Gieles, M., Balbinot, E., Hénault-Brunet, V., Sollima, A., Watkins, L. L., Claydon, I. 2019, MNRAS, 485, 4906")
        else:
            print("Harris, W.E. 1996 (2010 Edition), AJ, 112, 1487")

    return cluster


def _sample_limepy(g=1,model=None,**kwargs):
    """Get an Nbody realization of a LIMEPY model cluster

    Parameters
    ----------
    g : float
        model type for LIMEPY
    model : object
        LIMEPY model (default:None)
    phi0/W0 : float
        Central potential
    project : bool
        return projected values
    M : float
        cluster mass
    ro/rh/rv/rt : float
        radius used for scaling population from model units (scale,half-mass,virial,limitig radius)
    N : int
        number of stars to generate in model cluster

    Returns
    -------
        cluster : class
            StarCluster

    History
    -------
       2019 - Written - Webb (UofT)
    """

    if model is not None:
        lmodel=model
        M=lmodel._MS
    else:
        phi0=kwargs.get("W0",None)
        if phi0 is None:
            phi0 = float(kwargs.get("phi0"))
        else:
            phi0=float(phi0)

        project = bool(kwargs.get("projected", False))

        if "M" in kwargs:
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
            elif "rm" in kwargs:
                rh = float(kwargs.get("rm"))
                lmodel = limepy(phi0, g, M=M, rh=rh, project=project)
            elif "r0" in kwargs:
                r0 = float(kwargs.get("r0"))
                lmodel = limepy(phi0, g, M=M, r0=r0, project=project)
            elif "rc" in kwargs:
                r0 = float(kwargs.get("rc"))
                lmodel = limepy(phi0, g, M=M, r0=r0, project=project)
            else:
                lmodel = limepy(phi0, g, M=M, rv=1.0, project=project)
        else:
            M=1.
            lmodel = limepy(phi0, g, G=1, M=1, rv=1, project=project)

    mbar = kwargs.get("mbar", 1.)
    N = int(kwargs.get("N", M/mbar))

    ldata = sample(lmodel, N=N)

    cluster = StarCluster(units=None, origin="cluster")
    cluster.add_stars(
        ldata.x,
        ldata.y,
        ldata.z,
        ldata.vx,
        ldata.vy,
        ldata.vz,
        ldata.m,
        np.linspace(1, N, N, dtype=int),

    )
    cluster.find_centre()

    return cluster

def c_to_w0(c, invert=False):
    """ Convert King central concentration (c) values to central potential (W0) values 

    Parameters
    ----------
    c : float
        central concentration

    invert : bool
        convert from W0 to c instead, in which case the input c is the central potential (default: False)

    Returns
    -------
    W0 : float
        central potential (or c if invert==True)

    History
    -------
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
    """Convert central potential (W0) values to King central concentration (c)

    Parameters
    ----------
    W0 : float
        central potential

    Returns
    -------
    c : float
        central concentration

    History
    -------
       2019 - Written - Webb (UofT)

    """
    return c_to_w0(w0, invert=True)


def _get_cluster(gcname, source="default", mbar=1., params=False, **kwargs):
    """Generate a StarCluster based on a Galactic Globular Cluster

    By default, the functio looks for the best LIMEPY model fit to the cluster from de Boer et al. 2019. Since
    de Boer et al. 2019 does not fit every Galactic cluster, the function then looks to Harris 1996 (2010) edition
    for the best fit King (1966) model. Alternatively one can also specifiy de Boer et al. 2019 or Harris 1996 (2010).
    It is important to note that de Boer et al. 2019 find that their SPES models fit observed clusters better than LIMEPY models, 
    however a sampler for SPES models has not been implemented.

    Parameters
    ----------
    gcname : str
        name of cluster to be modelled.
    source : str
        source of model parameters to generate cluster. 
    mbar : float
        mean mass of stars in model (default: 1. Msun)
    params : bool
        return mass and size of clusters generate (default: False)

    Returns
    _______
    cluster : class
        StarCluster 
    if params == True:
        Mtot :float
            mass of cluster
        rm : fliat
            half-mass radius of cluster

    History
    -------
    2019 - Written - Webb (UofT)
    """
    data_dir=os.path.join(os.path.dirname(__file__), 'data/')

    ddata = np.loadtxt(data_dir+"deBoer2019.dat", str, skiprows=1
    )
    dname = ddata[:, 0]
    dmass = ddata[:, 7].astype(float)
    drad = ddata[:, 5].astype(float)

    hdata = np.loadtxt(data_dir+"harris2010.dat", str, skiprows=2
    )
    hname = hdata[:, 0]
    hname2 = hdata[:, 1]
    hmass = hdata[:, 2].astype(float)
    hrad = hdata[:, 4].astype(float)

    name_list = []
    mass_list = []
    rm_list = []

    gcname = gcname.upper()
    if (
        source == "default" or "deboer" in source or "deBoer" in source
    ) and gcname in dname:
        cluster = _get_deBoer_cluster(ddata, gcname, mbar, **kwargs)
    elif (source == "default" or "harris" in source or "Harris" in source) and (
        gcname in hname or gcname in hname2
    ):
        cluster = _get_harris_cluster(hdata, gcname, mbar, **kwargs)
    else:
        print('No match: ',source,gcname, gcname in dname, gcname in hname, gcname in hname2)
        print(dname)
        print(hname)
        print(hname2)
        
    return cluster

def _get_deBoer_cluster(data, gcname, mbar=0.4, **kwargs):
    """Generate a StarCluster instance based on a measurements of Galactic Globular Clusters by de Boer et al. 2019
       --> Cluster is also assigned an orbit based on Vasiliev 2019

    Parameters
    ----------
    data : float
        table of parameters from de Boer et al. 2019 (see ./data)

    gcname : str
        name of cluster to be modelled. 
    mbar : float
        mean mass of stars in model

    Returns
    -------
    cluster : class
        StarCluster

    History
    -------
    2019 - Written - Webb (UofT)
    """
    spes = False  # Not yet implemented into LIMEPY
    limepy = True

    name = data[:, 0]

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
        N = int(kwargs.get("N", M_pe / mbar))
        cluster = _get_spes(
            phi0=W_pe, B=B_pe, eta=eta_pe, fpe=fpe, M=M_pe, rt=rt_pe, N=N
        )

    else:
        W_lime = data[i_d, 1].astype(float)
        g_lime = data[i_d, 3].astype(float)
        rt_lime = data[i_d, 5].astype(float)
        M_lime = data[i_d, 7].astype(float)
        N = M_lime / mbar

        cluster = _sample_limepy(g=g_lime, phi0=W_lime, M=M_lime, rt=rt_lime, N=N)

    cluster.units='pckms'
    cluster.origin='cluster'

    cluster.orbit = _get_cluster_orbit(gcname)

    if cluster.orbit != -1:
        cluster.add_orbit(
            cluster.orbit.x(),
            cluster.orbit.y(),
            cluster.orbit.z(),
            cluster.orbit.vx(),
            cluster.orbit.vy(),
            cluster.orbit.vz(),
            ounits="kpckms",
        )
    else:
        cluster.orbit = None

    return cluster


def _get_harris_cluster(data, gcname, mbar=0.4, **kwargs):
    """Generate a StarCluster instance based on the Harris 1996 (2010 Edition) catalogue of Galactic Globular Clusters 
       --> Cluster is also assigned an orbit based on Vasiliev 2019

    Parameters
    ----------
    data : float
        table of parameters from de Boer et al. 2019 (see ./data)

    gcname : str
        name of cluster to be modelled. 
    mbar : float
        mean mass of stars in model

    Returns
    -------
    cluster : class
        StarCluster

    History
    -------
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
    N = int(kwargs.get("N", mgc / mbar))

    cluster = _sample_limepy(g=1.0, phi0=w0, M=mgc, rt=rl, N=N)

    cluster.units='pckms'
    cluster.origin='cluster'

    cluster.orbit = _get_cluster_orbit(gcname)

    if cluster.orbit != -1:
        cluster.add_orbit(
            cluster.orbit.x(),
            cluster.orbit.y(),
            cluster.orbit.z(),
            cluster.orbit.vx(),
            cluster.orbit.vy(),
            cluster.orbit.vz(),
            ounits="kpckms",
        )
    else:
        cluster.orbit = None

    return cluster

def _get_cluster_orbit(gcname,ro=solar_ro, vo=solar_vo):
    """Get the measured orbital parameters of a Galactic globular cluster
    - This is a simply wrapper for Orbit.from_name in galpy (Bovy 2015), which uses orbits measured by Vasiliev 2019 using Gaia DR2 via Galpy
    -- Bovy J., 2015, ApJS, 216, 29

    Parameters
    ----------
    gcname : str
        name of GC whose orbits is to be retrieved
    ro :float 
        galpy distance scale (Default: 8.)
    vo : float
        galpy velocity scale (Default: 220.)
    Returns
    -------
    orbit : class
        galpy orbit

    History
    -------
    2019 - Written - Webb (UofT)

    """
    return Orbit.from_name(gcname,ro=ro, vo=vo, solarmotion=solar_motion)

