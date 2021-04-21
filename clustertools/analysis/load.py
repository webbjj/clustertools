""" Read in cluster from Nbody simulations or generate an Nbody cluster

"""

__author__ = "Jeremy J Webb"
__all__ = [
    "load_cluster",
    "advance_cluster",
]

import numpy as np
from galpy.util import bovy_conversion
import os, struct
from .cluster import StarCluster
from .operations import *
from .orbit import initialize_orbit

# Try Importing AMUSE. Only necessary for _get_amuse_particles
try:
    import amuse.units.units as u
    from amuse.io import read_set_from_file
except:
    pass


def load_cluster(
    ctype="snapshot",
    units = "pckms",
    origin = "cluster",
    ofile=None,
    orbit=None,
    filename=None,
    particles=None,
    load_function=None,
    **kwargs,
):
    """Load a StarCluster snapshot from a generic code output

    Parameters
    __________
    ctype : str
        Type of file being loaded
        Currently supports:

            - nbody6
            - nbody6se
            - gyrfalcon
            - snaptrim
            - napauto
            - clustertools
            - snapshot
            - astropy_table

    units : str
        units of input data (default: kpckms)
    origin : str
        origin of input data (default: cluster)
    ofile : file
        an already opened file containing orbit information (default: None)
    orbit : class
        a galpy orbit to be used for the StarCluster's orbital information (default: None)
    filename : str 
        name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)
    particles : particles
        AMUSE particle dataset (default: None)
        or `~astropy.table.Table` instance if `ctype` is "astropy_table".
    load_function : function
        use a custom function to load data (default : None)

    Returns
    -------
    cluster : class
        StarCluster

    Other Parameters
    ----------------
    ofilename : str
        orbit filename if ofile or orbit is not given
    ounits : str
        units of orbital information (else assumed equal to StarCluster.units)
    nsnap : int
        if a specific snapshot is to be read in instead of starting from zero
    nzfill : int
        value for zfill when reading and writing snapshots (Default: 5)
    delimiter : str
        choice of delimiter when reading ascii/csv files (Default: ',')
    wdir : str
        working directory of snapshots if not current directory
    snapdir : str
        directory of snapshot (Default: './')
    snapbase : str
        base for snapshot filename (Default: '')
    snapend : str
        end for snapshot filename (Default: '')
    skiprows : int
        number of rows to skip when reading in snapshot (Default: 0)
    initialize : bool
        initialize a galpy orbit after reading in orbital information (default: False)
    projected : bool
        calculate projected values as well as 3D values (Default: True)
    sortstars : bool
        sort stars in order from closes to the origin to the farthest (default: True)
    column_mapper : dict
        see _get_astropy_table 
    verbose : bool
        print additional information to screen while loading (default : False)
    give : str
        set what parameters are read in from gyrfalcon (default: 'mxv')
        Currently only accepts 'mxvpqael' as an alternative.

    History
    _______
    2018 - Written - Webb (UofT)
    """
    wdir = kwargs.get("wdir", "./")
    initialize = kwargs.get("initialize", False)

    if "ofilename" in kwargs and ofile is None:
        ofile = open(wdir + kwargs["ofilename"], "r")

    if load_function is not None:
        ctype='custom'
        if particles is not None:
            cluster=load_function(ctype=ctype,units=units,origin=origin,ofile=ofile,orbit=orbit,particles=particles,**kwargs)
        elif filename is not None:
            cluster=load_function(ctype=ctype,units=units,origin=origin,ofile=ofile,orbit=orbit,filename=filename,**kwargs)
        else:
            cluster=load_function(ctype=ctype,units=units,origin=origin,ofile=ofile,orbit=orbit,**kwargs)

    elif ctype == "nbody6":

        # With stellar evolution turned ON, read in OUT3, OUT33, fort.82 and fort.83.
        if os.path.isfile("%sOUT3" % wdir):
            out3 = open("%sOUT3" % wdir, "rb")
        else:
            out3 = None

        if os.path.isfile("%sOUT33" % wdir):
            out33 = open("%sOUT33" % wdir, "rb")
        else:
            out33=None

        if os.path.isfile("%sfort.82" % wdir):
            fort82 = open("%sfort.82" % wdir, "r")
        else:
            fort82=None

        if os.path.isfile("%sfort.83" % wdir):
            fort83 = open("%sfort.83" % wdir, "r")
        else:
            fort83=None

        cluster = _get_nbody6(out3, out33, fort82=fort82, fort83=fort83, ofile=ofile, advance=False, **kwargs)

    elif ctype == "gyrfalcon":
        # Read in snapshot from gyrfalcon.
        filein = open(wdir + filename, "r")
        cluster = _get_gyrfalcon(filein, "WDunits", "galaxy", advance=False, **kwargs)

    elif ctype=='amuse':
        if filename is not None:
            filetype=kwargs.pop("filetype","hdf5")
            particles = read_set_from_file(filename, filetype)

        cluster = _get_amuse_particles(particles, units=units, origin=origin, ofile=ofile, **kwargs)

    elif ctype == "snapshot":
        # Read in standard generic snapshot
        col_names = kwargs.pop("col_names", ["m", "x", "y", "z", "vx", "vy", "vz"])
        col_nums = kwargs.pop("col_nums", [0, 1, 2, 3, 4, 5, 6])
        cluster = _get_snapshot(
            filename=filename,
            col_names=col_names,
            col_nums=col_nums,
            units=units,
            origin=origin,
            ofile=ofile,
            advance=False,
            **kwargs,
        )
    elif ctype.lower() == "astropy_table":
        column_mapper = kwargs.pop("column_mapper", None)

        # Read data from astropy table
        cluster = _get_astropy_table(
            particles,
            column_mapper=column_mapper,
            units=units,
            origin=origin,
            ofile=ofile,
            **kwargs
        )

    else:
        print("NO CTYPE GIVEN")
        return 0

    if ofile is not None:
        cluster.ofilename=ofile.name.split('/')[-1]

    # Add galpy orbit if given
    if orbit is not None:
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
            t = (cluster.tphys * cluster.tbar / 1000.0) / bovy_conversion.time_in_Gyr(
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

        units0, origin0, rorder0, rorder_origin0 = save_cluster(cluster)

        cluster.to_cluster(sortstars=False)
        cluster.find_centre()

        return_cluster(cluster, units0, origin0, rorder0, rorder_origin0)
    elif initialize:
        initialize_orbit(cluster)

    return cluster

def advance_cluster(
    cluster,
    load_function=None,
    ofile=None,
    orbit=None,
    filename=None,
    **kwargs,
):
    """Advance a loaded StarCluster snapshot to the next timestep

   - ofile or orbit need to be provded again, same as load_cluster.
   - Be sure that advance is set to True so next line of orbit file is read in
   - if last snapshot has been reached, returns an empty StarCluster

    Parameters
    ----------
    cluster - class 
        StarCluster to be advanced
    load_function : function
        use a custom function to load data (default : None)
    ofile : file
        an already opened file containing orbit information (default: None)
    orbit : class
        a galpy orbit to be used for the StarCluster's orbital information (default: None)
    filename :str
        name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)

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
    advance_kwargs = _get_advanced_kwargs(cluster, **kwargs)

    # Continue reading in cluster opened in _get_cluster()
    if load_function is not None:
        ctype='custom'
        if filename is not None:
            cluster=load_function(ctype=ctype,units=cluster.units_init,origin=cluster.origin_init,ofile=ofile,orbit=orbit,filename=filename,advance=True,**advance_kwargs,**kwargs)
        else:
            cluster=load_function(ctype=ctype,units=cluster.units_init,origin=cluster.origin_init,ofile=ofile,orbit=orbit,advance=True,**advance_kwargs,**kwargs)


    elif cluster.ctype == "nbody6":
        cluster = _get_nbody6(
            cluster.sfile, cluster.bfile, cluster.bsefile, cluster.ssefile, advance=True, **advance_kwargs
        )

    elif cluster.ctype == "gyrfalcon":

        cluster = _get_gyrfalcon(
            cluster.sfile,
            units="WDunits",
            origin="galaxy",
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )


    elif cluster.ctype == "snapshot":
        col_names = kwargs.pop("col_names", ["m", "x", "y", "z", "vx", "vy", "vz"])
        col_nums = kwargs.pop("col_nums", [0, 1, 2, 3, 4, 5, 6])

        cluster = _get_snapshot(
            filename=filename,
            col_names=col_names,
            col_nums=col_nums,
            units=cluster.units_init,
            origin=cluster.origin_init,
            ofile=ofile,
            advance=True,
            **advance_kwargs
        )
    else:
        cluster = StarCuster(ctype=cluster.ctype,units=cluster.units_init,origin=cluster.origin_init,**advance_kwargs)

    # Check for restart
    if cluster.ntot == 0.0:
        #print('NTOT = 0',cluster.wdir,advance_kwargs.get('wdir','./'))
        try:
            wdir = cluster.wdir + "cont/"
        except:
            print("WDIR NOT SET")
            wdir = "./cont/"

        try:
            ofilename = ofile.name.split('/')[-1]
            ofile=None
        except:
            ofile = None
            ofilename = None

        if os.path.exists(wdir):
            old_wdir=advance_kwargs.pop('wdir')
            cluster = load_cluster(
                ctype=cluster.ctype,units=cluster.units_init,origin=cluster.origin_init,orbit=orbit,filename=filename,load_function=load_function,wdir=wdir,ofilename=ofilename,**advance_kwargs
            )


    if cluster.ntot != 0.0:

        # Add galpy orbit if given
        if orbit != None:
            cluster.orbit = orbit
            if cluster.units == "pckms" or cluster.units == "kpckms":
                t = (cluster.tphys / 1000.0) / bovy_conversion.time_in_Gyr(
                    ro=8.0, vo=220.0
                )
            elif cluster.units == "nbody":
                t = (
                    cluster.tphys * cluster.tbar / 1000.0
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

    return cluster


def _get_advanced_kwargs(cluster, **kwargs):
    """get **kwargs from current cluster before advancing

    Parameters
    ----------
    cluster - class 
        StarCluster to be advanced

    Returns
    -------
    None

    Other Parameters
    ----------------
    Same as load_cluster

    History
    -------
    2018 - Written - Webb (UofT)
    """

    nsnap = np.maximum(int(kwargs.pop("nsnap", 0)), cluster.nsnap) + 1
    delimiter = kwargs.pop("delimiter", cluster.delimiter)
    wdir = kwargs.pop("wdir", cluster.wdir)
    nzfill = int(kwargs.pop("nzfill", cluster.nzfill))
    snapbase = kwargs.pop("snapbase", cluster.snapbase)
    snapend = kwargs.pop("snapend", cluster.snapend)
    snapdir = kwargs.pop("snapdir", cluster.snapdir)
    skiprows = kwargs.pop("skiprows", cluster.skiprows)

    projected = kwargs.pop("projected", cluster.projected)

    analyze = kwargs.pop("analyze", True)
    sortstars = kwargs.pop("sortstars", True)

    otime = kwargs.pop("otime", False)

    give = kwargs.pop('give','mxv')

    return {
        "nsnap": nsnap,
        "delimiter": delimiter,
        "wdir": wdir,
        "nzfill": nzfill,
        "snapbase": snapbase,
        "snapend": snapend,
        "snapdir": snapdir,
        "skiprows": skiprows,
        "projected": projected,
        "analyze": analyze,
        "sortstars": sortstars,
        "otime": otime,
        "give" : give,
    }  # ,"sfile":sfile,"bfile":bfile}


def _get_cluster_orbit(cluster, ofile, advance=False, **kwargs):
    """ Read in cluster oribit from an ascii file and apply it to StarCluster
    -Columns assumed to be time,x,y,z,vx,vy,vz.

    cluster - class 
        StarCluster to be advanced
    ofile : file
        an already opened file containing orbit information (default: None)
    advance : bool
        Is this a continuation from a previous timestep, in which case read next line (default: False)

    Returns
    -------
    cluster : class
        StarCluster

    Other Parameters
    ----------------
    nsnap : int 
        if nsnap is provided, read line # nsnap from the orbit file
    ounits : str
        if units are not the same as StarCluster units, provide them and they will be converted
    otime : bool
        use time in orbit file to set tphys (default:False)

    Same as load_cluster

    History
    -------
    2018 
    """
    nsnap = int(kwargs.get("nsnap", cluster.nsnap))
    ounits = kwargs.get("ounits", None)
    otime = kwargs.get("otime", False)


    # Read in orbital information from orbit
    if nsnap != 0 and not advance:
        for i in range(0, int(nsnap) + 1):
            data = ofile.readline().split()
    else:
        data = ofile.readline().split()


    tphys = float(data[0])
    xgc = float(data[1])
    ygc = float(data[2])
    zgc = float(data[3])
    vxgc = float(data[4])
    vygc = float(data[5])
    vzgc = float(data[6])

    print(cluster.tphys,otime)

    if cluster.tphys == 0.0 or otime:
        cluster.tphys = tphys

    cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc, ounits)

    return

# Get StarCluster from Gyrfalcon output
def _get_gyrfalcon(
    filein, units="WDunits", origin="galaxy", ofile=None, advance=False, **kwargs
):
    """Extract a single snapshot from an ascii file output from a gyrfalcon simulation

    Parameters
    ----------
    filein : file
        opened gyrfalcon file
    units : str
        units of data (default:'WDunits')
    ofile : file
        opened file containing orbital information
    advance : bool
        is this a snapshot that has been advanced to from initial  load_cluster? (default: False)

    kwargs
    ------

    give : str
        set what parameters are read in from gyrfalcon (default: 'mxv')
        Currently only accepts 'mxvpqael' as an alternative.

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

    if units == "WDunits":
        vcon = 220.0 / bovy_conversion.velocity_in_kpcGyr(220.0, 8.0)
        mcon = 222288.4543021174
        units = "kpckms"
    else:
        vcon = 1.0
        mcon = 1.0

    # Default **kwargs
    skiprows = kwargs.pop("skiprows", 13)
    give = kwargs.get('give','mxv')

    i_d = []
    m = []
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []

    if give == 'mxvpqael':
        gyrpot=[]
        gyrq=[]
        gyracc=[]
        gyreps=[]
        gyrlev=[]
    elif give =='mxve':
        gyreps=[]



    over_head = False
    ntot = 0
    tphys = 0.0

    for j in range(0, skiprows):
        data = filein.readline().split()
        if len(data) == 0:
            print("END OF FILE")
            return StarCluster(0.0,ctype="gyrfalcon",**kwargs)
        elif "#" not in data:
            over_head = True
            print("OVER HEAD")
            break
        if any("Ntot" in dat for dat in data):
            sntot = data[2]
            ntot = int(sntot[:-1])
        if any("time" in dat for dat in data):
            tphys = float(data[2]) * 1000.0

    cluster = StarCluster(
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


        if give == 'mxvpqael':
            gyrpot.append(float(data[7]))
            gyrq.append(float(data[8]))
            gyracc.append(float(data[9]))
            gyreps.append(float(data[10]))
            gyrlev.append(float(data[11]))
        elif give== 'mxve':
            gyreps.append(float(data[7]))



    if ntot > 0:

        cluster.add_stars(x, y, z, vx, vy, vz, m, i_d,sortstars=False)

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

        if give == 'mxvpqael':
            cluster.gyrpot=np.array(gyrpot)
            cluster.gyrq=np.array(gyrq)
            cluster.gyracc=np.array(gyracc)
            cluster.eps=np.array(gyreps)
            cluster.gyrlev=np.array(gyrlev)
        elif give== 'mxve':
            cluster.eps=np.array(gyreps)

    return cluster


def _get_nbody6(out3, out33=None, fort82=None, fort83=None, ofile=None, advance=False, **kwargs):
    """Extract a single snapshot from NBODY6 output

       - Called for Nbody6 simulations with or without stellar evolution

    Parameters
    ----------
    out3 : file
        opened OUT3 file
    out33 : file
        opened OUT33 file containing tail stars (default: None)
    fort82 : file
        opened fort.82 file containing BSE data (default: None)
    fort83 : file
        opened fort.83 file containing SSE data (default: None)
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
    2020 - Written - Webb (UofT)
    """
    
    initialize = kwargs.get("initialize", False)

    if out3 is not None:
    
        ntot,alist,x,y,z,vx,vy,vz,m,i_d=_get_nbody6_out3(out3,**kwargs)
        cluster = StarCluster(
            alist[0],
            units="nbody",
            origin="cluster",
            ctype="nbody6",
            sfile=out3,
        )
                
        if ntot > 0:
            cluster.add_nbody6(
            alist[13], alist[12], alist[2], alist[4], alist[6], alist[7], alist[8], alist[3], alist[11], alist[17], ntot, alist[1], ntot+alist[1]
        )
            cluster.add_stars(x, y, z, vx, vy, vz, m, i_d)


    if out33 is not None:
        cluster.bfile=out33

        ntot,alist,x,y,z,vx,vy,vz,m,i_d=_get_nbody6_out33(out33,**kwargs)
                
        if ntot > 0: 
            cluster.add_stars(x, y, z, vx, vy, vz, m, i_d)
            cluster.add_orbit(alist[0],alist[1],alist[2],alist[3],alist[4],alist[5])

    if fort82 is not None and fort83 is not None:
        cluster.ssefile=fort83
        cluster.bsefile=fort82
        i_d,kw,ri,m1,zl1,r1,te,i_d1,i_d2,kw1,kw2,kwb,rib,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,te1,te2=_get_nbody6se(fort82,fort83,**kwargs)
        cluster.add_sse(kw,zl1,r1)
        cluster.add_bse(i_d1,i_d2,kw1,kw2,kwb,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b)
    
    if kwargs.get("analyze", True) and cluster.ntot>0:
        sortstars=kwargs.get("sortstars", True)
        cluster.analyze(sortstars=sortstars)

        if ofile != None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)
            
    return cluster

def _get_nbody6_out3(f,**kwargs): 

    #Read in header
    try:
        start_header_block_size = struct.unpack('i',f.read(4))[0]
    except:
        return 0,np.zeros(20),0,0,0,0,0,0,0,0
    
    ntot = struct.unpack('i',f.read(4))[0] 
    model = struct.unpack('i',f.read(4))[0] 
    nrun =  struct.unpack('i',f.read(4))[0]
    nk = struct.unpack('i',f.read(4))[0]
    
    end_header_block_size = struct.unpack('i',f.read(4))[0] 

    if start_header_block_size != end_header_block_size:
        print('Error reading OUT3')
        return -1

    # Read in stellar data
    start_data_block_size = struct.unpack('i',f.read(4))[0] #begin data block size

    #Read in alist array from NBODY6
    alist = []
    for i in range(nk):
        alist.append(struct.unpack('f',f.read(4))[0]) #Sverre's 'as'

    #print(alist)
    #Read in masses, positions, velocities, and id's
    m=np.array([])
    x,y,z=np.array([]),np.array([]),np.array([])
    vx,vy,vz=np.array([]),np.array([]),np.array([])
    i_d=np.array([])
 
    for i in range(ntot):
        m=np.append(m,struct.unpack('f',f.read(4))[0])

    #print(m)
    
    for i in range(ntot):           
        x=np.append(x,struct.unpack('f',f.read(4))[0])
        y=np.append(y,struct.unpack('f',f.read(4))[0])
        z=np.append(z,struct.unpack('f',f.read(4))[0]) 

    for i in range(ntot):           
        vx=np.append(vx,struct.unpack('f',f.read(4))[0])
        vy=np.append(vy,struct.unpack('f',f.read(4))[0])
        vz=np.append(vz,struct.unpack('f',f.read(4))[0]) 

    for i in range(ntot):
        i_d=np.append(i_d,struct.unpack('i',f.read(4))[0])

    #print(i_d)
    
    end_data_block_size = struct.unpack('i',f.read(4))[0] #begin data block size

    if start_data_block_size != end_data_block_size:
        print('Error reading OUT3')
        return -1


    return ntot,alist,x,y,z,vx,vy,vz,m,i_d

def _get_nbody6_out33(f,**kwargs): 

    #Read in header
    try:
        start_header_block_size = struct.unpack('i',f.read(4))[0]
    except:
        return 0,np.zeros(20),0,0,0,0,0,0,0,0
    
    ntot = struct.unpack('i',f.read(4))[0] 
    model = struct.unpack('i',f.read(4))[0] 
    nk = struct.unpack('i',f.read(4))[0]
        
    end_header_block_size = struct.unpack('i',f.read(4))[0]

    if start_header_block_size != end_header_block_size:
        print('Error reading OUT33')
        return -1

    if ntot > 0:

        # Read in stellar data
        start_data_block_size = struct.unpack('i',f.read(4))[0] #begin data block size

        #Read in alist array from NBODY6
        alist = []
        for i in range(nk):
            alist.append(struct.unpack('f',f.read(4))[0]) #Sverre's 'as'

        #Read in masses, positions, velocities, and id's
        m=np.array([])
        x,y,z=np.array([]),np.array([]),np.array([])
        vx,vy,vz=np.array([]),np.array([]),np.array([])
        i_d=np.array([])
     
        for i in range(ntot):
            m=np.append(m,struct.unpack('f',f.read(4))[0])

        for i in range(ntot):           
            x=np.append(x,struct.unpack('f',f.read(4))[0])
            y=np.append(y,struct.unpack('f',f.read(4))[0])
            z=np.append(z,struct.unpack('f',f.read(4))[0]) 

        for i in range(ntot):           
            vx=np.append(vx,struct.unpack('f',f.read(4))[0])
            vy=np.append(vy,struct.unpack('f',f.read(4))[0])
            vz=np.append(vz,struct.unpack('f',f.read(4))[0]) 

        for i in range(ntot):
            i_d=np.append(i_d,struct.unpack('i',f.read(4))[0])

        end_data_block_size = struct.unpack('i',f.read(4))[0] #begin data block size
        
        if start_data_block_size != end_data_block_size:
            print('Error reading OUT33')
            return -1

        return ntot,alist,x,y,z,vx,vy,vz,m,i_d
    else:
        return 0,np.zeros(20),0,0,0,0,0,0,0,0

def _get_nbody6se(fort82, fort83, **kwargs):

    header=fort83.readline().split()
    ntot,tphys=int(header[2]),float(header[3])

    i_d=np.array([])
    kw=np.array([])
    ri=np.array([])
    m1=np.array([])
    zl1=np.array([])
    r1=np.array([])
    te=np.array([])

    for i in range(0,ntot):
        data=fort83.readline().split()
        if data[0]=='##': break

        i_d=np.append(i_d,int(data[0]))
        kw=np.append(kw,int(data[1]))
        ri=np.append(ri,float(data[2]))
        m1=np.append(m1,float(data[3]))

        if data[4]=='NaN':
            zl1=np.append(zl1,0.)
            r1=np.append(r1,0.)
            te=np.append(te,0.)
        else:
            zl1=np.append(zl1,float(data[4]))
            r1=np.append(r1,float(data[5]))
            te=np.append(te,float(data[6]))        

    if data[0]!='##':
        data=fort83.readline().split()
   
    header=fort82.readline().split()
    nb,tphys=int(header[2]),float(header[3])

    i_d1=np.array([])
    i_d2=np.array([])
    kw1=np.array([])
    kw2=np.array([])
    kwb=np.array([])
    rib=np.array([])
    ecc=np.array([])
    pb=np.array([])
    semi=np.array([])
    m1b=np.array([])
    m2b=np.array([])
    zl1b=np.array([])
    zl2b=np.array([])
    r1b=np.array([])
    r2b=np.array([])
    te1=np.array([])
    te2=np.array([])

    if nb>0:

        for i in range(0,nb):
            data=fort82.readline().split()
            if data[0]=='##': break

            i_d1=np.append(i_d1,int(data[0]))
            i_d2=np.append(i_d2,int(data[1]))
            kw1=np.append(kw1,int(data[2]))
            kw2=np.append(kw2,int(data[3]))
            kwb=np.append(kwb,int(data[4]))
            rib=np.append(rib,float(data[5]))
            ecc=np.append(ecc,float(data[6]))
            pb=np.append(pb,float(data[7]))
            semi=np.append(semi,float(data[8]))
            m1b=np.append(m1b,float(data[9]))
            m2b=np.append(m2b,float(data[10]))

            if data[11]=='NaN':
                zl1b=np.append(zl1b,0.)
                zl2b=np.append(zl2b,0.)
                r1b=np.append(r1b,0.)
                r2b=np.append(r2b,0.)
                te1=np.append(te1,0.)
                te2=np.append(te2,0.)
            else:
                zl1b=np.append(zl1b,float(data[11]))
                zl2b=np.append(zl2b,float(data[12]))
                r1b=np.append(r1b,float(data[13]))
                r2b=np.append(r2b,float(data[14]))
                te1=np.append(te1,float(data[15]))
                te2=np.append(te2,float(ata[16]))

        if data[0]!='##':
            data=fort82.readline().split()

    else:
        data=fort82.readline().split()

    return i_d,kw,ri,m1,zl1,r1,te,i_d1,i_d2,kw1,kw2,kwb,rib,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,te1,te2

def _get_snapshot(
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
    """Load a generic snapshot where column names and numbers can be manually assigned

    Parameters
    ----------
    ctype : str
        code used to generate data (nbody6/nbody6se/gyrfalcon/...)
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

    cluster = StarCluster(
        tphys, units=units, origin=origin, ctype=ctype, **kwargs
    )
    cluster.add_stars(x, y, z, vx, vy, vz, m, i_d,sortstars=False)
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
            sortstars=kwargs.get("sortstars", True)
            cluster.analyze(sortstars=sortstars)

        if ofile != None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)


    return cluster

def _get_amuse_particles(
    particles, units="kpckms", origin="galaxy", ofile=None, **kwargs
):
    """Convert AMUSE particle dataset to a StarCluster instance

    Parameters
    ----------
    particles : particles
        AMUSE particle dataset
    units : str
        units of input data (default: kpckms)
    origin : str
        origin of input data (default: cluster)
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

    cluster.add_stars(x, y, z, vx, vy, vz, m, i_d, sortstars=False)

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

    return cluster


def _get_astropy_table(
    table,
    column_mapper = None,
    units  = None,
    origin = None,
    ofile = None,
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
        units of input data (default: kpckms)
    origin : str
        origin of input data (default: cluster)
    ofile : file
        opened file containing orbital information

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
        0., units=units, origin=origin, ctype='table', **kwargs
    )

    cluster.add_stars(x, y, z, vx, vy, vz, m, ID, sortstars=False)

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

# /def
