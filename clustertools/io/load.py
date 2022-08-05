""" Read in cluster from Nbody simulations or generate an Nbody cluster

"""

__author__ = "Jeremy J Webb"
__all__ = [
    "load_cluster",
    "advance_cluster",
]

import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

import os, struct
from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from ..util.constants import *

#Import loaders for different code
from .gyrfalcon import _get_gyrfalcon,_get_new_gyrfalcon
from .nbody6pp import _get_nbody6pp
from .nbody6 import _get_nbody6
from .snapshot import _get_snapshot
from .amuse import _get_amuse_particles
from .astropy_table import _get_astropy_table
from .galpydf import _get_galpy_orbits
from .limepydf import _get_limepy

# Try Importing AMUSE. Only necessary for _get_amuse_particles
try:
    import amuse.units.units as u
    from amuse.io import read_set_from_file
except:
    pass

#Try importing hdf5. Only necessary with Nbody6++ and hdf5 output
try: 
    import h5py
except:
    pass

def load_cluster(
    ctype="snapshot",
    units = None,
    origin = None,
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

            - amuse
            - astropy_table
            - limepy
            - nbody6
            - nbody6se
            - nbody6pp or nbody6++
            - nemo or gyrfalcon
            - snapshot

    units : str
        units of input data (default: None)
    origin : str
        origin of input data (default: None)
    ofile : file
        an already opened file containing orbit information (default: None)
    orbit : class
        a galpy orbit to be used for the StarCluster's orbital information (default: None)
    filename : str 
        name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)
    particles : particles
        AMUSE particle dataset (default: None)
        or `~astropy.table.Table` instance if `ctype` is "astropy_table"
        or galpy orbits instance if 'ctype' is 'galpy'
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
        set what parameters are read in from nemo/gyrfalcon (default: 'mxv')
        Currently only accepts 'mxvpqael' as an alternative.
    dtout : integer
        number of nbody timesteps forward to advance to next Nbody6++ timestep (default = 1)

    History
    _______
    2018 - Written - Webb (UofT)
    """
    wdir = kwargs.get("wdir", "./")
    if wdir[-1] != '/':
        wdir+='/'

    #filename=_get_filename(filename,**kwargs)

    initialize = kwargs.get("initialize", False)

    if kwargs.get("ofilename", None) is not None and ofile is None:
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

    elif ctype == "nbody6pp" or ctype=='nbody6++':
        nsnap = kwargs.get("nsnap", 0)
        #Include deltat for legacy reasons
        deltat=kwargs.pop('deltat',1)
        dtout=kwargs.pop('dtout',deltat)
        hdf5=kwargs.pop('hdf5',False)

        if isinstance(dtout,float):
            nsnap=float(nsnap)

        if hdf5:
            if os.path.isfile("%sconf.3_%s" % (wdir,str(nsnap))):
                conf3 = open("%sconf.3_%s" % (wdir,str(nsnap)), "rb")
            else:
                conf3=None

            snap40 = h5py.File("%ssnap.40_%s.h5part" % (wdir,nsnap), "r")
            cluster = _get_nbody6pp(conf3, snap40=snap40, ofile=ofile, advance=False,dtout=dtout,**kwargs)
        else:

            if os.path.isfile("%sconf.3_%s" % (wdir,str(nsnap))):
                conf3 = open("%sconf.3_%s" % (wdir,str(nsnap)), "rb")
            else:
                conf3=None

            if os.path.isfile("%sbev.82_%s" % (wdir,str(nsnap))):
                bev82 = open("%sbev.82_%s" % (wdir,str(nsnap)), "r")
            else:
                bev82=None

            if os.path.isfile("%ssev.83_%s" % (wdir,str(nsnap))):
                sev83 = open("%ssev.83_%s" % (wdir,str(nsnap)), "r")
            else:
                sev83=None

            cluster = _get_nbody6pp(conf3, bev82=bev82, sev83=sev83, ofile=ofile, advance=False,dtout=dtout,**kwargs)


    elif ctype == "gyrfalcon" or ctype=='nemo':
        # Read in snapshot from gyrfalcon.
        filename=_get_filename(filename,**kwargs)
        filein = open(filename, "r")

        cluster = _get_gyrfalcon(filein, units=units, origin=origin, advance=False, **kwargs)

    elif ctype == "new_gyrfalcon" or ctype=='new_nemo':
        # Read in snapshot from gyrfalcon.
        filename=_get_filename(filename,**kwargs)

        cluster = _get_new_gyrfalcon(filename, units=units, origin=origin, advance=False, **kwargs)

    elif ctype=='amuse':
        filename=_get_filename(filename,**kwargs)

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

    elif ctype == 'galpy':

        cluster = _get_galpy_orbits(particles, units=units, origin=origin, ofile=ofile, **kwargs)

    elif ctype == 'limepy':
        cluster = _get_limepy(units=units, origin=origin, orbit=orbit, ofile=ofile, **kwargs)

    else:
        print("NO CTYPE GIVEN")
        return 0

    if ofile is not None:
        cluster.ofilename=ofile.name.split('/')[-1]
        cluster.ofile=ofile

    # Add galpy orbit if given
    if orbit is not None:
        cluster.orbit = orbit
        if cluster.units == "pckms":
            t = (cluster.tphys / 1000.0) / conversion.time_in_Gyr(ro=solar_ro, vo=solar_vo)
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
            t = (cluster.tphys * cluster.tbar / 1000.0) / conversion.time_in_Gyr(
                ro=solar_ro, vo=solar_vo
            )
            cluster.add_orbit(
                orbit.x(t) * 1000.0 / cluster.rbar,
                orbit.y(t) * 1000.0 / cluster.rbar,
                orbit.z(t) * 1000.0 / cluster.rbar,
                orbit.vx(t) / cluster.vbar,
                orbit.vy(t) / cluster.vbar,
                orbit.vz(t) / cluster.vbar,
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

        cluster.save_cluster()
        units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

        if cluster.origin is not None: 
            cluster.to_cluster(sortstars=False)
            cluster.find_centre()
            cluster.return_cluster(units0,origin0, rorder0, rorder_origin0 )
        else:
            cluster.find_centre()
    elif initialize:
        origin0=cluster.origin
        cluster.to_galaxy()
        solarmotion=kwargs.get('solarmotion',[-11.1, 24.0, 7.25])
        initialize_orbit(cluster,solarmotion=solarmotion)
        cluster.to_origin(origin0)
        
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
    Same as load_cluster except for:

    dtout : integer
        number of nbody timesteps forward to advance to next Nbody6++ timestep (default = 1)

    History
    -------
    2018 - Written - Webb (UofT)
    """
    nsnap=None
    advance_kwargs, kwargs = _get_advanced_kwargs(cluster, **kwargs)
    #filename=_get_filename(filename,**advance_kwargs)

    if ofile is None:
        ofile=cluster.ofile

    if kwargs.get("ofilename", None) is None:
        ofilename=cluster.ofilename

    wdir = advance_kwargs.get("wdir", "./")
    if wdir[-1] != '/':
        wdir+='/'

    # Continue reading in cluster opened in _get_cluster()
    if load_function is not None:
        ctype='custom'
        if filename is not None:
            cluster=load_function(ctype=ctype,units=cluster.units_init,origin=cluster.origin_init,ofile=ofile,orbit=orbit,filename=filename,advance=True,**advance_kwargs,**kwargs)
        else:
            cluster=load_function(ctype=ctype,units=cluster.units_init,origin=cluster.origin_init,ofile=ofile,orbit=orbit,advance=True,**advance_kwargs,**kwargs)


    elif cluster.ctype == "nbody6":
        cluster = _get_nbody6(
            cluster.sfile, cluster.bfile, cluster.bsefile, cluster.ssefile, ofile=ofile, advance=True, **advance_kwargs, **kwargs
        )

    elif cluster.ctype == "nbody6pp" or cluster.ctype == "nbody6++":
        

        hdf5=advance_kwargs.pop("hdf5")

        if hdf5:
            ngroups=advance_kwargs.get("ngroups")
            ngroup=advance_kwargs.get("ngroup")

            if ngroup < ngroups:

                nsnap = advance_kwargs.pop("nsnap") - 1

                if cluster.nc!=0. and cluster.n_p!=0.:

                    units0=cluster.units
                    if units0!='nbody': cluster.to_nbody()

                    nc,rc,rbar,rtide=cluster.nc,cluster.rc,cluster.rbar,cluster.rtide
                    xc,yc,zc=cluster.xc,cluster.yc,cluster.zc
                    zmbar,vbar,tbar,rscale=cluster.zmbar,cluster.vbar,cluster.tbar,cluster.rscale
                    ns,nb,n_p=cluster.ns,cluster.nb,cluster.n_p

                    nbody6list=[nc,rc,rbar,rtide,xc,yc,zc,zmbar,vbar,tbar,rscale,ns,nb,n_p]

                    if units0!='nbody': cluster.to_units(units0)

                else:
                    nbody6list=None

                conf3=None
                cluster = _get_nbody6pp(conf3, snap40=cluster.sfile, ofile=ofile, advance=True,nbody6list=nbody6list,nsnap=nsnap,**advance_kwargs, **kwargs)
                
            else:
                deltat=kwargs.pop('deltat',1)
                dtout=kwargs.pop('dtout',deltat)

                nsnap = advance_kwargs.pop("nsnap") + dtout - 1
                ngroup= advance_kwargs.pop("ngroup")

                if os.path.isfile("%sconf.3_%s" % (wdir,str(nsnap))):
                    conf3 = open("%sconf.3_%s" % (wdir,str(nsnap)), "rb")
                else:
                    conf3=None

                snap40 = h5py.File("%ssnap.40_%s.h5part" % (wdir,nsnap), "r")
                cluster = _get_nbody6pp(conf3, snap40=snap40, ofile=ofile, advance=True,nsnap=nsnap,dtout=dtout,**advance_kwargs, **kwargs)

        else:
            deltat=kwargs.pop('deltat',1)
            dtout=kwargs.pop('dtout',deltat)

            nsnap = advance_kwargs.pop("nsnap") + dtout - 1

            if os.path.isfile("%sconf.3_%s" % (wdir,str(nsnap))):
                conf3 = open("%sconf.3_%s" % (wdir,str(nsnap)), "rb")
            else:
                conf3=None

            if os.path.isfile("%sbev.82_%s" % (wdir,str(nsnap))):
                bev82 = open("%sbev.82_%s" % (wdir,str(nsnap)), "r")
            else:
                bev82=None

            if os.path.isfile("%ssev.83_%s" % (wdir,str(nsnap))):
                sev83 = open("%ssev.83_%s" % (wdir,str(nsnap)), "r")
            else:
                sev83=None

            cluster = _get_nbody6pp(conf3, bev82=bev82, sev83=sev83, ofile=ofile, advance=True,nsnap=nsnap,dtout=dtout,**advance_kwargs, **kwargs)

    elif cluster.ctype == 'nbody6':

        cluster = _get_nbody6(
            cluster.sfile, cluster.bfile, cluster.bsefile, cluster.ssefile, ofile=ofile, advance=True, **advance_kwargs, **kwargs
        )

    elif cluster.ctype == "gyrfalcon" or cluster.ctype=="nemo":

        filename=_get_filename(filename,**advance_kwargs)

        if filename is None:
            cluster = _get_gyrfalcon(
                cluster.sfile,
                units=cluster.units_init,
                origin=cluster.origin_init,
                ofile=ofile,
                advance=True,
                **advance_kwargs,
                **kwargs
            )

        else:
            filein = open(filename, "r")
            cluster = _get_gyrfalcon(
                filein,
                units=cluster.units_init,
                origin=cluster.origin_init,
                ofile=ofile,
                advance=True,
                **advance_kwargs,
                **kwargs
            )

    elif cluster.ctype == "new_gyrfalcon" or cluster.ctype=="new_nemo":

        filename=_get_filename(filename,**advance_kwargs)

        if filename is None:
            cluster = _get_new_gyrfalcon(
                cluster.sfile,
                units=cluster.units_init,
                origin=cluster.origin_init,
                ofile=ofile,
                advance=True,
                **advance_kwargs,
                **kwargs
            )

        else:
            cluster = _get_new_gyrfalcon(
                filename,
                units=cluster.units_init,
                origin=cluster.origin_init,
                ofile=ofile,
                advance=True,
                **advance_kwargs,
                **kwargs
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
            **advance_kwargs,
            **kwargs
        )
    elif ctype=='amuse':
        filename=_get_filename(filename,**advance_kwargs)

        if filename is not None:
            filetype=kwargs.pop("filetype","hdf5")
            particles = read_set_from_file(filename, filetype)
            cluster = _get_amuse_particles(particles, units=cluster.units_init, origin=cluster.origin_init, ofile=ofile,**advance_kwargs,**kwargs)

    else:
        cluster = StarCluster(ctype=cluster.ctype,units=cluster.units_init,origin=cluster.origin_init,ofile=ofile,**advance_kwargs, **kwargs)


    # Check for restart
    if cluster.ntot == 0.0:
        try:
            wdir = cluster.wdir + "cont/"
        except:
            wdir = "./cont/"


        if os.path.exists(wdir):

            if ofile is not None:
                ofilename = ofile.name.split('/')[-1]
                kwargs.pop('ofilename')
                ofile=None
            else:
                ofilename=kwargs.pop('ofilename',ofilename)

            old_wdir=advance_kwargs.pop('wdir')

            if nsnap is None:

                cluster = load_cluster(
                    ctype=cluster.ctype,units=cluster.units_init,origin=cluster.origin_init,orbit=orbit,filename=filename,load_function=load_function,wdir=wdir,ofilename=ofilename,**advance_kwargs, **kwargs,
                )

            else:
                   cluster = load_cluster(
                    ctype=cluster.ctype,units=cluster.units_init,origin=cluster.origin_init,orbit=orbit,filename=filename,nsnap=nsnap,load_function=load_function,wdir=wdir,ofilename=ofilename,**advance_kwargs, **kwargs,
                )         

    if cluster.ntot != 0.0:

        if ofile is not None:
            cluster.ofilename=ofile.name.split('/')[-1]
            cluster.ofile=ofile

        # Add galpy orbit if given
        if orbit != None:
            cluster.orbit = orbit
            if cluster.units == "pckms" or cluster.units == "kpckms":
                t = (cluster.tphys / 1000.0) / conversion.time_in_Gyr(
                    ro=solar_ro, vo=solar_vo
                )
            elif cluster.units == "nbody":
                t = (
                    cluster.tphys * cluster.tbar / 1000.0
                ) / conversion.time_in_Gyr(ro=solar_ro, vo=solar_vo)
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

def _get_filename(filename,**kwargs):
    """assemble filename from **kwargs

    Parameters
    ----------
    filename : str or None
       given filename to read in cluster data

    Returns
    -------
    filename : str

    Other Parameters
    ----------------
    Same as load_cluster

    History
    -------
    2021 - Written - Webb (UofT)
    """

    nzfill = int(kwargs.get("nzfill", 1))
    nsnap = int(kwargs.get("nsnap", "0"))
    wdir = kwargs.get("wdir", "./")
    snapdir = kwargs.get("snapdir", "snaps/")
    snapbase = kwargs.get("snapbase", "")
    snapend = kwargs.get("snapend", ".dat")

    if filename != None:
        if os.path.isfile(filename):
            pass
        elif os.path.isfile("%s%s%s" % (wdir, snapdir, filename)):
            filename="%s%s%s" % (wdir, snapdir, filename)
        elif os.path.isfile("%s%s" % (wdir, filename)):
            filename="%s%s" % (wdir, filename)

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
        filename = None

    return filename


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
    dt = kwargs.pop('dt',cluster.dt)
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

    hdf5 = kwargs.pop('hdf5',cluster.hdf5)
    ngroups = kwargs.pop('ngroups',cluster.ngroups)
    ngroup = np.maximum(int(kwargs.pop("ngroup", 0)), cluster.ngroup) + 1

    verbose = kwargs.pop('verbose',True)

    return {
        "nsnap": nsnap,
        "dt": dt,
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
        "hdf5" : hdf5,
        "ngroups" : ngroups,
        "ngroup" : ngroup,
        "verbose" : verbose
    }, kwargs

# /def
