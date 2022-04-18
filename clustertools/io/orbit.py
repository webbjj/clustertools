import numpy as np
import os, struct

def _get_cluster_orbit(cluster, ofile, advance=False, ocol_names=["t", "x", "y", "z", "vx", "vy", "vz"],ocol_nums=[0, 1, 2, 3, 4, 5, 6], **kwargs):
    """ Read in cluster oribit from an ascii file and apply it to StarCluster

    cluster - class 
        StarCluster to be advanced
    ofile : file
        an already opened file containing orbit information (default: None)
    advance : bool
        Is this a continuation from a previous timestep, in which case read next line (default: False)
    col_names : str
        names corresponding to time, position, and velocity
    col_nums : int
        column numbers corresponding to each column name
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


    #Testing 
    if True:
        tphys,xgc,ygc,zgc,vxgc,vygc,vzgc=0.,0.,0.,0.,0.,0.,0.

        for i in range(0,len(ocol_names)):
            if ocol_names[i]=="t":
                tphys=float(data[ocol_nums[i]])
            elif ocol_names[i]=="x":
                xgc=float(data[ocol_nums[i]])
            elif ocol_names[i]=="y":
                ygc=float(data[ocol_nums[i]])
            elif ocol_names[i]=="z":
                zgc=float(data[ocol_nums[i]])
            elif ocol_names[i]=="vx":
                vxgc=float(data[ocol_nums[i]])
            elif ocol_names[i]=="vy":
                vygc=float(data[ocol_nums[i]])
            elif ocol_names[i]=="vz":
                vzgc=float(data[ocol_nums[i]])
    else:
        tphys = float(data[0])
        xgc = float(data[1])
        ygc = float(data[2])
        zgc = float(data[3])
        vxgc = float(data[4])
        vygc = float(data[5])
        vzgc = float(data[6])

    if otime:
        cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc, ounits , tphys=tphys)
    else:
        cluster.add_orbit(xgc, ygc, zgc, vxgc, vygc, vzgc, ounits)

    return
