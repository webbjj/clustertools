import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from .orbit import _get_cluster_orbit
from ..util.constants import *

#Get StarCluster from Gyrfalcon output
def _get_gyrfalcon(
    filein, units=None, origin=None, ofile=None, advance=False, **kwargs
):
    """Extract a single snapshot from an ascii file output from a gyrfalcon simulation

    Parameters
    ----------
    filein : file
        opened nemo/gyrfalcon file
    units : str
        units of input data (default: None)
    origin : str
        origin of input data (default: None)
    ofile : file
        opened file containing orbital information
    advance : bool
        is this a snapshot that has been advanced to from initial  load_cluster? (default: False)

    kwargs
    ------

    give : str
        set what parameters are read in from nemo/gyrfalcon (default: 'mxv')
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
            return StarCluster(0.0,ctype="nemo",**kwargs)
        elif "#" not in data:
            over_head = True
            print("OVER HEAD")
            break
        if any("Ntot" in dat for dat in data):
            sntot = data[2]
            ntot = int(sntot[:-1])
        if any("time" in dat for dat in data):
            tphys = float(data[2])

    cluster = StarCluster(
        tphys,
        units=units,
        origin=origin,
        ctype="nemo",
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
        m.append(float(data[0]))
        x.append(float(data[1]))
        y.append(float(data[2]))
        z.append(float(data[3]))
        vx.append(float(data[4]))
        vy.append(float(data[5]))
        vz.append(float(data[6]))


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

        if ofile is not None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)
        else:
            cluster.find_centre()

        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            if cluster.origin is not None:
                cluster.to_cluster(sortstars=False)
                cluster.find_centre()
                cluster.to_centre(sortstars=sortstars)
                cluster.to_origin(origin)
            else:
                cluster.find_centre()
                cluster.analyze(sortstars=sortstars)

        if give == 'mxvpqael':
            cluster.gyrpot=np.array(gyrpot)
            cluster.gyrq=np.array(gyrq)
            cluster.gyracc=np.array(gyracc)
            cluster.eps=np.array(gyreps)
            cluster.gyrlev=np.array(gyrlev)
        elif give== 'mxve':
            cluster.eps=np.array(gyreps)


    return cluster

def _get_new_gyrfalcon(
    filename, units=None, origin=None, ofile=None, advance=False, **kwargs
):
    """Extract a single snapshot from an ascii file output from a gyrfalcon simulation

    Parameters
    ----------
    filename : file
        name of nemo/gyrfalcon file
    units : str
        units of input data (default: None)
    origin : str
        origin of input data (default: None)
    ofile : file
        opened file containing orbital information
    advance : bool
        is this a snapshot that has been advanced to from initial  load_cluster? (default: False)

    kwargs
    ------

    give : str
        set what parameters are read in from nemo/gyrfalcon (default: 'mxv')
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

    # Default **kwargs
    skiprows = kwargs.pop("skiprows", 13)
    give = kwargs.get('give','mxv')

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

    filein = open(filename, "r")

    for j in range(0, skiprows):
        data = filein.readline().split()
        if len(data) == 0:
            print("END OF FILE")
            return StarCluster(0.0,ctype="nemo",**kwargs)
        elif "#" not in data:
            over_head = True
            print("OVER HEAD")
            break
        if any("Ntot" in dat for dat in data):
            sntot = data[2]
            ntot = int(sntot[:-1])
        if any("time" in dat for dat in data):
            tphys = float(data[2])

    if advance and kwargs.get('dt',None) is None:

        for j in range(0,ntot):
            filein.readline()

        for j in range(0, skiprows):
            data = filein.readline().split()
            if len(data) == 0:
                print("END OF FILE")
                return StarCluster(0.0,ctype="nemo",**kwargs)
            elif "#" not in data:
                over_head = True
                print("OVER HEAD")
                break
            if any("time" in dat for dat in data):
                tphys2 = float(data[2])

        dtout=tphys2-tphys
    elif advance:
        dtout=kwargs.get('dt')

    filein.close()

    cluster = StarCluster(
        tphys,
        units=units,
        origin=origin,
        ctype="new_nemo",
        sfile=filename,
        bfile=None,
        skiprows=skiprows,
        **kwargs
    )

    if advance: 
        cluster.tphys+=dtout*cluster.nsnap
        cluster.dt=dtout


    if not advance:
        nskip=skiprows
    else:
        nskip=(skiprows+ntot)*cluster.nsnap+skiprows

    try:
        data=np.loadtxt(filename,skiprows=nskip,max_rows=ntot)
    except StopIteration:
        print("END OF FILE")
        return StarCluster(0.0,ctype="new_nemo",**kwargs)

    m=data[:,0].astype(float)
    x=data[:,1].astype(float)
    y=data[:,2].astype(float)
    z=data[:,3].astype(float)
    vx=data[:,4].astype(float)
    vy=data[:,5].astype(float)
    vz=data[:,6].astype(float)

    i_d=np.arange(0,len(m),1)

    if give == 'mxvpqael':
        gyrpot=data[:,7].astype(float)
        gyrq=data[:,8].astype(float)
        gyracc=data[:,9].astype(float)
        gyreps=data[:,10].astype(float)
        gyrlev=data[:,11].astype(float)
    elif give== 'mxve':
        gyreps=data[:,7].astype(float)



    if ntot > 0:

        cluster.add_stars(x, y, z, vx, vy, vz, m, i_d,sortstars=False)

        if ofile is not None:
            _get_cluster_orbit(cluster, ofile, advance=advance, **kwargs)
        else:
            cluster.find_centre()

        if kwargs.get("analyze", True):
            sortstars=kwargs.get("sortstars", True)
            if cluster.origin is not None:
                cluster.to_cluster(sortstars=False)
                cluster.find_centre()
                cluster.to_centre(sortstars=sortstars)
                cluster.to_origin(origin)
            else:
                cluster.find_centre()
                cluster.analyze(sortstars=sortstars)

        if give == 'mxvpqael':
            cluster.gyrpot=np.array(gyrpot)
            cluster.gyrq=np.array(gyrq)
            cluster.gyracc=np.array(gyracc)
            cluster.eps=np.array(gyreps)
            cluster.gyrlev=np.array(gyrlev)
        elif give== 'mxve':
            cluster.eps=np.array(gyreps)

    return cluster
