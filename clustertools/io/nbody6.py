import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

import os, struct
from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from .orbit import _get_cluster_orbit

def _get_nbody6(out3, out33=None, fort82=None, fort83=None, ofile=None, advance=False, **kwargs):
    """Extract a single snapshot from NBODY6 output

       - Called for Nbody6 simulations with or without stellar evolution
       - Note that binary centre of masses are loaded in the position, velocity and mass arrays


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

    ntot=0

    if out3 is not None:
    
        nstot,alist,xs,ys,zs,vxs,vys,vzs,ms,i_ds=_get_nbody6_out3(out3,**kwargs)

        cluster = StarCluster(
            alist[0],
            units="nbody",
            origin="cluster",
            ctype="nbody6",
            sfile=out3,
        )
                
        if nstot > 0:
            cluster.add_nbody6(
            alist[13], alist[12], alist[2], alist[4], alist[6], alist[7], alist[8], alist[3], alist[11],alist[10],alist[17], nstot, alist[1], nstot+alist[1]
        )
            ntot+=nstot

    #Add binaries first
    if out33 is not None:
        cluster.bfile=out33

        nbtot,alist,xb,yb,zb,vxb,vyb,vzb,mb,i_db=_get_nbody6_out33(out33,**kwargs)
                
        if nbtot > 0: 
            cluster.add_stars(xb, yb, zb, vxb, vyb, vzb, mb, i_db)
            cluster.add_orbit(alist[0],alist[1],alist[2],alist[3],alist[4],alist[5])
            ntot+=nbtot

    if nstot > 0:
        cluster.add_stars(xs, ys, zs, vxs, vys, vzs, ms, i_ds)

    if fort82 is not None and fort83 is not None:
        cluster.ssefile=fort83
        cluster.bsefile=fort82
        i_d,kw,ri,m1,zl1,r1,te,i_d1,i_d2,kw1,kw2,kwb,rib,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,te1,te2=_get_nbody6se(fort82,fort83,**kwargs)
        cluster.add_sse(kw,zl1,r1)

        pb=(10.0**pb)/cluster.tbar_days
        semi=(10.0**semi)/cluster.rbar_su
        m1b/=cluster.zmbar
        m2b/=cluster.zmbar

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

    #Single star arrays
    i_d=np.array([])
    kw=np.array([])
    ri=np.array([])
    m1=np.array([])
    zl1=np.array([])
    r1=np.array([])
    te=np.array([])

    #Binary star arrays
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

    #read binaries first
    header=fort82.readline().split()
    nb,tphys=int(header[2]),float(header[3])

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

            #Add select parametes to full array
            kw=np.append(kw,kwb[-1])
            lbtot=np.log10(10.0**np.asarray(zl1b)+10.0**np.asarray(zl2b))
            zl1=np.append(zl1,lbtot)
            r1=np.append(r1,np.maximum(r1b[-1],r2b[-1]))

        if data[0]!='##':
            data=fort82.readline().split()

    else:
        data=fort82.readline().split()


    #Read single stars second
    header=fort83.readline().split()
    ntot,tphys=int(header[2]),float(header[3])

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

    return i_d,kw,ri,m1,zl1,r1,te,i_d1,i_d2,kw1,kw2,kwb,rib,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,te1,te2


