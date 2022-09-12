import numpy as np
try:
    from galpy.util import conversion
except:
    import galpy.util.bovy_conversion as conversion

import os, struct
from ..cluster.cluster import StarCluster
from ..analysis.orbits import initialize_orbit
from .orbit import _get_cluster_orbit

#Try importing hdf5. Only necessary with Nbody6++ and hdf5 output
try: 
    import h5py
except:
    pass

def _get_nbody6pp(conf3, bev82=None, sev83=None, snap40=None, ofile=None, advance=False, nbody6list=None,  **kwargs):
    """Extract a single snapshot from NBODY6++ output

       - Note that for snap40=False, individual binary stars are loaded in the main position, velocity, and mass arrays
       - for snap40=True, binary centre of masses are loaded in the position, velocity and mass arrays as per nbody6

    Parameters
    ----------
    conf3 : file
        opened conf3 file
    bev82 : file
        opened bev82 file containing BSE data (default: None)
    sev83 : file
        opened sev83 file containing SSE data (default: None)
    snap40 : file
        opened snap40 file containing hdf5 data (default: None)
    ofile : file
        opened file containing orbital information
    advance : bool
        is this a snapshot that has been advanced to from initial  load_cluster? (default: False)
    nbody6list : float
        array of nbody6 parameters that need to get passed during an advance

    Returns
    -------
    cluster : class
        StarCluster

    Other Parameters
    ----------------
    Same as load_cluster

    History
    -------
    2021 - Written - Webb (UofT)
    """

    verbose=kwargs.get('verbose',True)


    initialize = kwargs.get("initialize", False)
    nsnap = kwargs.pop("nsnap", 0)
    wdir = kwargs.get("wdir", './')
    deltat=kwargs.get('deltat',1)
    dtout=kwargs.get('dtout',deltat)



    if snap40 is not None:
        ngroup=kwargs.pop('ngroup',0)
        tphys,ntot,x,y,z,vx,vy,vz,m,i_d,pot,kw,lum,rc,rs,te,binaries=_get_nbody6pp_hdf5(snap40,ngroup=ngroup,**kwargs)

        if binaries:
            bdata=_get_nbody6pp_hdf5_binaries(snap40,ngroup=ngroup,**kwargs)
            semi,ecc,gb,kw1,kw2,kwb,zl1b,zl2b,m1b,m2b,mc1,mc2,i_d1,i_d2,idc,pb,potb,rc1,rc2,r1b,r2b,te1,te2,vc1,vc2,vc3,vr1,vr2,vr3,xc1,xc2,xc3,xr1,xr2,xr3=bdata
            mbtot=np.asarray(m1b)+np.asarray(m2b)
            lbtot=np.log10(10.0**np.asarray(zl1b)+10.0**np.asarray(zl2b))
            nb=len(semi)
        else:
            nb=0

        cluster = StarCluster(
            tphys,
            units="pckms",
            origin="cluster",
            ctype="nbody6++",
            sfile=snap40,
            nsnap=nsnap,
            wdir=wdir,
        )

        cluster.hdf5=True
        cluster.ngroups=len(snap40)
        cluster.ngroup=ngroup
        
        if conf3 is not None:
            alist=_get_nbody6pp_conf3(conf3,nsnap=nsnap,return_alist_only=True,**kwargs)
            cluster.add_nbody6(
            alist[13], alist[12], alist[2], alist[4], alist[6], alist[7], alist[8], alist[3], alist[11],alist[10],alist[17], ntot, nb, ntot+alist[1])
        elif nbody6list is not None:
            cluster.add_nbody6(nbody6list)

        if binaries:
            cluster.add_stars(xc1,xc2,xc3,vc1,vc2,vc3,mbtot,i_d1,analyze=False)
            pb=np.array(pb)
            semi=np.array(semi)
            m1b=np.array(m1b)
            m2b=np.array(m2b)

            cluster.add_bse(i_d1,i_d2,kw1,kw2,kwb,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,ep1=te1,ep2=te2)
            cluster.bunits='audays'

            cluster.add_sse(kw1,lbtot,np.maximum(r1b,r2b))

            bm1, bm2 = m1b, m2b
            bm = bm1 + bm2
            bmr1 =  bm2 / bm
            bmr2 = -bm1 / bm
            bmr1a = bmr1 / 206265.0
            bmr2a = bmr2 / 206265.0
            #xb1, yb1, zb1 = xc1 + bmr1a * xr1, xc2 + bmr1a * xr2, xc3 + bmr1a * xr3
            #vxb1, vyb1, vzb1 = vc1 + bmr1 * vr1,  vc2 + bmr1 * vr2,  vc3 + bmr1 * vr3
            #xb2, yb2, zb2 = xc1 + bmr2a * xr1, xc2 + bmr2a * xr2, xc3 + bmr2a * xr3
            #vxb2, vyb2, vzb2  = vc1 + bmr2 * vr1,  vc2 + bmr2 * vr2,  vc3 + bmr2 * vr3

            xb1, yb1, zb1 =  bmr1 * xr1, bmr1 * xr2, bmr1 * xr3
            vxb1, vyb1, vzb1 = bmr1 * vr1,  bmr1 * vr2, bmr1 * vr3
            xb2, yb2, zb2 = bmr2 * xr1, bmr2 * xr2, bmr2 * xr3
            vxb2, vyb2, vzb2  = bmr2 * vr1,  bmr2 * vr2,  bmr2 * vr3

            cluster.add_binary_stars(xb1, yb1, zb1, vxb1, vyb1, vzb1, xb2, yb2, zb2, vxb2, vyb2, vzb2,set_com=False)


        cluster.add_stars(x, y, z, vx, vy, vz, m, i_d, analyze=False)
        cluster.add_sse(kw,lum,rs)

        if binaries: 
            cluster.pot=np.append(potb,pot)
        else:
            cluster.pot=pot

        if conf3 is not None or nbody6list is not None:
            cluster.xc*=cluster.rbar
            cluster.yc*=cluster.rbar
            cluster.zc*=cluster.rbar
            cluster.tphys*=cluster.tbar

        elif not advance:
            cluster.reset_nbody_scale(rvirial=True)
            cluster.xc*=cluster.rbar
            cluster.yc*=cluster.rbar
            cluster.zc*=cluster.rbar
            cluster.tphys*=cluster.tbar

        cluster.to_nbody(analyze=False)

        if binaries: cluster.nb = len(semi)


    else:
        ntot,alist,x,y,z,vx,vy,vz,m,i_d,rhos,xns,pot=_get_nbody6pp_conf3(conf3,nsnap=nsnap,**kwargs)

        cluster = StarCluster(
            alist[0],
            units="nbody",
            origin="cluster",
            ctype="nbody6++",
            sfile=conf3,
            nsnap=nsnap,
            wdir=wdir,
        )

        if ntot > 0:
            cluster.add_nbody6(
            alist[13], alist[12], alist[2], alist[4], alist[6], alist[7], alist[8], alist[3], alist[11],alist[10],alist[17], ntot, alist[1], ntot+alist[1]
        )
            nb=int(alist[1])
            cluster.add_stars(x, y, z, vx, vy, vz, m, i_d,nb=nb,analyze=False)
            cluster.rhos=np.zeros(ntot-nb)
            pots=np.zeros(ntot-nb)

            if nb>0:
                binargs=np.arange(0,nb,1)
                binargs1=np.arange(0,int(2*nb),2)
                binargs2=binargs1+1
                cluster.rhos[binargs]=rhos[binargs1]+rhos[binargs2]
                pots[binargs]=pot[binargs1]+pot[binargs2]

            cluster.rhos[nb:]=rhos[2*nb:]
            pots[nb:]=pot[2*nb:]

            v=np.sqrt(cluster.vx**2.+cluster.vy**2.+cluster.vz**2.)
            ek=0.5*cluster.m*v**2.
            cluster.add_energies(ek,pots)

        if bev82 is not None and sev83 is not None:

            kws=np.ones(ntot-nb)*-10.0
            zl1s=np.ones(ntot-nb)*-10.0
            r1s=np.ones(ntot-nb)*-10.0

            arg,i_d,kw,ri,m1,zl1,r1,te,i_d1,i_d2,kw1,kw2,kwb,rib,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,te1,te2=_get_nbody6pp_ev(bev82,sev83,nsnap=nsnap,**kwargs)
            #Convert from fortran array address to python

            if nb>0:
                kws[binargs]=np.maximum(kw1,kw2)
                zl1s[binargs]=zl1b+zl2b
                r1s[binargs]=r1b+r2b
                args=arg[2*nb:]-nb-1
            else:
                args=arg-1

            args=args.astype(int)

            kws[args]=kw[2*nb:]
            zl1s[args]=zl1[2*nb:]
            r1s[args]=r1[2*nb:]

            cluster.add_sse(kws,zl1s,r1s)

            pb=(10.0**pb)/cluster.tbar_days
            semi=(10.0**semi)/cluster.rbar_su
            m1b/=cluster.zmbar
            m2b/=cluster.zmbar

            cluster.add_bse(i_d1,i_d2,kw1,kw2,kwb,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b)


    if np.sum(cluster.kw<0) >0:
        if verbose: print('SSE/BSE NBODY6++ ERROR',np.sum(cluster.kw<0))

    if kwargs.get("analyze", True) and cluster.ntot>0:
        sortstars=kwargs.get("sortstars", True)
        cluster.analyze(sortstars=sortstars)

    if ofile != None:
        _get_cluster_orbit(cluster, ofile, advance=advance, nsnap=int(nsnap/dtout),**kwargs)
           

    return cluster

def _get_nbody6pp_conf3(f,return_alist_only=False,**kwargs): 

    verbose=kwargs.get('verbose',True)

    #Read in header
    try:
        start_header_block_size = struct.unpack('i',f.read(4))[0]
    except:
        return 0,np.zeros(20),0,0,0,0,0,0,0,0,0,0,0
        
    ntot = struct.unpack('i',f.read(4))[0] 
    model = struct.unpack('i',f.read(4))[0] 
    nrun = struct.unpack('i',f.read(4))[0]
    nk = struct.unpack('i',f.read(4))[0]
             
    end_header_block_size = struct.unpack('i',f.read(4))[0]
    
    if start_header_block_size != end_header_block_size:
        if verbose: print('Error reading CONF3')
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
        rhos=np.array([])
        xns=np.array([])
        x,y,z=np.array([]),np.array([]),np.array([])
        vx,vy,vz=np.array([]),np.array([]),np.array([])
        phi=np.array([])
        i_d=np.array([])
     
        for i in range(ntot):
            m=np.append(m,struct.unpack('f',f.read(4))[0])

        for i in range(ntot):
            rhos=np.append(rhos,struct.unpack('f',f.read(4))[0])
            
        for i in range(ntot):
            xns=np.append(xns,struct.unpack('f',f.read(4))[0])

        for i in range(ntot):           
            x=np.append(x,struct.unpack('f',f.read(4))[0])
            y=np.append(y,struct.unpack('f',f.read(4))[0])
            z=np.append(z,struct.unpack('f',f.read(4))[0]) 

        for i in range(ntot):           
            vx=np.append(vx,struct.unpack('f',f.read(4))[0])
            vy=np.append(vy,struct.unpack('f',f.read(4))[0])
            vz=np.append(vz,struct.unpack('f',f.read(4))[0]) 

        for i in range(ntot):
            phi=np.append(phi,struct.unpack('i',f.read(4))[0])            

        for i in range(ntot):
            i_d=np.append(i_d,struct.unpack('i',f.read(4))[0])

        end_data_block_size = struct.unpack('i',f.read(4))[0] #begin data block size
        
        if start_data_block_size != end_data_block_size:
            if verbose: print('Error reading CONF3')
            return -1

        if return_alist_only:
            return alist
        else:
            return ntot,alist,x,y,z,vx,vy,vz,m,i_d,rhos,xns,phi
    else:
        return 0,np.zeros(20),0,0,0,0,0,0,0,0,0,0,0

def _get_nbody6pp_ev(bev, sev, **kwargs):
    
    verbose=kwargs.get('verbose',True)

    arg=np.array([])
    i_d=np.array([])
    kw=np.array([])
    ri=np.array([])
    m1=np.array([])
    zl1=np.array([])
    r1=np.array([])
    te=np.array([])


    #Read in binary data first 
  
    header=bev.readline().split()
    nb,tphys=int(header[0]),float(header[1])

    if nb>0:

        data=np.loadtxt(bev.name,skiprows=1)

        if nb==1:
            data=data.reshape(1,len(data))

        arg1=data[:,1].astype(int)
        arg2=data[:,2].astype(int)
        i_d1=data[:,3].astype(int)
        i_d2=data[:,4].astype(int)
        kw1=data[:,5].astype(int)
        kw2=data[:,6].astype(int)
        kwb=data[:,7].astype(int)
        rib=data[:,8].astype(float)
        ecc=data[:,9].astype(float)
        pb=data[:,10].astype(float)
        semi=data[:,11].astype(float)
        m1b=data[:,12].astype(float)
        m2b=data[:,13].astype(float)

        if len(data)>=14:
            zl1b=data[:,14].astype(float)
            zl2b=data[:,15].astype(float)
            r1b=data[:,16].astype(float)
            r2b=data[:,17].astype(float)
            te1=data[:,18].astype(float)
            te2=data[:,19].astype(float)
        else:
            zl1b=np.zeros(len(i_d1))
            zl2b=np.zeros(len(i_d1))
            r1b=np.zeros(len(i_d1))
            r2b=np.zeros(len(i_d1))
            te1=np.zeros(len(i_d1))
            te2=np.zeros(len(i_d1))

        argbs=np.array([])
        idbs=np.array([])
        kwbs=np.array([])
        ribs=np.array([])
        m1bs=np.array([])
        zl1bs=np.array([])
        r1bs=np.array([])
        tebs=np.array([])

        for i in range(0,len(arg1)):
            argbs=np.append(argbs,arg1[i])
            argbs=np.append(argbs,arg2[i])
            idbs=np.append(idbs,i_d1[i])
            idbs=np.append(idbs,i_d2[i])
            kwbs=np.append(kwbs,kw1[i])
            kwbs=np.append(kwbs,kw2[i])
            ribs=np.append(ribs,rib[i])
            ribs=np.append(ribs,rib[i])
            m1bs=np.append(m1bs,m1b[i])
            m1bs=np.append(m1bs,m2b[i])
            zl1bs=np.append(zl1bs,zl1b[i])
            zl1bs=np.append(zl1bs,zl2b[i])
            r1bs=np.append(r1bs,r1b[i])
            r1bs=np.append(r1bs,r2b[i])
            tebs=np.append(tebs,te1[i])
            tebs=np.append(tebs,te2[i])

    else:
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

    header=sev.readline().split()
    ntot,tphys=int(header[0]),float(header[1])

    data=np.loadtxt(sev.name,skiprows=1)

    arg=data[:,1].astype(int)
    i_d=data[:,2].astype(int)
    kw=data[:,3].astype(int)
    ri=data[:,4].astype(float)
    m1=data[:,5].astype(float)
    zl1=data[:,6].astype(int)
    r1=data[:,7].astype(float)
    te=data[:,8].astype(float)

    np.nan_to_num(zl1,copy=False)
    np.nan_to_num(r1,copy=False)
    np.nan_to_num(te,copy=False)

    #Add select parameters to single star array
    if len(i_d1)>0:
        arg=np.append(argbs,arg)
        i_d=np.append(idbs,i_d)
        kw=np.append(kwbs,kw)
        ri=np.append(ribs,ri)
        m1=np.append(m1bs,m1)
        zl1=np.append(zl1bs,zl1)
        r1=np.append(r1bs,r1)
        te=np.append(tebs,te)

    return arg,i_d,kw,ri,m1,zl1,r1,te,i_d1,i_d2,kw1,kw2,kwb,rib,ecc,pb,semi,m1b,m2b,zl1b,zl2b,r1b,r2b,te1,te2


def _get_nbody6pp_hdf5(f,ngroup=0,**kwargs):
        
    #datakeys=['NAM', 'X1', 'X2', 'X3', 'V1', 'V2', 'V3', 'A1', 'A2', 'A3', 'J1', 'J2', 'J3', 'M']       
    snapshot=f['/Step#%d' % ngroup]

    ntot=snapshot.attrs['N_SINGLE']
    tphys=snapshot.attrs['Time']
    
    i_d=snapshot['NAM']
    x,y,z=snapshot['X1'],snapshot['X2'],snapshot['X3']
    vx,vy,vz=snapshot['V1'],snapshot['V2'],snapshot['V3']
    m=snapshot['M']
    
    kw,lum,rc,rs,te=snapshot['KW'],np.log10(snapshot['L']),snapshot['RC'],snapshot['RS'],snapshot['TE']
    pot=snapshot['POT']

    if 'Binaries' in snapshot:
        binaries=True
    else:
        binaries=False

    return tphys,ntot,x,y,z,vx,vy,vz,m,i_d,pot,kw,lum,rc,rs,te,binaries
    
def _get_nbody6pp_hdf5_binaries(f,ngroup=0,**kwargs):
        
    #datakeys=['A', 'ECC', 'G', 'KW1', 'KW2', 'KWC', 'L1', 'L2', 'M1', 'M2', 'MC1', 'MC2', 'NAM1', 'NAM2', 'NAMC', 'P', 'POT', 'RC1', 'RC2', 'RS1', 'RS2', 'TE1', 'TE2', 'VC1', 'VC2', 'VC3', 'VR1', 'VR2', 'VR3', 'XC1', 'XC2', 'XC3', 'XR1', 'XR2', 'XR3']      
    snapshot=f['/Step#%d/Binaries' % ngroup]

    a,ecc,gb=snapshot['A'],snapshot['ECC'],snapshot['G']
    kw1,kw2,kwc=snapshot['KW1'],snapshot['KW2'],snapshot['KWC']
    l1,l2,m1,m2,mc1,mc2=np.log10(snapshot['L1']),np.log10(snapshot['L2']),snapshot['M1'],snapshot['M2'],snapshot['MC1'],snapshot['MC2']
    id1,id2,idc=snapshot['NAM1'],snapshot['NAM2'],snapshot['NAMC']
    pb,pot,rc1,rc2,rs1,rs2,te1,te2=snapshot['P'],snapshot['POT'],snapshot['RC1'],snapshot['RC2'],snapshot['RS1'],snapshot['RS2'],snapshot['TE1'],snapshot['TE2']
    vc1,vc2,vc3,vr1,vr2,vr3=snapshot['VC1'],snapshot['VC2'],snapshot['VC3'],snapshot['VR1'],snapshot['VR2'],snapshot['VR3']
    xc1,xc2,xc3,xr1,xr2,xr3=snapshot['XC1'],snapshot['XC2'],snapshot['XC3'],snapshot['XR1'],snapshot['XR2'],snapshot['XR3']

    return a,ecc,gb,kw1,kw2,kwc,l1,l2,m1,m2,mc1,mc2,id1,id2,idc,pb,pot,rc1,rc2,rs1,rs2,te1,te2,vc1,vc2,vc3,vr1,vr2,vr3,xc1,xc2,xc3,xr1,xr2,xr3
