import math
from cluster import StarCluster
import numpy as np
from galpy.util import bovy_conversion
import os

#get_cluster is the main function for loading clusters from all types of codes (e.g. NBODY6, GYRFALCON, etc.). Individual functions are then called that have been customised to call files from a given software package and any necessary **kwargs for those files are then required. Due to the common practice of editing output files in commonly used and publically availble codes, get_cluster is written to be modular in the sense that it is trivial to add your own get_mycode() function to this list by creating a new ctype.

def get_cluster(ctype,units0='realpc',origin0='cluster',**kwargs):
    #ctype = fort for Nbody6/Nbody6tt with Jarrod Hurley's version of hrplot.f
    #ctype = out is for Nbody6/Nbody6tt with Jeremy Webb's version of OUT9 and OUT34
    #ctype= snapauto is for Nbody6/Nbody6tt/Nbody6++ with Jongsuk Hong's code snap6++.f for conversion OUT3/conf_ binary files to ascii snapshots
    #ctype = gyrfalcon is for a gyrfalcon simulation file that has been converted to ascii via the s2a command.
    #ctype = snapshot is for snapshot files created with nbodypy
    
    #Default units are set to realpc and default origin is set to cluster. These are really only used when reading in nbodypy snapshots, as other codes have their own default units and origin
    
    #kwargs:
    #orbit - allows for the file containing cluster orbital information to be passed through for get_cluster_orbit()
    #nsnap - if a specific snapshot is to be read in instead of starting from zero
    #filename - when ctype=gyrfalcon, the file name needs to read in as its a customizable option in gyrfalcon
    #nzfill - value for zfill when reading and writing snapshots (Default: 5)
    #delimiter - choice of delimiter when reading ascii/csv files (Default: ',')
    
    if kwargs['orbit']:
        ofile=kwargs.get('orbit')
    else:
        ofile=None
    
    if ctype=='fort':
        fort82=open('fort.82','r')
        fort83=open('fort.83','r')
        cluster=get_nbody6_jarrod(fort82,fort83,ofile)
    elif ctype=='out':
        if os.path.isfile('OUT9'):
            out9=open('OUT9','r')
        else:
            out9=None
        out34=open('OUT34','r')
        cluster=get_nbody6_out(out9,out34)
    #Read in snapshot produced from snapauto.f which reads binary files from either NBODY6 or NBODY6++
    elif ctype=='snapauto':
        nsnap=kwargs.get('nsnap')
        
        if kwargs['nzfill']:
            nzfill=kwargs.get('nzfill')
        else:
            nzfill=5

        snapfile=('./snapshot/snap.dat.%s' % str(nsnap).zfill(nzfill))
        cluster=get_nbody6_snapauto(snapfile,units0,origin0,nsnap,nzfill,ofile)
    #Read in snapshot from gyrfalcon. Note option for swapping y-z coordinates (necessary when using a logarithmic potential)
    elif ctype=='gyrfalcon':
        filein=kwargs.get('filename')
        cluster=get_gyrfalcon(filein,'WDunits',ofile)
    #Read in standard nbodypy snapshot
    elif ctype=='snapshot':
        nsnap=kwargs.get('nsnap')

        if 'nzfill' in kwargs:
            nzfill=kwargs.get('nzfill')
        else:
            nzfill=5

        if 'delimiter' in kwargs:
            delimiter=kwargs.get('delimiter')
        else:
            delimiter=','

        cluster=get_nbodypy_snapshot(nsnap,units0,origin0,delimiter,nzfill,ofile)

    return cluster

def advance_cluster(cluster,ofile):
    print(cluster.ctype)
    #Either read in fort.82/fort83, OUT34, or OUT3 and OUT34 from NBODY6. Note these output files have been customised to match my personal versions
    if cluster.ctype=='fort':
        cluster=get_nbody6_jarrod(cluster.bfile,cluster.sfile,ofile)
    elif cluster.ctype=='out':
        cluster=get_nbody6_out(cluster.bfile,cluster.sfile,ofile)
    #Read in snapshot produced from snapauto.f which reads binary files from either NBODY6 or NBODY6++
    elif cluster.ctype=='snapauto':
        nsnap=int(cluster.nsnap)+1
        cluster=get_nbody6_snapauto(nsnap,cluster.units,cluster.origin,ofile,nsnap=cluster.nsnap+1,nzfill=cluster.nzfill,ocontinue=True)
    #Read in snapshot from gyrfalcon. Note option for swapping y-z coordinates (necessary when using a logarithmic potential)
    elif cluster.ctype=='gyrfalcon':
        cluster=get_gyrfalcon(cluster.sfile,'WDunits',ofile)
    elif cluster.ctype=='snapshot':
        nsnap=int(cluster.nsnap)+1
        cluster=get_nbodypy_snapshot(nsnap,cluster.units,cluster.origin,cluster.delimiter,cluster.nzfill,ofile,ocontinue=True)

    return cluster

def get_cluster_orbit(orbit,cluster,nsnap=None,ocontinue=False):
    if nsnap!=None and not ocontinue:
        for i in range(0,int(nsnap)+1):
            data=orbit.readline().split()
    else:
        data=orbit.readline().split()
    
    if orbit.name=='gc_orbit.dat':
        #Saved orbit from doing a grep of NBODY6 logfile

        xgc=float(data[8])
        ygc=float(data[9])
        zgc=float(data[10])
        vxgc=float(data[11])
        vygc=float(data[12])
        vzgc=float(data[13])
    else:
        #Input logfile, assumed units is kpc and km/s
        xgc=float(data[1])
        ygc=float(data[2])
        zgc=float(data[3])
        vxgc=float(data[4])
        vygc=float(data[5])
        vzgc=float(data[6])
    
    if cluster.units=='nbody':
        xgc*=(1000.0/cluster.rbar)
        ygc*=(1000.0/cluster.rbar)
        zgc*=(1000.0/cluster.rbar)
        vxgc/=(cluster.vstar)
        vygc/=(cluster.vstar)
        vzgc/=(cluster.vstar)
    elif cluster.units=='realpc':
        xgc*=1000.0
        ygc*=1000.0
        zgc*=1000.0
    elif cluster.units=='galpy':
        xgc/=8.
        ygc/=8.
        zgc/=8.
        vxgc/=220.
        vygc/=220.
        vzgc/=220.
    return xgc,ygc,zgc,vxgc,vygc,vzgc


#Get StarCluster from Gyrfalcon output
def get_gyrfalcon(filein,units='WDunits',ofile=None):
    if units=='WDunits':
        vcon=220.0/bovy_conversion.velocity_in_kpcGyr(220.0,8.0)
        mcon=222288.4543021174
    else:
        vcon=1.
        mcon=1.

    nhead=13
    id=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]

    over_head=False

    for j in range(0,nhead):
        data=filein.readline().split()
        if '#' not in data:
            over_head=True
            print('OVER HEAD')
            break
        if (len(data)==0):
            print('END OF FILE')
            return StarCluster(0,0.0)
        if any ('Ntot' in dat for dat in data):
                sntot=data[2]
                ntot=int(sntot[:-1])
        if any ('time' in dat for dat in data):
                tphys=float(data[2])*1000.0        

    cluster=StarCluster(ntot,tphys,units='realkpc',origin='galaxy',ctype='gyrfalcon',sfile=filein,bfile=None)
            
    for j in range(ntot):
        if over_head:
            over_head=False
        else:
            data=filein.readline().split()
        if '#' in data:
            print('REACHED HEADER')
            break

        id.append(j+1)
        m.append(float(data[0])*mcon)
        x.append(float(data[1]))
        y.append(float(data[2]))
        z.append(float(data[3]))
        vx.append(float(data[4])*vcon)
        vy.append(float(data[5])*vcon)
        vz.append(float(data[6])*vcon)

    kw=np.zeros(ntot)

    cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw=kw)

    if ofile!=None:
        xgc,ygc,zgc,vxgc,vygc,vzgc=get_cluster_orbit(ofile,cluster)
        xc=np.median(cluster.x)
        yc=np.median(cluster.y)
        zc=np.median(cluster.z)
        cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xc,yc,zc)
    else:
        #Estimate center of distribution using median function
        xgc=np.median(cluster.x)
        ygc=np.median(cluster.y)
        zgc=np.median(cluster.z)
        vxgc=np.median(cluster.vx)
        vygc=np.median(cluster.vy)
        vzgc=np.median(cluster.vz)
        xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.find_center(xgc,ygc,zgc)

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

#Get StarCluster from NBODY6 using Jarrod Hurley's version of hrplot.f
def get_nbody6_jarrod(fort82,fort83,ofile=None):
    
    line1=fort83.readline().split()
    if (len(line1)==0):
        print('END OF FILE')
        return StarCluster(0,0.0)

    line2=fort83.readline().split()
    line3=fort83.readline().split()
    line1b=fort82.readline().split()


    ns=int(line1[0])
    tphys=float(line1[1])
    nc=int(line2[0])
    rc=max(float(line2[1]),0.01)
    rbar=float(line2[2])
    rtide=float(line2[3])
    xc=float(line2[4])
    yc=float(line2[5])
    zc=float(line2[6])
    zmbar=float(line3[0])
    vstar=0.06557*math.sqrt(zmbar/rbar)
    rscale=float(line3[2])
    nb=int(line1b[0])
    ntot=ns+nb

    nsbnd=0
    nbbnd=0
    nbnd=0
    

    id=[]
    id1=[]
    id2=[]
    kw=[]
    kw1=[]
    kw2=[]
    kcm=[]
    ecc=[]
    pb=[]
    semi=[]
    m1=[]
    m2=[]
    m=[]
    logl1=[]
    logl2=[]
    logl=[]
    logr1=[]
    logr2=[]
    logr=[]
    x=[]
    y=[]
    z=[]
    rxy=[]
    r=[]
    vx=[]
    vy=[]
    vz=[]
    v=[]
    ep=[]
    ep1=[]
    ep2=[]
    ospin=[]
    ospin1=[]
    ospin2=[]
    kin=[]
    pot=[]
    etot=[]

    data=fort82.readline().split()

    while int(data[0]) > 0 and len(data)>0:
        id1.append(int(data[0]))
        id2.append(int(data[1]))
        id.append(id1[-1])
        kw1.append(int(data[2]))
        kw2.append(int(data[3]))
        kw.append(max(kw1[-1],kw2[-1]))
        kcm.append(float(data[4]))
        ecc.append(float(data[5]))
        pb.append(float(data[6]))
        semi.append(float(data[7]))
        m1.append(float(data[8])/zmbar)
        m2.append(float(data[9])/zmbar)
        m.append(m1[-1]+m2[-1])
        logl1.append(float(data[10]))
        logl2.append(float(data[11]))
        logl.append(max(logl1[-1],logl2[-1]))
        logr1.append(float(data[12]))
        logr2.append(float(data[13]))
        logr.append(max(logr1,logr2))
        x.append(float(data[14]))
        y.append(float(data[15]))
        z.append(float(data[16]))
        vx.append(float(data[17]))
        vy.append(float(data[18]))
        vz.append(float(data[19]))

        if 'bnd' in fort82.name or 'esc' in fort82.name:
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

        nbbnd+=1
        data=fort82.readline().split()


    data=fort83.readline().split()
    while int(data[0]) > 0 and len(data)>0:
        id.append(int(data[0]))
        kw.append(int(data[1]))
        m.append(float(data[2])/zmbar)
        logl.append(float(data[3]))
        logr.append(float(data[4]))
        x.append(float(data[5]))
        y.append(float(data[6]))
        z.append(float(data[7]))
        vx.append(float(data[8]))
        vy.append(float(data[9]))
        vz.append(float(data[10]))

        if 'bnd' in fort83.name or 'esc' in fort83.name:
            kin.append(float(data[11]))
            pot.append(float(data[12]))
            etot.append(float(data[14]))
        else:
            kin.append(0.0)
            pot.append(0.0)
            etot.append(0.0)
            ep.append(float(data[11]))
            ospin.append(float(data[12]))

        nsbnd+=1
        data=fort83.readline().split()

    nbnd=nsbnd+nbbnd

    if nbnd > 0:
        cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',ctype='fort',sfile=fort83,bfile=fort82)
        cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd)
        cluster.add_stars(id,m,x,y,z,vx,vy,vz)
        cluster.add_se(kw,logl,logr,ep,ospin)
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2,ep1,ep2,ospin1,ospin2)
        cluster.add_energies(kin,pot,etot)
    
        if ofile!=None:
            xgc,ygc,zgc,vxgc,vygc,vzgc=get_cluster_orbit(ofile,cluster)
            cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

        #Estimate center of distribution using median function
        xc=np.median(cluster.x)
        yc=np.median(cluster.y)
        zc=np.median(cluster.z)
        cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xc,yc,zc)
    
    else:
        cluster=StarCluster(0,tphys)

    return cluster

#Get StarCluster from custom version of OUT34 and OUT9 (if binaries)
def get_nbody6_out(out9,out34):
    
    line1=out34.readline().split()
    if (len(line1)==0):
        print('END OF FILE')
        return StarCluster(0,0.0)

    line2=out34.readline().split()
    line3=out34.readline().split()

    ns=int(line1[0])
    tphys=float(line1[1])
    n_p=int(line1[4])
    nb=int(float(line1[11]))

    if out9!=None:
        line1b=out9.readline().split()
        line2b=out9.readline().split()
        line3b=out9.readline().split()

        if nb!=int(line1b[0]):
            print('ERROR: NUMBER OF BINARIES DO NOT MATCH')

    nc=int(line2[0])
    rc=max(float(line2[1]),0.01)
    rbar=float(line2[2])
    rtide=float(line2[3])
    xc=float(line2[4])
    yc=float(line2[5])
    zc=float(line2[6])
    zmbar=float(line3[0])
    vstar=0.06557*math.sqrt(zmbar/rbar)
    rscale=float(line3[2])
    ntot=ns+nb

    #Put Orbital Properties in Real Units
    xgc=float(line1[5])
    ygc=float(line1[6])
    zgc=float(line1[7])
    vxgc=float(line1[8])
    vygc=float(line1[9])
    vzgc=float(line1[10])

    nsbnd=0
    nbbnd=0
    nbnd=0

    id=[]
    kw=[]
    m=[]
    logl=[]
    logr=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    kin=[]
    pot=[]
    etot=[]

    if out9!=None:
        
        yrs = (rbar*1296000./(2.0*np.pi))**1.5/np.sqrt(zmbar)
        days = 365.25*yrs
        
        id1=[]
        kw1=[]
        m1=[]
        logl1=[]
        logr1=[]
        id2=[]
        kw2=[]
        m2=[]
        logl2=[]
        logr2=[]
        pb=[]
        kcm=[]
        ecc=[]
        semi=[]
        
        for i in range(0,nb):
            data=out9.readline().split()
            nbbnd+=1

            ecc.append(float(data[1]))
            m1.append(float(data[4])/zmbar)
            m2.append(float(data[5])/zmbar)
            pb.append(float(data[6])/days)
            id1.append(int(float(data[7])))
            id2.append(int(float(data[8])))
            kw1.append(int(data[9]))
            kw2.append(int(data[10]))
            kcm.append(int(data[11]))
            
            
            logl1.append(1.0)
            logl2.append(1.0)
            logr1.append(1.0)
            logr2.append(1.0)

            x1=float(data[12])
            y1=float(data[13])
            z1=float(data[14])
            vx1=float(data[15])
            vy1=float(data[16])
            vz1=float(data[17])
            x2=float(data[18])
            y2=float(data[19])
            z2=float(data[20])
            vx2=float(data[21])
            vy2=float(data[22])
            vz2=float(data[23])

            x.append((x1*m1[-1]+x2*m2[-1])/(m1[-1]+m2[-1]))
            y.append((y1*m1[-1]+y2*m2[-1])/(m1[-1]+m2[-1]))
            z.append((z1*m1[-1]+z2*m2[-1])/(m1[-1]+m2[-1]))
            vx.append((vx1*m1[-1]+vx2*m2[-1])/(m1[-1]+m2[-1]))
            vy.append((vy1*m1[-1]+vy2*m2[-1])/(m1[-1]+m2[-1]))
            vz.append((vz1*m1[-1]+vz2*m2[-1])/(m1[-1]+m2[-1]))
            m.append(m1[-1]+m2[-1])
            id.append(id1[-1])
            kw.append(max(kw1[-1],kw2[-1]))
            logl.append(1.0)
            logr.append(1.0)

            r1=np.sqrt((x1-x[-1])**2.0+(y1-y[-1])**2.0+(z1-z[-1])**2.0)
            r2=np.sqrt((x2-x[-1])**2.0+(y2-y[-1])**2.0+(z2-z[-1])**2.0)
            
            semi.append((abs((pb[-1]**2.0)*(m[-1]/2.)))**(1./3.))

    data=out34.readline().split()
    while int(float(data[0])) >= -999 and len(data)>0:
        # IGNORE GHOST PARTICLES
        if float(data[2])==0.0:
            ns-=1
            ntot-=1
        else:
            id.append(int(float(data[0])))
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

            if len(data)>14:
                kin.append(float(data[13]))
                pot.append(float(data[14]))
                etot.append(float(data[15]))
            else:
                kin.append(0.0)
                pot.append(0.0)
                etot.append(0.0)
        
            nsbnd+=1
        data=out34.readline().split()

    nbnd=nsbnd+nbbnd

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',ctype='out',sfile=out34,bfile=out9)
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,n_p)
    #Add back on the center of mass which has been substracted off by NBODY6
    cluster.add_stars(id,m,x+xc,y+yc,z+zc,vx,vy,vz)
    cluster.add_se(kw,logl,logr)
    cluster.add_energies(kin,pot,etot)
    if out9!=None:
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2)

    #Estimate center of distribution using median function
    xc=np.median(cluster.x)
    yc=np.median(cluster.y)
    zc=np.median(cluster.z)
    cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xc,yc,zc)

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

#Get StarCluster from custom version of OUT34
def get_nbody6_out34(out34):
    
    line1=out34.readline().split()
    if (len(line1)==0):
        print('END OF FILE')
        return StarCluster(0,0.0)

    line2=out34.readline().split()
    line3=out34.readline().split()

    ns=int(line1[0])
    tphys=float(line1[1])
    n_p=int(line1[4])

    nc=int(line2[0])
    rc=max(float(line2[1]),0.01)
    rbar=float(line2[2])
    rtide=float(line2[3])
    xc=float(line2[4])
    yc=float(line2[5])
    zc=float(line2[6])
    zmbar=float(line3[0])
    vstar=0.06557*math.sqrt(zmbar/rbar)
    rscale=float(line3[2])
    nb=0
    ntot=ns+nb

    #Orbit Properties
    xgc=float(line1[5])
    ygc=float(line1[6])
    zgc=float(line1[7])
    vxgc=float(line1[8])
    vygc=float(line1[9])
    vzgc=float(line1[10])
    
    nsbnd=0
    nbbnd=0
    nbnd=0
    
    id=[]
    kw=[]
    m=[]
    logl=[]
    logr=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
    kin=[]
    pot=[]
    etot=[]

    data=out34.readline().split()
    while int(float(data[0])) >= -999 and len(data)>0:
        # IGNORE GHOST PARTICLES
        if float(data[2])==0.0:
            ns-=1
            ntot-=1
        else:
            id.append(int(float(data[0])))
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

            if len(data)>14:
                kin.append(float(data[13]))
                pot.append(float(data[14]))
                etot.append(float(data[15]))
            else:
                kin.append(0.0)
                pot.append(0.0)
                etot.append(0.0)
        
            nsbnd+=1
        data=out34.readline().split()

    nbnd=nsbnd+nbbnd

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',ctype='out',sfile=out34,bfile=None)
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,n_p)
    #Add back on the center of mass which has been substracted off by NBODY6
    cluster.add_stars(id,m,x+xc,y+yc,z+zc,vx,vy,vz)
    cluster.add_se(kw,logl,logr,np.zeros(nbnd),np.zeros(nbnd))
    cluster.add_energies(kin,pot,etot)
    
    #Estimate center of distribution using median function
    xc=np.median(cluster.x)
    yc=np.median(cluster.y)
    zc=np.median(cluster.z)
    cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xc,yc,zc)

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)


    return cluster

#Get StarCluster snapshots produces by snapauto
def get_nbody6_snapauto(snap,units0='realpc',origin0='cluster',ofile=None,**kwargs):
    
    nsnap=kwargs.get('nsnap')
    nzfill=kwargs.get('nzfill')
    if 'ocontinue' in kwargs:
        ocontinue=kwargs.get('ocontinue')
    else:
        ocontinue=False
    
    data=snap.readline().split()
    tphys=float(data[0])
    ntot=float(data[1])
    nb=int(float(data[2]))
    xc=float(data[3])
    yc=float(data[4])
    zc=float(data[5])
    
    id=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]

    data=snap.readline().split()
    nbnd=0
    while len(data)>0:
            if int(data[7]) <=0:
                break

            m.append(float(data[0]))
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))
            vx.append(float(data[4]))
            vy.append(float(data[5]))
            vz.append(float(data[6]))
            id.append(int(float(data[7])))
        
            nbnd+=1

            data=snap.readline().split()

    kw=np.zeros(int(ntot)) 
    cluster=StarCluster(nbnd,tphys,units='realpc',origin='cluster',ctype='snapauto',nsnap=nsnap,nzfill=nzfill)
    cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw=kw)
    
    #Estimate center of distribution using median function
    xc=np.median(cluster.x)
    yc=np.median(cluster.y)
    zc=np.median(cluster.z)
    cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xc,yc,zc)

    if ofile!=None:
        xgc,ygc,zgc,vxgc,vygc,vzgc=get_cluster_orbit(ofile,cluster,nsnap,ocontinue=ocontinue)
        cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

#Get StarCluster from standard nbodypy snapshot
def get_nbodypy_snapshot(nsnap,units0='realpc',origin0='cluster',delimiter=',',nzfill=5,ofile=None,**kwargs):
    
    if 'ocontinue' in kwargs:
        ocontinue=kwargs.get('ocontinue')
    else:
        ocontinue=False
   
    if os.path.isfile('./snaps/snap.dat.%s' % str(nsnap).zfill(nzfill)):
        nsnapfile=('./snaps/snap.dat.%s' % str(nsnap).zfill(nzfill))
        data=np.loadtxt(nsnapfile,delimiter=delimiter)
    elif os.path.isfile('./snaps/fort_%s.dat' % str(nsnap).zfill(nzfill)):
        nsnapfile=('./snaps/fort_%s.dat' % str(nsnap).zfill(nzfill))
        data=np.loadtxt(nsnapfile,delimiter=delimiter)
    else:
        print('NO SNAPSHOT FOUND')
        return 0

    m=data[:,0]
    x=data[:,1]
    y=data[:,2]
    z=data[:,3]
    vx=data[:,4]
    vy=data[:,5]
    vz=data[:,6]
    id=data[:,7]
    kw=data[:,8]
    
    nbnd=len(m)

    cluster=StarCluster(nbnd,0.0,units=units0,origin=origin0,ctype='snapshot',sfile=nsnapfile,bfile=None,delimiter=delimiter,nsnap=nsnap,nzfill=nzfill)
    cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw=kw)

    if origin0=='galaxy':

        if ofile==None:
            xgc,ygc,zgc,vxgc,vygc,vzgc=cluster.find_center(xgc,ygc,zgc)
        else:
            xgc,ygc,zgc,vxgc,vygc,vzgc=get_cluster_orbit(ofile,cluster,nsnap,ocontinue=ocontinue)

        cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)
        #Estimate center of distribution using median function
        cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xgc,ygc,zgc)
    elif origin0=='cluster':
        #Estimate center of distribution using median function
        xc=np.median(cluster.x)
        yc=np.median(cluster.y)
        zc=np.median(cluster.z)
        cluster.xc,cluster.yc,cluster.zc,cluster.vxc,cluster.vyc,cluster.vzc=cluster.find_center(xc,yc,zc)

        if ofile!=None:
            xgc,ygc,zgc,vxgc,vygc,vzgc=get_cluster_orbit(ofile,cluster,nsnap,ocontinue=ocontinue)
            cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

