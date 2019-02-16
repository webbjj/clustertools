import math
from cluster import StarCluster
import numpy as np

#To Do:
#Add get functions for public versions of OUT34, fort.82 and fort.83

#Get StarCluster from Gyrfalcon output
def get_gyrfalcon(filein,r0=8.0,v0=220.0,vcon=1.0,mcon=1.0,yzswap=False,find_center=False):
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

    cluster=StarCluster(ntot,tphys,units='realkpc',origin='galaxy')
            
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

    if yzswap:
        cluster.add_stars(id,m,x,z,y,vx,vz,vy,kw=kw)
    else:
        cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw=kw)

    #Estimate center of distribution using median function
    xgc=np.median(cluster.x)
    ygc=np.median(cluster.y)
    zgc=np.median(cluster.z)
    vxgc=np.median(cluster.vx)
    vygc=np.median(cluster.vy)
    vzgc=np.median(cluster.vz)

    if find_center:
        xgc,ygc,zgc,vxgc,vygz,vzgc=cluster.find_center(xgc,ygc,zgc)

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

#Get StarCluster from NBODY6 using Jarrod Hurley's version of hrplot.f
def get_nbody6_jarrod(fort82,fort83,do_keyparams=True):
    
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
        cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster')
        cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd)
        cluster.add_stars(id,m,x,y,z,vx,vy,vz)
        cluster.add_se(kw,logl,logr,ep,ospin)
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2,ep1,ep2,ospin1,ospin2)
        cluster.add_energies(kin,pot,etot)
    else:
        cluster=StarCluster(0,tphys)

    return cluster

#Get StarCluster from custom version of OUT34 and OUT9 (if binaries)
def get_nbody6_out(out9,out34,debug=False):
    
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
    xgc=float(line1[5])*rbar/1000.0
    ygc=float(line1[6])*rbar/1000.0
    zgc=float(line1[7])*rbar/1000.0
    vxgc=float(line1[8])*vstar
    vygc=float(line1[9])*vstar
    vzgc=float(line1[10])*vstar

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
    if debug:
        debug_column=[]

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

            if debug:
                debug_column.append(float(data[13]))

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

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster')
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,n_p)
    cluster.add_stars(id,m,x,y,z,vx,vy,vz)
    cluster.add_se(kw,logl,logr)
    cluster.add_energies(kin,pot,etot)
    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)
    if out9!=None:
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2)

    if debug:
        return cluster,debug_column
    else:
        return cluster

#Get StarCluster from custom version of OUT34
def get_nbody6_out34(out34,debug=False):
    
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

    #Put Orbital Properties in Real Units
    xgc=float(line1[5])*rbar/1000.0
    ygc=float(line1[6])*rbar/1000.0
    zgc=float(line1[7])*rbar/1000.0
    vxgc=float(line1[8])*vstar
    vygc=float(line1[9])*vstar
    vzgc=float(line1[10])*vstar

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
    if debug:
        debug_column=[]

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

            if debug:
                debug_column.append(float(data[13]))

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

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster')
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,n_p)
    cluster.add_stars(id,m,x,y,z,vx,vy,vz)
    cluster.add_se(kw,logl,logr,np.zeros(nbnd),np.zeros(nbnd))
    cluster.add_energies(kin,pot,etot)
    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    if debug:
        return cluster,debug_column
    else:
        return cluster

#Get StarCluster snapshots produces by snapauto
def get_nbody6_snapauto(snap,debug=False):
    
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
    if debug:
        debug_column=[]

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

            if debug:
                debug_column.append(float(data[13]))
        
            nbnd+=1

            data=snap.readline().split()

    kw=np.zeros(int(ntot)) 
    cluster=StarCluster(nbnd,tphys,units='realpc',origin='cluster')
    cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw=kw)

    if debug:
        return cluster,debug_column
    else:
        return cluster

#Get StarCluster from standard nbodypy snapshot
def get_nbodypy_snapshot(snap,delimiter=',',units0='realpc',origin0='cluster',find_center=True,debug=False):
    
    data=np.loadtxt(snap,delimiter=delimiter)

    m=data[:,0]
    x=data[:,1]
    y=data[:,2]
    z=data[:,3]
    vx=data[:,4]
    vy=data[:,5]
    vz=data[:,6]
    id=data[:,7]
    kw=data[:,8]


    if debug:
        debug_column=data[:,9]
        
    nbnd=len(m)

    cluster=StarCluster(nbnd,0.0,units=units0,origin=origin0)
    cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw=kw)

    #Estimate center of distribution using median function
    xgc=np.median(cluster.x)
    ygc=np.median(cluster.y)
    zgc=np.median(cluster.z)
    vxgc=np.median(cluster.vx)
    vygc=np.median(cluster.vy)
    vzgc=np.median(cluster.vz)
    
    if find_center:
        xgc,ygc,zgc,vxgc,vygz,vzgc=cluster.find_center(xgc,ygc,zgc)

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    if debug:
        return cluster,debug_column
    else:
        return cluster

