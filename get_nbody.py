import math
from cluster import StarCluster
from star import Star
import numpy as np

#To Do:
#Add get functions for public versions of OUT34, fort.82 and fort.83

#Get StarCluster from Gyrfalcon output
def get_gyrfalcon(filein,r0=8.0,v0=220.0,vcon=1.0,do_keyparams=True):
    nhead=13
    id=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]
   
    for j in range(0,nhead):
        data=filein.readline().split()
        if (len(data)==0):
            print('END OF FILE')
            return StarCluster(0,0.0)
        if any ('Ntot' in dat for dat in data):
                sntot=data[2]
                ntot=int(sntot[:-1])
        if any ('time' in dat for dat in data):
                tphys=float(data[2])*1000.0
        
    cluster=StarCluster(ntot,tphys,units='realkpc',origin='galaxy',keyparams=do_keyparams)
            
    for j in range(ntot):
        data=filein.readline().split()
        id.append(j)
        m.append(float(data[0]))
        x.append(float(data[1]))
        y.append(float(data[2]))
        z.append(float(data[3]))
        vx.append(float(data[4])*vcon)
        vy.append(float(data[5])*vcon)
        vz.append(float(data[6])*vcon)

    cluster.add_stars(id,m,x,y,z,vx,vy,vz)

    xgc=np.mean(x)
    ygc=np.mean(y)
    zgc=np.mean(z)
    vxgc=np.mean(vx)
    vygc=np.mean(vy)
    vzgc=np.mean(vz)

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
        m.append(float(data[8])+float(data[9])/zmbar)
        logl1.append(float(data[10]))
        logl2.append(float(data[11]))
        logl.append(max(logl1,logl2))
        logr1.append(float(data[12]))
        logr2.append(float(data[13]))
        logr.append(max(logr1,logr2))
        x.append(float(data[14]))
        y.append(float(data[15]))
        z.append(float(data[16]))
        vx.append(float(data[17]))
        vy.append(float(data[18]))
        vz.append(float(data[19]))
 
        if len(data) > 22:
            kin.append(float(data[20]))
            pot.append(float(data[21]))
            etot.append(float(data[23]))
        else:
            kin.append(0.0)
            pot.append(0.0)
            etot.append(0.0)

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

        if len(data)>13:
            kin.append(float(data[11]))
            pot.append(float(data[12]))
            etot.append(float(data[14]))
        else:
            kin.append(0.0)
            pot.append(0.0)
            etot.append(0.0)
        
        nsbnd+=1
        data=fort83.readline().split()

    nbnd=nsbnd+nbbnd

    if nbnd > 0:
        cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',keyparams=do_keyparams)
        cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd)
        cluster.add_stars(id,m,x,y,z,vx,vy,vz)
        cluster.add_se(kw,logl,logr)
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2)
        cluster.add_energies(kin,pot,etot)
    else:
        cluster=StarCluster(0,tphys)

    return cluster

#Get StarCluster from custom version of OUT34
def get_nbody6_out34(out34,do_keyparams=True):
    
    line1=out34.readline().split()
    if (len(line1)==0):
        print('END OF FILE')
        return StarCluster(0,0.0)

    line2=out34.readline().split()
    line3=out34.readline().split()

    ns=int(line1[0])
    tphys=float(line1[1])
    np=int(line1[4])

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

    data=out34.readline().split()
    while int(data[0]) > 0 and len(data)>0:

        # IGNORE GHOST PARTICLES
        if float(data[2])==0.0:
            ns-=1
            ntot-=1
        else:
            id.append(int(data[0]))
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

            if len(data)>13:
                kin.append(float(data[11]))
                pot.append(float(data[12]))
                etot.append(float(data[14]))
            else:
                kin.append(0.0)
                pot.append(0.0)
                etot.append(0.0)
        
            nsbnd+=1
        data=out34.readline().split()

    nbnd=nsbnd+nbbnd

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',keyparams=do_keyparams)
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,np)
    cluster.add_stars(id,m,x,y,z,vx,vy,vz)
    cluster.add_se(kw,logl,logr)
    cluster.add_energies(kin,pot,etot)
    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

#Get StarCluster from array of Stars
def get_starcluster(stars,do_keyparams=True):
    id=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]

    for i in range(0,len(stars)):
        id.append(stars[i].id)
        m.append(stars[i].m)
        x.append(stars[i].x)
        y.append(stars[i].y)
        z.append(stars[i].z)
        vx.append(stars[i].vx)
        vy.append(stars[i].vy)
        vz.append(stars[i].vz)
    cluster=StarCluster(len(stars),stars[0].tphys,id,m,x,y,z,vx,vy,vz,units=stars[0].units,origin=stars[0].origin,keyparams=do_keyparams)

    return cluster

#Get Star from an individual star
def get_star(cluster,id):
    indx=cluster.id.index(id)
    star=Star(cluster.id[indx],cluster.m[indx],cluster.x[indx],cluster.y[indx],cluster.z[indx],cluster.vx[indx],cluster.vy[indx], cluster.vx[indx],tphys=cluster.tphys,units=cluster.units,origin=cluster.origin)
    return star

