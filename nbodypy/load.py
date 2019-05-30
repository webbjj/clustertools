#Read in cluster from Nbody simulations or generate an Nbody cluster

import numpy as np
from galpy.util import bovy_conversion
import os
from .cluster import StarCluster
from .operations import *
from .orbit import initialize_orbit

def load_cluster(ctype='snapshot',units='realpc',origin='cluster',ofile=None,orbit=None,filename=None,**kwargs):
    """
    NAME:

       load_cluster

    PURPOSE:

       Load a StarCluster snapshot from a generic code output
       --> Note user can easily create their own def my_code() for new code output types

    INPUT:

       ctype - Type of file being loaded (Currently supports nbody6, nbody6se, gyrfalcon, snapauto, nbodypy, snapshot)

       units - units of input data (default: realkpc)

       origin - origin of input data (default: cluster)

       ofile - an already opened file containing orbit information (default: None)

       orbit - a galpy orbit to be used for the StarCluster's orbital information (default: None)

       filename - name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)

    KWARGS:

        ofilename - orbit filename if ofile is not given

        ounits - units of orbital information (else assumed equal to StarCluster.units)

        nsnap - if a specific snapshot is to be read in instead of starting from zero
    
        nzfill - value for zfill when reading and writing snapshots (Default: 5)
    
        delimiter - choice of delimiter when reading ascii/csv files (Default: ',')
    
        wdir - working directory of snapshots if not current directory

        intialize - initialize a galpy orbit after reading in orbital information (default: False)

        kwfile - open file containing stellar evolution type (kw) of individual stars (used if snapshot file does not contain kw and kw is needed)

        ocontinue - if True read next line (or line # nsnap) of oribital information file instead of starting from beginning

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    wdir=kwargs.get('wdir','./')
    initialize=kwargs.get('initialize',False)

    if 'ofilename' in kwargs and ofile==None:
        ofile=open(wdir+kwargs.get('ofilename'),'r')

    if ctype=='nbody6se':
        #When stellar evolution is turned on, read in fort.82 and fort.83 and if possible gc_orbit.dat
        fort82=open('%sfort.82' % wdir,'r')
        fort83=open('%sfort.83' % wdir,'r')
        cluster=get_nbody6_jarrod(fort82,fort83,ofile,**kwargs)
    elif ctype=='nbody6':
        #With stellar evolution turned off, read in OUT9 and OUT34. Orbit data already in OUT34
        if os.path.isfile('%sOUT9'):
            out9=open('%sOUT9' % wdir,'r')
        else:
            out9=None
        out34=open('%sOUT34' % wdir,'r')
        cluster=get_nbody6_out(out9,out34,**kwargs)
    elif ctype=='snapauto':
        #Read in snapshot produced from snapauto.f which reads binary files from either NBODY6 or NBODY6++
        cluster=get_nbody6_snapauto(filename,units,origin,ofile,**kwargs)
    elif ctype=='gyrfalcon':
        #Read in snapshot from gyrfalcon.
        filein=open(wdir+filename,'r')
        cluster=get_gyrfalcon(filein,'WDunits',ofile,**kwargs)
    elif ctype=='nbodypy':
        #Read in standard nbodypy snapshot
        cluster=get_nbodypy_snapshot(filename,units,origin,ofile,**kwargs)
    elif ctype=='snapshot':
         #Read in standard generic snapshot
        col_names=kwargs.pop('col_names',['m','x','y','z','vx','vy','vz'])
        col_nums=kwargs.pop('col_nums',[0,1,2,3,4,5,6])
        cluster=get_snapshot(filename,col_names,col_nums,units,origin,ofile,**kwargs)       
    elif ctype=='mycode':
        #Read in new cluster type
        cluster=get_mycode()
    else:
        print('NO CTYPE GIVEN')
        return 0

    kwfile=kwargs.get('kwfile',None)
    if kwfile!=None:
        cluster.kw[cluster.id-1]=get_kwtype(cluster,kwfile)

   #Add galpy orbit if given
    if orbit!=None:
        cluster.orbit=orbit
        if cluster.units=='realpc':
            t=(cluster.tphys/1000.)/bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
            cluster.add_orbit(orbit.x(t)*1000.,orbit.y(t)*1000.,orbit.z(t)*1000.,orbit.vx(t),orbit.vy(t),orbit.vz(t))
        if cluster.units=='realkpc':
            cluster.add_orbit(orbit.x(t),orbit.y(t),orbit.z(t),orbit.vx(t),orbit.vy(t),orbit.vz(t))
        elif cluster.units=='nbody':
            t=(cluster.tphys*cluster.tstar/1000.)/bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
            cluster.add_orbit(orbit.x(t)*1000./cluster.rbar,orbit.y(t)*1000./cluster.rbar,orbit.z(t)*1000./cluster.rbar,orbit.vx(t)/cluster.vstar,orbit.vy(t)/cluster.vstar,orbit.vz(t)/cluster.vstar)
        elif cluster.units=='galpy':
            t=cluster.tphys
            cluster.add_orbit(orbit.x(t)/8.,orbit.y(t)/8.,orbit.z(t)/8.,orbit.vx(t)/220.,orbit.vy(t)/220.,orbit.vz(t)/220.)

        units0,origin0=save_cluster(cluster)

        cluster.to_cluster()
        cluster.find_centre()

        return_cluster(cluster,units0,origin0)
    elif initialize:
        initialize_orbit(cluster)

    cluster.key_params()


    return cluster

def advance_cluster(cluster,ofile=None,orbit=None,filename=None,**kwargs):
    """
    NAME:

       advance_cluster

    PURPOSE:

       Advance a loaded StarCluster snapshot to the next timestep
       --> ofile or orbit need to be provded again, same as load_cluster.
       --> Be sure that ocontinue is set to True so next line of orbit file is read in

    INPUT:

       filename - name of file to be opened (optional - necessary if no defaults assigned to ctype) (default: None)

    KWARGS:

        same as load_cluster

    OUTPUT:

       StarCluster instance

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    if 'kwfile' in kwargs:
        kw0=cluster.kw

    #Continue reading in cluster opened in get_cluster()
    if cluster.ctype=='nbody6se':
        cluster=get_nbody6_jarrod(cluster.bfile,cluster.sfile,ofile,ocontinue=True,**kwargs)
    elif cluster.ctype=='nbody6':
        cluster=get_nbody6_out(cluster.bfile,cluster.sfile,ocontinue=True,**kwargs)
    elif cluster.ctype=='snapauto':
        nsnap=np.maximum(int(kwargs.pop('nsnap','1')),cluster.nsnap)+1
        cluster=get_nbody6_snapauto(filename,cluster.units,cluster.origin,ofile,ocontinue=True,nsnap=nsnap,**kwargs)
    elif cluster.ctype=='gyrfalcon':
        cluster=get_gyrfalcon(cluster.sfile,'WDunits',ofile,**kwargs)
    elif cluster.ctype=='nbodypy':
        nsnap=np.maximum(int(kwargs.pop('nsnap','0')),cluster.nsnap)+1
        cluster=get_nbodypy_snapshot(filename,cluster.units,cluster.origin,ofile,ocontinue=True,nsnap=nsnap,**kwargs)
    elif ctype=='snapshot':
        col_names=kwargs.pop('col_names',['m','x','y','z','vx','vy','vz'])
        col_nums=kwargs.pop('col_nums',[0,1,2,3,4,5,6])
        cluster=get_snapshot(filename,col_names,col_nums,cluster.units,cluster.origin,ofile,ocontinue=True,**kwargs)  
    elif ctype=='mycode':
        cluster=get_mycode()
    else:
        cluster=StarCuster()

    #Check for restart
    if cluster.ntot==0.0:
        try:
            wdir=cluster.wdir+'cont/'
        except:
            print('WDIR NOT SET')
            wdir='./cont/'
        
        try:
            ofilename=ofile.name
        except:
            print('OFILE NOT SET')
            ofile=None

        if os.path.exists(wdir): cluster=npy.get_cluster(cluster.ctype,ofile=ofile,kwfile=kwfile,wdir=wdir)

    if cluster.ntot!=0.0:
        if 'kwfile' in kwargs:
            cluster.kw[cluster.id-1]=kw0

        #Add galpy orbit if given
        if orbit!=None:
            cluster.orbit-orbit
            if cluster.units=='realpc' or cluster.units=='realkpc':
                t=(cluster.tphys/1000.)/bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
            elif cluster.units=='nbody':
                t=(cluster.tphys*cluster.tstar/1000.)/bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
            elif cluster.units=='galpy':
                t=cluster.tphys

            cluster.add_orbit(orbit.x(t),orbit.y(t),orbit.z(t),orbit.vx(t),orbit.vy(t),orbit.vz(t))

        cluster.key_params()

    return cluster

def get_cluster_orbit(cluster,ofile,**kwargs):
    
    ocontinue=kwargs.get('ocontinue',False)
    nsnap=int(kwargs.get('nsnap',cluster.nsnap))
    ounits=kwargs.get('ounits',None)
    
    #Read in orbital information from orbit
    if nsnap!=0 and not ocontinue:
        for i in range(0,int(nsnap)+1):
            data=ofile.readline().split()
    else:
        data=ofile.readline().split()
    
    if 'gc_orbit.dat' in ofile.name:
        #Saved orbit from doing a grep of NBODY6 or NBODY6++ logfile

        if len(data)==18:
            xgc=float(data[9])
            ygc=float(data[10])
            zgc=float(data[11])
            vxgc=float(data[12])
            vygc=float(data[13])
            vzgc=float(data[14])
        else:
            xgc=float(data[8])
            ygc=float(data[9])
            zgc=float(data[10])
            vxgc=float(data[11])
            vygc=float(data[12])
            vzgc=float(data[13])
    else:
        xgc=float(data[1])
        ygc=float(data[2])
        zgc=float(data[3])
        vxgc=float(data[4])
        vygc=float(data[5])
        vzgc=float(data[6])


    if ounits==None and 'gc_orbit.dat' in ofile.name:
        ounits='realkpc'

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc,ounits)
                
    return

def get_kwtype(cluster,kwfile):
    kw= np.loadtxt(kwfile)[:,1]
    indx=cluster.id-1
    kw0=kw[indx]
    return kw0

#Get StarCluster from Gyrfalcon output
def get_gyrfalcon(filein,units='WDunits',ofile=None,**kwargs):
    if units=='WDunits':
        vcon=220.0/bovy_conversion.velocity_in_kpcGyr(220.0,8.0)
        mcon=222288.4543021174
    else:
        vcon=1.
        mcon=1.

    nhead=13
    ntot=0
    i_d=[]
    m=[]
    x=[]
    y=[]
    z=[]
    vx=[]
    vy=[]
    vz=[]

    over_head=False
    ntot=0
    tphys=0.

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

        i_d.append(j+1)
        m.append(float(data[0])*mcon)
        x.append(float(data[1]))
        y.append(float(data[2]))
        z.append(float(data[3]))
        vx.append(float(data[4])*vcon)
        vy.append(float(data[5])*vcon)
        vz.append(float(data[6])*vcon)

    if ntot>0:

        cluster.add_stars(i_d,m,x,y,z,vx,vy,vz)

        if ofile==None:
            cluster.find_centre()
        else:
            get_cluster_orbit(cluster,ofile,**kwargs)

        cluster.to_cluster()
        cluster.find_centre()
        cluster.to_galaxy()

    return cluster

#Get StarCluster from NBODY6 using Jarrod Hurley's version of hrplot.f
def get_nbody6_jarrod(fort82,fort83,ofile=None,**kwargs):
    
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
    vstar=0.06557*np.sqrt(zmbar/rbar)
    rscale=float(line3[2])
    nb=int(line1b[0])
    ntot=ns+nb

    nsbnd=0
    nbbnd=0
    nbnd=0
    

    i_d=[]
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
        i_d.append(id1[-1])
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
        i_d.append(int(data[0]))
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
        cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',ctype='nbody6se',sfile=fort83,bfile=fort82)
        cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd)
        cluster.add_stars(i_d,m,x,y,z,vx,vy,vz)
        cluster.add_se(kw,logl,logr,ep,ospin)
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2,ep1,ep2,ospin1,ospin2)
        cluster.add_energies(kin,pot,etot)
    
        if ofile!=None:
            get_cluster_orbit(cluster,ofile,**kwargs)

        #Estimate centre
        cluster.find_centre()
    
    else:
        cluster=StarCluster(0,tphys)

    return cluster

#Get StarCluster from custom version of OUT34 and OUT9 (if binaries)
def get_nbody6_out(out9,out34,**kwargs):
    
    line1=out34.readline().split()
    if (len(line1)==0):
        print('END OF FILE')
        return StarCluster(0,0.0)

    line2=out34.readline().split()
    line3=out34.readline().split()

    ns=int(line1[0])
    tphys=float(line1[1])
    n_p=int(line1[4])
    if len(line1)>11:
        nb=int(float(line1[11]))
    else:
        nb=0

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
    vstar=0.06557*np.sqrt(zmbar/rbar)
    rscale=float(line3[2])
    ntot=ns+nb

    #Orbital Properties 
    xgc=float(line1[5])
    ygc=float(line1[6])
    zgc=float(line1[7])
    vxgc=float(line1[8])
    vygc=float(line1[9])
    vzgc=float(line1[10])

    nsbnd=0
    nbbnd=0
    nbnd=0

    i_d=[]
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

            x.append((x1*m1[-1]+x2*m2[-1])/(m1[-1]+m2[-1])+xc)
            y.append((y1*m1[-1]+y2*m2[-1])/(m1[-1]+m2[-1])+yc)
            z.append((z1*m1[-1]+z2*m2[-1])/(m1[-1]+m2[-1])+zc)
            vx.append((vx1*m1[-1]+vx2*m2[-1])/(m1[-1]+m2[-1]))
            vy.append((vy1*m1[-1]+vy2*m2[-1])/(m1[-1]+m2[-1]))
            vz.append((vz1*m1[-1]+vz2*m2[-1])/(m1[-1]+m2[-1]))
            m.append(m1[-1]+m2[-1])
            i_d.append(id1[-1])
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
            i_d.append(int(float(data[0])))
            kw.append(int(data[1]))
            m.append(float(data[2]))
            logl.append(float(data[3]))
            logr.append(float(data[4]))
            x.append(float(data[5])+xc)
            y.append(float(data[6])+yc)
            z.append(float(data[7])+zc)
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

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',ctype='nbody6',sfile=out34,bfile=out9)
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,n_p)
    #Add back on the centre of mass which has been substracted off by NBODY6
    cluster.add_stars(i_d,m,x,y,z,vx,vy,vz)
    cluster.add_se(kw,logl,logr,np.zeros(nbnd),np.zeros(nbnd))
    cluster.add_energies(kin,pot,etot)
    if out9!=None:
        cluster.add_bse(id1,id2,kw1,kw2,kcm,ecc,pb,semi,m1,m2,logl1,logl2,logr1,logr2)

    #Estimate centre of distribution
    cluster.find_centre()

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)

    return cluster

#Get StarCluster from custom version of OUT34
def get_nbody6_out34(out34,**kwargs):
    
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
    vstar=0.06557*np.sqrt(zmbar/rbar)
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
    
    i_d=[]
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
            i_d.append(int(float(data[0])))
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

    cluster=StarCluster(nbnd,tphys,units='nbody',origin='cluster',ctype='nbody6',sfile=out34,bfile=None)
    cluster.add_nbody6(nc,rc,rbar,rtide,xc,yc,zc,zmbar,vstar,rscale,nsbnd,nbbnd,n_p)
    #Add back on the centre of mass which has been substracted off by NBODY6
    cluster.add_stars(i_d,m,x+xc,y+yc,z+zc,vx,vy,vz)
    cluster.add_se(kw,logl,logr,np.zeros(nbnd),np.zeros(nbnd))
    cluster.add_energies(kin,pot,etot)
    
    #Estimate centre of distribution using median function
    cluster.find_centre()

    cluster.add_orbit(xgc,ygc,zgc,vxgc,vygc,vzgc)


    return cluster

#Get StarCluster snapshots produces by snapauto
def get_nbody6_snapauto(filename=None,units='realpc',origin='cluster',ofile=None,**kwargs):
    
    nsnap=kwargs.get('nsnap','1')
    nzfill=int(kwargs.get('nzfill',5))   
    ocontinue=kwargs.get('ocontinue',False)
    wdir=kwargs.get('wdir','./')
    nskip=int(kwargs.get('nskip',1))
    snapdir=kwargs.get('snapdir','snapshot/')
    snapbase=kwargs.get('snapbase','snap.dat.')
    snapend=kwargs.get('snapend','')

    if filename!=None:
        if os.path.isfile('%s%s%s' % (wdir,snapdir,filename)):
            snap=np.loadtxt('%s%s%s' %(wdir,snapdir,filename),skiprows=nskip)
            data=open('%s%s%s' %(wdir,snapdir,filename),'r').readline().split()
        else:
            print('NO FILE FOUND')
            cluster=StarCluster()
            print(cluster.ntot)
            return cluster

    elif os.path.isfile('%s%s%s%s%s' % (wdir,snapdir,snapbase,str(nsnap).zfill(nzfill),snapend)):
        snap=np.loadtxt(('%s%s%s%s%s' % (wdir,snapdir,snapbase,str(nsnap).zfill(nzfill),snapend)),skiprows=nskip)
        data=open('%s%s%s%s%s' % (wdir,snapdir,snapbase,str(nsnap).zfill(nzfill),snapend),'r').readline().split()
    else:
        print('NO FILE FOUND')
        cluster=StarCluster()
        print(cluster.ntot)
        return cluster

    if len(data)<=5:
        tphys=float(data[0])
        ntot=float(data[1])
        xgc=int(float(data[2]))
        ygc=float(data[3])
        zgc=float(data[4])
    else:
        tphys=float(data[0])
        ntot=float(data[1])
        nb=int(float(data[2]))
        xc=float(data[3])
        yc=float(data[4])
        zc=float(data[5])

    m=snap[:,0]
    x=snap[:,1]
    y=snap[:,2]
    z=snap[:,3]
    vx=snap[:,4]
    vy=snap[:,5]
    vz=snap[:,6]
    i_d=snap[:,7]
        
    nbnd=len(m)

    cluster=StarCluster(nbnd,tphys,units='realpc',origin='cluster',ctype='snapauto',nsnap=nsnap,nzfill=nzfill,wdir=wdir)
    cluster.add_stars(i_d,m,x,y,z,vx,vy,vz)
    
    #Estimate centre of distribution
    cluster.find_centre()

    if ofile!=None:
        get_cluster_orbit(cluster,ofile,**kwargs)

    return cluster

#Get StarCluster from standard nbodypy snapshot
def get_nbodypy_snapshot(filename=None,units='realpc',origin='cluster',ofile=None,**kwargs):
   
    nsnap=int(kwargs.get('nsnap','0'))
 
    nzfill=int(kwargs.get('nzfill',5))
    delimiter=kwargs.get('delimiter',None)
    wdir=kwargs.get('wdir','./')
    ocontinue=kwargs.get('ocontinue',False)
    snapdir=kwargs.get('snapdir','snaps/')
    snapbase=kwargs.get('snapbase','')
    snapend=kwargs.get('snapend','.dat')
    nskip=int(kwargs.get('nskip','0'))

    if filename!=None:
        if os.path.isfile('%s%s%s' % (wdir,snapdir,filename)):
            data=np.loadtxt('%s%s%s' %(wdir,snapdir,filename),delimiter=delimiter,skiprows=nskip)
        elif os.path.isfile('%s%s' % (wdir,filename)):
            data=np.loadtxt('%s%s' %(wdir,filename),delimiter=delimiter,skiprows=nskip)
        else:
            print('NO FILE FOUND')
            cluster=StarCluster()
            print(cluster.ntot)
            return cluster
    elif os.path.isfile('%s%s%s%s%s' % (wdir,snapdir,snapbase,str(nsnap).zfill(nzfill),snapend)):
        data=np.loadtxt(('%s%s%s%s%s' % (wdir,snapdir,snapbase,str(nsnap).zfill(nzfill),snapend)),delimiter=delimiter,skiprows=nskip)
    elif os.path.isfile('%s%s%s%s' % (wdir,snapbase,str(nsnap).zfill(nzfill),snapend)):
        data=np.loadtxt(('%s%s%s%s' % (wdir,snapbase,str(nsnap).zfill(nzfill),snapend)),delimiter=delimiter,skiprows=nskip)
    else:
        print('NO FILE FOUND - %s%s%s%s%s' % (wdir,snapdir,snapbase,str(nsnap).zfill(nzfill),snapend) )
        cluster=StarCluster()
        print(cluster.ntot)
        return cluster

    m=data[:,0]
    x=data[:,1]
    y=data[:,2]
    z=data[:,3]
    vx=data[:,4]
    vy=data[:,5]
    vz=data[:,6]
    i_d=data[:,7]
    kw=data[:,8]
    
    nbnd=len(m)

    cluster=StarCluster(nbnd,0.0,units=units,origin=origin,ctype='nbodypy',sfile=filename,bfile=None,delimiter=delimiter,nsnap=nsnap,nzfill=nzfill,wdir=wdir)
    cluster.add_stars(i_d,m,x,y,z,vx,vy,vz)
    cluster.kw=kw

    if origin=='galaxy':
        if ofile==None:
            cluster.find_centre()
        else:
            get_cluster_orbit(cluster,ofile,**kwargs)

        cluster.to_cluster()
        cluster.find_centre()
        cluster.to_galaxy()

    elif origin=='cluster':
        #Estimate centre of distribution
        cluster.find_centre()

        if ofile!=None:
            get_cluster_orbit(cluster,ofile,**kwargs)

    return cluster

#Get generic from standard nbodypy snapshot
def get_snapshot(filename,col_names=['m','x','y','z','vx','vy','vz'],col_nums=[0,1,2,3,4,5,6],units='realpc',origin='cluster',ofile=None,**kwargs):
    
    delimiter=kwargs.get('delimiter',None)
    wdir=kwargs.get('wdir','./')
    ocontinue=kwargs.get('ocontinue',False)
    skiprows=kwargs.get('skiprows',0)
    col_names=np.array(col_names)
    col_nums=np.array(col_nums)

    if units=='WDunits':
        vcon=220.0/bovy_conversion.velocity_in_kpcGyr(220.0,8.0)
        mcon=222288.4543021174
        units='realkpc'
    else:
        vcon=1.
        mcon=1.


    snapdir=kwargs.get('snapdir','')
   
    if os.path.isfile(('%s%s%s') % (wdir,snapdir,filename)):
        data=np.loadtxt('%s%s%s' % (wdir,snapdir,filename),delimiter=delimiter,skiprows=skiprows)
    else:
        print('NO SNAPSHOT FOUND')
        return 0

    mindx=np.argwhere(col_names=='m')
    m=data[:,col_nums[mindx]]*mcon

    xindx=np.argwhere(col_names=='x')
    x=data[:,col_nums[xindx]]
    yindx=np.argwhere(col_names=='y')
    y=data[:,col_nums[yindx]]
    zindx=np.argwhere(col_names=='z')
    z=data[:,col_nums[zindx]]
    vxindx=np.argwhere(col_names=='vx')
    vx=data[:,col_nums[vxindx]]*vcon
    vyindx=np.argwhere(col_names=='vy')
    vy=data[:,col_nums[vyindx]]*vcon
    vzindx=np.argwhere(col_names=='vz')
    vz=data[:,col_nums[vzindx]]*vcon

    if 'id' in col_names:
        idindx=np.argwhere(col_names=='id')
        i_d=data[:,col_nums[idindx]]
    else:
        i_d=np.linspace(0,len(x),len(x))+1

    if 'kw' in col_names:
        kwindx=np.argwhere(col_names=='kd')
        kw=data[:,col_nums[kwindx]]
    else:
        kw=np.zeros(len(x))   
    
    nbnd=len(m)

    cluster=StarCluster(nbnd,0.0,units=units,origin=origin,ctype='snapshot',sfile=filename,bfile=None,delimiter=delimiter,wdir=wdir)
    cluster.add_stars(i_d,m,x,y,z,vx,vy,vz)
    cluster.kw=kw

    if origin=='galaxy':
        if ofile==None:
            cluster.find_centre()
        else:
            get_cluster_orbit(cluster,ofile,**kwargs)

        cluster.to_cluster()
        cluster.find_centre()
        cluster.to_galaxy()

    elif origin=='cluster':
        #Estimate centre of distribution
        cluster.find_centre()

        if ofile!=None:
            get_cluster_orbit(cluster,ofile,**kwargs)

    return cluster

def get_mycode():
    cluster=StarCluster(0.0)
    return cluster

