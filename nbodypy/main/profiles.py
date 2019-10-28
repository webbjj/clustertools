#Determine radial profiles of key properties

import numpy as np
from galpy.util import bovy_coords

from ..util.constants import *
from ..util.recipes import *
from .operations import *
from ..util.plots import *
from ..util.coordinates import rect_to_sphere

def rho_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
    """
    NAME:

       rho_prof

    PURPOSE:

       Measure the density profile of the cluster

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

        Same as for ..util.plot.nplot

    OUTPUT:

        rprof,pprof,nprof (radius, density, number of stars)

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()


    rprof=np.array([])
    pprof=np.array([])
    nprof=np.array([])

    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r <= rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(r[indx],nrad)

    for i in range(0,len(r_mean)):
        rindx=indx * (r >= r_lower[i]) * (r <= r_upper[i])
        rprof=np.append(rprof,r_mean[i])
        if projected:
            vol=np.pi*(r_upper[i]**2-r_lower[i]**2.)
        else:
            vol=(4./3.)*np.pi*(r_upper[i]**3-r_lower[i]**3.)

        pprof=np.append(pprof,np.sum(cluster.m[rindx]/vol))
        nprof=np.append(nprof,np.sum(rindx))


    if plot:
        filename=kwargs.pop('filename',None)   
        overplot=kwargs.pop('overplot',False)        
     
        if cluster.units=='nbody':
            xunits=' (NBODY)'
            yunits=' (NBODY)'
        elif cluster.units=='realpc':
            xunits=' (pc)'
            if projected:
                yunits=' Msun/pc^2'
            else:
                yunits=' Msun/pc^3'
        elif cluster.units=='realkpc':
            xunits=' (kpc)'
            if projected:
                yunits=' Msun/kpc^2'
            else:
                yunits=' Msun/kpc^3'
        elif cluster.units=='galpy':
            xunits=' (GALPY)'
            yunits=' (GALPY)'

        else:
            xunits=''
            yunits=''

        x,y,n=rprof,pprof,nprof
        nlplot(x,y,xlabel=r'$R %s$' % xunits,ylabel=r'$\rho %s$' % yunits,title='Time = %f' % cluster.tphys,log=True,overplot=overplot,filename=filename)

        if filename!=None:
            plt.savefig(filename)

    return_cluster(cluster,units0,origin0)

    return rprof,pprof,nprof

def m_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,indx=None,projected=False,cumulative=False,obs_cut=None,plot=False,**kwargs):
    """
    NAME:

       m_prof

    PURPOSE:

       Measure the mass profile of the cluster

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       cumalitive - determine the cumulative mass profile instead (Default: False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    OUTPUT:

        rprof,mprof,nprof (radius, mass, number of stars)

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 
    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    rprof=[]
    mprof=[]
    nprof=[]

    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(r[indx],nrad)

    for i in range(0,len(r_mean)):
        if cumulative:
            rindx=indx * (r <= r_upper[i])
        else:
            rindx=indx * (r >= r_lower[i]) * (r <= r_upper[i])
        rprof.append(r_mean[i])

        mprof.append(np.sum(cluster.m[rindx]))
        nprof.append(np.sum(rindx))

    if plot:
        filename=kwargs.pop('filename',None)   
        overplot=kwargs.pop('overplot',False)        
     
        if cluster.units=='nbody':
            xunits=' (NBODY)'
            yunits=' (NBODY)'
        elif cluster.units=='realpc':
            xunits=' (pc)'
            yunits=' Msun'
        elif cluster.units=='realkpc':
            xunits=' (kpc)'
            yunits=' Msun'
        elif cluster.units=='galpy':
            xunits=' (GALPY)'
            yunits=' (GALPY)'
        else:
            xunits=''
            yunits=''

        x,y,n=rprof,mprof,nprof
        nlplot(x,y,xlabel=r'$R %s $' % xunits,ylabel=r'$M %s $' % yunits, title='Time = %f' % cluster.tphys,log=True,overplot=overplot,filename=filename)

        if filename!=None:
            plt.savefig(filename)


    return_cluster(cluster,units0,origin0)

    return rprof,mprof,nprof

def alpha_prof(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,mcorr=None,projected=False,obs_cut=None,plot=False,**kwargs):
    """
    NAME:

       alpha_prof

    PURPOSE:

       Measure the radial variation in the mass function

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       nmass - number of mass bins to calculate slope of mass function

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       mcorr - correction function for masses

       projected - use projected values and constraints (Default:False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    OUTPUT:

        lrprofn - natural log of radius (normalized by half-mass radius)

        aprof - slope of the mass function

        dalpha - delta_alpha = d(alpha)/d(ln(r/rm) 

        edalpha - error in dalpha

        ydalpha,eydalpha - y-intercept and error in fit to alpha vs ln(r/rm)

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    if mcorr is None:
        mcorr=np.ones(cluster.ntot)

    lrprofn=[]
    aprof=[]
    
    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(r[indx],nrad)

    for i in range(0,len(r_mean)):
        rindx=indx * (r >= r_lower[i]) * (r <= r_upper[i])

        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=dx_function(cluster.m[rindx],nmass,mcorr[rindx])

        if alpha > -100:
            if projected:
                lrprofn.append(np.log(r_mean[i]/cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i]/cluster.rm))

            aprof.append(alpha)

    if len(lrprofn)>3:
        (dalpha,ydalpha),V=np.polyfit(lrprofn,aprof,1,cov=True)
        edalpha=np.sqrt(V[0][0])
        eydalpha=np.sqrt(V[1][1])
    else:
        dalpha=-100.0
        ydalpha=0.0
        edalpha=0.0
        eydalpha=0.0

    if plot:
        filename=kwargs.pop('filename',None)
        overplot=kwargs.pop('overplot',False)        

        nplot(lrprofn,aprof,xlabel=r'$\ln(r/r_m)$',ylabel=r'$\alpha$',overplot=overplot,**kwargs)
        rfit=np.linspace(np.min(lrprofn),np.max(lrprofn),nrad)
        afit=dalpha*rfit+ydalpha
        nlplot(rfit,afit,overplot=True,label=(r'd$\alpha$ = %f' % dalpha))
        plt.legend()

        if filename!=None:
            plt.savefig(filename)

    cluster.dalpha=dalpha

    return_cluster(cluster,units0,origin0)

    return lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha

def sigv_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,projected=False,obs_cut=None):
    """
    NAME:

       sigv_prof

    PURPOSE:

       Measure the radial variation in the velocity dispersion

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    OUTPUT:

        lrprofn - natural log of radius (normalized by half-mass radius)

        sigvprof - velocity dispersion

        betaprof - anisotropy parameter 

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()


    lrprofn=[]
    sigvprof=[]
    betaprof=[]

    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    #Build subcluster containing only stars in the full radial and mass range:
    indx=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)
   
    #Convert to cylindrical or spherical coordinates:
    if projected:
        r,theta,z=bovy_coords.rect_to_cyl(cluster.x,cluster.y,cluster.z)
        vr,vtheta,vz=bovy_coords.rect_to_cyl_vec(cluster.vx,cluster.vy,cluster.vz,cluster.x,cluster.y,cluster.z)
    else:
        r,theta,phi,vr,vt,vp=rect_to_sphere(cluster)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(r[indx],nrad)

    for i in range(0,len(r_mean)):
        rindx=indx * (r >= r_lower[i]) * (r <= r_upper[i])

        if np.sum(rindx) > 3.:

            sigr=np.std(vr[rindx])
            sigt=np.std(vt[rindx])

            if projected:
                sigp=np.zeros(len(vr))
                beta=sigt/sigr-1.
            else:
                sigp=np.std(vp[rindx])
                beta=1.0-(sigt**2.0+sigp**2.0)/(2.*(sigr**2.0))
            
            sigv=np.sqrt(sigr**2.0+sigt**2.0+sigp**2.0)

            if projected:
                lrprofn.append(np.log(r_mean[i]/cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i]/cluster.rm))

            sigvprof.append(sigv)
            betaprof.append(beta)
           

    return_cluster(cluster,units0,origin0)
     
    return lrprofn,sigvprof,betaprof

def v_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,indx=None,projected=False,obs_cut=None):
    """
    NAME:

       v_prof

    PURPOSE:

       Measure the radial variation in the mean velocity 

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    OUTPUT:

        lrprofn - natural log of radius (normalized by half-mass radius)

        vprof - mean velocity

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    lrprofn=[]
    sigvprof=[]

    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)
  
    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    #Convert to cylindrical or spherical coordinates:
    if projected:
        r,theta,z=bovy_coords.rect_to_cyl(cluster.x,cluster.y,cluster.z)
        vr,vtheta,vz=bovy_coords.rect_to_cyl_vec(cluster.vx,cluster.vy,cluster.vz,cluster.x,cluster.y,cluster.z)
    else:
        r,theta,phi,vr,vt,vp=rect_to_sphere(cluster)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(r[indx],nrad)

    for i in range(0,len(r_mean)):
        rindx=indx * (r >= r_lower[i]) * (r <= r_upper[i])

        if np.sum(rindx) > 3.:

            vrmean=np.mean(vr[rindx])
            vtmean=np.mean(vt[rindx])

            if projected:
                vpmean=np.zeros(len(vr))
            else:
                vpmean=np.mean(vp[rindx])
            
            vmean=np.sqrt(vrmean**2.0+vtmean**2.0+vpmean**2.0)

            if projected:
                lrprofn.append(np.log(r_mean[i]/cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i]/cluster.rm))

            vprof.append(vmean)

    return_cluster(cluster,units0,origin0)

    return lrprofn,vprof

def eta_prof(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,projected=False,obs_cut=None):
    """
    NAME:

       eta_prof

    PURPOSE:

       Measure the radial variation in eta

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       nmass - number of mass bins over which to measure eta

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

       Same as for ..util.plot.nplot

    OUTPUT:

        lrprofn - natural log of radius (normalized by half-mass radius)

        eprof - slope of the sigma_v-mass function

        deta - deta = d(eta)/d(ln(r/rm) 

        edeta - error in deta

        ydeta,eydeta - y-intercept and error in fit to eta vs ln(r/rm)

    HISTORY:

       2018 - Written - Webb (UofT)

    """ 

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    lrprofn=[]
    eprof=[]
    
    if projected:
        r=cluster.rpro
        v=cluster.vpro
    else:
        r=cluster.r
        v=cluster.v

    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r<=rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)
  
    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    r_lower,r_mean,r_upper,r_hist=nbinmaker(cluster.r[indx],nrad)

    for i in range(0,len(r_mean)):
        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=mmin,mmax=mmax,nmass=nmass,rmin=r_lower[i],rmax=r_upper[i],vmin=vmin,vmax=vmax,kwmin=kwmin,kwmax=kwmax,projected=projected,obs_cut=obs_cut)

        if alpha > -100:
            if projected:
                lrprofn.append(np.log(r_mean[i]/cluster.rmpro))
            else:
                lrprofn.append(np.log(r_mean[i]/cluster.rm))

            eprof.append(eta)

    if len(lrprofn)>3:
        (deta,ydeta),V=np.polyfit(lrprofn,eprof,1,cov=True)
        edeta=np.sqrt(V[0][0])
        eydeta=np.sqrt(V[1][1])
    else:
        deta=-100.0
        ydeta=0.0
        edeta=0.0
        eydeta=0.0

    return_cluster(cluster,units0,origin0)

    return lrprofn,eprof,deta,edeta,ydeta,eydeta

def vcirc_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=15,indx=None,projected=False,obs_cut=None,plot=False,**kwargs):
    """
    NAME:

       vcirc_prof

    PURPOSE:

       Measure the circulr velocity profile of the cluster

    INPUT:

       cluster - StarCluster instance

       mmin/mmax - minimum and maximum stellar mass

       rmin/rmax - minimum and maximum stellar radii

       nrad - number of radial bins

       vmin/vmax - minimum and maximum stellar velocity

       emin/emax - minimum and maximum stellar energy

       kwmin/kwmax - minimum and maximum stellar type (kw)

       indx - user defined boolean array from which to extract the subset

       projected - use projected values and constraints (Default:False)

       obs_cut - apply an observational mask to the dataset (Default: False)

       plot - plot the density profile (Default: False)

    KWARGS:

        Same as for ..util.plot.nplot

    OUTPUT:

        rprof,vprof,rvmax,vmax (radius, circular velocity, radius of maximum virc, maximum vcirc)

    HISTORY:

       2019 - Written - Webb (UofT)

    """ 

    units0,origin0=save_cluster(cluster)
    cluster.to_centre()

    if cluster.units=='nbody':
        grav=1.0
    elif cluster.units=='realpc':
        #G has units of pc (km/s)^2 / Msun
        grav=4.302e-3
    elif cluster.units=='realkpc':
        #G has units of kpc (km/s)^2 / Msun
        grav=4.302e-6
    else:
        grav=1.0



    rprof=np.array([])
    vcprof=np.array([])

    if projected:
        if cluster.rproorder is None:
            cluster.key_params(do_order=True)
        r=cluster.rpro[cluster.rproorder]
        v=cluster.vpro[cluster.rproorder]
        m=cluster.m[cluster.rproorder]
    else:
        if cluster.rorder is None:
            cluster.key_params(do_order=True)
        r=cluster.r[cluster.rorder]
        v=cluster.v[cluster.rorder]
        m=cluster.m[cluster.rorder]


    if rmin==None: rmin=np.min(r)
    if rmax==None: rmax=np.max(r)
    if vmin==None: vmin=np.min(v)
    if vmax==None: vmax=np.max(v)
    if mmin==None: mmin=np.min(cluster.m)
    if mmax==None: mmax=np.max(cluster.m)

    if indx is None:
        indx=(cluster.id > -1)

    #Build subcluster containing only stars in the full radial and mass range:
    indx*=(r >= rmin) * (r <= rmax) * (cluster.m >= mmin) * (cluster.m <= mmax) * (v >=vmin) * (v <=vmax) * (cluster.kw >=kwmin) * (cluster.kw <=kwmax)

    if emin!=None:
        indx*=(cluster.etot >= emin)
    if emin!=None:
        indx*=(cluster.etot <= emax)

    r=r[indx]
    v=v[indx]
    m=m[indx]

    msum=np.cumsum(m)
    vcirc=np.sqrt(grav*msum/r)
    vmax=np.amax(vcirc)
    rvmax=r[np.argmax(vcirc)]

    rprof=r
    vcprof=vcirc

    if plot:
        filename=kwargs.pop('filename',None)   
        overplot=kwargs.pop('overplot',False)        
     
        if cluster.units=='nbody':
            xunits=' (NBODY)'
            yunits=' (NBODY)'
        elif cluster.units=='realpc':
            xunits=' (pc)'
            yunits=' km/s'
        elif cluster.units=='realkpc':
            xunits=' (kpc)'
            yunits=' km/s'

        elif cluster.units=='galpy':
            xunits=' (GALPY)'
            yunits=' (GALPY)'

        else:
            xunits=''
            yunits=''

        x,y=rprof,vcprof
        nlplot(x,y,xlabel=r'$R %s$' % xunits,ylabel=r'$vc %s $' % yunits,title='Time = %f' % cluster.tphys,log=True,overplot=overplot,filename=filename)
        nlplot([rvmax,rvmax],[np.amin(y),np.amax(y)],'--',overplot=True)
        nlplot([np.amin(x),np.amax(x)],[vmax,vmax],'--',overplot=True)

        if filename!=None:
            plt.savefig(filename)

    return_cluster(cluster,units0,origin0)

    return rprof,vcprof,rvmax,vmax

