#Routines for analysing Nbody models as if they were Observations
import numpy as np
from ..main.operations import save_cluster,return_cluster
from ..util.recipes import *
from ..util.plots import *

def obsrbinmaker(r,rm,obs_mask):
    """
    NAME:

       obsrbinmaker (WORK IN PROGRESS)

    PURPOSE:

       Radially bin data based on a set of observations
       --> Want to include spatial and magnitude/mass cuts 

    INPUT:

       r - stellar radii

       rm - cluster half-mass radius

       obs_mask - name of observational mask to be placed on cluster (Only current working mask is M30)

    OUTPUT:

       r_lower,r_mean,r_upper,r_hist

    HISTORY:

       2018 - Written - Webb (UofT)
    """

    if obs_mask=='M30':
        rh=61.800000000000004
        #In arcseconds:
        r_lower=np.array([10.0,20.0,40.0,200.0,250.0,350.0,650.0])
        r_upper=np.array([20.0,40.0,100.0,250.0,350.0,650.0,1000.0])

        r_lower=rm*r_lower/rh
        r_upper=rm*r_upper/rh

        r_hist=np.zeros(len(r_lower))
        r_sum=np.zeros(len(r_lower))

        for j in range(0,len(r_lower)):
            indx=(r>=r_lower[j]) * (r<=r_upper[j])
            r_hist[j]=len(r[indx])
            r_sum[j]=np.sum(r[indx])

        r_mean=[]
        for i in range(0,len(r_lower)):
            if r_hist[i]>0:
                r_mean=np.append(r_mean,r_sum[i]/r_hist[i])
            else:
                r_mean=np.append(r_mean,(r_lower[i]+r_upper[i])/2.0)

    return r_lower,r_mean,r_upper,r_hist

def obs_mass_function(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,mcorr=None,projected=False,obs_cut=None,plot=False,**kwargs):
    """
    NAME:

       obs_mass_function

    PURPOSE:

       Find mass function over a given mass range using nmass bins containing an equal number of stars.
       obs_mass_function differs from mass_function as it incorporates a mass completeness correction term and returns the completeness correction per bin

    INPUT:

       cluster - StarCluster instance
       mmin/mmax - specific mass range
       nmass - number of mass bins used to calculate alpha
       rmin/rmax - specific radial range
       vmin/vmax - specific velocity range
       emin/emax - specific energy range
       kwmin/kwmax - specific stellar evolution type range
       indx - specific subset of stars
       mcorr - correction for masses
       projected - use projected values
       obs_cut - trim data to match an observed dataset (WIP)
       plot - plot the mass function
       **kwargs - key words for plotting

    OUTPUT:

       m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha

    HISTORY:

       2018 - Written - Webb (UofT)
    """

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

    if np.sum(indx) >= nmass:

        m_lower,m_mean,m_upper,m_hist=nbinmaker(cluster.m[indx],nmass)       

        m_corr_hist=np.zeros(nmass)
        for i in range(0,len(m_hist)):
          mindx=(cluster.m>=m_lower[i]) * (cluster.m<=m_upper[i]) * indx
          m_corr_hist[i]=np.sum(1.0/mcorr[mindx])

        mbincorr=m_hist/m_corr_hist

        lm_mean=np.log10(m_mean)
        dm=m_hist/(m_upper-m_lower)
        ldm=np.log10(dm)

        (alpha,yalpha),V=np.polyfit(lm_mean,ldm,1,cov=True)
        ealpha=np.sqrt(V[0][0])
        eyalpha=np.sqrt(V[1][1])

        if plot:
            filename=kwargs.get('filename',None)
            nplot(m_mean,np.log10(dm),xlabel='M',ylabel='LOG(dN/dM)',**kwargs)
            mfit=np.linspace(np.min(m_mean),np.max(m_mean),nmass)
            dmfit=10.0**(alpha*np.log10(mfit)+yalpha)
            nlplot(mfit,np.log10(dmfit),overplot=True,label=(r'$\alpha$ = %f' % alpha))

            plt.legend()

            if filename!=None:
                plt.savefig(filename)


        return m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha,mbincorr
    else:
        print('NOT ENOUGH STARS TO ESTIMATE MASS FUNCTION')
        return np.zeros(nmass),np.zeros(nmass),np.zeros(nmass),-1000.,-1000.,-1000.,-1000.

def obs_alpha_prof(cluster,mmin=None,mmax=None,nmass=10,rmin=None,rmax=None,nrad=20,vmin=None,vmax=None,emin=None,emax=None,kwmin=0,kwmax=1,indx=None,mcorr=None,projected=False,obs_cut=None,plot=False,**kwargs):
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
    cluster.to_centre(do_order=True,do_key_params=True)

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

        m_mean,m_hist,dm,alpha,ealpha,yalpha,eyalpha=dx_corr_function(cluster.m[rindx],nmass,mcorr[rindx])

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

    return_cluster(cluster,units0,origin0,do_order=True,do_key_params=True)

    return lrprofn,aprof,dalpha,edalpha,ydalpha,eydalpha

def dx_corr_function(x, nx=10, xcorr=None, bintype='num'):
  """
  NAME:

     dx_function

  PURPOSE:

     Find distribution function using nx bins containing an equal number of points
  INPUT:

     x - input array

     nx - number of bins

     xcorr - correction function for x (default:None)

     bintype - bin with equal number of stars per bin (bin) or evenly in x (fix) (default: num)

  OUTPUT:

     x_mean,x_hist,dx,alpha,ealpha,yalpha,eyalpha


  HISTORY:

     2018 - Written - Webb (UofT)

  """

  if bintype == 'num':
    x_lower, x_mean, x_upper, x_hist = nbinmaker(x, nx)
  else:
    x_lower, x_mean, x_upper, x_hist = binmaker(x, nx)

  if xcorr is not None:
    x_hist = np.zeros(nx)
    for i in range(0, len(x_hist)):
      indx = (x >= x_lower[i]) * (x <= x_upper[i])
      x_hist[i] = np.sum(1.0 / xcorr[indx])

  lx_mean = np.log10(x_mean)
  dx = x_hist / (x_upper - x_lower)
  ldx = np.log10(dx)

  (alpha, yalpha), V = np.polyfit(lx_mean, ldx, 1, cov=True)
  ealpha = np.sqrt(V[0][0])
  eyalpha = np.sqrt(V[1][1])

  return x_mean, x_hist, dx, alpha, ealpha, yalpha, eyalpha
