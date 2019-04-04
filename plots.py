#Make key plots given StarCluster
#import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from galpy.util import bovy_plot
import seaborn as sns
import os
from scipy.ndimage import gaussian_filter
import matplotlib.colors as colors

bovy_plot.bovy_print(axes_labelsize=18.,xtick_labelsize=14.,ytick_labelsize=14.)
current_palette = sns.color_palette()


# @mpl_decorator(fig='new', )
# def _nplot(x, y, alpha=.1, **kw):
#     splt.plot(x, y, **kw)

#     if xlim!=None:
#         plt.xlim(xlim)
#         if ylim==None and scale:
#             xindx=(x>=xlim[0]) * (x <=xlim[1])
#             ymin=np.min(y[xindx])
#             ymax=np.max(y[xindx])
#             plt.ylim((ymin,ymax))

#     if ylim!=None:
#         plt.ylim(ylim)
#         if xlim==None and scale:
#             yindx=(y>=ylim[0]) * (y <=ylim[1])
#             xmin=np.min(x[yindx])
#             xmax=np.max(x[yindx])
#             plt.xlim((xmin,xmax))

def nscatter(x,y,z=None,xlabel='',ylabel='',legend=False,title='',xlim=None,ylim=None,scale=True,filename=None,overplot=False,**kwargs):

    x=np.asarray(x)
    y=np.asarray(y)
  
    if not overplot:
        plt.figure()

    alpha=kwargs.pop('alpha',1.0)

    if z!=None:
        z=np.asarray(z)
        out=plt.scatter(x,y,c=z,alpha=alpha,**kwargs)
        plt.colorbar(out)
    else:
        out=plt.scatter(x,y,alpha=alpha,**kwargs)

    if overplot:
        return out
    else:
        if xlim!=None:
            plt.xlim(xlim)
            if ylim==None and scale:
                xindx=(x>=xlim[0]) * (x <=xlim[1])
                ymin=np.min(y[xindx])
                ymax=np.max(y[xindx])
                plt.ylim((ymin,ymax))

        if ylim!=None:
            plt.ylim(ylim)
            if xlim==None and scale:
                yindx=(y>=ylim[0]) * (y <=ylim[1])
                xmin=np.min(x[yindx])
                xmax=np.max(x[yindx])
                plt.xlim((xmin,xmax))

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        if legend:
            plt.legend()

        plt.tight_layout()
        
        if filename!=None:
            plt.savefig(filename)
            return 0
        else:
            return 0

def nplot(x,y,ptype='.',xlabel='',ylabel='',legend=False,title='',xlim=None,ylim=None,scale=True,log=False,logx=False,logy=False,filename=None,overplot=False,**kwargs):

    if not overplot:
        plt.figure()

    alpha=kwargs.pop('alpha',1.0)

    if log or (logx and logy):
        out=plt.loglog(x,y,ptype,alpha=alpha,**kwargs)
    elif logx:
        out=plt.semilogx(x,y,ptype,alpha=alpha,**kwargs)
    elif logy:
        out=plt.semilogy(x,y,ptype,alpha=alpha,**kwargs)
    else:
        out=plt.plot(x,y,ptype,alpha=alpha,**kwargs)

    if overplot:
        return out
    else:
        if xlim!=None:
            plt.xlim(xlim)
            if ylim==None and scale:
                xindx=(x>=xlim[0]) * (x <=xlim[1])
                ymin=np.min(y[xindx])
                ymax=np.max(y[xindx])
                plt.ylim((ymin,ymax))

        if ylim!=None:
            plt.ylim(ylim)
            if xlim==None and scale:
                yindx=(y>=ylim[0]) * (y <=ylim[1])
                xmin=np.min(x[yindx])
                xmax=np.max(x[yindx])
                plt.xlim((xmin,xmax))

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)

        if legend:
            plt.legend()


        plt.tight_layout()
        
        if filename!=None:
            plt.savefig(filename)
            return 0
        else:
            return 0

def nlplot(x,y,ltype='-',xlabel='',ylabel='',legend=False,title='',xlim=None,ylim=None,scale=True,log=False,logx=False,logy=False,filename=None,overplot=False,**kwargs):

    return nplot(x,y,ptype=ltype,xlabel=xlabel,ylabel=ylabel,legend=legend,title=title,xlim=xlim,ylim=ylim,scale=scale,log=log,logx=logx,logy=logy,filename=filename,overplot=overplot,**kwargs)


def nhist(x,nbin=10,xlabel='',legend=False,title='',xlim=None,fill=False,filename=None,overplot=False,**kwargs):
    
    nbin=kwargs.pop('bins',nbin)

    if not overplot:
        plt.figure()

    if fill:
        out=plt.hist(x,bins=nbin,**kwargs)
    else:
        histtype=kwargs.pop('histtype','step')
        out=plt.hist(x,bins=nbin,histtype=histtype,**kwargs)

    if overplot:
        return out
    else:
        if xlim!=None:
            plt.xlim(xlim)

        plt.xlabel(xlabel)
        plt.title(title)

        if legend:
            plt.legend()

        if filename!=None:
            plt.savefig(filename)
            return 0
        else:
            return 0

def nhist2d(x,y,nbin=10,xlabel='',ylabel='',title='',xlim=None,ylim=None,filename=None,**kwargs):
    

    if xlim==None:
        xlim=[np.min(x),np.max(x)]
    if ylim==None:
        ylim=[np.min(y),np.max(y)]
   
    cmap=kwargs.pop('cmap',plt.cm.jet)

    plt.hist2d(x,y,bins=nbin,range=[xlim,ylim],cmap=cmap,**kwargs)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    plt.title(title)
    plt.colorbar()

    if filename!=None:
        plt.savefig(filename)
        return 0
    else:
        return 0

def ndens(x,y,z=[None],normed=False,xstep=1.,ystep=1.,sigx=4.,sigy=4.,order=0,extent=None,interpolation='nearest',cmap='Spectral'):

    xbins = np.arange(x.min(), x.max()+xstep, xstep)
    ybins = np.arange(y.min(), y.max()+ystep, ystep)

    if None in z:
        H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins]) 
    else:
        H, xedges, yedges = np.histogram2d(x, y, bins=[xbins, ybins],normed=normed,weights=z) 

    fimg = gaussian_filter(H, sigma=(sigx, sigy), order=order)

    if extent==None:
        extent=[x.min(),x.max(),y.min(),y.max()]

    plt.imshow(fimg, norm=colors.LogNorm(), interpolation=interpolation, cmap=cmap,extent=extent)
    plt.colorbar()

    return 0

def starplot(cluster,coords='xyz',xlim=None,ylim=None,legend=False,filename=None,overplot=False,do_center=False,**kwargs):

    alpha=kwargs.pop('alpha',0.1)

    if coords!='xyz':

        if not overplot:
            plt.figure()

        if coords=='xz':
            x=cluster.x
            y=cluster.z
            if do_center:
                xgc=cluster.xgc
                ygc=cluster.zgc
        elif coords=='yz':
            x=cluster.y
            y=cluster.z
            if do_center:
                xgc=cluster.ygc
                ygc=cluster.zgc
        else:
            x=cluster.x
            y=cluster.y
            if do_center:
                xgc=cluster.xgc
                ygc=cluster.ygc

        out=plt.plot(x,y,'.',alpha=alpha,**kwargs)

        if overplot:
            return out
        else:
            if do_center:
                plt.plot(xgc,ygc,'.',alpha=1.,label='COM',**kwargs)       
            if cluster.units=='nbody':
                units='(NBODY)'
            elif cluster.units=='realpc':
                units='(pc)'
            elif cluster.units=='realkpc':
                units='(kpc)'
            elif cluster.units=='galpy':
                units='(GALPY)'
            else:
                units=''

            if coords=='xz':
                x=cluster.x
                y=cluster.z
                plt.xlabel('X '+units)
                plt.ylabel('Z '+units)
            elif coords=='yz':
                x=cluster.y
                y=cluster.z
                plt.xlabel('Y '+units)
                plt.ylabel('Z '+units)
            else:
                x=cluster.x
                y=cluster.y
                plt.xlabel('X '+units)
                plt.ylabel('Y '+units)

            plt.title('Time = %f' % cluster.tphys)

            if xlim!=None:
                plt.xlim(xlim)
            if ylim!=None:
                plt.ylim(ylim)

            plt.tight_layout()

            if filename!=None:
                plt.savefig(filename)
                return 0
            else:
                return 0

    else:
        
        plt.subplot(1,2,1)
        
        plt.plot(cluster.x, cluster.z,'.',alpha=alpha,**kwargs)
        
        if cluster.origin=='galaxy':
            plt.plot(cluster.xgc,cluster.zgc,'o',label='Center')
            plt.plot(cluster.xgc+cluster.xc,cluster.zgc+cluster.zc,'o',label='COM')
        else:
            plt.plot(cluster.xc,cluster.zc,'o',label='COM')


        if xlim!=None:
            plt.xlim(xlim)
        if ylim!=None:
            plt.ylim(ylim)
        plt.xlabel('X (kpc)')
        plt.ylabel('Z (kpc)')

        plt.subplot(1,2,2)

        plt.plot(cluster.x, cluster.y,'.',alpha=alpha,**kwargs)

        if cluster.origin=='galaxy':
            plt.plot(cluster.xgc,cluster.ygc,'o',label='Center')
            plt.plot(cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,'o',label='COM')
        else:
            plt.plot(cluster.xc,cluster.yc,'o',label='COM')

        if xlim!=None:
            plt.xlim(xlim)
        if ylim!=None:
            plt.ylim(ylim)
        plt.xlabel('X (kpc)')
        plt.ylabel('Y (kpc)')

        if legend:
            plt.legend()

        plt.tight_layout()

        if filename!=None:
            plt.savefig(filename)
            return 0
        else:
            return 0

    return 0

def skyplot(cluster,coords='radec',xlim=None,ylim=None,legend=False,filename=None,overplot=False,do_center=False,**kwargs):

    cluster.to_sky()
    alpha=kwargs.pop('alpha',0.1)

    if not overplot:
        plt.figure()

    if coords=='radec':
        x=cluster.ra
        y=cluster.dec
    elif coords=='pm':
        x=cluster.pmra
        y=cluster.pmdec

    out=plt.plot(x,y,'.',alpha=alpha,**kwargs)

    if overplot:
        return out
    else:
        if do_center:
            xgc,ygc=cluster.ragc,cluster.decgc
            plt.plot(xgc,ygc,'.',alpha=1.,label='COM',**kwargs)

        if coords=='radec':
            plt.xlabel('Ra (degree)')
            plt.ylabel('Dec (degree)')
        elif coords=='pm':
            plt.xlabel(r'$\mu_{Ra}$ (mas $yr^{-1}$)')
            plt.ylabel(r'$\mu_{Dec}$ (mas $yr^{-1}$)')

        plt.title('Time = %f' % cluster.tphys)

        if xlim!=None:
            plt.xlim(xlim)
        if ylim!=None:
            plt.ylim(ylim)

        plt.tight_layout()

        if filename!=None:
            plt.savefig(filename)
            return 0
        else:
            return 0

def mplot(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=10,kwmin=0,kwmax=15,projected=False,cumulative=False,obs_cut=None,log=False, filename=None,overplot=False,**kwargs):

    if not overplot:
        plt.figure()

    if cluster.units=='nbody':
        xunits='(NBODY)'
        yunits='(NBODY)'
    elif cluster.units=='realpc':
        xunits='(pc)'
        yunits='Msun'

    elif cluster.units=='realkpc':
        xunits='(kpc)'
        yunits='Msun'
    elif cluster.units=='galpy':
        xunits='(GALPY)'
        yunits='(GALPY)'

    else:
        xunits=''
        yunits=''

    x,y,n=m_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=10,kwmin=0,kwmax=15,projected=False,cumulative=cumulative,obs_cut=None)

    if log:
        out=plt.loglog(x,y,**kwargs)
    else:
        out=plt.plot(x,y,**kwargs)

    if overplot:
        return out
    else:
        if projected:
            plt.xlabel('R'+xunits)
        else:
            plt.xlabel('r'+xunits)

        plt.ylabel('M'+yunits)

        plt.title('Time = %f' % cluster.tphys)

        if filename!=None:
            plt.savefig(filename)
            return 0
        else:
            return 0

def explot(data,prefix='',**kwargs):
    #Plot data saved in extrct.npy

    t=data[:,2]
    m=data[:,4]
    rm=data[:,9]
    
    plt.subplot(2,1,1)

    plt.plot(t,m,**kwargs)
    plt.xlabel('Time (Myr)')
    plt.ylabel('Mass (Msun)')

    plt.subplot(2,1,2)

    plt.plot(t,rm,**kwargs)
    plt.xlabel('Time (Myr)')
    plt.ylabel('rm (pc)')

    plt.savefig(prefix+'explot.png')
    
    plt.close()

def dvplot(data,nsnap=0,tsnap=None,prefix='',nrad=10,filename=None,**kwargs):
    #Plot data saved in dvprof.npy

    t=data[:,0]
    m=data[:,1]
    rm=data[:,-1]
    nstart=0.
    nstop=0.

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])
    
    rbin=data[nsnap,2:2+nrad]
    sig=data[nsnap,2+nrad:2+nrad+nrad]

    plt.plot(rbin,sig,**kwargs)
    plt.xlabel(r'ln(r/$r_m$)')
    plt.ylabel(r'$\sigma_v$ (km/s)')

    if filename!=None:
        plt.savefig('%sdvplot%s.png' % (prefix,str(i)))
        return 0
    else:
        return 0

def bplot(data,nsnap=0,tsnap=None,prefix='',nrad=10,filename=None,**kwargs):
    #Plot data saved in dvprof.npy

    t=data[:,0]
    m=data[:,1]
    rm=data[:,-1]
    nstart=0.
    nstop=0.

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    rbin=data[nsnap,2:2+nrad]
    beta=data[nsnap,2+nrad+nrad:2+nrad+nrad+nrad]

    plt.plot(rbin,beta,**kwargs)
    plt.xlabel(r'ln(r/$r_m$)')
    plt.ylabel(r'$\beta$')

    if filename!=None:
        plt.savefig('%sbplot%s.png' % (prefix,str(i)))
        return 0
    else:
        return 0        