import matplotlib.pyplot as plt
import numpy as np
from galpy.util import bovy_plot
import seaborn as sns
import matplotlib.colors as colors
from ..util.plots import *
from ..nbodypy.profiles import m_prof

bovy_plot.bovy_print(axes_labelsize=18.,xtick_labelsize=14.,ytick_labelsize=14.)
current_palette = sns.color_palette()

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

    out=plt.plot(t,rm,**kwargs)
    plt.xlabel('Time (Myr)')
    plt.ylabel('rm (pc)')

    plt.savefig(prefix+'explot.png')
    
    plt.close()

    return out

def dvplot(data,nsnap=0,tsnap=None,prefix='',nrad=20,filename=None,**kwargs):
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

    out=plt.plot(rbin,sig,**kwargs)
    plt.xlabel(r'ln(r/$r_m$)')
    plt.ylabel(r'$\sigma_v$ (km/s)')

    if filename!=None:
        plt.savefig('%sdvplot%s.png' % (prefix,str(i)))

    return out

def bplot(data,nsnap=0,tsnap=None,prefix='',nrad=20,filename=None,**kwargs):
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

    out=plt.plot(rbin,beta,**kwargs)
    plt.xlabel(r'ln(r/$r_m$)')
    plt.ylabel(r'$\beta$')

    if filename!=None:
        plt.savefig('%sbplot%s.png' % (prefix,str(nsnap)))

    return out       

def aplot(data,nsnap=0,tsnap=None,prefix='',nmass=10,nrad=20,filename=None,**kwargs):
    #Plot data saved in alpha_prof.npy

    t=data[:,0]
    alpha=data[:,1]
    ealpha=data[:,2]
    yalpha=data[:,3]
    eyalpha=data[:,4]

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    nstart,nstop=5,5+nmass
    mmean=data[nsnap,nstart:nstop]
  
    nstart=nstop
    nstop=nstart+nmass
    dm=data[nsnap,nstart:nstop]

    nstart=nstop
    nstop=nstart+nrad
    lrprof=data[nsnap,nstart:nstop]

    nstart=nstop
    nstop=nstart+nrad
    aprof=data[nsnap,nstart:nstop]

    nstart=nstop
    dalpha=data[nsnap,nstart]
    edalpha=data[nsnap,nstart+1]
    ydalpha=data[nsnap,nstart+2]
    eydalpha=data[nsnap,nstart+3]

    rfrac=kwargs.pop('rfrac',False)
    if rfrac:
        lrprof=np.linspace(0.,1.,len(aprof))
        xlabel=r'$r_n$'
    else:
        xlabel=r'ln(r/$r_m$)'

    if filename!=None:
        filename=('%saplot%s.png' % (prefix,str(nsnap)))
    

    out=nlplot(lrprof,aprof,xlabel=xlabel,ylabel=r'$\alpha$',filename=filename,**kwargs)



    return out 

def mplot(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,kwmin=0,kwmax=15,projected=False,cumulative=False,obs_cut=None,log=False, filename=None,overplot=False,**kwargs):

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

    x,y,n=m_prof(cluster,mmin=None,mmax=None,rmin=None,rmax=None,nrad=20,kwmin=0,kwmax=15,projected=False,cumulative=cumulative,obs_cut=None)

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

        return out   
