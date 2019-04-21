"""
ANIMATE

Routines to easily generate the same figures at different timesteps

"""

import matplotlib.pyplot as plt
import numpy as np
from galpy.util import bovy_plot
import os


from plots import dvplot, bplot

def dvanimate(data,nsnap=0,tsnap=None,prefix='',nrad=20,save=True,**kwargs):
    """
    NAME:

       dvanimate

    PURPOSE:

       Plot velocity dispersion profiles using dvplot over a range of timesteps
       To Do:
        - use data from StarCluster as opposed to a written file

    INPUT:

       data - array from output.sigv_out (t,m,rm,r[0:nrad],sig[0:nrad],beta[0:nrad])

       nsnap - starting snapshot (default:0)

       tsnap - starting time step (default: 0 Myr, overwrites nsnap)

       prefix - string noting the begining of directory to save images

       nrad - number of radial bins in profile

       save - option to save images or simply just show them (default: True)

    KWARGS:

       kwargs - for passing to plots.dvplot

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    t=data[:,0]
    m=data[:,1]
    rm=data[:,-1]

    
    if not os.path.exists('./%sdvmovie/' % prefix):
        os.makedirs('./%sdvmovie/' % prefix)

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    for i in range(nsnap,len(t)+1):

        if save:
            filename='./%sdvmovie/%sdvplot%s.png' % (prefix,prefix,str(i))
        else:
            filename=None

        npy.dvplot(data,nsnap,tsnap,prefix='',nrad=nrad,nsnap=i,filename=filename,**kwargs)

    return 0

def banimate(data,nsnap=0,tsnap=None,prefix='',nrad=20,save=True,**kwargs):
    """
    NAME:

       banimate

    PURPOSE:

       Plot anisotropy profiles using bplot over a range of timesteps
       To Do:
        - use data from StarCluster as opposed to a written file

    INPUT:

       data - array from output.sigv_out (t,m,rm,r[0:nrad],sig[0:nrad],beta[0:nrad])

       nsnap - starting snapshot (default:0)

       tsnap - starting time step (default: 0 Myr, overwrites nsnap)

       prefix - string noting the begining of directory to save images

       nrad - number of radial bins in profile

       save - option to save images or simply just show them (default: True)

    KWARGS:

       kwargs - for passing to plots.bplot

    OUTPUT:

       None

    HISTORY:

       2019 - Written - Webb (UofT)
    """
    t=data[:,0]
    m=data[:,1]
    rm=data[:,-1]

    
    if not os.path.exists('./%sbmovie/' % prefix):
        os.makedirs('./%sbmovie/' % prefix)

    if tsnap!=None:
        nsnap=int(np.argwhere(t>=tsnap)[0])

    for i in range(nsnap,len(t)+1):

        if save:
            filename='./%sbmovie/%sbplot%s.png' % (prefix,prefix,str(i))
        else:
            filename=None

        npy.bplot(data,nsnap,tsnap,prefix='',nrad=nrad,nsnap=i,filename=filename,**kwargs)

    return 0