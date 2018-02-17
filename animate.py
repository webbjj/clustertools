#Create a time sequence of plots in order to animate the evolution
#of various parameters
import matplotlib.pyplot as plt
import os
from get_nbody import *

def animate_nbody6_jarrod(fort82,fort83,xparam,yparam,xlim=None,ylim=None,savedir='./',units=None):

    if not os.path.isdir(savedir):
        print('SAVE DIRECTORY DOES NOT EXIST')
        return -1

    filenum=0
    cluster=get_nbody6_jarrod(fort82,fort83)

    while cluster.ntot>0:
        filename=str(filenum).zfill(5)+'.png'        
        makeplot(cluster,xparam,yparam,xlim,ylim,savedir,filename,units)
        cluster=get_nbody6_jarrod(fort82,fort83)
        filenum+=1

def animate_nbody6_out34(out34,xparam,yparam,xlim=None,ylim=None,savedir='./',units=None):

    if not os.path.isdir(savedir):
        print('SAVE DIRECTORY DOES NOT EXIST')
        return -1

    filenum=0
    cluster=get_nbody6_out34(out34)

    while cluster.ntot>0:
        filename=str(filenum).zfill(5)+'.png'
        makeplot(cluster,xparam,yparam,xlim,ylim,savedir,filename,units)
        cluster=get_nbody6_out34(out34)
        filenum+=1

def makeplot(cluster,xparam,yparam,xlim=None,ylim=None,savedir='./',filename=None,units=None):

    #Set Units
    if units!=cluster.units:
        if units=='realpc' and cluster.units=='nbody':
            nbody_to_realpc(cluster)
        elif units=='realkpc' and cluster.units=='nbody':
            nbody_to_realkpc(cluster)
        elif units=='galpy' and cluster.units=='nbody':
            nbody_to_galpy(cluster)
            

    #Find which parameters to plot (this list is coninually being added to as I need to make plots)

    if xparam=='x':
        plotx=cluster.x
        plt.xlabel('X')
    elif xparam=='y':
        plotx=cluster.y
        plt.xlabel('Y')
    elif xparam=='z':
        plotx=cluster.z
        plt.xlabel('Z')

    if yparam=='x':
        ploty=cluster.x
        plt.ylabel('X')
    elif yparam=='y':
        ploty=cluster.y
        plt.ylabel('Y')
    elif yparam=='z':
        ploty=cluster.z
        plt.ylabel('Z')

    if filename==None:
        filename=xparam + yparam + '_plot.png'

    if xlim!=None:
        plt.xlim(xlim)
    if ylim!=None:
        plt.ylim(ylim)

    plt.plot(plotx,ploty,'k.')
    plt.title('Time = %f' % cluster.tphys)
    plt.close() 
