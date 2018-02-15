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

def makeplot(object,xparam,yparam,xlim=None,ylim=None,savedir='./',filename=None,units=None):

    #Set Units
    if units!=object.units:
        if units=='realpc' and object.units=='nbody':
            nbody_to_realpc(object)
        elif units=='realkpc' and object.units=='nbody':
            nbody_to_realkpc(object)
        elif units=='galpy' and object.units=='nbody':
            nbody_to_galpy(object)
            

    #Find which parameters to plot (this list is coninually being added to)
    if xparam=='x':
        plotx=object.x
        plt.xlabel('X')
    elif xparam=='y':
        plotx=object.y
        plt.xlabel('Y')
    elif xparam=='z':
        plotx=object.z
        plt.xlabel('Z')

    if yparam=='x':
        ploty=object.x
        plt.ylabel('X')
    elif yparam=='y':
        ploty=object.y
        plt.ylabel('Y')
    elif yparam=='z':
        ploty=object.z
        plt.ylabel('Z')

    if filename==None:
        filename=xparam + yparam + '_plot.png'

    if xlim!=None:
        plt.xlim(xlim)
    if ylim!=None:
        plt.ylim(ylim)

    plt.plot(plotx,ploty,'k.')
    plt.title(xparam+' VS '+yparam)
    plt.savefig(savedir+filename)   
    plt.close() 
