#Make key plots given StarCluster
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from galpy.util import bovy_plot
import seaborn as sns

bovy_plot.bovy_print(axes_labelsize=18.,xtick_labelsize=14.,ytick_labelsize=14.)
current_palette = sns.color_palette()

def posplot(cluster,coords='xy',filename=None):

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
    plt.plot(x,y,'.',alpha=0.5)

    if filename!=None:
        plt.savefig(filename)
        plt.close()   
    else:
        return plt

def double_posplot(cluster,xlim=(-100,100),ylim=(-100,100),filename=None):
    
    plt.subplot(1,2,1)
    
    plt.plot(cluster.x, cluster.z,'k.',alpha=0.1)
    
    if cluster.units=='galaxy':
        plt.plot(cluster.xgc,cluster.zgc,'ro')
        plt.plot(cluster.xgc+cluster.xc,cluster.zgc+cluster.zc,'bo')
    else:
        plt.plot(cluster.xc,cluster.zc,'bo')

    
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel('X (kpc)')
    plt.ylabel('Z (kpc)')

    plt.subplot(1,2,2)

    plt.plot(cluster.x, cluster.y,'k.',alpha=0.1)

    if cluster.units=='galaxy':
        plt.plot(cluster.xgc,cluster.ygc,'ro')
        plt.plot(cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,'bo')
    else:
        plt.plot(cluster.xc,cluster.yc,'bo')

    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.xlabel('X (kpc)')
    plt.ylabel('Z (kpc)')

    if filename!=None:
        plt.savefig(filename)
    plt.close()
