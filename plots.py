#Make key plots given StarCluster
import matplotlib.pyplot as plt

def xyplot(cluster,filename=None):

    x=cluster.x
    y=cluster.y

    if cluster.units=='nbody':
        units='(NBODY)'
    elif cluster.units=='realpc':
        units='(pc)'
    elif cluster.units=='realkpc':
        units='(kpc)'
    else:
        units=''

    plt.xlabel('X '+units)
    plt.ylabel('Y '+units)
    plt.title('Time = %f' % cluster.tphys)
    plt.plot(x,y,'k.')

    if filename!=None:
        plt.savefig(filename)   
    else:
        return plt

def xzplot(cluster,filename=None):

    x=cluster.x
    y=cluster.z

    if cluster.units=='nbody':
        units='(NBODY)'
    elif cluster.units=='realpc':
        units='(pc)'
    elif cluster.units=='realkpc':
        units='(kpc)'
    else:
        units=''

    plt.xlabel('X '+units)
    plt.ylabel('Z '+units)
    plt.title('Time = %f' % cluster.tphys)
    plt.plot(x,y,'k.')

    if filename!=None:
        plt.savefig(filename)
    else:
        return plt
 

def yzplot(cluster,filename=None):

    x=cluster.y
    y=cluster.z

    if cluster.units=='nbody':
        units='(NBODY)'
    elif cluster.units=='realpc':
        units='(pc)'
    elif cluster.units=='realkpc':
        units='(kpc)'
    else:
        units=''

    plt.xlabel('Y '+units)
    plt.ylabel('Z '+units)
    plt.title('Time = %f' % cluster.tphys)
    plt.plot(x,y,'k.')

    if filename!=None:
        plt.savefig(filename)
    else:
        return plt
