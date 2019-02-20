#Make key plots given StarCluster
import matplotlib.pyplot as plt

from galpy.util import bovy_plot
import seaborn as sns

bovy_plot.bovy_print(axes_labelsize=18.,xtick_labelsize=14.,ytick_labelsize=14.)
current_palette = sns.color_palette()

def posplot(cluster,filename=None,coords='xy'):

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
    else:
        return plt
