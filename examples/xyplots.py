import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import clustertools as cts
import sys

def simplot(ctype,**kwargs):

    nsnap=0

    print('CTYPE = ',ctype)
    print(kwargs)
    wdir=kwargs.get('wdir','./')
    units=kwargs.pop('units',None)
    origin=kwargs.pop('origin',None)
    filename=kwargs.pop('filename',None)


    dtout=int(kwargs.pop('dtout',1))

    coords=kwargs.pop('coords','xyz')
    xlower=kwargs.pop('xlower',None)
    xupper=kwargs.pop('xupper',None)
    if xlower is None or xupper is None:
        xlim=None
    else:
        xlim=(float(xlower),float(xupper))

    ylower=kwargs.pop('ylower',None)
    yupper=kwargs.pop('yupper',None)
    if ylower is None or yupper is None:
        ylim=None
    else:
        ylim=(float(ylower),float(yupper))

    plotunits=kwargs.pop('plotunits',units)
    plotorigin=kwargs.pop('plotorigin',origin)


    #See if orbit file is provided
    if 'ofilename' in kwargs:
        ofile=open(kwargs.pop('ofilename'),'r')
    else:
        ofile=None

    cluster=cts.load_cluster(ctype,wdir=wdir,units=units,origin=origin,ofile=ofile,filename=filename,**kwargs)


    while cluster.ntot > 0:
        filename='%s.png' % str(nsnap).zfill(5)

        print(nsnap,filename,cluster.ntot)


        cluster.to_units(plotunits)
        cluster.to_origin(plotorigin)
        cts.starplot(cluster,coords=coords,xlim=xlim,ylim=ylim,filename=filename)
        plt.close()
        cluster=cts.advance_cluster(cluster,dtout=dtout,**kwargs)
        nsnap+=1


if __name__ in '__main__':

    ctype=sys.argv[1]

    if len(sys.argv) < 2:
        raise SyntaxError("Insufficient arguments")
    elif len(sys.argv)==2:
        simplot(ctype)
    else:
        # If there are keyword arguments
        simplot(ctype,**dict(arg.split('=') for arg in sys.argv[2:]))
