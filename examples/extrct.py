import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import clustertools as cts
import sys
import numpy as np

def simplot(ctype,**kwargs):

    nsnap=0

    print('CTYPE = ',ctype)
    print(kwargs)
    wdir=kwargs.get('wdir','./')
    units=kwargs.pop('units',None)
    origin=kwargs.pop('origin',None)
    filename=kwargs.pop('filename',None)

    dtout=int(kwargs.pop('dtout',1))

    exunits=kwargs.pop('exunits',units)
    exorigin=kwargs.pop('exorigin',origin)

    exfile=open('extrct.dat','w')

    #See if orbit file is provided
    if 'ofilename' in kwargs:
        ofile=open(kwargs.pop('ofilename'),'r')
    else:
        ofile=None

    cluster=cts.load_cluster(ctype,wdir=wdir,units=units,origin=origin,ofile=ofile,filename=filename,**kwargs)


    while cluster.ntot > 0:
        print(nsnap,cluster.ntot)

        cluster.to_units(exunits)
        cluster.to_origin(exorigin)


        cluster.half_mass_relaxation_time()
        cluster.rlagrange()


        exfile.write('%f %f %f %f %f' % (cluster.ntot,cluster.nb,cluster.tphys,cluster.trh,cluster.mtot))

        for rn in cluster.rn:
            exfile.write(' %f' % rn)

        exfile.write('\n')

        cluster=cts.advance_cluster(cluster,dtout=dtout,**kwargs)
        
        nsnap+=1

    exfile.close()

if __name__ in '__main__':

    ctype=sys.argv[1]

    if len(sys.argv) < 2:
        raise SyntaxError("Insufficient arguments")
    elif len(sys.argv)==2:
        simplot(ctype)
    else:
        # If there are keyword arguments
        simplot(ctype,**dict(arg.split('=') for arg in sys.argv[2:]))
