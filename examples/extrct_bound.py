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

    #Extract bound stars only
    exfile=open('extrct_bound.dat','w')

    #Look for cluster centre using clustertools
    findcentre=kwargs.pop('findcentre',False)
    if findcentre=='True':
        findcentre=True

    #See if orbit file is provided
    if 'ofilename' in kwargs:
        ofile=open(kwargs.pop('ofilename'),'r')
    else:
        ofile=None

    cluster=cts.load_cluster(ctype,wdir=wdir,units=units,origin=origin,ofile=ofile,filename=filename,**kwargs)

    #Assume that initially all stars are bound
    bndindx=np.ones(cluster.ntot,dtype=bool)
    bndids=cluster.id[bndindx]

    while cluster.ntot > 0:
        print(nsnap,cluster.ntot)

        if findcentre: cluster.find_centre(indx=bndindx)

        cluster.to_centre()
        cluster.energies()
        bndindx=(cluster.etot<0)
        bndids=cluster.id[bndindx]

        cluster.to_units(exunits)
        cluster.to_origin(exorigin)

        subcluster=cts.sub_cluster(cluster,indx=bndindx)
        subcluster.half_mass_relaxation_time()
        subcluster.rlagrange()

        exfile.write('%f %f %f %f %f' % (subcluster.ntot,subcluster.nb,subcluster.tphys,subcluster.trh,subcluster.mtot))

        for rn in subcluster.rn:
            exfile.write(' %f' % rn)

        exfile.write('\n')

        cluster=cts.advance_cluster(cluster,dtout=dtout,**kwargs)
        bndindx=np.in1d(cluster.id,bndids)
        
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
