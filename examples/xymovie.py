import matplotlib
matplotlib.use('Agg')
import nbodypy as npy
import matplotlib.pyplot as plt
import os

filepath=''
f82=open('fort.82','r')
f83=open('fort.83','r')
savedir='./xymovie/'

cluster=npy.get_nbody6_jarrod(f82,f83)

nfig=0

while cluster.ntot>0:
    npy.nbody_to_realpc(cluster)
    print(cluster.tphys)    
    plt.plot(cluster.x,cluster.y,'k.')

    plt.xlabel('X (pc)')
    plt.ylabel('Y (pc)')
    plt.xlim(-100,100)
    plt.ylim(-100,100)
    plt.title('Time = %f' % cluster.tphys)
    filename=(savedir+'%s.png' % str(nfig).zfill(5))
    plt.savefig(filename)
    plt.close()
    nfig+=1

    cluster=npy.get_nbody6_jarrod(f82,f83)

    if cluster.ntot==0:
        f82.close()
        f83.close()
        filepath+='cont/'
            if os.path.isfile(filepath+sfilename):
                f82=open(filepath+bfilename,'r')
                f83=open(filepath+sfilename,'r')
                cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)





