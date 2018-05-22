#Plot x and y position of all stars, grouped into subsets

import matplotlib
matplotlib.use('Agg')
import nbodypy as npy
import matplotlib.pyplot as plt
import os

filepath=''
f82=open('fort.82','r')
f83=open('fort.83','r')
savedir='./xymovie/'
colours=['k.','r.','b.','g.','m.']

cluster=npy.get_nbody6_jarrod(f82,f83)
r0=np.zeros(cluster.ntot)
for i in range(0,cluster.ntot):
    r0[cluster.id[i]]=cluster.r[i]*cluster.rbar

rlimit=[]
rlimit.append(cluster.rnfind(0.2)*cluster.rbar)
rlimit.append(cluster.rnfind(0.4)*cluster.rbar)
rlimit.append(cluster.rnfind(0.6)*cluster.rbar)
rlimit.append(cluster.rnfind(0.8)*cluster.rbar)
rlimit.append(cluster.rnfind(1.0)*cluster.rbar)

nfig=0

while cluster.ntot>0:
    npy.nbody_to_realpc(cluster)
    print(cluster.tphys)
    
    for i in range(0,len(rlimit)):
        x=[]
        y=[]
        if i==0:
            rlower=0.0
            rupper=rlimit[0]
        else:
            rlower=rlimit[i-1]
            rupper=rlimit[i]
        
        for j in range(0,cluster.ntot):
             if r0[cluster.id[j]]>=rlower and r0[cluster.id[j]]<=rupper:
                 x.append(cluster.x[j])
                 y.append(cluster.y[j])

        plt.plot(x,y,colours[i])

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





