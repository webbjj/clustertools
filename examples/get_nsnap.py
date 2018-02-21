#Write out id,m,x,y,z,vx,vy,vz for all stars at a given time
import nbodypy as npy
import sys,os

nsnap=float(sys.argv[1])

filepath=''
bfilename='fort.82'
sfilename='fort.83'

if not os.path.isfile(sfilename):
    filepath+='./start/'

f82=open(filepath+bfilename,'r')
f83=open(filepath+sfilename,'r')

cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)


for i in range(0,nsnap):
    print(cluster.tphys,cluster.ntot)
    cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)

    if cluster.ntot==0:
        f82.close()
        f83.close()
        filepath+='cont/'
        if os.path.isfile(filepath+sfilename):
            f82=open(filepath+bfilename,'r')
            f83=open(filepath+sfilename,'r')
            cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)
  
npy.nbody_to_realpc(cluster)
npy.snapout(cluster,'nsnap_%s.dat' % str(nsnap))
