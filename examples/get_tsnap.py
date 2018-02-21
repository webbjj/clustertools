#Write out id,m,x,y,z,vx,vy,vz for all stars at a given time
import nbodypy as npy
import sys,os

tsnap=float(sys.argv[1])

filepath=''
bfilename='fort.82'
sfilename='fort.83'

if not os.path.isfile(sfilename):
    filepath+='./start/'

f82=open(filepath+bfilename,'r')
f83=open(filepath+sfilename,'r')

cluster=npy.get_nbody6_jarrod(f82,f83,do_keyparams=False)

while cluster.tphys<tsnap:
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
  

if cluster.tphys>=tsnap:
    npy.nbody_to_realpc(cluster)
    npy.snapout(cluster,'tsnap_%s.dat' % str(tsnap))
else:
    print('Time stamp out of range')
