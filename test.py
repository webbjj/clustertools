import nbodypy as npy
import math

x,y,z=[5.0,-1.0,6.0]
vx,vy,vz=[1.0,8.0,2.0]
r=math.sqrt(x*x+y*y+z*z)
v=math.sqrt(vx*vx+vy*vy+vz*vz)

print(x,y,z,r,vx,vy,vz,v)

print(npy.cart_to_sphere(x,y,z))
