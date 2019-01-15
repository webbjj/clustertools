#Routines for converting into different coordinate systems:
import numpy as np

#Scale positions and velocity by a specific value
def xvscale(cluster,rscale=1.0,vscale=1.0):

    cluster.x*=rscale
    cluster.y*=rscale
    cluster.z*=rscale
    cluster.vx*=vscale
    cluster.vy*=vscale
    cluster.vz*=vscale
    
    if cluster.keyparams:
        cluster.key_params()


#Shift cluster to/from cluster-centric and galactocentric coordinates
def xvshift(cluster,x,y,z,vx,vy,vz,origin=None):
    cluster.x+=x
    cluster.y+=y
    cluster.z+=z
    cluster.vx+=vx
    cluster.vy+=vy
    cluster.vz+=vz

    if origin!=None:
        cluster.origin=origin

    if cluster.keyparams:   
        cluster.key_params()


def rect_to_sphere(x,y,z):
    r=np.sqrt(x*x+y*y+z*z)
    rxy=np.sqrt(x*x+y*y)
    phi=np.arctan2(y,x)
    theta=np.arctan2(rxy,z)

    return r,phi,theta

def rect_to_sphere_vel(x,y,z,vx,vy,vz):
    r,phi,theta=rect_to_sphere(x,y,z)
    v=np.sqrt(vx*vx+vy*vy+vz*vz)

    #v DOT r_hat,phi_hat,and theta_hat
    
    vr=vx*np.cos(phi)*np.sin(theta)+vy*np.sin(phi)*np.sin(theta)+vz*np.cos(theta)
    vp=-1.0*vx*np.sin(phi)+vy*np.cos(phi)
    vt=vx*np.cos(phi)*np.cos(theta)+vy*np.sin(phi)*np.cos(theta)-vz*np.sin(theta)

    return vr,vp,vt



