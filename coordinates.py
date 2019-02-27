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


    cluster.key_params()



