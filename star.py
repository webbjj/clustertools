#Define the Star class and add properties to the star
import numpy
import math

class Star(object):
    def __init__(self,id,m,x,y,z,vx,vy,vz,tphys='0.0',units='nbody',origin='cluster'):
        self.ntot=1
        self.id=id
        self.m=m
        self.x=x
        self.y=y
        self.z=z
        self.vx=vx
        self.vy=vy
        self.vz=vz
        self.tphys=tphys
        self.units=units
        self.origin=origin
    
        #Calculate Key Parameters
        self.v=math.sqrt(self.vx**2.0+self.vy**2.0+self.vz**2.0)
        self.rxy=math.sqrt(self.x**2.0+self.y**2.0)
        self.r=math.sqrt(self.x**2.0+self.y**2.0+self.z**2.0)

    def add_se(self,kw,logl,logr):
        self.kw=kw
        self.logl=logl
        self.logr=logr

    def add_energies(self,kin,pot,etot):
        self.kin=kin
        self.pot=pot
        self.etot=etot

    def add_orbit(self,xgc,ygc,zgc,vxgc,vygc,vzgc):
        self.xgc=xgc
        self.ygc=ygc
        self.zgc=zgc
        self.rgc=math.sqrt(xgc**2.0+ygc**2.0+zgc**2.0)
        self.vxgc=vxgc
        self.vygc=vygc
        self.vzgc=vzgc
