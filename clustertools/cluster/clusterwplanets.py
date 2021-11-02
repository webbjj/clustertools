""" The PlanetarySystem class and key internal functions

"""

__author__ = "Jeremy J Webb"

__all__ = [
    "PlanetarySystem",
    "_get_lonelyplanets",
]

import numpy as np
import h5py

def _get_lonelyplanets(npsys=1,p_sys_id=None,pwdir=None,psnapdir='',p_sys_ic='planetary_systems.h5ic',psnapbase='p_sys_',p_sys_id_sub='',psnapend='.hdf5',**kwargs):

    psys=np.array([])
    pid=np.array([])

    if p_sys_id is None:
        if npsys==1:
            p_sys_id=np.array([0])
        else:
            p_sys_id=np.linspace(0,npsys-1,npsys)

    p_sys_id=p_sys_id.astype(int)

    icfilename='%s%s%s' % (pwdir,psnapdir,p_sys_ic)
    icdata=h5py.File(icfilename,'r')


    for i in range(0,len(p_sys_id)):
        pid=np.append(pid,icdata['p_sys_%s' % str(p_sys_id[i])].attrs['host_star_id'])

        filename='%s%s%s%s%s%s' % (pwdir,psnapdir,psnapbase,str(p_sys_id[i]),p_sys_id_sub,psnapend)
        data=h5py.File(filename,'r')

        time=data['time'][:]
        x=np.swapaxes(data['x'],0,1)[1:]
        y=np.swapaxes(data['y'],0,1)[1:]
        z=np.swapaxes(data['z'],0,1)[1:]
        vx=np.swapaxes(data['vx'],0,1)[1:]
        vy=np.swapaxes(data['vy'],0,1)[1:]
        vz=np.swapaxes(data['vz'],0,1)[1:]
        ax=np.swapaxes(data['ax'],0,1)[1:]
        ay=np.swapaxes(data['ay'],0,1)[1:]
        az=np.swapaxes(data['az'],0,1)[1:]
        ecc=np.swapaxes(data['ecc'],0,1)
        inc=np.swapaxes(data['inc'],0,1)
        a=np.swapaxes(data['a'],0,1)
        energy=data['energy'][:]

        m_pert=np.swapaxes(data['m_pert'],0,1)
        pert_dist=np.swapaxes(data['pert_dist'],0,1)


        data.close()

        psys=np.append(psys,PlanetarySystem(p_sys_id[i],pwdir,psnapdir,p_sys_ic,psnapbase,p_sys_id_sub,psnapend))
        psys[-1].add_planet_histories(time,x,y,z,vx,vy,vz,ax,ay,az,ecc,inc,a,energy)
        psys[-1].add_perturber_histories(m_pert,pert_dist)

    icdata.close()

    return pid,p_sys_id,psys


class PlanetarySystem(object):

    def __init__(self,p_sys_id=None,pwdir=None,psnapdir='',p_sys_ic='planetary_systems.h5ic',psnapbase='p_sys_',p_sys_id_sub='',psnapend='.hdf5'):

        self.p_sys_id=p_sys_id
        self.pwdir=pwdir
        self.psnapdir=psnapdir
        self.psnapbase=psnapbase
        self.psnapend=psnapend

        
    def add_planet_histories(self,time,x=None,y=None,z=None,vx=None,vy=None,vz=None,ax=None,ay=None,az=None,ecc=None,inc=None,a=None,energy=None):

        self.time=time
        self.xt=x
        self.yt=y
        self.zt=z
        self.vxt=vx
        self.vyt=vy
        self.vzt=vz
        self.axt=ax
        self.ayt=ay
        self.azt=az
        self.ecct=ecc
        self.inct=inc
        self.at=a
        self.energyt=energy

        self.update_planets(0.)

    def update_planets(self,tphys):
        tindx=self.time>=(tphys)

        if np.sum(tindx)>0:
            if self.xt is not None: self.x=self.xt[:,tindx][:,0]
            if self.yt is not None: self.y=self.yt[:,tindx][:,0]
            if self.zt is not None: self.z=self.zt[:,tindx][:,0]
            if self.vxt is not None: self.vx=self.vxt[:,tindx][:,0]
            if self.vyt is not None: self.vy=self.vyt[:,tindx][:,0]
            if self.vzt is not None: self.vz=self.vzt[:,tindx][:,0]
            if self.axt is not None:self.ax=self.axt[:,tindx][:,0]
            if self.ayt is not None:self.ay=self.ayt[:,tindx][:,0]
            if self.azt is not None: self.az=self.azt[:,tindx][:,0]
            if self.ecct is not None: self.ecc=self.ecct[:,tindx][:,0]
            if self.inct is not None: self.inc=self.inct[:,tindx][:,0]
            if self.at is not None: self.a=self.at[:,tindx][:,0]
            if self.energyt is not None: self.energy=self.energyt[tindx][0]

        else:
            print('PLANET DATA DOES NOT EXIST AT THIS TIMESTEP')

    def add_perturber_histories(self,m_pert,pert_dist):
        self.m_pert=m_pert
        self.pert_dist=pert_dist