""" The PlanetarySystem class and key internal functions

"""

__author__ = "Jeremy J Webb"

__all__ = [
    "PlanetarySystem",
]

import numpy as np
import h5py
from ..cluster.cluster import StarCluster

class StarClusterwPlanets(StarCluster):
    """ A class that represents a star cluster population with planets that ooperations and functions can be performed on
    
    Parameters
    ----------
    tphys : float
        Time (units not necessary) associated with the population (default: 0)
    units : str
        Units of stellar positions and velocties. Options include 'pckms',
        'kpckms','radec','nbody',and 'galpy'. For 'pckms' and 'kpckms', 
        stellar velocities are assumed to be km/s. (default: None)
    origin : str
        Origin of coordinate system within which stellar positions and velocities are defined. 
        Options include 'centre', 'cluster', 'galaxy' and 'sky'. Note that 'centre' corresponds
        to the systems centre of density or mass (as calculated by clustertools) while 'cluster' corresponds
        to the orbital position of the cluster's initial centre of mass given 'tphys'. In some cases
        these two coordinate systems may be the same. (default: None)
    ctype : str
        Code type used to generate the cluster population. Current options include 'snapshot',
        'nbody6','nbody6se','nemo','gyrfalcon','snaptrim', amuse','clustertools','snapauto'. The parameter
        informs clustertools how to load the stellar popuation and advance to the next snapshot.
        (default: 'snapshot')
    projected : bool
        return projected values instead of 3D value. (default: False)

    Other Parameters
    ----------------
    sfile : str
        file containing single star data
    bfile : str
        file contain binary star data
    ssefile : str
        file containing single star evolution data
    bsefile : str
        file contain binary star evolution data
    tfile : str
        name of file contain tail star data
    ofile : str
        orbit file
    ofilename : str
        orbit filename if ofile is not given
    orbit : str
        galpy orbit instance
    ounits : str
        {'pckms','kpckms','radec','nbody','galpy'} units of orbital information (else assumed equal to StarCluster.units)
    nsnap : int
        if a specific snapshot is to be read in instead of starting from zero
    nzfill : int
        value for zfill when reading and writing snapshots (default: 5)
    snapbase : str
        string of characters in filename before nsnap (default: '')
    snapend : str
        string of character in filename after nsnap (default: '.dat')
    snapdir : str
        string name of directory of snapshots if different than wdir (Default: '')
    delimiter : str
        choice of delimiter when reading ascii/csv files (default: ',')
    wdir : str
        working directory of snapshots if not current directory
    intialize : bool
        initialize a galpy orbit after reading in orbital information (default: False)
    advance : bool
        set to True if this a snapshot that has been advanced to from an initial one? (default: False)
    centre_method: str
        {None,'orthographic','VandeVen'} method to convert to clustercentric coordinates when units are in radec (Default: None)
    give : str
        set what parameters are read in from nemo/gyrfalcon (default: 'mxv')
        Currently only accepts 'mxvpqael' as an alternative.

    History
    -------
    2021 - Written - Webb (UofT)

    """

    def __init__(
        self,
        tphys=0.0,
        units=None,
        origin=None,
        ctype="snapshot",
        projected=False,
        **kwargs,
    ):

        super().__init__(tphys,units,origin,ctype,projected,**kwargs)

        self.planets=True
        self.npsys=0
        self.pid=np.array([])
        self.p_sys_id=np.array([])
        self.psys=np.array([])
        self.parg=np.array([])

    def add_planets(self,pctype='LonelyPlanets',npsys=1,p_sys_id=None,pwdir=None,psnapdir='',p_sys_ic='planetary_systems.h5ic',psnapbase='p_sys_',p_sys_id_sub='',psnapend='.hdf5',**kwargs):

        self.pctype=pctype
        if pwdir is None:
            pwdir=self.wdir

        if self.pctype=='LonelyPlanets':
            pid,psysid,psys=_get_lonelyplanets(npsys,p_sys_id,pwdir,psnapdir,p_sys_ic,psnapbase,p_sys_id_sub,psnapend,**kwargs)
        
        self.pid=np.append(self.pid,pid)
        self.p_sys_id=np.append(self.p_sys_id,psysid)
        self.psys=np.append(self.psys,psys)


        for i in range(0,len(psys)):
            self.parg=np.append(self.parg,np.argmin(np.fabs(self.id-self.pid[i])))

        self.parg=self.parg.astype(int)

        self.npsys=len(self.psys)

    def update_planets(self):

        self.save_cluster()
        self.to_pckms()

        for i in range(0,len(self.psys)):
            print(self.psys[i])
            self.psys[i].update_planets(self.tphys*1000000.)

        self.return_cluster()


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

def sub_clusterwplanets(
    cluster,
    rmin=None,
    rmax=None,
    mmin=None,
    mmax=None,
    vmin=None,
    vmax=None,
    emin=None,
    emax=None,
    kwmin=0,
    kwmax=15,
    indx=[None],
    projected=False,
    sortstars=True,
    reset_centre=False,
    reset_nbody=False,
    reset_nbody_mass=False,
    reset_nbody_radii=False,
    reset_rvirial=False,
    reset_projected=False,
    **kwargs
):
    """Extract a sub population of stars from a StarCluster

    - automatically moves cluster to centre of mass, so all constraints are in clustercentric coordinates and current StarCluster.units

    Parameters
    ----------
    rmin/rmax : float
        minimum and maximum stellar radii
    mmin/mmax : float
        minimum and maximum stellar mass
    vmin/vmax : float
        minimum and maximum stellar velocity
    emin/emax : float
        minimum and maximum stellar energy
    kwmin/kwmax : int
        minimum and maximum stellar type (kw)
    indx : bool
        user defined boolean array from which to extract the subset
    projected : bool 
        use projected values and constraints (default:False)
    sortstars: bool
        order stars by radius (default: True)
    reset_centre : bool
        re-calculate cluster centre after extraction (default:False)
    reset_nbody : bool
        reset nbody scaling factors (default:False)
    reset_nbody_mass : bool 
        find new nbody scaling mass (default:False)
    reset_nbody_radii : bool
        find new nbody scaling radius (default:False)
    reset_rvirial : bool
        use virial radius to find nbody scaling radius (default: True)
    reset_projected : bool 
        use projected radii to find nbody scaling radius (default: False)

    Returns
    -------
    StarCluster

    History
    -------
    2018 - Written - Webb (UofT)

    """
    cluster.save_cluster()
    units0,origin0, rorder0, rorder_origin0 = cluster.units0,cluster.origin0, cluster.rorder0, cluster.rorder_origin0

    if projected:
        r = cluster.rpro
        v = cluster.vpro
    else:
        r = cluster.r
        v = cluster.v

    if rmin == None:
        rmin = np.amin(r)
    if rmax == None:
        rmax = np.amax(r)
    if vmin == None:
        vmin = np.amin(v)
    if vmax == None:
        vmax = np.amax(v)
    if mmin == None:
        mmin = np.amin(cluster.m)
    if mmax == None:
        mmax = np.amax(cluster.m)

    if emin == None and emax != None:
        eindx = cluster.etot <= emax
    elif emin != None and emax == None:
        eindx = cluster.etot >= emin
    elif emin != None and emax != None:
        eindx = (cluster.etot <= emax) * (cluster.etot >= emin)
    else:
        eindx = cluster.id > -1

    if None in indx:
        indx = cluster.id > -1

    indx *= (
        (r >= rmin)
        * (r <= rmax)
        * (cluster.m >= mmin)
        * (cluster.m <= mmax)
        * (v >= vmin)
        * (v <= vmax)
        * eindx
    )

    if len(cluster.kw) > 0:
        indx *= (cluster.kw >= kwmin) * (cluster.kw <= kwmax)

    if np.sum(indx) > 0:

        if cluster.planets:
            subcluster = StarClusterwPlanets(
                cluster.tphys,
                units=cluster.units,
                origin=cluster.origin,
                ctype=cluster.ctype,
            )
        else:
            subcluster = StarCluster(
                cluster.tphys,
                units=cluster.units,
                origin=cluster.origin,
                ctype=cluster.ctype,
            )

        subcluster.add_stars(
            cluster.x[indx],
            cluster.y[indx],
            cluster.z[indx],
            cluster.vx[indx],
            cluster.vy[indx],
            cluster.vz[indx],
            cluster.m[indx],
            cluster.id[indx],
            sortstars=sortstars,
        )

        if len(cluster.ra)==len(cluster.x):
            subcluster.ra, subcluster.dec, subcluster.dist = (
                cluster.ra[indx],
                cluster.dec[indx],
                cluster.dist[indx],
            )
            subcluster.pmra, subcluster.pmdec, subcluster.vlos = (
                cluster.pmra[indx],
                cluster.pmdec[indx],
                cluster.vlos[indx],
            )

        subcluster.add_nbody6(cluster.nc,cluster.rc,cluster.rbar,
            cluster.rtide,cluster.xc,cluster.yc,cluster.zc,
            cluster.zmbar,cluster.vbar,cluster.tbar,cluster.rscale,
            cluster.ns,cluster.nb,cluster.n_p)

        subcluster.projected = cluster.projected
        subcluster.centre_method = cluster.centre_method

        if len(cluster.logl) > 0:
            if len(cluster.ep) !=0 and len(cluster.ospin) != 0:
                subcluster.add_sse(
                    cluster.kw[indx],
                    cluster.logl[indx],
                    cluster.logr[indx],
                    cluster.ep[indx],
                    cluster.ospin[indx],
                )
            else:
                subcluster.add_sse(
                        cluster.kw[indx],
                        cluster.logl[indx],
                        cluster.logr[indx],
                    )

        else:
            subcluster.kw = cluster.kw[indx]


        if len(cluster.id2) > 0:
            bindx1 = np.in1d(cluster.id1, cluster.id[indx])
            bindx2 = np.in1d(cluster.id2, cluster.id[indx])
            bindx=np.logical_or(bindx1,bindx2)


            if len(cluster.ep1) !=0 and len(cluster.ospin1) != 0:

                subcluster.add_bse(
                    cluster.id1[bindx],
                    cluster.id2[bindx],
                    cluster.kw1[bindx],
                    cluster.kw2[bindx],
                    cluster.kcm[bindx],
                    cluster.ecc[bindx],
                    cluster.pb[bindx],
                    cluster.semi[bindx],
                    cluster.m1[bindx],
                    cluster.m2[bindx],
                    cluster.logl1[bindx],
                    cluster.logl2[bindx],
                    cluster.logr1[bindx],
                    cluster.logr2[bindx],
                    cluster.ep1[bindx],
                    cluster.ep2[bindx],
                    cluster.ospin1[bindx],
                    cluster.ospin2[bindx],
                )
            else:
                subcluster.add_bse(
                    cluster.id1[bindx],
                    cluster.id2[bindx],
                    cluster.kw1[bindx],
                    cluster.kw2[bindx],
                    cluster.kcm[bindx],
                    cluster.ecc[bindx],
                    cluster.pb[bindx],
                    cluster.semi[bindx],
                    cluster.m1[bindx],
                    cluster.m2[bindx],
                    cluster.logl1[bindx],
                    cluster.logl2[bindx],
                    cluster.logr1[bindx],
                    cluster.logr2[bindx],
                )

        if len(cluster.etot) > 0:
            subcluster.add_energies(
                cluster.kin[indx], cluster.pot[indx],
            )

        if cluster.give == 'mxvpqael':
            subcluster.give=cluster.give
            subcluster.gyrpot=cluster.gyrpot[indx]
            subcluster.gyrq=cluster.gyrq[indx]
            subcluster.gyracc=cluster.gyracc[indx]
            subcluster.eps=cluster.eps[indx]
            subcluster.gyrlev=cluster.gyrlev[indx]
        elif cluster.give =='mxve':
            subcluster.give=cluster.give
            subcluster.eps=cluster.eps[indx]


        if reset_centre:
            subcluster.add_orbit(
                cluster.xgc,
                cluster.ygc,
                cluster.zgc,
                cluster.vxgc,
                cluster.vygc,
                cluster.vzgc,
            )

            if cluster.origin=='centre' or cluster.origin=='cluster':
                subcluster.find_centre(0.0, 0.0, 0.0, reset_centre=reset_centre)
            else:
                subcluster.find_centre(reset_centre=reset_centre)

        else:
            subcluster.add_orbit(
                cluster.xgc,
                cluster.ygc,
                cluster.zgc,
                cluster.vxgc,
                cluster.vygc,
                cluster.vzgc,
            )
            subcluster.xc, subcluster.yc, subcluster.zc = (
                cluster.xc,
                cluster.yc,
                cluster.zc,
            )                                                                                                                          
            subcluster.vxc, subcluster.vyc, subcluster.vzc = (
                cluster.vxc,
                cluster.vyc,
                cluster.vzc,
            )

            subcluster.ra_gc, subcluster.dec_gc, subcluster.dist_gc - cluster.ra_gc, cluster.dec_gc, cluster.dist_gc
            subcluster.pmra_gc, subcluster.pmdec_gc, subcluster.vlos_gc = (
                cluster.pmra_gc,
                cluster.pmdec_gc,
                cluster.vlos_gc,
            )

        if reset_nbody:
            subcluster.to_pckms()
            subcluster.analyze()
            subcluster.reset_nbody_scale(mass=True,radius=True,rvirial=reset_rvirial,projected=reset_projected,**kwargs)
        elif reset_nbody_mass or reset_nbody_radii:
            subcluster.to_pckms()
            subcluster.analyze()
            subcluster.reset_nbody_scale(mass=reset_nbody_mass,radius=reset_nbody_radii,rvirial=reset_rvirial,projected=reset_projected,**kwargs)

        #Copy planet over to subcluster
        if cluster.planets:
            for i in range(0,cluster.npsys):
                if cluster.pid[i] in subcluster.id:
                    subcluster.npsys+=1
                    subcluster.psys=np.append(subcluster.psys,cluster.psys[i])

                    subcluster.parg=np.append(subcluster.parg,np.argmin(np.fabs(subcluster.id-cluster.pid[i])))
                    subcluster.p_sys_id=np.append(subcluster.p_sys_id,cluster.p_sys_id[i])

            subcluster.parg=subcluster.parg.astype(int)

    else:
        subcluster = StarCluster(cluster.tphys)

    if subcluster.ntot > 0:
        subcluster.analyze(sortstars=sortstars)

    cluster.return_cluster(units0,origin0, rorder0, rorder_origin0)
    subcluster.return_cluster(units0,origin0, rorder0, rorder_origin0)

    return subcluster
