import numpy as np
from cluster import StarCluster


#Routine for joining cluster timesteps
def join_clusters(clusters):

    base_cluster=clusters[0]

    id=np.asarray([])
    m=np.asarray([])
    x=np.asarray([])
    y=np.asarray([])
    z=np.asarray([])
    vx=np.asarray([])
    vy=np.asarray([])
    vz=np.asarray([])
    kw=np.asarray([])
    
    for i in range(1,len(clusters)):
        id=np.append(id,clusters[i].id)
        m=np.append(m,clusters[i].m)
        x=np.append(x,clusters[i].x)
        y=np.append(y,clusters[i].y)
        z=np.append(z,clusters[i].z)
        x=np.append(vx,clusters[i].vx)
        y=np.append(vy,clusters[i].vy)
        z=np.append(vz,clusters[i].vz)
        kw=np.append(kw,clusters[i].kw)

    base_cluster.add_stars(id,m,x,y,z,vx,vy,vz,kw)

    base_cluster.key_params()


    return base_cluster

# Extract a sub population of stars from a cluster
# Extraction criteria include radius, mass and stellar evolution (KW) type (default all stars)
def sub_cluster(cluster,rmin=None,rmax=None,mmin=None,mmax=None,se_min=0,se_max=100,e_min=None,e_max=None,projected=False,new_center=False):

    origin0=cluster.origin
    center0=cluster.center
    if origin0!='cluster':
        cluster.to_cluster()
    if not center0:
        cluster.to_center()

    if rmin==None and projected:
        rmin=np.min(cluster.rpro)
    elif rmin==None and not projected:
        rmin=np.min(cluster.r)
       
    if rmax==None and projected:
        rmax=np.max(cluster.rpro)
    elif rmax==None and not projected:
        rmax=np.max(cluster.r)

    if mmin==None:
        mmin=np.min(cluster.m)
    if mmax==None:
        mmax=np.max(cluster.m)

    if e_min==None and e_max!=None:
        eindx=cluster.etot<=e_max
    elif e_min!=None and e_max==None:
        eindx=cluster.etot>=e_min
    elif e_min!=None and e_max!=None:
        eindx=(cluster.etot<=e_max) * (cluster.etot>=e_min)
    else:
        eindx=cluster.id > -1

    indx=(cluster.r>=rmin) * (cluster.r<=rmax) * (cluster.m>=mmin) * (cluster.m<=mmax) * (cluster.kw>=se_min) * (cluster.kw<=se_max)


    indx=np.logical_and(indx,eindx)

    if len(cluster.id[indx])>0:
        subcluster=StarCluster(len(cluster.id[indx]),cluster.tphys,units=cluster.units,origin=cluster.origin,center=cluster.center)
        subcluster.add_stars(cluster.id[indx],cluster.m[indx],cluster.x[indx],cluster.y[indx],cluster.z[indx],cluster.vx[indx],cluster.vy[indx],cluster.vz[indx],cluster.kw[indx])

        if len(cluster.logl)>0:
            subcluster.add_se(cluster.kw[indx],cluster.logl[indx],cluster.logr[indx],cluster.ep[indx],cluster.ospin[indx])
        if len(cluster.id2)>0:
            bindx=np.in1d(cluster.id1,cluster.id[indx])
            subcluster.add_bse(cluster.id1[bindx],cluster.id2[bindx],cluster.kw1[bindx],cluster.kw2[bindx],cluster.kcm[bindx],cluster.ecc[bindx],cluster.pb[bindx],cluster.semi[bindx],cluster.m1[bindx],cluster.m2[bindx],cluster.logl1[bindx],cluster.logl2[bindx],cluster.logr1[bindx],cluster.logr2[bindx],cluster.ep1[bindx],cluster.ep2[bindx],cluster.ospin1[bindx],cluster.ospin2[bindx])
        if len(cluster.etot)>0:
            subcluster.add_energies(cluster.kin[indx],cluster.pot[indx],cluster.etot[indx])

        if new_center:
            subcluster.add_orbit(cluster.xgc+cluster.xc,cluster.ygc+cluster.yc,cluster.zgc+cluster.zc,cluster.vxgc+cluster.vxc,cluster.vygc+cluster.vyc,cluster.vzgc+cluster.vzc)
            subcluster.xc,subcluster.yc,subcluster.zc=0.0,0.0,0.0
            subcluster.vxc,subcluster.vyc,subcluster.vzc=0.0,0.0,0.0
            subcluster.xc,subcluster.yc,subcluster.zc,subcluster.vxc,subcluster.vyc,subcluster.vzc=subcluster.find_center(0.0,0.0,0.0)
        else:
            subcluster.add_orbit(cluster.xgc,cluster.ygc,cluster.zgc,cluster.vxgc,cluster.vygc,cluster.vzgc)
            subcluster.xc,subcluster.yc,subcluster.zc=cluster.xc,cluster.yc,cluster.zc
            subcluster.vxc,subcluster.vyc,subcluster.vzc=cluster.vxc,cluster.vyc,cluster.vzc

        subcluster.key_params()

        if not center0:
            subcluster.from_center()
        if origin0!='cluster':
            cluster.to_galaxy()

    else:
        subcluster=StarCluster(0,cluster.tphys)

    return subcluster

