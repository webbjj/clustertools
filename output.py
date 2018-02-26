#Output files containing key cluster properties 
from functions import *

#Output dvprof.dat (WIP)
def eta_out(cluster,fileout):
    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    m_lower=[0.1,0.1,0.1]
    m_upper=[50.0,0.5,1.8]


    yeta_all=[]

    rn=rlagrange(cluster,nlagrange=10)
    
    for i in range(0,len(m_lower)):
        for j in range(0,len(rn)):
            if j==0:
                m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10,rmin=0.0,rmax=rn[i])
            else:
                m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10,rmin=rn[i-1],rmax=rn[i])
            fileout.write("%f " % eta)
            yeta_all.append(yeta)
        
        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10,rmin=0.0,rmax=cluster.rm)
        fileout.write("%f " % eta)
        yeta_all.append(yeta)

        m_mean,sigvm,eta,eeta,yeta,eyeta=eta_function(cluster,mmin=m_lower[i],mmax=m_upper[i],nmass=10)
        fileout.write("%f " % eta)
        yeta_all.append(yeta)

    for yeta in yeta_all:
        fileout.write("%f " % yeta)

    fileout.write("\n")


#Output dalpha_prof.dat (TODO - NEED TO ADD PROJECTED VALUES)
def dalpha_out(cluster,fileout):

    fileout.write("%f %f " % (cluster.tphys,cluster.mtot))

    m_lower=[0.1,0.3,0.5,0.1,0.3,0.5]
    m_upper=[0.5,0.8,0.8,0.5,0.8,0.8]
    r_lower=[0.0,0.0,0.0,None,None,None]
    r_upper=[cluster.rm,cluster.rm,cluster.rm,None,None,None]

    for i in range(0,len(m_lower)):

        m_mean,m_hist,alpha,ealpha,yalpha,eyalpha=mass_function(cluster,m_lower[i],m_upper[i],10,rmin=r_lower[i],rmax=r_upper[i])
        alphag.append(alpha)

        lrprofn,aprof,da,eda,yda,eyda=npy.alpha_prof(cluster,m_lower[i],m_upper[i],10,rmin=r_lower[i],rmax=r_upper[i],nrad=10)
        dalpha.append(da)
        ydalpha.append(yda)

    for ag in alphag:
        fileout.write("%f " % ag)
    for da in dalpha:
        fileout.write("%f " % da)
    for yda in ydalpha:
        fileout.write("%f " % yda)
    fileout.write("%f\n" % cluster.rm)

#Output a snapshot
def snapout(cluster,filename):

    fileout=open(filename,'w')

    for i in range(0,cluster.ntot):
        fileout.write("%i %f %f %f %f %f %f %f" % (cluster.id[i],cluster.x[i],cluster.y[i],cluster.z[i],cluster.vx[i],cluster.vy[i],cluster.vz[i]))

    fileout.close()
