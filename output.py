#Output files containing key cluster properties 
#INCOMPLETE

#Output a snapshot
def snapout(cluster,filename):

    fileout=open(filename,'w')

    for i in range(0,cluster.ntot):
        fileout.write("%i %f %f %f %f %f %f %f" % (cluster.id[i],cluster.x[i],cluster.y[i],cluster.z[i],cluster.vx[i],cluster.vy[i],cluster.vz[i]))

    fileout.close()
