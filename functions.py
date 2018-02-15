import math
import numpy as np

#TO DO Add References
#Half-Mass Relaxation Time (Units = Myr)
def Half_Mass_Relaxation_Time(cluster,full=False):
    print('DEBUG: ',cluster.r50,cluster.mtot,cluster.ntot)

    p50=3.0*(0.5*cluster.mtot)/(4.0*math.pi*(cluster.r50**3.0))
    
    #If full use caluclation of trh that does not assume virial equilibrium
    if full:
        trh=0.34*(cluster.sigv50**3.0)/(((4.302e-3)**2.0)*(1.0/3.0)*p50*math.log(0.4*cluster.mtot/(1.0/3.0)))
        trh=trh*3.086e13/(3600.0*24.0*365.0*1000000.0)
    else:
        trh = float(cluster.ntot)/math.log(0.4*float(cluster.ntot),10)
        trh=0.858*trh*math.sqrt(math.pow(cluster.r50,3.0)/cluster.mtot)

    cluster.trh=trh

def Inter_Stellar_Distances(cluster):
    s2sd=np.zeros((cluster.ntot,cluster.ntot))
    for i in range(0,cluster.ntot):
        for j in range(i,cluster.ntot):
            d=math.sqrt((cluster.x[i]-cluster.x[j])**2.0+(cluster.y[i]-cluster.y[j])**2.0+(cluster.z[i]-cluster.z[j])**2.0)
            s2sd[i][j]=d
            s2sd[j][i]=d

    return s2sd
