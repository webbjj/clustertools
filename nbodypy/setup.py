import numpy as np
from galpy.util import bovy_conversion
import os
import limepy
from limepy import limepy
from limepy import sample

from .cluster import StarCluster
from .profiles import m_prof
from .orbit import get_cluster_orbit

def setup_cluster(ctype,units0='realpc',origin0='cluster',orbit=None,gcname='list',**kwargs):
    #setup_cluster is the main function for setting up clusters, built around LIMEPY
   
    #ctype = limepy - use limepy to setup initial cluster conditions (limepy parameters read in via **kwargs)

    #Default units are set to realpc and default origin is set to cluster. These are really only used when reading in nbodypy snapshots, as other codes have their own default units and origin

    #orbit - allows for a galpy orbit to be used to set the cluster's orbit

    #kwargs:
    #Limepy parameters
   

    #Generate cluster using plummer of limepy
    if ctype=='plummer':
        M=float(kwargs.get('M'))
        rh=float(kwargs.get('rh'))
        cluster=get_plummer(M,rh)
    elif ctype=='limepy':
        g=kwargs.pop('g')
        cluster=get_limepy(g=g,**kwargs)
    elif ctype=='woolley':
        g=kwargs.pop('g',0)
        cluster=get_limepy(g=g,**kwargs)        
    elif ctype=='king':
        g=kwargs.pop('g',1)
        cluster=get_limepy(g=g,**kwargs)
    elif ctype=='wilson':
        g=kwargs.pop('g',2)
        cluster=get_limepy(g=g,**kwargs)
    else:
        source=kwargs.pop('source','default')
        mbar=kwargs.pop('mbar',0.4)
        names=kwargs.pop('names',False)
        cluster=get_cluster(gcname,source,mbar,names)


   #Add galpy orbit if given
    if orbit!=None:
        cluster.orbit=orbit
        t=(cluster.tphys/1000.)/bovy_conversion.time_in_Gyr(ro=8.,vo=220.)
        cluster.add_orbit(orbit.x(t),orbit.y(t),orbit.z(t),orbit.vx(t),orbit.vy(t),orbit.vz(t),'realkpc')

    cluster.key_params()

    return cluster

def get_limepy(g=1,**kwargs):
    phi0=float(kwargs.get('phi0'))
    project=bool(kwargs.get('project',False))

    if 'M' in kwargs:
        units='realpc'
        M=float(kwargs.get('M'))
        if 'rt' in kwargs:
            rt=float(kwargs.get('rt'))  
            lmodel=limepy(phi0,g,M=M,rt=rt,project=project)
        elif 'rv' in kwargs:
            rv=float(kwargs.get('rv'))      
            lmodel=limepy(phi0,g,M=M,rv=rv,project=project)
        elif 'rh' in kwargs:
            rh=float(kwargs.get('rh'))      
            lmodel=limepy(phi0,g,M=M,rh=rh,project=project)
        elif 'r0' in kwargs:
            r0=float(kwargs.get('r0'))      
            lmodel=limepy(phi0,g,M=M,r0=r0,project=project)
        else:
            lmodel=limepy(phi0,g,M=M,r0=1.,project=project)
    else:
        units='nbody'
        lmodel=limepy(phi0,g,G=1, M=1, rv=1,project=project)

    N=int(kwargs.get('N',1000))

    ldata=sample(lmodel,N=N)

    cluster=StarCluster(N,units=units,origin='cluster')
    cluster.ctype='limepy'
    cluster.add_stars(np.linspace(1,N,N,dtype=int),ldata.m,ldata.x,ldata.y,ldata.z,ldata.vx,ldata.vy,ldata.vz)
    cluster.find_centre()
    cluster.key_params()

    return cluster

def get_plummer(M,rm):

    cluster=StarCluster()
    
    return cluster

def get_imf(mlimits = (0.1, 50.0), alphas = (-1.35), ntot=1, mtot=None, do_random=True):
    """
    NAME:

       get_imf

    PURPOSE:

       Generate stellar masses from a defined IMF
       Notes:
        -- This has been heaviliy borrowed (stolen) from AMUSE (amusecode.org), specifically
           the routine /src/amuse/ic/brokenimf.py

    INPUT:

       mlimits - mass limits (default: 0.1-50.0)

       alphas - slopes of the mass function between mass limits (default: -1.35 - Salpeter)

       ntot - number of stars to be generated (default: 1)

       mtot - total mass of system (overrides ntot) (default: None)

       do_random - use randomly generated masses or evenly distributed masses (default: True)

    OUTPUT:

       trelax

    HISTORY:

       2019 - Written - Webb (UofT)

    """ 
    
    mlimits = np.array(mlimits)
    alphas = np.array(alphas)

    #If mtot is given, initially assume all stars are equal to the lower mass limit
    if mtot!=None:
        ntot=int(mtot/mlimits[0])
        print(ntot)

    if do_random:
        random = np.random.random(ntot)
    else:
        random = np.linspace(0.0, 1.0, ntot)

    nbins = len(alphas)
    
    #Calculate fraction per bin
    nbin = []

    for i in range(0,nbins):
        alpha=alphas[i]
        if alpha == -1:
            factor = np.log(mlimits[i+1] / mlimits[i])
        else:
            factor = (mlimits[i+1]**(alpha+1) - mlimits[i]**(alpha+1)) / (alpha+1)

        for j in range(0,nbins - i - 1):
            factor *= mlimits[-j-2]**(alphas[-j-1]-alphas[-j-2])

        nbin.append(factor)
    total = sum(nbin, 0.0)
    fbin=np.array(nbin)/total 

    #Calculate cumultive fraction per bin, factors, and inverse alphas 
    cumulative_fractions = np.array([sum(fbin[:i]) for i in range(nbins+1)])
    cumulative_fractions[-1] = 1.0     # In case of round-off errors
    factors = pow(mlimits[1:] / mlimits[:-1], alphas + 1.0) - 1.0
    inv_alpha1s = np.array([np.inf if alpha==-1 else (1.0 / (alpha + 1.0)) for alpha in alphas])

    #Calculate masses
    indices = np.searchsorted(cumulative_fractions[:-1], random, 'right') - 1
    scaled = ((random - cumulative_fractions[indices]) / 
        (cumulative_fractions[indices+1] - cumulative_fractions[indices]))

    result = np.empty_like(random)
    zerodiv = alphas[indices]==-1
    normal = np.logical_not(zerodiv)
    result[zerodiv] = pow(mlimits[1:][indices[zerodiv]] / mlimits[:-1][indices[zerodiv]], scaled[zerodiv])
    result[normal] = pow(1.0 + factors[indices[normal]] * scaled[normal], inv_alpha1s[indices[normal]])

    masses=mlimits[:-1][indices] * result

    #If mtot is given, only return enough stars so that total mass ~ mtot
    if mtot!=None:
        ntot=int(mtot/np.mean(masses))
        masses=masses[0:ntot]

    return np.array(masses)

def c_to_w0(c,invert=False):
    #From gridfit (McLaughlin)
    w0=np.array([0.300000,0.400000,0.500000,0.600000,0.700000,0.800000,0.900000,1.000000, \
        1.100000,1.200000,1.300000,1.400000,1.500000,1.600000,1.700000,1.800000,1.900000, \
        2.000000,2.100000,2.200000,2.300000,2.400000,2.500000,2.600000,2.700000,2.800000, \
        2.900000,3.000000,3.100000,3.200000,3.300000,3.400000,3.500000,3.600000,3.700000, \
        3.800000,3.900000,4.000000,4.100000,4.200000,4.300000,4.400000,4.500000,4.600000, \
        4.700000,4.800000,4.900000,5.000000,5.100000,5.200000,5.300000,5.400000,5.500000, \
        5.600000,5.700000,5.800000,5.900000,6.000000,6.100000,6.200000,6.300000,6.400000, \
        6.500000,6.600000,6.700000,6.800000,6.900000,7.000000,7.100000,7.200000,7.300000, \
        7.400000,7.500000,7.600000,7.700000,7.800000,7.900000,8.000000,8.100000,8.200000, \
        8.300000,8.400000,8.500000,8.600000,8.700000,8.800000,8.900000,9.000000,9.100000 \
        ,9.200000,9.300000,9.400000,9.500000,9.600000,9.700000,9.800000,9.900000,10.000000, \
        10.100000,10.200000,10.300000,10.400000,10.500000,10.600000,10.700000,10.800000, \
        10.900000,11.000000,11.100000,11.200000,11.300000,11.400000,11.500000,11.600000, \
        11.700000,11.800000,11.900000,12.000000,12.100000,12.200000,12.300000,12.400000, \
        12.500000,12.600000,12.700000,12.800000,12.900000,13.000000,13.100000,13.200000, \
        13.300000,13.400000,13.500000,13.600000,13.700000,13.800000,13.900000,14.000000, \
        14.100000,14.200000,14.300000,14.400000,14.500000,14.600000,14.700000,14.800000, \
        14.900000,15.000000,15.100000,15.200000,15.300000,15.400000,15.500000,15.600000, \
        15.700000,15.800000,15.900000,16.000000,16.100000,16.200000,16.300000,16.400000, \
        16.500000,16.600000,16.700000,16.800000,16.900000,17.000000,17.100000,17.200000, \
        17.300000,17.400000,17.500000,17.600000,17.700000,17.800000,17.900000,18.000000, \
        18.100000,18.200000,18.300000,18.400000,18.500000,18.600000,18.700000,18.800000, \
        18.900000,19.000000,19.100000,19.200000,19.300000,19.400000,19.500000])
    
    conc=np.array([1.004710,1.171350,1.322660,1.463770,1.597760,1.726690,1.851980,1.974730, \
        2.095780,2.215830,2.335470,2.455190,2.575460,2.696690,2.819260,2.943540,3.069880, \
        3.198640,3.330170,3.464800,3.602910,3.744850,3.891000,4.041760,4.197530,4.358750, \
        4.525890,4.699410,4.879840,5.067720,5.263660,5.468270,5.682240,5.906290,6.141230, \
        6.387900,6.647220,6.920200,7.207920,7.511580,7.832460,8.171960,8.531630,8.913140, \
        9.318330,9.749200,10.208000,10.697100,11.219100,11.777000,12.374000,13.013600, \
        13.699700,14.436500,15.228700,16.081400,17.000300,17.991600,19.062000,20.218900, \
        21.470400,22.825300,24.293000,25.883700,27.608600,29.479300,31.508300,33.708600, \
        36.093800,38.678000,41.475300,44.499900,47.765800,51.286400,55.074300,59.140700, \
        63.495500,68.146800,73.100800,78.361400,83.930800,89.808800,95.993700,102.482000, \
        109.270000,116.352000,123.724000,131.381000,139.319000,147.537000,156.034000, \
        164.810000,173.871000,183.222000,192.871000,202.828000,213.106000,223.721000, \
        234.690000,246.032000,257.769000,269.924000,282.522000,295.592000,309.162000, \
        323.263000,337.930000,353.196000,369.099000,385.679000,402.977000,421.038000, \
        439.907000,459.634000,480.270000,501.871000,524.493000,548.199000,573.053000, \
        599.122000,626.479000,655.201000,685.368000,717.064000,750.382000,785.415000, \
        822.265000,861.038000,901.847000,944.812000,990.059000,1037.720000,1087.940000, \
        1140.860000,1196.650000,1255.460000,1317.470000,1382.870000,1451.860000,1524.620000, \
        1601.400000,1682.400000,1767.870000,1858.070000,1953.250000,2053.690000,2159.690000, \
        2271.550000,2389.600000,2514.160000,2645.600000,2784.280000,2930.590000,3084.940000, \
        3247.740000,3419.440000,3600.510000,3791.430000,3992.700000,4204.840000,4428.420000, \
        4664.010000,4912.200000,5173.620000,5448.940000,5738.830000,6044.000000,6365.210000, \
        6703.230000,7058.880000,7433.020000,7826.530000,8240.350000,8675.460000,9132.880000, \
        9613.680000,10119.000000,10650.000000,11207.900000,11794.100000,12409.800000,13056.600000, \
        13735.900000,14449.300000,15198.400000,15985.100000,16811.200000,17678.500000,18589.200000, \
        19545.300000,20549.200000,21603.100000,22709.700000])

    if invert:

        w=c
        indx=np.argmin(abs(w0-w))

        if w0[indx] < w:
            m=(conc[indx+1]-conc[indx])/(w0[indx+1]-w0[indx])
        else:
            m=(conc[indx]-conc[indx-1])/(w0[indx]-w0[indx-1])

        b=conc[indx]-m*w0[indx]

        return np.log10(m*w+b)

    else:

        c=10.0**c

        indx=np.argmin(abs(conc-c))

        if conc[indx] < c:
            m=(w0[indx+1]-w0[indx])/(conc[indx+1]-conc[indx])
        else:
            m=(w0[indx]-w0[indx-1])/(conc[indx]-conc[indx-1])

        b=w0[indx]-m*conc[indx]

        return m*c+b

def w0_to_c(w0):
    return c_to_w0(w0,invert=True)

def get_cluster(gcname='list',source='default',mbar=0.4,names=False):

    ddata=np.loadtxt('/Users/webbjj/Codes/nbodypy/tables/deBoer2019.dat',str,skiprows=1)
    dname=ddata[:,0]

    hdata=np.loadtxt('/Users/webbjj/Codes/nbodypy/tables/harris2010.dat',str,skiprows=2)
    hname=hdata[:,0]
    hname2=hdata[:,1]

    name=[]

    if isinstance(gcname,str):

        if (gcname == 'list' or gcname == 'all') and ('deboer' in source or 'deBoer' in source):
            for i in range(0,len(dname)):
                print(dname[i])
                name.append(dname[i])

        elif (gcname == 'list' or gcname == 'all')  and ('harris' in source or 'Harris' in source):
            for i in range(0,len(hname)):
                print(hname[i],hname2[i])
                name.append(hname[i])
        elif gcname == 'list' or gcname == 'all':
            for i in range(0,len(dname)):
                print('DEBOER: ', dname[i])
                name.append(dname[i])

            indx=np.in1d(hname,dname,invert=True) * np.in1d(hname2,dname,invert=True)
            for i in range(0,np.sum(indx)):
                print('HARRIS: ',hname[indx][i])
                name.append(hname[indx][i])

        else:
            gcname=gcname.upper()
            if (source=='default' or 'deboer' in source or 'deBoer' in source) and gcname in dname:
                cluster=get_deBoer_cluster(ddata,gcname,mbar,names)
            elif (source=='default' or 'harris' in source or 'Harris' in source) and (gcname in hname or gcname in hname2):
                cluster=get_harris_cluster(hdata,gcname,mbar,names)

            if names:
                return cluster,gcname
            else:
                return cluster
    else:
        name=gcname

    if len(name) > 0 and 'list' not in gcname:
        cluster=[]
        for i in range(0,len(name)):
            name[i]=name[i].upper()
            if (source=='default' or 'deboer' in source or 'deBoer' in source) and name[i] in dname:
                cluster.append(get_deBoer_cluster(ddata,name[i],mbar,names))
            elif (source=='default' or 'harris' in source or 'Harris' in source) and (name[i] in hname or name[i] in hname2):
                cluster.append(get_harris_cluster(hdata,name[i],mbar,names))
            else:
                print('COULD NOT FIND CLUSTER %s' % name[i])
            cluster[-1].ctype=name[i]

        if names:
            return cluster,name
        else:
            return cluster
    elif 'list' in gcname and names:
        return name
    else:
        return
    

def get_deBoer_cluster(data,gcname,mbar=0.4,names=False):
    name=data[:,0]

    gcname=gcname.upper()
    indx=(name==gcname)
    i_d=np.argwhere(indx==True)[0]
    print('GETTING: ',name[i_d])
    W_lime=data[i_d,1].astype(float)
    g_lime=data[i_d,3].astype(float)
    rt_lime=data[i_d,5].astype(float)
    M_lime=data[i_d,7].astype(float)
    N=M_lime/mbar
    cluster=get_limepy(g=g_lime,phi0=W_lime,M=M_lime,rt=rt_lime,N=N)
    cluster.orbit=get_cluster_orbit(name[i_d])
    if cluster.orbit!=-1:
        cluster.add_orbit(cluster.orbit.x(),cluster.orbit.y(),cluster.orbit.z(),cluster.orbit.vx(),cluster.orbit.vy(),cluster.orbit.vz(),ounits='realkpc')
    else:
        cluster.orbit=None

    return cluster

def get_harris_cluster(data,gcname,mbar=0.4,names=False):

    name=data[:,0]
    name2=data[:,1]

    indx=np.logical_or(np.in1d(name,gcname),np.in1d(name2,gcname))
    i_d=np.argwhere(indx==True)[0]
    print('GETTING: ',name[i_d])
    mgc=data[i_d,2].astype(float)
    rc=data[i_d,3].astype(float)
    rh=data[i_d,4].astype(float)
    rl=data[i_d,5].astype(float)
    c=np.log10(rl/rc)
    w0=c_to_w0(c)
    N=mgc/mbar
    cluster=get_limepy(g=1.,phi0=w0,M=mgc,rt=rl,N=N)
    cluster.orbit=get_cluster_orbit(name[i_d])
    if cluster.orbit!=-1:
        cluster.add_orbit(cluster.orbit.x(),cluster.orbit.y(),cluster.orbit.z(),cluster.orbit.vx(),cluster.orbit.vy(),cluster.orbit.vz(),ounits='realkpc')
    else:
        cluster.orbit=None

    return cluster

