"""For changing units

"""
__author__ = "Jeremy J Webb"


import numpy as np
from ..util.constants import *

def _convert_length(x,units,cluster):
    """Convert x from units to cluster.units

    Parameters
    ----------
    x : float
      measure of distance in units
    units : str
      units associated with x
    cluster : class
      StarCluster

    Returns
    -------
    x : float
      measure of distance with units the same as cluster.units

    History
    -------
    2022 - Written - Webb (UofT)
    """


    if units=='nbody':
        if cluster.units=='pckms' or cluster.units=='pcmyr':
            x*=cluster.rbar
        elif cluster.units=='kpckms' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
            x*=cluster.rbar/1000.0
        elif cluster.units=='galpy':
            x*=cluster.rbar/1000.0/cluster._ro
        elif cluster.units=='radec':
            x*=cluster.rbar/1000.0
            x=np.degrees(np.arctan2(x,cluster.zgc))

    elif units=='pckms' or units=='pcmyr':
        if cluster.units=='nbody':
            x/=cluster.rbar
        elif cluster.units=='kpckms' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
            x/=1000.0
        elif cluster.units=='galpy':
            x/=(1000.0*cluster._ro)
        elif cluster.units=='radec':
            x/=1000.0
            x=np.degrees(np.arctan2(x,cluster.zgc))

    elif units=='kpckms' or units=='kpcgyr' or units=='WDunits':
        if cluster.units=='nbody':
            x*=1000.0/cluster.rbar
        elif cluster.units=='pckms' or cluster.units=='pcmyr':
            x*=1000.0
        elif cluster.units=='galpy':
            x/=cluster._ro
        elif cluster.units=='radec':
            x=np.degrees(np.arctan2(x,cluster.zgc))

    elif units=='galpy':
        if cluster.units=='pckms' or cluster.units=='pcmyr':
            x*=(1000.0*cluster._ro)
        elif cluster.units=='kpckms' or cluster.units=='kpcgyr' or cluster.units=='WDunits':
            x*=cluster._ro
        elif cluster.units=='nbody':
            x*=(cluster._ro*1000.0/cluster.rbar)
        elif cluster.units=='radec':
            x*=cluster._ro
            x=np.degrees(np.arctan2(x,cluster.zgc))

    elif units=='radec':
        if cluster.units!='radec':
            dist=np.sqrt(cluster.xgc**2.+cluster.ygc**2.0+cluster.zgc**2.0)
            x=np.tan(np.radians(x))*dist

            if cluster.units=='nbody':
                x*=1000.0/cluster.rbar
            elif cluster.units=='pckms' or cluster.units=='pcmyr':
                x*=1000.0
            elif cluster.units=='galpy':
                x/=cluster._ro

    return x