""" Generic recipes for making key calculations

"""
__author__ = "Jeremy J Webb"
__all__ = [
    'nbinmaker',
    "binmaker",
    'roaming_nbinmaker',
    "roaming_binmaker",
    "power_law_distribution_function",
    "dx_function",
    "tapered_dx_function",
    "x_hist",
    "mean_prof",
    "smooth",
    "interpolate",
    "minimum_distance",
    "distance",
]

import numpy as np
import numba
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from ..util.plots import _plot,_lplot


def nbinmaker(x, nbin=10, nsum=False):
    """Split an array into bins with equal numbers of elements

    Parameters
    ----------
    x : float
      input array
    nbin : int
      number of bins
    nsum : bool
      return number of points in each bin (default: False)

    Returns
    -------
    x_lower : float
      lower bin values
    x_mid : float
      mean value in each bin
    x_upper : float
      upper bin values
    x_hist : 
      number of points in bin

    if nsum==True:
      x_sum : float
        sum of point values in each bin

    History
    -------
    2018 - Written - Webb (UofT)
    """
    x = np.asarray(x)

    xorder = np.argsort(x)

    x_lower = np.array([])
    x_upper = np.array([])
    x_hist = np.array([])
    x_sum = np.array([])
    x_mid = np.array([])

    for i in range(0, nbin):
        indx = int(float(i) * float(len(x)) / float(nbin))
        x_lower = np.append(x_lower, x[xorder[indx]])

    x_upper=x_lower[1:]
    x_upper=np.append(x_upper,np.amax(x))

    indx = x_lower != x_upper
    x_lower = x_lower[indx]
    x_upper = x_upper[indx]

    for i in range(0, np.sum(indx)):
        if i<np.sum(indx)-1:
            xindx = (x >= x_lower[i]) * (x < x_upper[i])
        else:
            xindx = (x >= x_lower[i])

        x_hist = np.append(x_hist, np.sum(xindx))
        x_sum = np.append(x_sum, np.sum(x[xindx]))
        x_mid = np.append(x_mid, x_sum[i] / x_hist[i])

    if nsum:
        return x_lower, x_mid, x_upper, x_hist, x_sum
    else:
        return x_lower, x_mid, x_upper, x_hist


def binmaker(x, nbin=10, nsum=False, steptype="linear"):
    """Split an array into bins of equal width

    Parameters
    ----------
    x : float
      input array
    nbin : int
      number of bins
    nsum : bool
      return number of points in each bin (default: False)
    steptype : str
      linear or logarithmic steps (default: linear)

    Returns
    -------
    x_lower : float
      lower bin values
    x_mid : float
      mean value in each bin
    x_upper : float
      upper bin values
    x_hist : 
      number of points in bin

    if nsum==True:
      x_sum : float
        sum of point values in each bin

    History
    -------
    2018 - Written - Webb (UofT)
    """

    x_hist = np.zeros(nbin)
    x_sum = np.zeros(nbin)
    x = np.array(x)

    if steptype == "linear":
        steps = np.linspace(np.amin(x), np.amax(x), nbin + 1)
    else:
        steps = np.logspace(np.log10(np.amin(x)), np.log10(np.amax(x)), nbin + 1)

    x_lower = steps[:-1]
    x_upper = steps[1:]

    x_mid = (x_upper + x_lower) / 2.0

    for j in range(0, nbin):
        if j<nbin-1:
            indx = (x >= x_lower[j]) * (x < x_upper[j])
        else:
            indx = (x >= x_lower[j]) * (x <= x_upper[j])

        x_hist[j] = len(x[indx])
        x_sum[j] = np.sum(x[indx])

    if nsum:
        return x_lower, x_mid, x_upper, x_hist, x_sum
    else:
        return x_lower, x_mid, x_upper, x_hist

def roaming_nbinmaker(x, nbin=10, ntot=20, nsum=False):

    """Split an array into bins with equal numbers of elements

    Parameters
    ----------
    x : float
      input array
    nbin : int
      number of bins to set bin fraction
    ntot : int
      number of total bins
    nsum : bool
      return number of points in each bin (default: False)

    Returns
    -------
    x_lower : float
      lower bin values
    x_mid : float
      mean value in each bin
    x_upper : float
      upper bin values
    x_hist : 
      number of points in bin

    if nsum==True:
      x_sum : float
        sum of point values in each bin

    History
    -------
    2018 - Written - Webb (UofT)
    """
    xs = np.sort(x)

    dx=1.0/float(nbin)
    
    xfrac_lower=np.linspace(0.,1-dx,ntot)
    xfrac_upper=xfrac_lower+dx
        
    x_lower = xs[(xfrac_lower*(len(x)-1)).astype(int)]
    x_upper = xs[(xfrac_upper*(len(x)-1)).astype(int)]

    x_hist = np.array([])
    x_sum = np.array([])
    x_mid = np.array([])

    indx = x_lower != x_upper
    x_lower = x_lower[indx]
    x_upper = x_upper[indx]

    for i in range(0, np.sum(indx)):
        if i<np.sum(indx)-1:
            xindx = (x >= x_lower[i]) * (x < x_upper[i])
        else:
            xindx = (x >= x_lower[i])

        x_hist = np.append(x_hist, np.sum(xindx))
        x_sum = np.append(x_sum, np.sum(x[xindx]))
        x_mid = np.append(x_mid, x_sum[i] / x_hist[i])

    if nsum:
        return x_lower, x_mid, x_upper, x_hist, x_sum
    else:
        return x_lower, x_mid, x_upper, x_hist
        
def roaming_binmaker(x, nbin=10, ntot=20, nsum=False, steptype="linear"):
    """Split an array into bins of equal width with a roaming average

    Parameters
    ----------
    x : float
      input array
    nbin : int
      number of bins to set bin width
    ntot : int
      number of total bins
    nsum : bool
      return number of points in each bin (default: False)
    steptype : str
      linear or logarithmic steps (default: linear)

    Returns
    -------
    x_lower : float
      lower bin values
    x_mid : float
      mean value in each bin
    x_upper : float
      upper bin values
    x_hist : 
      number of points in bin

    if nsum==True:
      x_sum : float
        sum of point values in each bin

    History
    -------
    2018 - Written - Webb (UofT)
    """

    if steptype=='linear':
    
        xmin=np.amin(x)
        xmax=np.amax(x)
        dx=np.fabs(xmax-xmin)/nbin

        x_lower=np.linspace(xmin,xmax-dx,ntot)
        x_upper=x_lower+dx
        x_mid=(x_upper+x_lower)/2.
        
    else:
        xmin=np.amin(np.log10(x))
        xmax=np.amax(np.log10(x))
        dx=np.fabs(xmax-xmin)/nbin
        
        x_lower=np.logspace(xmin,xmax-dx,ntot)
        x_upper=np.logspace(xmin+dx,xmax,ntot)
                
        x_mid=10.0**((np.log10(x_upper)+np.log10(x_lower))/2.)

    x_hist=np.array([])
    x_sum=np.array([])

    for j in range(0, len(x_lower)):
        if j<len(x_lower)-1:
            indx = (x >= x_lower[j]) * (x < x_upper[j])
        else:
            indx = (x >= x_lower[j]) * (x <= x_upper[j])

        x_hist = np.append(x_hist,np.sum(indx))
        x_sum = np.append(x_sum,np.sum(x[indx]))

    if nsum:
        return x_lower, x_mid, x_upper, x_hist, x_sum
    else:
        return x_lower, x_mid, x_upper, x_hist

def power_law_distribution_function(n, alpha, xmin, xmax):
    """Generate points from a power-law distribution function

    Parameters
    ----------
    n : int
      number of points
    alpha : float
      power-law slope of distribution function
    xmin,xmax : float
      minimum and maximum values of distribution

    Returns
    -------
    x : float
      array of values drawn from distribution

    History
    -------
    2019 - Written - Webb (UofT)
    """

    eta = alpha + 1.0

    if xmin == xmax:
        x = xmin
    elif alpha == 0:
        x = xmin + np.random.random(n) * (xmax - xmin)
    elif alpha > 0:
        x = xmin + np.random.power(eta, n) * (xmax - xmin)
    elif alpha < 0 and alpha != -1.0:
        x = (xmin ** eta + (xmax ** eta - xmin ** eta) * np.random.rand(n)) ** (
            1.0 / eta
        )
    elif alpha == -1:
        x = np.log10(xmin) + np.random.random(n) * (np.log10(xmax) - np.log10(xmin))
        x = 10.0 ** x

    return np.array(x)


def dx_function(x, nx=10, bintype="num", x_lower=None, x_mean=None,x_upper=None, plot=False, **kwargs):
    """Find distribution function using nx bins

    Parameters
    ----------
    x : float
      input arrayå
    nx : int
      number of bins (default : 10)
    bintype : str
      bin with equal number of stars per bin (num) or evenly in x (fix) (default: num)
    x_lower,x_mean,x_upper : float
      preset lower limit, mean value, and upper limit bins

    Returns
    -------
    x_mean : float
      mean value in each bin
    x_hist : float
      number of stars in each bin
    dx : float
      number of stars in each bin divided by width of bin
    alpha : float
      power law slope fit to dx vs x_mean
    ealpha : float
      error in alpha
    yalpha : float
      y-intercept of fit to log(dx) vs lod(x_mean)
    eyalpha : float
      error in yalpha

    History
    -------
    2018 - Written - Webb (UofT)
    """
    if x_lower is None:
        if bintype == "num":
            x_lower, x_mean, x_upper, x_hist = nbinmaker(x, nx)
        else:
            x_lower, x_mean, x_upper, x_hist = binmaker(x, nx)
    else:
        x_hist=np.array([])
        for i in range(0, len(x_lower)):
            indx = (x >= x_lower[i]) * (x < x_upper[i])
            x_hist = np.append(x_hist, np.sum(indx))

    dx = x_hist / (x_upper - x_lower)
    indx=dx>0

    lx_mean = np.log10(x_mean[indx])
    ldx = np.log10(dx[indx])

    (alpha, yalpha), V = np.polyfit(lx_mean, ldx, 1, cov=True)
    ealpha = np.sqrt(V[0][0])
    eyalpha = np.sqrt(V[1][1])

    if plot:
            filename = kwargs.get("filename", None)
            _plot(x_mean[indx], np.log10(dx[indx]), xlabel="x", ylabel="LOG(dN/dx)", **kwargs)
            xfit = np.linspace(np.min(x_mean), np.max(x_mean), nx)
            dxfit = 10.0 ** (alpha * np.log10(xfit) + yalpha)
            kwargs.pop("overplot",None)

            _lplot(
                xfit, np.log10(dxfit), overplot=True, label=(r"$\alpha$ = %f" % alpha), **kwargs
            )

            plt.legend()

            if filename != None:
                plt.savefig(filename)

    return x_mean[indx], x_hist[indx], dx[indx], alpha, ealpha, yalpha, eyalpha

def tapered_dx_function(x, nx=10, bintype="num", x_lower=None, x_mean=None,x_upper=None, plot=False, **kwargs):
    """Find distribution function using nx bins

    Parameters
    ----------
    x : float
      input arrayå
    nx : int
      number of bins (default : 10)
    bintype : str
      bin with equal number of stars per bin (num) or evenly in x (fix) (default: num)
    x_lower,x_mean,x_upper : float
      preset lower limit, mean value, and upper limit bins

    Returns
    -------
    x_mean : float
      mean value in each bin
    x_hist : float
      number of stars in each bin
    dx : float
      number of stars in each bin divided by width of bin
    alpha : float
      power law slope fit to dx vs x_mean
    ealpha : float
      error in alpha
    yalpha : float
      y-intercept of fit to log(dx) vs lod(x_mean)
    eyalpha : float
      error in yalpha

    History
    -------
    2018 - Written - Webb (UofT)
    """
    if x_lower is None:
        if bintype == "num":
            x_lower, x_mean, x_upper, x_hist = nbinmaker(x, nx)
        else:
            x_lower, x_mean, x_upper, x_hist = binmaker(x, nx)
    else:
        x_hist=np.array([])
        for i in range(0, len(x_lower)):
            indx = (x >= x_lower[i]) * (x < x_upper[i])
            x_hist = np.append(x_hist, np.sum(indx))

    dx = x_hist / (x_upper - x_lower)
    indx=dx>0

    lx_mean = np.log10(x_mean[indx])
    ldx = np.log10(dx[indx])

    (A, alpha, xc, beta), V=curve_fit(tapered_func,10.0**np.array(lx_mean),10.0**np.array(ldx) ,bounds=([0.,-1.*np.inf,np.amin(x),-1.*np.inf],[np.inf,np.inf,np.amax(x),np.inf]))

    eA = np.sqrt(V[0][0])
    ealpha = np.sqrt(V[1][1])
    exc = np.sqrt(V[2][2])
    ebeta = np.sqrt(V[3][3])

    if plot:
            filename = kwargs.get("filename", None)
            _plot(x_mean[indx], np.log10(dx[indx]), xlabel="x", ylabel="LOG(dN/dx)", **kwargs)
            xfit = np.linspace(np.min(x_mean), np.max(x_mean), nx)
            dxfit = tapered_func(xfit,A,alpha,xc,beta)

            kwargs.pop("overplot",None)

            _lplot(
                xfit, np.log10(dxfit), overplot=True, label=(r"$\alpha$ = %f" % alpha), **kwargs
            )

            plt.legend()

            if filename != None:
                plt.savefig(filename)

    return x_mean[indx], x_hist[indx], dx[indx], A, eA, alpha, ealpha, xc, exc, beta, ebeta

def tapered_func(x,A,alpha,xc,beta):

    dx=A*(x**alpha)*(1.0-np.exp(-1.*(x/xc)**beta))

    return dx

def x_hist(x, nx=10, bintype="num", bins=False, x_lower=None, x_mean=None,x_upper=None):
    """Find histogram data using nx bins

    Parameters
    ----------
    x : float
      input arrayå
    nx : int
      number of bins (default : 10)
    bintype : str
      bin with equal number of stars per bin (num) or evenly in x (fix) (default: num)
    bins : bool
      return bins
    x_lower,x_mean,x_upper : float
      preset lower limit, mean value, and upper limit bins

    Returns
    -------
    if bins:
        x_lower,x_mean,x_upper,x_hist
    else:
        x_mean : float
          mean value in each bin
        x_hist : float
          number of stars in each bin

    History
    -------
    2019 - Written - Webb (UofT)
    """

    if x_lower is None:
        if bintype == "num":
            x_lower, x_mean, x_upper, x_hist = nbinmaker(x, nx)
        else:
            x_lower, x_mean, x_upper, x_hist = binmaker(x, nx)
    else:
        x_hist=np.array([])
        for i in range(0, len(x_lower)):
            indx = (x >= x_lower[i]) * (x < x_upper[i])
            x_hist = np.append(x_hist, np.sum(indx))

    if bins:
        return x_lower,x_mean,x_upper,x_hist
    else:
        return x_mean,x_hist

def mean_prof(x, y, nbin=10, bintype="num", steptype="linear", median=False, x_lower=None, x_mean=None,x_upper=None):
    """ Calculate mean profile of parameter y that depends on x

    Parameters
    ----------
    x,y : float
      coordinates from which to measure the mean profile
    nbin : int
      number of bins
    bintype : str
      can be bins of fixed size ('fix') or equal number of stars ('num') (default: num)
    steptype : str
      for fixed size arrays, set step size to 'linear' or 'log'
    median : bool
      find median instead of mean (Default: False)
    x_lower,x_mean,x_upper : float
      preset lower limit, mean value, and upper limit bins

    Returns
    -------
    x_bin : float
      x values of mean profile
    y_bin : float
      y values of mean profile
    y_sig : float
      dispersion about the mean profile

    History
    -------
    2018 - Written - Webb (UofT)

    """
    if x_lower is None:
        if bintype == "num":
            x_lower, x_mid, x_upper, x_hist = nbinmaker(x, nbin)
        else:
            x_lower, x_mid, x_upper, x_hist = binmaker(x, nbin, steptype=steptype)
    else:
        x_mid=x_mean
        x_hist=np.array([])
        for i in range(0, len(x_lower)):
            indx = (x >= x_lower[i]) * (x < x_upper[i])
            x_hist = np.append(x_hist, np.sum(indx))

    y_bin = []
    y_sig = []
    x_bin = []

    for i in range(0, len(x_lower)):
        indx = (x >= x_lower[i]) * (x < x_upper[i])

        if True in indx:
            x_bin = np.append(x_bin, x_mid[i])
            if x_hist[i] > 1:
                if median:
                    y_bin = np.append(y_bin, np.median(y[indx]))
                else:
                    y_bin = np.append(y_bin, np.mean(y[indx]))

                y_sig = np.append(y_sig, np.std(y[indx]))
            elif x_hist[i] == 1:
                y_bin = np.append(y_bin, y[indx])
                y_sig = np.append(y_sig, 0.0)
            else:
                y_bin = np.append(y_bin, 0.0)
                y_sig = np.append(y_sig, 0.0)

    return np.array(x_bin), np.array(y_bin), np.array(y_sig)


def smooth(x, y, dx, bintype="num", median=False):
    """Smooth a profile

    Parameters
    ----------
    x,y : float
      coordinates from which to measure the mean profile
    dx : float
      width of x smoothening bin
    bintype : str
      can be bins of fixed size ('fix') or equal number of stars ('num') (default: num)
    steptype : str
      for fixed size arrays, set step size to 'linear' or 'log'
    median : bool
      find median instead of mean (Default: False)

    Returns
    -------
    x_bin : float
      x values of mean profile
    y_bin : float
      y values of mean profile
    y_sig : float
      dispersion about the mean profile

    Other Parameters
    ----------------
    dx : float
      width of smoothening bin

    History
    -------
    2018 - Written - Webb (UofT)

    """
    x = np.array(x)
    y = np.array(y)

    # Smooth by intervals in dx
    nbin = int((np.max(x) - np.min(x)) / dx)

    if bintype == "num":
        x_lower, x_mid, x_upper, x_hist = nbinmaker(x, nbin)
    else:
        x_lower, x_mid, x_upper, x_hist = binmaker(x, nbin)

    y_bin = []
    y_sig = []
    x_bin = []
    y_max = []
    y_min = []

    for i in range(0, nbin):
        indx = (x >= x_lower[i]) * (x < x_upper[i])

        if True in indx:
            x_bin = np.append(x_bin, x_mid[i])
            if x_hist[i] > 1:
                if median:
                    y_bin = np.append(y_bin, np.median(y[indx]))
                else:
                    y_bin = np.append(y_bin, np.mean(y[indx]))
                y_sig = np.append(y_sig, np.std(y[indx]))
                y_min = np.append(y_min, np.min(y[indx]))
                y_max = np.append(y_max, np.max(y[indx]))
            elif x_hist[i] == 1:
                y_bin = np.append(y_bin, y[indx])
                y_sig = np.append(y_sig, 0.0)
                y_min = np.append(y_min, y[indx])
                y_max = np.append(y_max, y[indx])
            else:
                y_bin = np.append(y_bin, 0.0)
                y_sig = np.append(y_sig, 0.0)
                y_min = np.append(y_min, 0.0)
                y_max = np.append(y_max, 0.0)

    return x_bin, y_bin, y_sig, y_min, y_max


def interpolate(r1, r2, x=None, y=None):
    """Perform simple linear interpolation between two points in 2D
      - one of x or y must be defined

    Parameters
    ----------

    r1,r2 : float
      x,y coordinates from which to interpolate

    x : float
      x-value from which to interpolate y (default: None)

    y : float
      y-value from which to interpolate x (default: None)

    Returns
    -------

    val : float
      interpolated value

    History

    2019 - Written - Webb (UofT)

    """
    x1, y1 = r1
    x2, y2 = r2

    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1

    if x != None:
        val= m * x + b
    elif y != None:
        val=(y - b) / m
    else:
        print("NO INTERMEDIATE COORDINATE GIVEN")
        val=0

    return val

@numba.njit
def minimum_distance(x):
    """Find distance to each point's nearest neighbour

    Parameters
    ----------
    x : float (array)
      3D position of each point of the form [x,y,z].Transpose

    Returns
    -------
    min_distance : float
      distance to each points nearest neighbour

    History
    -------
    2019 - Written - Webb (UofT)
    """
    min_distance = [-1.0] * len(x)
    for i in range(len(x) - 1):
        for j in range(i + 1, len(x)):
            r = distance(x[i], x[j])
            if min_distance[i] < 0:
                min_distance[i] = r
            else:
                min_distance[i] = np.minimum(min_distance[i], r)

            if min_distance[j] < 0:
                min_distance[j] = r
            else:
                min_distance[j] = np.minimum(min_distance[j], r)

    return min_distance


@numba.njit
def distance(x1, x2):
    """Find distance between two points (made for use with numba)

    Parameters
    ----------
    x1 : float
      3D position of first point of the form [x,y,z]
    x2 : float
      3D position of first point of the form [x,y,z]

    Returns
    -------
    distance : float
      distance between the two points

    History
    -------
    2019 - Written - Webb (UofT)
    """
    dx = x2[0] - x1[0]
    dy = x2[1] - x1[1]
    dz = x2[2] - x1[2]

    r = (dx * dx + dy * dy + dz * dz) ** 0.5

    return r
