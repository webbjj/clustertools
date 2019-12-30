# Generic recipes for making key calculations

import numpy as np
import numba
import matplotlib.pyplot as plt
from ..util.plots import *


def nbinmaker(x, nbin=10, nsum=False):
    """
  NAME:

     nbinmaker

  PURPOSE:

     Split an array into nbin bin's of equal number elements

  INPUT:

     x - input array

     nbin - number of bins

     nsum - return number of points in each bin

  OUTPUT:

     x_lower,x_mid,x_upper,x_hist (if nsum==False)

     x_lower,x_mid,x_upper,x_hist,x_sum (if nsum==True)

  HISTORY:

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
        indx = int(float(i + 1) * float(len(x)) / float(nbin)) - 1
        x_upper = np.append(x_upper, x[xorder[indx]])

    indx = x_lower != x_upper
    x_lower = x_lower[indx]
    x_upper = x_upper[indx]

    for i in range(0, np.sum(indx)):
        indx = (x >= x_lower[i]) * (x < x_upper[i])
        x_hist = np.append(x_hist, np.sum(indx))
        x_sum = np.append(x_sum, np.sum(x[indx]))
        x_mid = np.append(x_mid, x_sum[i] / x_hist[i])

    if nsum:
        return x_lower, x_mid, x_upper, x_hist, x_sum
    else:
        return x_lower, x_mid, x_upper, x_hist


def binmaker(x, nbin=10, nsum=False, steptype="linear"):
    """
  NAME:

     binmaker

  PURPOSE:

     Split an array into nbin bin's of equal size

  INPUT:

     x - input array

     nbin - number of bins

     nsum - return number of points in each bin

     steptype - linear or logarithmic steps

  OUTPUT:

     x_lower,x_mid,x_upper,x_hist (if nsum==False)

     x_lower,x_mid,x_upper,x_hist,x_sum (if nsum==True)

  HISTORY:

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
        indx = (x >= x_lower[j]) * (x < x_upper[j])
        x_hist[j] = len(x[indx])
        x_sum[j] = np.sum(x[indx])

    if nsum:
        return x_lower, x_mid, x_upper, x_hist, x_sum
    else:
        return x_lower, x_mid, x_upper, x_hist


def power_law_distribution_function(n, alpha, xmin, xmax):
    """
  NAME:

     power_law_distribution_function

  PURPOSE:

     Generate points from a power-law distribution function
  INPUT:

     n - number of points

     alpha - power-law slope of distribution function

     xmin,xmax - minimum and maximum values of distribution

  OUTPUT:

     c

  HISTORY:

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


def distribution_function(n, x=None, xfunc=None, galpy=False, galpy_param=None):
    """
  rad=np.array([])
  n=10000
  nhalo=0

  #Amp - total mass of sub-halos
  amp=pot.mass(100./8.,0.,ro=8.,vo=220.)

  #Assume they all have the same mass
  m=np.ones(n)*amp/float(n)

  while nhalo<n:
      rtemp=np.random.rand()*100.0
      
      if np.random.rand()*amp > pot.mass(rtemp/8.,0.,ro=8.,vo=220.):
          rad=np.append(rad,rtemp)
          nhalo+=1
  """
    return 0


def dx_function(x, nx=10, bintype="num", x_lower=None, x_mean=None,x_upper=None, plot=False, **kwargs):
    """
  NAME:

     dx_function

  PURPOSE:

     Find distribution function using nx bins containing an equal number of points
  INPUT:

     x - input array

     nx - number of bins

     bintype - bin with equal number of stars per bin (bin) or evenly in x (fix) (default: num)

     x_lower,x_mean,x_upper - preset bins

  OUTPUT:

     x_mean,x_hist,dx,alpha,ealpha,yalpha,eyalpha


  HISTORY:

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

    lx_mean = np.log10(x_mean)
    dx = x_hist / (x_upper - x_lower)
    ldx = np.log10(dx)

    (alpha, yalpha), V = np.polyfit(lx_mean, ldx, 1, cov=True)
    ealpha = np.sqrt(V[0][0])
    eyalpha = np.sqrt(V[1][1])

    if plot:
            filename = kwargs.get("filename", None)
            nplot(x_mean, np.log10(dx), xlabel="M", ylabel="LOG(dN/dM)", **kwargs)
            xfit = np.linspace(np.min(x_mean), np.max(x_mean), nx)
            dxfit = 10.0 ** (alpha * np.log10(xfit) + yalpha)
            nlplot(
                xfit, np.log10(dxfit), overplot=True, label=(r"$\alpha$ = %f" % alpha)
            )

            plt.legend()

            if filename != None:
                plt.savefig(filename)

    return x_mean, x_hist, dx, alpha, ealpha, yalpha, eyalpha

def x_hist(x, nx=10, bintype="num", x_lower=None, x_mean=None,x_upper=None):
    """
  NAME:

     x_hist

  PURPOSE:

     Find histogram
  INPUT:

     x - input array

     nx - number of bins

     bintype - bin with equal number of stars per bin (bin) or evenly in x (fix) (default: num)

     x_lower,x_mean,x_upper - preset bins

  OUTPUT:

     x_mean,x_hist


  HISTORY:

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

    return x_mean,x_hist

def mean_prof(x, y, nbin=10, bintype="fix", steptype="linear", median=False):
    """
  NAME:

     mean_prof

  PURPOSE:

     Calculate mean profile of parameter y that depends on x  
  INPUT:

     x,y coordinates from which to measure the mean profile

     nbin - number of bins

     bintype - can be bins of fixed size ('fix') or equal number of stars ('num')

     steptype - for fixed size arrays, set step size to 'linear' or 'log'

     median - find median instead of mean (Default: False)

  OUTPUT:

     x_bin,y_bin,y_sig

  HISTORY:

     2018 - Written - Webb (UofT)

  """
    if bintype == "num":
        x_lower, x_mid, x_upper, x_hist = nbinmaker(x, nbin)
    else:
        x_lower, x_mid, x_upper, x_hist = binmaker(x, nbin, steptype=steptype)

    y_bin = []
    y_sig = []
    x_bin = []

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
            elif x_hist[i] == 1:
                y_bin = np.append(y_bin, y[indx])
                y_sig = np.append(y_sig, 0.0)
            else:
                y_bin = np.append(y_bin, 0.0)
                y_sig = np.append(y_sig, 0.0)

    return np.array(x_bin), np.array(y_bin), np.array(y_sig)


def smooth(x, y, nbin=10, bintype="fix", median=False, **kwargs):
    """
  NAME:

     smooth

  PURPOSE:

     Smooth a profile

  INPUT:

     x,y coordinates from which to measure the mean profile

     nbin - number of bins to smooth over

     bintype - can be bins of fixed size ('fix') or equal number of stars ('num')

     median - find median instead of mean (Default: False)

  KWARGS:

     dx - width of smoothening bin

  OUTPUT:

     x_bin,y_bin,y_sig

  HISTORY:

     2018 - Written - Webb (UofT)

  """
    x = np.array(x)
    y = np.array(y)

    # Smooth by intervals in dx
    if "dx" in kwargs:
        dx = float(kwargs.get("dx"))
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
    """
  NAME:

     interpolate

  PURPOSE:

     Perform simple linear interpolate between two points in 2D

  INPUT:

     r1,r2 - x,y coordinates from which to interpolate

     x - x-value from which to interpolate y

     y - y-value from which to interpolate x

  OUTPUT:

     interpolated value

  HISTORY:

     2019 - Written - Webb (UofT)

  """
    x1, y1 = r1
    x2, y2 = r2

    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1

    if x != None:
        return m * x + b
    elif y != None:
        return (y - b) / m
    else:
        print("NO INTERMEDIATE COORDINATE GIVEN")
        return 0


def rotate(x, y, z, thetax=0.0, thetay=0.0, thetaz=0.0):
    """
  NAME:

     rotate
     --> This function is a work in progress. It is meant to eventually be a function for rotating streams to phi1/phi2 coordinates

  PURPOSE:

     Rotate StarCluster

  INPUT:

     x,y,z - coordinates of stars

     thetax,thetay,thetaz - angle about corresponding axis to rotate coordinates

  OUTPUT:

     rotated values of x,y,z

  HISTORY:

     2019 - Written - Webb (UofT)

  """
    # rotate about x axis:
    y2 = y * np.cos(thetax) - z * np.sin(thetax)
    z2 = -y * np.sin(thetax) + z * np.cos(thetax)

    y = y2
    z = z2

    # rotate about y axis:
    x2 = x * np.cos(thetay) - z * np.sin(thetay)
    z2 = -x * np.sin(thetay) + z * np.cos(thetay)

    x = x2
    z = z2

    # rotate about z axis:
    x2 = x * np.cos(thetaz) - y * np.sin(thetaz)
    y2 = -x * np.sin(thetaz) + y * np.cos(thetaz)

    x = x2
    y = y2

    return np.array(x), np.array(y), np.array(z)


def area_enclosed(
    x, y, thresh=None, nrand=1000, method="sphere", full=False, plot=False
):
    """
    NAME:

     area_enclosed

    PURPOSE:

     calculate area enclosed by a collection of points by finding what fraction of a random distribution overlaps with points

    INPUT:

     x,y - x and y coorindates
     
     thresh - threshold for overlap between random distribution and points (Default: None - use maximum nearest neighbour)

     nrand - number of random points to be generated in uniform distribution (Default: 1000)

     method - generate spherical or rectangular distribution of random points (Default: sphere)
     
     plot - plot overlap (Default: False)

    OUTPUT:

     area

    HISTORY:

     2019 - Written - Webb (UofT)

    """

    if thresh is None:
        z = np.zeros(len(x))
        rmin = minimum_distance(np.array([x, y, z]).T)
        thresh = np.amax(rmin)

    if method == "rectangle":
        xmin, xmax = np.amin(x), np.amax(x)
        ymin, ymax = np.amin(y), np.amax(y)

        dx = np.fabs(xmax - xmin)
        dy = np.fabs(ymax - ymin)

        atot = dx * dy

        # Generate random particles in area
        xrand = np.random.rand(nrand) * dx + xmin
        yrand = np.random.rand(nrand) * dy + ymin

    elif method == "sphere":
        r = np.sqrt(x ** 2.0 + y ** 2.0)
        rmin, rmax = np.amin(r), np.amax(r)

        phi = np.arctan2(y, x)
        phi[phi < 0] += 2.0 * np.pi
        phimin, phimax = np.amin(phi), np.amax(phi)

        if phimin < 0:
            phimin += 2.0 * np.pi

        dr = np.amax(r) - np.amin(r)

        atot = np.pi * (rmax ** 2.0 - rmin ** 2.0)
        rrand = np.random.rand(nrand) * dr + rmin
        phirand = np.random.rand(nrand) * np.fabs(phimax - phimin) + phimin
        xrand = rrand * np.cos(phirand)
        yrand = rrand * np.sin(phirand)

    drmin = np.array([])
    if full:
        xrfull = np.tile(xrand, len(x))
        yrfull = np.tile(yrand, len(x))

        xfull = np.repeat(x, nrand)
        yfull = np.repeat(y, nrand)

        dr = np.sqrt((xrfull - xfull) ** 2.0 + (yrfull - yfull) ** 2.0)
        dr = dr.reshape(len(x), nrand)
        drmin = np.amin(dr, 0)
    else:
        for i in range(0, nrand):
            dr = np.sqrt((xrand[i] - x) ** 2.0 + (yrand[i] - y) ** 2.0)
            drmin = np.append(drmin, np.amin(dr))

    indx = drmin < thresh

    area = atot * np.sum(indx) / float(nrand)

    if plot:
        plt.plot(x, y, ".")
        plt.plot(xrand[indx], yrand[indx], "o", alpha=0.1)
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    return area


@numba.njit
def minimum_distance(x):
    """
    NAME:

       minimum_distance

    PURPOSE:

       Find distance to each point's nearest neighbour

    INPUT:

       x=[x,y,z].T

    OUTPUT:

       distance to nearest neighbour

    HISTORY:

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
    """
    NAME:

       distance

    PURPOSE:

       Find distance between two points (made for use with numba)

    INPUT:

       x1=[x,y,z]
       x2=[x,y,z]

    OUTPUT:

       distance

    HISTORY:

       2019 - Written - Webb (UofT)
    """

    dx = x2[0] - x1[0]
    dy = x2[1] - x1[1]
    dz = x2[2] - x1[2]

    r = (dx * dx + dy * dy + dz * dz) ** 0.5

    return r
