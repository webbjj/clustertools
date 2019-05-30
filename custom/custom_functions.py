import numpy as np

def def_dc(m, x, v=None, r_min=0.0045):
    """
    Calculates the center of density. Translated to Python from phiGRAPE.

    Parameters
    ----------
    m, x, v : array_like
        Particle parameters. v is optional.
    r_min : scalar
        For star clusters, should be 0.1 pc in N-body units.

    Returns
    -------
    xdc, vdc : ndarrays
    """
    calc_vdc = not v is None
    x_ = x.copy()
    xdc = np.zeros(3)
    if calc_vdc:
        v_ = v.copy()
        vdc = np.zeros(3)
    else:
      v_ = None
    r_lim = np.sqrt(np.max(np.sum(x**2, axis=1)))
    num_iter = 0

    while (r_lim > r_min) and (num_iter < 100):
        ncm, mcm, xcm, vcm = cenmas_lim(m, x_, v_, r_lim)
        if((mcm > 0) and (ncm > 100)):
            x_ -= xcm
            xdc += xcm
            if calc_vdc:
              v_ -= vcm
              vdc += vcm
        else:
            break
        r_lim *= 0.8
        num_iter += 1
    if calc_vdc:
        return xdc, vdc
    else:
        return xdc

def cenmas_lim(m, x, v, r_lim):
    r2 = np.sum(x**2, axis=1)
    cond = r2 < r_lim**2
    ncm = np.sum(cond)
    mcm = np.sum(m[cond])
    if mcm == 0: return ncm, 0., np.zeros(3), np.zeros(3)
    xcm = np.sum(m[cond,None] * x[cond,:], axis=0) / mcm
    if not v is None: vcm = np.sum(m[cond,None] * v[cond,:], axis=0) / mcm
    else: vcm = None
    return ncm, mcm, xcm, vcm
