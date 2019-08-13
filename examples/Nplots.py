import nbodypy as npy
import numpy as np


explot=False
rtexplot=False

dvplot=False
rtdvplot=False

dvmovie=False
rtdvmovie=True


if explot:
    exdata=np.loadtxt('extrct.npy')
    npy.explot(exdata)

if rtexplot:
    exdata=np.loadtxt('rtextrct.npy')
    npy.explot(exdata,'rt')

if dvplot:
    dvdata=np.loadtxt('dvprof.npy')
    npy.dvplot(dvdata,nsnap=0)
    npy.dvplot(dvdata,tsnap=12000.0)

if rtdvplot:
    dvdata=np.loadtxt('rtdvprof.npy')
    npy.dvplot(dvdata,'rt',nsnap=0)
    npy.dvplot(dvdata,'rt',tsnap=12000.0)

if dvmovie:
    dvdata=np.loadtxt('dvprof.npy')
    npy.dvplot(dvdata,animate=True)
    
if rtdvmovie:
    dvdata=np.loadtxt('rtdvprof.npy')
    npy.dvplot(dvdata,'rt',animate=True)
