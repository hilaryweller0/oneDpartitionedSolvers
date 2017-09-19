from __future__ import absolute_import, division, print_function

from pylab import *

# Fucntion for calculating ddx using 1st order upwind differencing
# and centred, staggered differencing

def ddxUp(u, dir, dx):
    "Calculate ddx(u) using upwind differencing with direction defined in dir."
    "dir is co-located with u"
    
    dudx = zeros_like(u)
    nx = len(u)
    
    # Boundary values (dependent on dir)
    if dir[0] >= 0:
        dudx[0] = (u[0] - u[-1])/dx
    else:
        dudx[0] = (u[1] - u[0])/dx
    
    if dir[-1] >=0:
        dudx[-1] = (u[-1] - u[-2])/dx
    else:
        dudx[-1] = (u[0] - u[-1])/dx

    # Internal values (dependent of dir)
    for i in range(1,nx-1):
        if dir[i] >= 0:
            dudx[i] = (u[i] - u[i-1])/dx
        else:
            dudx[i] = (u[i+1] - u[i])/dx

    return dudx

def dudx(u, dx):
    "Calculate ddx(u) using centred, compact diffrencing."
    "Result is at mid-points"
    
    return (roll(u,-1) - u)/dx

def dpdx(p, dx):
    "Calculate ddx(p) using centred, compact diffrencing."
    "Result is at mid-points"
    
    return (p - roll(p,1))/dx

def pAtu_upwind(p, u):
    "p is assumed to be offset by 0.5dx from u. pAtu_upwind is p on the upwind"
    "side"
    
    pAtU = zeros_like(u)
    nx = len(u)
    for i in range(1,nx-1):
        if u[i] >= 0:
            pAtU[i] = p[i-1]
        else:
            pAtU[i] = p[i]
    
    if u[0] >= 0:
        pAtU[0] = p[-1]
    else:
        pAtU[0] = p[0]
    
    if u[-1] >= 0:
        pAtU[-1] = p[-2]
    else:
        pAtU[-1] = p[-1]
    
    return pAtU

