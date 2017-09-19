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

def ddx(phi, dx):
    "Calculates ddx(phi) using centred, leapfrog differencing."
    "Co-located results"
    
    return 0.5*(roll(phi,-1) - roll(phi,1))/dx

def dudx(u, dx):
    "Calculate ddx(u) using centred, compact diffrencing."
    "Result is at mid-points"
    
    return (roll(u,-1) - u)/dx

def dhdx(h, dx):
    "Calculate ddx(h) using centred, compact diffrencing."
    "Result is at mid-points"
    
    return (h - roll(h,1))/dx

def hAtu_upwind(h, u):
    "h is assumed to be offset by 0.5dx from u. pAtu_upwind is p on the upwind"
    "side"
    
    hAtU = zeros_like(u)
    nx = len(u)
    for i in range(1,nx-1):
        if u[i] >= 0:
            hAtU[i] = h[i-1]
        else:
            hAtU[i] = h[i]
    
    if u[0] >= 0:
        hAtU[0] = h[-1]
    else:
        hAtU[0] = h[0]
    
    if u[-1] >= 0:
        hAtU[-1] = h[-2]
    else:
        hAtU[-1] = h[-1]
    
    return hAtU

def hAtu_centred(h):
    return 0.5*(h + roll(h,1))

def uAth(u):
    return 0.5*(u + roll(u,-1))

def KE(h, u, dx):
    return 0.5*sum(h*uAth(u)**2)*dx

def PE(h1, h2, dx):
    return 0.5*sum((h1+h2)**2)*dx

