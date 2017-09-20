#!/usr/bin/python

# Code to solve the conditionally averaged 1D SWE:
#
# ddt(u1) + 0.5 ddx(u1**2) = -ddx(h)
# ddt(u2) + 0.5 ddx(u2**2) = -ddx(h)
# ddt(h1) + ddx(h1 u1) = 0
# ddt(h2) + ddx(h2 u2) = 0
# h = h1+h2
# sigma = h1/h
#
# All advection terms use upwind differencing
# Periodic boundary conditions
# h1 and h2 are offset by half dx from u1 and u2 (C-grid staggering)

from __future__ import absolute_import, division, print_function
from pylab import *
import os
import sys

execfile(os.path.join(sys.path[0],"operators.py"))
execfile(os.path.join(sys.path[0],"initialConditions.py"))
execfile(os.path.join(sys.path[0],"plots.py"))

def main():
    # Parameters
    nx = 40
    nt = 200
    dt = 0.01
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    dx = 1./nx

    # Space for the u and p grids
    xu = arange(0.,1.,dx)
    xh = arange(0.5*dx, 1, dx)

    # Initial conditions
    h1 = 0.1 + 0.8*squareWave(xh, squareWaveMin, squareWaveMax)
    h2 = 1-h1
    u1 = np.ones(len(xu))
    u2 = u1.copy()
    
    h1old = h1.copy()
    h2old = h2.copy()
    u1old = u1.copy()
    u2old = u2.copy()
    
    # Plot the initial conditions
    plotSolution(xh, h1, h2, xu, u1, u2, 0)
    
    energy = zeros(nt+1)
    KE1 = zeros(nt+1)
    KE2 = zeros(nt+1)
    KE1[0] = KE(h1, u1, dx)
    KE2[0] = KE(h2, u2, dx)
    energy[0] = KE1[0] + KE2[0] + PE(h1, h2, dx)

    # Loop over time
    for it in xrange(int(nt)):
        print("time step ", it+1)
        
        # Outer loop
        for iter in range(2):
        
            h1Atu = hAtu_upwind(h1,u1)
            h2Atu = hAtu_upwind(h2,u2)
#            h1Atu = hAtu_centred(h1)
#            h2Atu = hAtu_centred(h2)
        
            # Advect h1 and h2
            h1 = h1old - dt*dudx(h1Atu*u1, dx)
            h2 = h2old - dt*dudx(h2Atu*u2, dx)
        
            # Update u1 and u2 with upwind advection
            gradh = dhdx(h1+h2, dx)
            for innerIter in range(2):
                u1 = u1old - dt*(u1*ddxUp(u1, u1, dx) + gradh)
                u2 = u2old - dt*(u2*ddxUp(u2, u2, dx) + gradh)
            
#            # Update u1 and u2 vector invariant form
#            u1 = u1old - dt*dhdx(0.5*uAth(u1)**2 + h1+h2, dx)
#            u2 = u2old - dt*dhdx(0.5*uAth(u2)**2 + h1+h2, dx)
        
        h1old = h1.copy()
        h2old = h2.copy()
        u1old = u1.copy()
        u2old = u2.copy()

        # Energy diagnostics
        KE1[it+1] = KE(h1, u1, dx)
        KE2[it+1] = KE(h2, u2, dx)
        energy[it+1] = KE1[it+1] + KE2[it+1] + PE(h1, h2, dx)
        
        print("Minimum h1 = ", min(h1))
        print("Minimum h2 = ", min(h2))
        print("Maximum u1 Courant number = ", max(abs(u1)*dt/dx))
        print("Maximum u2 Courant number = ", max(abs(u2)*dt/dx))
        print("Normalised change in energy = ", (energy[it+1]-energy[0])/energy[0])

        # Plot the solution
        plotSolution(xh, h1, h2, xu, u1, u2, it+1)
        if it%10 == 0:
            plotEnergy(energy, KE1, KE2, dt)
    
    plotEnergy(energy, KE1, KE2, dt)
main()

