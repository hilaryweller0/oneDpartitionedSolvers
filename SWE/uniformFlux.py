#!/usr/bin/python

# Code to solve the conditionally averaged 1D SWE in flux form:
#
# ddt(hiui) + ddx(hi ui ui) = -hi ddx(h)
# ddt(hi) + ddx(hi ui) = 0
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

SMALL = 1e-9

def main():
    # Parameters
    nx = 40
    nt = 200
    dt = 0.005
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
    u2 = np.ones(len(xu))
    h1u1 = hAtu_centred(h1)*u1
    h2u2 = hAtu_centred(h2)*u2
    
    h1old = h1.copy()
    h2old = h2.copy()
    h1u1old = h1u1.copy()
    h2u2old = h2u2.copy()
    
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
        for iter in range(1):
        
            h1Atu = hAtu_upwind(h1,u1)
            h2Atu = hAtu_upwind(h2,u2)
        
            # Advect h1 and h2
            h1 = h1old - dt*dudx(h1Atu*u1, dx)
            h2 = h2old - dt*dudx(h2Atu*u2, dx)
        
            # Update h1u1 and h2u2 with upwind advection
            gradh = dhdx(h1+h2, dx)
            for innerIter in range(1):
                h1u1 = h1u1old - dt*(ddxUp(h1u1*u1, u1, dx) + h1*gradh)
                h2u2 = h2u2old - dt*(ddxUp(h2u2*u2, u2, dx) + h2*gradh)
                u1 = h1u1/hAtu_centred(h1+SMALL)
                u2 = h2u2/hAtu_centred(h2+SMALL)
        
        h1old = h1.copy()
        h2old = h2.copy()
        h1u1old = h1u1.copy()
        h2u2old = h2u2.copy()

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

