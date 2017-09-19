#!/usr/bin/python

# Code to solve the 1D SWE:
#
# ddt(u) + 0.5 ddx(u**2) = -ddx(h)
# ddt(h) + ddx(h u) = 0
#
# Periodic boundary conditions
# h is offset by half dx from u (C-grid staggering)

from __future__ import absolute_import, division, print_function
from pylab import *

execfile("operators.py")
execfile("initialConditions.py")
execfile("plots.py")

def main():
    # Parameters
    nx = 40
    nt = 10
    dt = 0.01
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    dx = 1./nx

    # Space for the u and p grids
    xu = arange(0.,1.,dx)
    xh = arange(0.5*dx, 1, dx)

    # Initial conditions
    h = 0.2 + 0.7*squareWave(xh, squareWaveMin, squareWaveMax)
    u = 0.5*(1+sin(2*pi*xu))
    
    hold = h.copy()
    uold = u.copy()
    
    # Plot the initial conditions
    plotSolution(xh, h, 0*h, xu, u, 0*u, 0)
    
    energy = zeros(nt+1)
    KE1 = zeros(nt+1)
    KE1[0] = KE(h, u, dx)
    energy[0] = KE1[0] + PE(h, 0*h, dx)

    # Loop over time
    for it in xrange(int(nt)):
        print("time step ", it+1)
        
        # Outer loop
        for iter in range(2):
        
            hAtu = hAtu_centred(h)
        
            # Advect h
            h = hold - dt*dudx(hAtu*u, dx)
        
            # Update u vector invariant form
#            u = uold - dt*dhdx(0.5*uAth(u)**2 + h, dx)
        
            # Update u advective form
            u = uold - dt*(u*ddxUp(u,u,dx) + dhdx(h,dx))

        hold = h.copy()
        uold = u.copy()

        # Energy diagnostics
        KE1[it+1] = KE(h, u, dx)
        energy[it+1] = KE1[it+1] + PE(h, 0*h, dx)
        
        print("Minimum h = ", min(h))
        print("Maximum u Courant number = ", max(abs(u)*dt/dx))
        print("Normalised energy change = ", (energy[it+1]-energy[0])/energy[0])

        # Plot the solution
        plotSolution(xh, h, 0*h, xu, u, 0*u, it+1)
    
    plotEnergy(energy, KE1, 0*KE1, dt)
main()

