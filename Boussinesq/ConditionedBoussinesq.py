#!/usr/bin/python

# Code to solve the conditionally averaged 1D Boussinesq equations:
# ddt(sigma) + ddx(sigma u1) = 0               1
# ddt(u1) + u1 ddx(u1) = -ddx(Phi)             2
# ddt(u2) + u2 ddx(u2) = -ddx(Phi)             3
# Substituting 2 and 3 into the total continuity eqn gives the Possion eqn:
# d2dx2(Phi) = -(sigma u1 ddx(u1) + (1-sigma) u2 ddx(u2))
# All advection terms use upwind differencing
# Periodic boundary conditions
# sigma and Phi are offset by half dx from u1 and u2 (C-grid staggering)

from __future__ import absolute_import, division, print_function
from pylab import *

# The linear algebra package for solving the matrix equation
import scipy.linalg as la

execfile("operators.py")
execfile("initialConditions.py")
execfile("plots.py")

def main():
    # Parameters
    nx = 40
    nt = 100
    dt = 0.01
    squareWaveMin = 0.4
    squareWaveMax = 0.6
    dx = 1./nx

    # Space
    x = arange(0.,1.,dx)

    # Initial conditions
    sigma = 0.1 + 0.8*squareWave(x, squareWaveMin, squareWaveMax)
    u1 = 0.1*sigma
    u2 = 1-sigma*u1/(1-sigma)
    Phi = zeros(nx)
    
    # Plot the initial conditions
    plotSolution(x, sigma, u1, u2, Phi, 0)

    # Fixed matrix to solve for the Poisson equation
    laplacianMatrix = zeros([nx,nx])
    for i in xrange(nx):
        laplacianMatrix[i,i] = -2/dx**2
        laplacianMatrix[i,(i-1)%nx] = 1/dx**2
        laplacianMatrix[i,(i+1)%nx] = 1/dx**2

    # Loop over time
    for it in xrange(int(nt)):
        time = (it+1)*dt
        print("time = ", time)
        
        # Calculate the RHS of the Poisson equation
        sigmaAtu = pAtu_upwind(sigma, u1)
        PoissonRHS = -dudx \
        (
            sigmaAtu*u1*ddxUp(u1, u1, dx)
          - (1-sigmaAtu)*u2*ddxUp(u2, u2, dx),
            dx
        )

        # Solve the Poisson equation
        Phi = la.solve(laplacianMatrix, PoissonRHS)

        # Update u1 and u2
        u1 -= dt*(u1*ddxUp(u1, u1, dx) - dpdx(Phi, dx))
        u2 -= dt*(u2*ddxUp(u2, u2, dx) - dpdx(Phi, dx))

        print("Maximum u1 Courant number = ", max(u1*dt/dx))
        print("Maximum u2 Courant number = ", max(u2*dt/dx))

        # Advect sigma
        sigma -= dt*dudx(sigmaAtu*u1, dx)

        print("sigma goes from ", min(sigma), " to ", max(sigma))
        
        # Plot the solution
        plotSolution(x, sigma, u1, u2, Phi, 0) # time)

main()

