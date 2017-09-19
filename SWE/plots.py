# Plot the solution of the conditionally averaged Boussinesq equations

def plotSolution(xp, p, q, xu, u, v, time):
    "Plots p,q,u,v as a function of xp and xu and labels with time"
    
    figure(1)
    font = {'size'   : 18}
    rc('font', **font)
    clf()
#    ion()
    plot(xp, p/(p+q), label='$\sigma$', color='black')
    plot(xp, p+q, label='$h$', color='red')

    plot(xu, u, label='$u_1$', color='cyan')
    plot(xu, v, label='$u_2$', color='magenta')

    axhline(0, linestyle=':', color='black')
    axhline(1, linestyle=':', color='black')
    legend(loc=2,ncol=4)
    xlabel('$x$')
    ylim([-1.5,2])
    title("Time "+str(time))
    savefig("plots/pquv_nx_" + str(len(xu)) + "_time_"+str(time)+".pdf")

def plotEnergy(energy, KE1, KE2, dt):

    figure(2)
    font = {'size'   : 18}
    rc('font', **font)
    clf()
#    ion()
    nt = len(energy)
    t = linspace(0, (nt-1)*dt, nt)
    fig,ax1 = subplots()
    ax1.plot(t, KE1, label='KE$_1$', color='cyan')
    ax1.plot(t, KE2, label='KE$_2$', color='magenta')
    xlabel('time')
    ax2 = ax1.twinx()
    ax2.plot(t, (energy-energy[0])/energy[0], label='Normalised energy change',
             color='black')
    #legend(loc=2,ncol=4)
    savefig("plots/energy.pdf")

