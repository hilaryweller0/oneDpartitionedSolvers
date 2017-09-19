# Plot the solution of the conditionally averaged Boussinesq equations

def plotSolution(x, sigma, u1, u2, Phi, time):
    "Plots sigma, u1 and u2 as a function of x and labels with time"
    
    font = {'size'   : 20}
    rc('font', **font)
    clf()
    ion()
    plot(x, sigma, label='$\sigma$', color='black')
    plot(x, u1, label='$u_1$', color='blue')
    plot(x, u2, label='$u_2$', color='red')
    plot(x, Phi, label='$\Phi$', color='grey')
    axhline(0, linestyle=':', color='black')
    axhline(1, linestyle=':', color='black')
    legend(bbox_to_anchor=(1.1, 1))
    xlabel('$x$')
    ylim([-0.2,1.2])
    title("Time "+str(time))
    savefig("plots/simga_u1_u2_nx_" + str(len(x)) + "_time_"+str(time)+".pdf")

