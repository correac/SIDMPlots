import numpy as np
from pylab import *
import matplotlib.pyplot as plt


def read_fits(file):
    data = np.loadtxt(file)
    r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso = data[0,0:6]
    r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse = data[1,0:6]
    return r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso,\
           r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse



if __name__ == "__main__":

    DSphs = ['Carina','Draco','Fornax','LeoI','Sculptor','Sextans','UMi']
    # Log10 of M200c masses taken from Read et al. (2019)
    M200c = np.array([0.8,1.8,21.9,5.6,5.7,2.0,2.8]) * 1e9
    error_M200c = np.array([0.3,0.7,7.4,2.2,2.3,0.8,1.1]) * 1e9

    # Plot parameters
    params = {
        "font.size": 10,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.25,
        "figure.subplot.hspace": 0.25,
        "lines.markersize": 6,
        "lines.linewidth": 1.5,
        "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    grid(True)

    i = 0
    color=['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'black']
    for dwarf in DSphs:

        file = "DSh_"+dwarf+"_fit.txt"
        r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso,\
            r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse= read_fits(file)

        yerr = np.array([(10**rho0_iso-10**low_rho0_iso, 10**high_rho0_iso-10**rho0_iso)]).T
        xerr = error_M200c[i]
        plt.errorbar([M200c[i]], [10**rho0_iso], xerr=xerr, yerr=yerr, fmt='',color=color[i])
        plt.plot([M200c[i]], [10**rho0_iso],'o',ms=4,color='white')
        plt.plot([M200c[i]], [10**rho0_iso],'o',ms=3,color=color[i],label=dwarf)

        yerr = np.array([(10**rho0_pse-10**low_rho0_pse, 10**high_rho0_pse-10**rho0_pse)]).T
        plt.errorbar([M200c[i]], [10**rho0_pse], xerr=xerr, yerr=yerr, fmt='',color=color[i])
        plt.plot([M200c[i]], [10**rho0_pse],'v',ms=4,color='white')
        plt.plot([M200c[i]], [10**rho0_pse],'v',ms=3,color=color[i])

        i+=1

    yscale('log')
    xscale('log')
    xlabel(r'$M_{200}$ [M$_{\odot}$]')
    ylabel(r'$\rho_{0}$ [M$_{\odot}$/kpc$^{3}$]')
    axis([5e8, 5e11, 1e7, 1e10])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig("rho0_mass_relation.png", dpi=200)
    plt.close()


    figure()
    ax = plt.subplot(1, 1, 1)
    grid(True)

    i = 0
    color=['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'black']
    for dwarf in DSphs:

        file = "DSh_"+dwarf+"_fit.txt"
        r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso,\
            r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse= read_fits(file)

        xerr = error_M200c[i]
        yerr = np.array([(r0_iso-low_r0_iso, high_r0_iso-r0_iso)]).T
        plt.errorbar([M200c[i]], [r0_iso], xerr=xerr, yerr=yerr, fmt='',color=color[i])
        plt.plot([M200c[i]], [r0_iso],'o',ms=4,color='white')
        plt.plot([M200c[i]], [r0_iso],'o',ms=3,color=color[i],label=dwarf)

        yerr = np.array([(r0_pse-low_r0_pse, high_r0_pse-r0_pse)]).T
        plt.errorbar([M200c[i]], [r0_pse], xerr=xerr, yerr=yerr, fmt='',color=color[i])
        plt.plot([M200c[i]], [r0_pse],'v',ms=4,color='white')
        plt.plot([M200c[i]], [r0_pse],'v',ms=3,color=color[i])

        i+=1

    yscale('log')
    xscale('log')
    xlabel(r'$M_{200}$ [M$_{\odot}$]')
    ylabel(r'$r_{0}$ [kpc]')
    axis([5e8, 5e11, 1e-2, 1e0])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig("r0_mass_relation.png", dpi=200)
    plt.close()