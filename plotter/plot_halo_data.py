import numpy as np
from pylab import *
import matplotlib.pyplot as plt

class make_halo_data:

    def __init__(self, mass, conc, vmax, rs, counter):

        self.mass = mass
        self.concentration = conc
        self.vmax = vmax
        self.rs = rs
        self.counter = counter

    def add_data(self, mass, conc, vmax, rs):

        self.mass = np.append(self.mass,mass)
        self.concentration = np.append(self.concentration,conc)
        self.vmax = np.append(self.vmax,vmax)
        self.rs = np.append(self.rs, rs)

def c_M_relation(log10_M0):
    """
    Concentration-mass relation from Correa et al. (2015).
    This relation is most suitable for Planck cosmology.
    """
    z = 0
    # Best-fit params:
    alpha = 1.7543 - 0.2766 * (1. + z) + 0.02039 * (1. + z) ** 2
    beta = 0.2753 + 0.0035 * (1. + z) - 0.3038 * (1. + z) ** 0.0269
    gamma = -0.01537 + 0.02102 * (1. + z) ** (-0.1475)

    log_10_c200 = alpha + beta * log10_M0 * (1. + gamma * log10_M0 ** 2)
    c200 = 10 ** log_10_c200
    return c200

def median_relations(x, y, xrange):

    xvalues = np.ones(len(xrange) - 1) * (-10)
    yvalues = np.zeros(len(xrange) - 1)
    yvalues_err_down = np.zeros(len(xrange) - 1)
    yvalues_err_up = np.zeros(len(xrange) - 1)

    perc = [16, 84]

    for i in range(0, len(xrange) - 2):
        mask = (x > xrange[i]) & (x < xrange[i + 1])
        if len(x[mask]) > 4:
            xvalues[i] = np.median(x[mask])
            yvalues[i] = np.median(y[mask])
            yvalues_err_down[i], yvalues_err_up[i] = np.transpose(np.percentile(y[mask], perc))

    mask = xvalues>-10
    xvalues = xvalues[mask]
    yvalues = yvalues[mask]
    yvalues_err_down = yvalues_err_down[mask]
    yvalues_err_up = yvalues_err_up[mask]

    return xvalues, yvalues, yvalues_err_down, yvalues_err_up


def store_halo_data(sim_info, halo_data):

    mass = sim_info.halo_data.log10_halo_mass
    conc = sim_info.halo_data.concentration
    vmax = sim_info.halo_data.vmax
    rs = sim_info.halo_data.scale_radius * 1e3 #kpc
    counter = len(mass)
    if halo_data == None:
        halo_data = make_halo_data(mass, conc, vmax, rs, counter)
    else:
        halo_data.add_data(mass, conc, vmax, rs)

    return halo_data



def plot_relations(halo_data, sim_info, name_list):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.45,
        "figure.subplot.hspace": 0.35,
        "lines.markersize": 2,
        "lines.linewidth": 2,
        "axes.axisbelow": False,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    counter = 0
    for name in name_list:

        mass = 10**halo_data.mass[counter:counter+halo_data.counter]
        c200 = halo_data.concentration[counter:counter+halo_data.counter]

        xrange = np.arange(6,15,0.2)
        xrange = 10**xrange
        xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass, c200, xrange)

        plt.plot(mass, c200, 'o', color=color[i],alpha=0.2)
        plt.plot(xvalues, yvalues, '-', color='white',zorder=10)
        plt.plot(xvalues, yvalues, '-', lw=1.5, color=color[i], label=name,zorder=10)
        plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color=color[i])
        i += 1
        counter += halo_data.counter

    plt.plot(xrange, c_M_relation(xrange), lw=1, color='black',label='Correa et al. (2015)')

    plt.axis([1e9,1e15,1,50])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("c$_{200}$")
    plt.xscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/cM_relation"+ name_list[0] + "_" + name_list[1] +".png", dpi=200)
    plt.close()

    ###########

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    counter = 0
    for name in name_list:

        mass = 10**halo_data.mass[counter:counter+halo_data.counter]
        rs = halo_data.rs[counter:counter+halo_data.counter]

        xrange = np.arange(6,15,0.2)
        xrange = 10**xrange
        xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass, rs, xrange)

        plt.plot(mass, rs, 'o', color=color[i],alpha=0.2)
        plt.plot(xvalues, yvalues, '-', color='white',zorder=10)
        plt.plot(xvalues, yvalues, '-', lw=1.5, color=color[i], label=name,zorder=10)
        plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color=color[i])
        i += 1
        counter += halo_data.counter

    c200 = c_M_relation(xrange)
    z = sim_info.z
    Om = sim_info.Omega_m
    Ol = sim_info.Omega_l
    rho_crit = sim_info.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**xrange / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc
    rs = R200/c200
    plt.plot(xrange, rs, lw=1, color='black',label='Correa et al. (2015)')


    plt.axis([1e9,1e15,1,60])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("Scale radius [kpc]")
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/rsM_relation"+ name_list[0] + "_" + name_list[1] +".png", dpi=200)
    plt.close()


    ###########
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    counter = 0
    for name in name_list:

        mass = 10**halo_data.mass[counter:counter+halo_data.counter]
        vmax = halo_data.vmax[counter:counter+halo_data.counter]

        xrange = np.arange(6, 15, 0.2)
        xrange = 10 ** xrange
        xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass, vmax, xrange)

        plt.plot(mass, vmax, 'o', color=color[i], alpha=0.2)
        plt.plot(xvalues, yvalues, '-', color='white',zorder=10)
        plt.plot(xvalues, yvalues, '-', lw=1.5, color=color[i], label=name,zorder=10)
        plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color=color[i])
        i += 1
        counter += halo_data.counter

    plt.axis([1e9, 1e15, 1, 1e3])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("V$_{\mathrm{max}}$ [km/s]")
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/VmaxM_relation"+ name_list[0] + "_" + name_list[1] +".png", dpi=200)
    plt.close()