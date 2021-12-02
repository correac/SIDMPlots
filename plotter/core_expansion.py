from pylab import *
import matplotlib.pyplot as plt
import numpy as np


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


def plot_t0(sim_info):

    M = np.arange(8, 16, 0.2)

    z = sim_info.z
    Om = sim_info.Omega_m
    Ol = sim_info.Omega_l
    rho_crit = sim_info.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    c200 = c_M_relation(M)
    rs = R200 / c200 #kpc

    sigma = 1.
    t0_sigma_1 = 14.04 * sigma**(-1.) * (rs/10)**(7./2.) * (10**M/1e10)**(-3./2.)

    sigma = 10.
    t0_sigma_10 = 14.04 * sigma**(-1.) * (rs/10)**(7./2.) * (10**M/1e10)**(-3./2.)

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
    }
    rcParams.update(params)
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(10**M, 25 * t0_sigma_1, lw=1, color='tab:blue',label=r'$\sigma=1$ cm$^{2}$/g')
    plt.plot(10**M, 25 * t0_sigma_10, lw=1, color='tab:orange',label=r'$\sigma=10$ cm$^{2}$/g')

    # plt.ylim(0, 80)
    plt.xlim(1e8, 1e16)
    plt.xscale("log")
    plt.xlabel("M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Core expansion [Gyr]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

    plt.savefig(f"{sim_info.output_path}/core_expansion_time.png", dpi=200)
    plt.close()

