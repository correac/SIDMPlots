from pylab import *
import matplotlib.pyplot as plt
from scipy.special import spence
import numpy as np


def sigma_1D(x, log10_M, a, siminfo):
    z = siminfo.z
    Om = siminfo.Omega_m
    Ol = siminfo.Omega_l
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    G = 4.3e-6  # kpc km^2 Msun^-1 s^-2
    ff = np.zeros(len(x))
    for i in range(0, len(x)):
        s = x[i] / R200
        c = R200 / a
        u = 1. + c * s
        gc = 1. / (np.log(1. + c) - c / (1. + c))
        f = 0.5 * c ** 2 * gc * s * u ** 2 * G * 10**log10_M / R200
        u = 1. + c * s
        Li = spence(u)
        f *= (np.pi ** 2 - np.log(c * s) - 1. / (c * s) - 1. / u ** 2 - 6. / u + (
                    1. + 1. / (c ** 2 * s ** 2) - 4. / (c * s) - 2. / u) * np.log(u) + 3 * (np.log(u)) ** 2 + 6. * Li)
        ff[i] = f
    return np.sqrt(ff)


def calc_density(x, log10_M, siminfo):
    z = siminfo.z
    Om = siminfo.Omega_m
    Ol = 1. - Om
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    # NFW profile #
    c = c_M_relation(log10_M)
    a = R200 / c
    delta = 200. / 3.
    delta *= c ** 3 / (np.log(1. + c) - c / (1. + c))
    f = np.zeros(len(x))
    for i in range(0, len(x)): f[i] = rho_crit * delta * 1e-9 / ((x[i] / a) * (1. + x[i] / a) ** 2)
    return f

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



def plot_profiles(profile_data, sim_info, output_name_list):

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (6, 3),
        "figure.subplot.left": 0.12,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.4,
        "figure.subplot.hspace": 0.4,
        "lines.markersize": 6,
        "lines.linewidth": 1.0,
        "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    ######################
    figure()

    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    colors = ['tab:blue','tab:orange']
    counter = 0
    k = 0
    for name in output_name_list:
        centers = profile_data.centers
        density = profile_data.density[counter:counter+len(centers)]
        sig_density = profile_data.sig_density[counter:counter + len(centers)]

        plt.plot(centers, density, lw=2, color=colors[k],label=name)
        plt.fill_between(centers, density - sig_density / 2,
                         density + sig_density / 2, alpha=0.4,
                         color=colors[k])

        k += 1
        counter += len(centers)

    log10_M200 = profile_data.log10_M200
    NFWrho = calc_density(centers, log10_M200, sim_info)  # Msun/kpc^3
    plt.plot(centers, NFWrho, lw=1, color='black', label="NFW profile")

    if profile_data.structure_type == 10:
        plt.text(0.05, 0.05, 'Central haloes of $10^{%0.1f' % log10_M200 + '}$M$_{\odot}$', transform=ax.transAxes)
    else:
        plt.text(0.05, 0.05, 'Satellite haloes of $10^{%0.1f' % log10_M200 + '}$M$_{\odot}$', transform=ax.transAxes)

    plt.plot([sim_info.softening, sim_info.softening], [1e3, 1e9], '--', lw=1, color='grey')

    plt.ylim(1e4, 1e9)
    plt.xlim(0.1, 1e3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Density [M$_{\odot}$/kpc$^{3}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

    ######################
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    counter = 0
    k = 0
    for name in output_name_list:
        velocity = profile_data.velocity[counter:counter+len(centers)]
        sig_velocity = profile_data.sig_velocity[counter:counter + len(centers)]

        plt.plot(centers, velocity, lw=2, color=colors[k],label=name)
        plt.fill_between(centers, velocity - sig_velocity / 2,
                         velocity + sig_velocity / 2,
                         alpha=0.4, color=colors[k])

        k += 1
        counter += len(centers)


    rs = profile_data.rs
    NFWsig1D = sigma_1D(centers, log10_M200, rs, sim_info)  # km/s
    plt.plot(centers, NFWsig1D, lw=1, color='black')

    # plt.ylim(0, 80)
    plt.xlim(0.1, 1e3)
    plt.xscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Velocity dispersion [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    if profile_data.structure_type == 10:
        plt.savefig(f"{sim_info.output_path}/Density_%0.1f"%log10_M200+"_halos_"+
                    output_name_list[0]+"_"+
                    output_name_list[1]+".png", dpi=200)
    else:
        plt.savefig(f"{sim_info.output_path}/Density_%0.1f"%log10_M200+"_subhalos_"+
                    output_name_list[0]+"_"+
                    output_name_list[1]+".png", dpi=200)
    plt.close()

    ######################

    ######################
    figure()

    ax = plt.subplot(1, 2, 1)
    plt.grid("True")

    colors = ['tab:blue', 'tab:orange']
    counter = 0
    k = 0
    for name in output_name_list:
        centers = profile_data.centers
        density = profile_data.density[counter:counter + len(centers)]
        sig_density = profile_data.sig_density[counter:counter + len(centers)]

        plt.plot(centers, density, lw=2, color=colors[k], label=name)
        plt.fill_between(centers, density - sig_density / 2,
                         density + sig_density / 2, alpha=0.4,
                         color=colors[k])

        k += 1
        counter += len(centers)

    log10_M200 = profile_data.log10_M200
    NFWrho = calc_density(centers, log10_M200, sim_info)  # Msun/kpc^3
    plt.plot(centers, NFWrho, lw=1, color='black', label="NFW profile")

    if profile_data.structure_type == 10:
        plt.text(0.05, 0.05, 'Central haloes of $10^{%0.1f' % log10_M200 + '}$M$_{\odot}$', transform=ax.transAxes)
    else:
        plt.text(0.05, 0.05, 'Satellite haloes of $10^{%0.1f' % log10_M200 + '}$M$_{\odot}$', transform=ax.transAxes)

    plt.plot([sim_info.softening, sim_info.softening], [1e3, 1e9], '--', lw=1, color='grey')

    plt.ylim(1e4, 1e9)
    plt.xlim(0.1, 1e3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Density [M$_{\odot}$/kpc$^{3}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

    ######################
    ax = plt.subplot(1, 2, 2)
    plt.grid("True")

    counter = 0
    k = 0
    for name in output_name_list:
        sigma = profile_data.sigma[counter:counter + len(centers)]
        sig_sigma = profile_data.sig_sigma[counter:counter + len(centers)]

        plt.plot(centers, sigma, lw=2, color=colors[k], label=name)
        plt.fill_between(centers, sigma - sig_sigma / 2,
                         sigma + sig_sigma / 2,
                         alpha=0.4, color=colors[k])

        k += 1
        counter += len(centers)

    # plt.ylim(0, 80)
    plt.xlim(0.1, 1e3)
    plt.xscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Cross section [cm$^{2}$/g]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    if profile_data.structure_type == 10:
        plt.savefig(f"{sim_info.output_path}/Cross_section_profile_%0.1f" % log10_M200 + "_halos_" +
                    output_name_list[0] + "_" +
                    output_name_list[1] + ".png", dpi=200)
    else:
        plt.savefig(f"{sim_info.output_path}/Cross_section_profile_%0.1f" % log10_M200 + "_subhalos_" +
                    output_name_list[0] + "_" +
                    output_name_list[1] + ".png", dpi=200)
    plt.close()