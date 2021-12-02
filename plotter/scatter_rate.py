import h5py
import numpy as np
import glob
import os
from pylab import *
import matplotlib.pyplot as plt
from scipy.special import spence
from colossus.cosmology import cosmology
cosmology.setCosmology('planck13');
from colossus.lss import mass_function

class make_cosmic_scatter_rate:

    def __init__(self, redshift, scatter_rate):
        self.redshift = redshift
        self.scatter_rate = scatter_rate

    def add_data(self, scatter_rate):
        self.scatter_rate = np.append(self.scatter_rate, scatter_rate)

def calc_density(x, log10_M, a, siminfo):
    z = siminfo.z
    Om = siminfo.Omega_m
    Ol = siminfo.Omega_l
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    # NFW profile #
    c = R200 / a
    delta = 200. / 3.
    delta *= c ** 3 / (np.log(1. + c) - c / (1. + c))
    f = np.zeros(len(x))
    for i in range(0, len(x)): f[i] = rho_crit * delta * 1e-9 / ((x[i] / a) * (1. + x[i] / a) ** 2)
    return f


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


def gamma(log10_M, siminfo):

    z = siminfo.z
    Om = siminfo.Omega_m
    Ol = siminfo.Omega_l
    sigma = siminfo.cross_section
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    # Duffy et al. (2008) relation
    c200 = 5.72 * (10**log10_M / 1e14) ** (-0.081) * (1. + z) ** (-0.71)

    a = R200 / c200  # r200/ r200/rs [kpc]

    radial_bins = np.arange(-5, np.log10(R200), 0.1)
    x = 10 ** radial_bins
    deltax = x[1:] - x[:1]
    deltax = np.append(deltax, deltax[-1])
    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21

    sig1D = sigma_1D(x, log10_M, a, siminfo)  # km/s
    sig1D[np.isnan(sig1D)] = 0
    sig1D *= 1e5  # cm/s

    rho = calc_density(x, log10_M, a, siminfo)  # Msun/kpc^3
    rho *= Msun_in_cgs / kpc_in_cgs ** 3  # g/cm^2

    f = sig1D * rho * (4. / np.sqrt(np.pi)) * sigma  # 1/s
    Gyr = 1e9 * 365.25 * 24 * 3600

    f *= Gyr  # Gyr^-1 particle^-1

    rho = calc_density(x, log10_M, a, siminfo)  # Msun/kpc^3
    vol = 4 * np.pi * x ** 2 * deltax * rho  # Msun
    Total = np.sum(f * vol) / 10**log10_M
    return Total


def analytic_scatter(siminfo, z, Mmin, Mmax):

    h = siminfo.h
    rho_crit = siminfo.rhocrit0

    M = np.arange(Mmin, Mmax, 0.25)
    M = 10 ** M
    mfunc = mass_function.massFunction(M / h, z, mdef='200c', model='tinker08', q_out='dndlnM')
    mfunc /= h**3
    mfunc *= M
    mfunc /= rho_crit
    f = mfunc

    G = np.zeros(len(M))
    for i in range(0, len(M)): G[i] = gamma(M[i], z, siminfo)  # [Gyr^-1 particle^-1]

    n = f * G
    n = np.sum(n * 0.25 / np.log(10.))
    return n  # [Gyr^-1 particle^-1]


def compute_scatter_rate(siminfo):

    snapshot = load(f"{siminfo.directory}/{siminfo.snapshot_base_name}"+"_%04i.hdf5"%0)
    old_time = snapshot.metadata.cosmology.hubble_time.value # Gyr

    nparts = siminfo.num_dm_particles
    scatter_rate = np.zeros(siminfo.n_snapshots+1)
    redshift = np.zeros(siminfo.n_snapshots+1)
    redshift[0] = snapshot.metadata.redshift
    previous_events = 0

    for i in range(1,siminfo.n_snapshots+1):

        snapshot_file = f"{siminfo.directory}/{siminfo.snapshot_base_name}" + "_%04i.hdf5" % i
        with h5py.File(snapshot_file, "r") as sim:
            n_sidm_events = np.sum(sim["/PartType1/SIDM_events"][:])

        snapshot = load(f"{siminfo.directory}/{siminfo.snapshot_base_name}" + "_%04i.hdf5" %i)
        redshift[i] = snapshot.metadata.redshift
        time = snapshot.metadata.cosmology.hubble_time.value # Gyr
        delta_time = time - old_time # Gyr
        old_time = time # Gyr

        n_collisions = n_sidm_events-previous_events
        previous_events = n_sidm_events
        scatter_rate[i] = n_collisions / ( delta_time * nparts )


    if scatter_rate == None:
        scatter_rate = make_cosmic_scatter_rate(redshift, scatter_rate)
    else:
        scatter_rate.add_data(scatter_rate)

    return scatter_rate


def plot_cosmic_scatter_rate(cosmic_scatter_rate, siminfo, output_name_list):

    #######################
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
        "lines.markersize": 6,
        "lines.linewidth": 2,
    }
    rcParams.update(params)

    #######################
    # Plot the interesting quantities
    figure()
    ax = plt.subplot(1, 1, 1)
    grid(True)

    z = np.arange(0, 20, 0.2)
    analytic_solution = np.zeros(len(z))
    for i in range(0, len(z)): analytic_solution[i] = analytic_scatter(siminfo, z[i], 9.0, 14)
    plot(1 + z, analytic_solution, '-', lw=1.5, color='tab:blue',
         label=r'Analytic (m$_{\mathrm{halo}}{>}10^{9}$M$_{\odot}$)')

    analytic_solution = np.zeros(len(z))
    for i in range(0, len(z)): analytic_solution[i] = analytic_scatter(siminfo, z[i], 10.0, 14)
    plot(1 + z, analytic_solution, '-', lw=1.5, color='tab:green',
         label=r'Analytic (m$_{\mathrm{halo}}{>}10^{10}$M$_{\odot}$)')

    analytic_solution = np.zeros(len(z))
    for i in range(0, len(z)): analytic_solution[i] = analytic_scatter(siminfo, z[i], 11.0, 14)
    plot(1 + z, analytic_solution, '-', lw=1.5, color='tab:red',
         label=r'Analytic (m$_{\mathrm{halo}}{>}10^{11}$M$_{\odot}$)')

    k = 0
    counter = 0
    color = ['black', 'tab:grey']
    for name in output_name_list:
        redshift = cosmic_scatter_rate.redshift
        scatter = cosmic_scatter_rate.scatter_rate[counter:counter+len(redshift)]

        plot(1 + redshift, scatter, '-', color=color[k], label=name)
        k += 1
        counter += len(redshift)


    xscale('log')
    xlabel('$1+z$')
    ylabel(r'$\Gamma(r)$ [Gyr$^{-1}$ particle$^{-1}$]')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc=[0.45, 0.65], labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{siminfo.output_path}/Cosmic_scatter_rate.png", dpi=200)
    plt.clf()