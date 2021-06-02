import h5py
import numpy as np
import glob
import os.path
from pylab import *
import matplotlib.pyplot as plt
from scipy.special import spence
from colossus.cosmology import cosmology
cosmology.setCosmology('planck13');
from colossus.lss import mass_function


def calc_density(x, M, a, z, siminfo):
    Om = siminfo.Om
    Ol = 1. - Om
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    # NFW profile #
    c = R200 / a
    delta = 200. / 3.
    delta *= c ** 3 / (np.log(1. + c) - c / (1. + c))
    f = np.zeros(len(x))
    for i in range(0, len(x)): f[i] = rho_crit * delta * 1e-9 / ((x[i] / a) * (1. + x[i] / a) ** 2)
    return f


def sigma_1D(x, M, a, z, siminfo):
    Om = siminfo.Om
    Ol = 1. - Om
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    G = 4.3e-6  # kpc km^2 Msun^-1 s^-2
    ff = np.zeros(len(x))
    for i in range(0, len(x)):
        s = x[i] / R200
        c = R200 / a
        u = 1. + c * s
        gc = 1. / (np.log(1. + c) - c / (1. + c))
        f = 0.5 * c ** 2 * gc * s * u ** 2 * G * M / R200
        u = 1. + c * s
        Li = spence(u)
        f *= (np.pi ** 2 - np.log(c * s) - 1. / (c * s) - 1. / u ** 2 - 6. / u + (
                    1. + 1. / (c ** 2 * s ** 2) - 4. / (c * s) - 2. / u) * np.log(u) + 3 * (np.log(u)) ** 2 + 6. * Li)
        ff[i] = f
    return np.sqrt(ff)


def gamma(M, z, siminfo):

    sigma = siminfo.sigma
    Om = siminfo.Om
    Ol = 1. - Om
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    # output = commah.run('Planck13',zi=z,Mi=M,z=z)
    # c200 = output['c']
    # c200 = 5.72 * (M/(1e14/h))**(-0.081) * (1.+z)**(-0.71)
    c200 = 5.72 * (M / 1e14) ** (-0.081) * (1. + z) ** (-0.71)

    a = R200 / c200  # r200/ r200/rs [kpc]

    radial_bins = np.arange(-5, np.log10(R200), 0.1)
    x = 10 ** radial_bins
    deltax = x[1:] - x[:1]
    deltax = np.append(deltax, deltax[-1])
    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21

    sig1D = sigma_1D(x, M, a, z, siminfo)  # km/s
    sig1D[np.isnan(sig1D)] = 0
    sig1D *= 1e5  # cm/s

    rho = calc_density(x, M, a, z, siminfo)  # Msun/kpc^3
    rho *= Msun_in_cgs / kpc_in_cgs ** 3  # g/cm^2

    f = sig1D * rho * (4. / np.sqrt(np.pi)) * sigma  # 1/s
    Gyr = 1e9 * 365.25 * 24 * 3600

    f *= Gyr  # Gyr^-1 particle^-1

    rho = calc_density(x, M, a, z, siminfo)  # Msun/kpc^3
    vol = 4 * np.pi * x ** 2 * deltax * rho  # Msun
    Total = np.sum(f * vol) / M
    return Total


def analytic_scatter(siminfo, z, Mmin, Mmax):
    h = siminfo.h
    Om = siminfo.Om
    rho_crit0 = siminfo.rhocrit0
    rho_crit = rho_crit0 #* ( Om * (1.+ z)**3 + Ol)

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



def cosmic_scatter_rate(siminfo):

    unit = 3.085677e19 / (3600 * 24 * 365.25 * 1e9) # Gyr
    sim = h5py.File(siminfo.snapshot(0), "r")
    old_time = sim["/Cosmology"].attrs["Universe age [internal units]"] * unit
    n_sidm_events = sim["/PartType1/SIDM_events"][:]
    nparts = len(n_sidm_events)
    scatter_rate = np.zeros(siminfo.n_snapshots+1)
    redshift = np.zeros(siminfo.n_snapshots+1)
    redshift[0] = sim["/Header"].attrs["Redshift"]
    previous_events = 0

    for i in range(1,siminfo.n_snapshots+1):

        sim = h5py.File(siminfo.snapshot(i), "r")
        n_sidm_events = sim["/PartType1/SIDM_events"][:]
        redshift[i] = sim["/Header"].attrs["Redshift"]
        time = sim["/Cosmology"].attrs["Universe age [internal units]"] * unit
        delta_time = time - old_time
        old_time = time

        n_collisions = np.sum(n_sidm_events)-previous_events
        previous_events = np.sum(n_sidm_events)
        scatter_rate[i] = n_collisions / ( delta_time * nparts )

    return redshift, scatter_rate


def output_cosmic_scatter_rate(siminfo):

    redshift, scatter_rate = cosmic_scatter_rate(siminfo)

    np.savetxt(f"{siminfo.output_path}/Cosmic_scatter_rate_" + siminfo.name + ".txt",
               np.transpose([redshift, scatter_rate]))


def plot_cosmic_scatter_rate(siminfo, name_list):

    #######################
    # Plot parameters
    params = {
        "font.size": 12,
    #    "font.family": "Times",
    #    "text.usetex": True,
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
    color = ['black', 'tab:grey']
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Cosmic_scatter_rate_" + name + ".txt")
        redshift = data[:,0]
        scatter_rate = data[:,1]
        plot(1 + redshift, scatter_rate, '-', color=color[k], label=name)
        k += 1

    xscale('log')
    xlabel('$1+z$')
    ylabel(r'$\Gamma(r)$ [Gyr$^{-1}$ particle$^{-1}$]')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc=[0.45, 0.65], labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    # axis([1,50,0,0.1])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{siminfo.output_path}/Cosmic_scatter_rate.png", dpi=200)
    plt.clf()
