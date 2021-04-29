
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter

from colossus.cosmology import cosmology
cosmology.setCosmology('planck13');
from colossus.lss import mass_function

def colossus_hmf(log_m, siminfo):
    z = np.around(siminfo.z)[0]  # Redshift
    M = 10 ** log_m # Msun/h
    f = mass_function.massFunction(M, z, mdef = '200c', model = 'tinker08', q_out='dndlnM')
    return f

def HMF(siminfo,halo):

    dlogm = 0.2
    bins = 10 ** (np.arange(6, 15, dlogm))
    V = siminfo.boxSize ** 3 # Mpc

    # Load the data
    g = h5py.File(siminfo.halo_properties, "r")
    mass = g["Mass_200crit"][:] * 1e10  # convert to Msun

    if halo == 'subhalo':
        subtype = g["Structuretype"][:]
        subhalo = subtype > 10
        mass = mass[subhalo]

    binnedmass, massrange = np.histogram(mass, bins=bins)

    massnlarger = np.zeros(len(binnedmass))
    for i in range(0, len(massnlarger)):
        massnlarger[i] = np.sum(binnedmass[i:])

    # Determine the HMF
    errormassn = massnlarger ** 0.5
    numbden = massnlarger / V
    numbdenerr = errormassn / V
    massplot = (massrange[:-1] + massrange[1:]) / 2
    dernumbden = -np.diff(numbden) / np.diff(np.log10(massplot))
    dererr = 2 ** 0.5 / dlogm * (numbdenerr[:-1] + numbdenerr[1:]) / 2

    M = (massplot[:-1] + massplot[1:]) / 2.
    dndM = dernumbden
    error = dererr

    return M, dndM, error


def HMF_output_files(siminfo):

    M, dndM, error = HMF(siminfo, 'halo')
    np.savetxt(f"{siminfo.output_path}/HMF_" + siminfo.name + ".txt",
               np.transpose([M, dndM, error]))

    M, dndM, error = HMF(siminfo, 'subhalo')
    np.savetxt(f"{siminfo.output_path}/subHMF_" + siminfo.name + ".txt",
               np.transpose([M, dndM, error]))


def make_HMF(siminfo, name_list):

    limit = 100 * siminfo.dmpart_mass #Msun

    # Plot parameters
    params = {
        "font.size": 14,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 4),
        "figure.subplot.left": 0.17,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.15,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.3,
        "figure.subplot.hspace": 0.3,
        "lines.markersize": 6,
        "lines.linewidth": 1.0,
        "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1,1,1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/HMF_" + name + ".txt")

        M = data[:,0]
        dndM = data[:,1]
        error = data[:,2]
        plt.plot(M, dndM, label=name, lw=2, color=color[i])
        plt.fill_between(M, dndM - error, dndM + error, alpha=0.4, color=color[i])
        i += 1

    M = np.arange(6, 16, 0.2)
    mfunc = colossus_hmf(M, siminfo) / (siminfo.h)**3
    plt.plot(10 ** M * siminfo.h, mfunc, '-', label = 'Tinker et al. (2008)',color='black',zorder=1)

    plt.xscale("log")
    plt.ylim(1e-5, 1e2)
    plt.xlim(10 ** 6, 10 ** 15)
    plt.axvline(x=limit, linestyle="--", lw=1, color="grey")

    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("dn/d($\log$10(M${}_{200}$) (Mpc$^{-3}$)")
    ax.tick_params(direction='in',axis='both',which='both',pad=4.5)
    plt.yscale("log")
    plt.legend(loc="upper right",labelspacing=0.2,handlelength=1.5,handletextpad=0.4,frameon=False)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{siminfo.output_path}/HMF_%04d.png" % siminfo.n_snapshots, dpi=200)
    plt.close()


    figure()
    ax = plt.subplot(1,1,1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/subHMF_" + name + ".txt")

        M = data[:, 0]
        dndM = data[:, 1]
        error = data[:, 2]
        plt.plot(M, dndM, label=name, lw=2, color=color[i])
        plt.fill_between(M, dndM - error, dndM + error, alpha=0.4, color=color[i])
        i += 1

    plt.xscale("log")
    plt.ylim(1e-5, 1e2)
    plt.xlim(10 ** 6, 10 ** 15)
    plt.axvline(x=limit, linestyle="--", lw=1, color="grey")

    M = np.arange(6, 16, 0.2)
    mfunc = colossus_hmf(M, siminfo) / (siminfo.h)**3
    plt.plot(10 ** M * siminfo.h, mfunc, '-', label='Tinker et al. (2008)', color='black', zorder=1)

    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("(Subhalo) dn/d($\log$10(M${}_{200}$) (Mpc$^{-3}$)")
    ax.tick_params(direction='in',axis='both',which='both',pad=4.5)
    plt.yscale("log")
    plt.legend(loc="upper right",labelspacing=0.2,handlelength=1.5,handletextpad=0.4,frameon=False)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{siminfo.output_path}/SubHMF_%04d.png" % siminfo.n_snapshots, dpi=200)
    plt.close()
