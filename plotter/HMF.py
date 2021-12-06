
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter

from colossus.cosmology import cosmology
cosmology.setCosmology('planck13');
from colossus.lss import mass_function

from velociraptor import load


class make_HMF_data:

    def __init__(self, hmf_M, hmf_dndM, hmf_error,
                 sub_hmf_M, sub_hmf_dndM, sub_hmf_error,
                 hmf_counter, sub_hmf_counter):

        self.hmf_counter = hmf_counter
        self.sub_hmf_counter = sub_hmf_counter
        self.hmf_M = hmf_M
        self.hmf_dndM = hmf_dndM
        self.hmf_error = hmf_error
        self.sub_hmf_M = sub_hmf_M
        self.sub_hmf_dndM = sub_hmf_dndM
        self.sub_hmf_error = sub_hmf_error

    def add_data(self, hmf_M, hmf_dndM, hmf_error,
                 sub_hmf_M, sub_hmf_dndM, sub_hmf_error):

        self.hmf_M = np.append(self.hmf_M,hmf_M)
        self.hmf_dndM = np.append(self.hmf_dndM,hmf_dndM)
        self.hmf_error = np.append(self.hmf_error,hmf_error)
        self.sub_hmf_M = np.append(self.sub_hmf_M,sub_hmf_M)
        self.sub_hmf_dndM = np.append(self.sub_hmf_dndM,sub_hmf_dndM)
        self.sub_hmf_error = np.append(self.sub_hmf_error,sub_hmf_error)

def colossus_hmf(log_m, sim_info):
    z = sim_info.z # Redshift
    M = 10 ** log_m # Msun/h
    f = mass_function.massFunction(M, z, mdef = '200c', model = 'tinker08', q_out='dndlnM')
    return f

def HMF(mass, bins, dlogm, V):
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

def make_HMF(sim_info, HMF_data):

    dlogm = 0.4
    bins = 10 ** (np.arange(6, 15, dlogm))
    V = (sim_info.boxSize * 1e-3) ** 3 # Mpc

    # Load the data
    catalogue = load(sim_info.halo_data.path_to_catalogue)
    mass = catalogue.masses.mass_200crit.to("Msun").value
    type = catalogue.structure_type.structuretype
    subhalo = type > 10

    hmf_M, hmf_dndM, hmf_error = HMF(mass, bins, dlogm, V)
    sub_hmf_M, sub_hmf_dndM, sub_hmf_error = HMF(mass[subhalo], bins, dlogm, V)

    hmf_counter = len(hmf_M)
    sub_hmf_counter = len(sub_hmf_M)

    if HMF_data == None:
        HMF_data = make_HMF_data(hmf_M, hmf_dndM, hmf_error,
                                 sub_hmf_M, sub_hmf_dndM, sub_hmf_error, hmf_counter, sub_hmf_counter)
    else:
        HMF_data.add_data(hmf_M, hmf_dndM, hmf_error,
                          sub_hmf_M, sub_hmf_dndM, sub_hmf_error)

    return HMF_data


def plot_HMF(HMF_data, sim_info, name_list):

    limit = 100 * sim_info.dm_particle_mass #Msun

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
    counter = 0
    for name in name_list:

        M = HMF_data.hmf_M[counter:counter+HMF_data.hmf_counter]
        dndM = HMF_data.hmf_dndM[counter:counter + HMF_data.hmf_counter]
        error = HMF_data.hmf_error[counter:counter + HMF_data.hmf_counter]

        plt.plot(M, dndM, label=name, lw=2, color=color[i])
        plt.fill_between(M, dndM - error, dndM + error, alpha=0.4, color=color[i])
        i += 1
        counter += HMF_data.hmf_counter

    M = np.arange(6, 16, 0.2)
    mfunc = colossus_hmf(M, sim_info) / sim_info.h ** 3
    plt.plot(10 ** M * sim_info.h, mfunc, '-', label = 'Tinker et al. (2008)',color='black',zorder=1)

    plt.xscale("log")
    plt.ylim(1e-5, 1e2)
    plt.xlim(10 ** 6, 10 ** 15)
    plt.axvline(x=limit, linestyle="--", lw=1, color="grey")

    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("dn/d$\log$10(M${}_{200}$) (Mpc$^{-3}$)")

    ax.tick_params(direction='in',axis='both',which='both',pad=4.5)
    plt.yscale("log")
    plt.legend(loc="upper right",labelspacing=0.2,handlelength=1.5,handletextpad=0.4,frameon=False)

    output_file = f"{sim_info.output_path}/HMF_"+ name_list[0]
    if len(name_list) > 1:
        output_file += "_" + name_list[1] +".png"
    else:
        output_file += ".png"

    plt.savefig(output_file, dpi=200)
    plt.close()


    figure()
    ax = plt.subplot(1,1,1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    counter = 0
    for name in name_list:

        M = HMF_data.sub_hmf_M[counter:counter+HMF_data.sub_hmf_counter]
        dndM = HMF_data.sub_hmf_dndM[counter:counter + HMF_data.sub_hmf_counter]
        error = HMF_data.sub_hmf_error[counter:counter + HMF_data.sub_hmf_counter]

        plt.plot(M, dndM, label=name, lw=2, color=color[i])
        plt.fill_between(M, dndM - error, dndM + error, alpha=0.4, color=color[i])
        i += 1
        counter += HMF_data.sub_hmf_counter

    plt.xscale("log")
    plt.ylim(1e-5, 1e2)
    plt.xlim(10 ** 6, 10 ** 15)
    plt.axvline(x=limit, linestyle="--", lw=1, color="grey")

    M = np.arange(6, 16, 0.2)
    mfunc = colossus_hmf(M, sim_info) / (sim_info.h)**3
    plt.plot(10 ** M * sim_info.h, mfunc, '-', label='Tinker et al. (2008)', color='black', zorder=1)

    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("(Subhalo) dn/d($\log$10(M${}_{200}$) (Mpc$^{-3}$)")
    ax.tick_params(direction='in',axis='both',which='both',pad=4.5)
    plt.yscale("log")
    plt.legend(loc="upper right",labelspacing=0.2,handlelength=1.5,handletextpad=0.4,frameon=False)
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    output_file = f"{sim_info.output_path}/SubHMF_"+ name_list[0]
    if len(name_list) > 1:
        output_file += "_" + name_list[1] +".png"
    else:
        output_file += ".png"

    plt.savefig(output_file, dpi=200)
    plt.close()
