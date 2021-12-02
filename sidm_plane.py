"""
Description here
"""
from argumentparser import ArgumentParser
from object import simulation_data
from tqdm import tqdm
import h5py
import numpy as np
import scipy.stats as stat
from scipy.optimize import curve_fit
from pylab import *
import matplotlib.pyplot as plt
from fit_profile.fit_profile import calculate_profile_params

def fit_model(x, a, b):
    f = a + b * x
    return f


def plot_data(sim_info, output_name_list):


    # Plot parameters
    params = {
    "font.size": 11,
    "font.family": "Times",
    "text.usetex": True,
    "figure.figsize": (4, 3),
    "figure.subplot.left": 0.18,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "lines.markersize": 1,
    "lines.linewidth": 1.,
    "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    grid(True)

    x = []
    y = []

    for name in output_name_list:

        filename = f"{sim_info.output_path}/Profile_params_" + name + ".hdf5"
        with h5py.File(filename, "r") as file:
            log10_M200c = file["Data/M200c"][:]
            Type = file["Data/StructureType"][:]
            r0 = file["Data/r0_isothermal_fit"][:]
            rho0 = file["Data/rho0_isothermal_fit"][:]
            v0 = file["Data/velocity_dispersion_0"][:]

        sample = np.where((log10_M200c >= 9) & (log10_M200c <= 11) & (r0 > 0.01))[0]
        cen = np.where(Type[sample] == 10)[0]
        cen = sample[cen]

        r0 = r0[cen]
        rho0 = rho0[cen]
        v0 = v0[cen]

        mass_bins = np.arange(9, 11, 0.1)
        median_rho0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=rho0, statistic="median", bins=mass_bins, )
        median_r0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=r0, statistic="median", bins=mass_bins, )
        median_v0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=v0, statistic="median", bins=mass_bins, )

        # #if name == 'DML025N752SigmaConstant00':
        #     # plot(0.001 * median_r0, median_v0 * 10**median_rho0,'o',label=name)
        if name == 'DML025N752SigmaConstant01':
            x_sim = np.log10(v0)
            y_sim = rho0
            plot(median_v0, 10 ** median_rho0, 'o', label=name)
        if name == 'DML025N752SigmaConstant10':
            x_sim = np.log10(10 * v0)
            y_sim = rho0
            plot(10*median_v0, 10**median_rho0,'o',label=name)

        #plot(10**x_sim, 10**y_sim,'o',label=name)
        select = np.where(mass_bins>10.6)[0]

        x = np.append(x, x_sim[select])
        y = np.append(y, y_sim[select])

    popt, pcov = curve_fit(fit_model, x, y, p0=[1, 10])
    print(popt, 'x=sigma0 x r0, y=rho0 (Full-ISO)')
    #xrange = np.arange(-1, 2, 0.1)
    #yrange = fit_model(xrange, *popt)
    #plt.plot(10 ** xrange, 10 ** yrange, ':',label=name)
        #yrange = fit_model(xrange, popt[0], -0.66)
        #plt.plot(10 ** xrange, 10 ** yrange, '--')
    # plt.text(0.1, 0.04, "a={0:.2f},".format(popt[0]) + " b={0:.2f}".format(popt[1]) + " (Full-ISO)",
    #          transform=ax.transAxes)

    xscale('log')
    yscale('log')
    xlabel(r'($\sigma/m$)$\times$ r$_{0}$ [(cm$^{2}$/g)$\times$ kpc]')
    ylabel(r'$\rho_{0}$ [M$_{\odot}$/kpc$^{3}$]')
    # axis([0, 50, 5, 10])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="lower left", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/SIDM_plane.png", dpi=200)
    plt.close()



if __name__ == "__main__":

    config = ArgumentParser()
    output_name_list = []

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        # Fetch relevant input parameters from lists
        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]
        output = config.output_directory

        # Load all data and save it in SimInfo class
        sim_info = simulation_data.SimInfo(
            directory=directory,
            snapshot=snapshot,
            catalogue=catalogue,
            name=sim_name,
            output=output,
        )

        output_name_list.append(sim_info.simulation_name)

        # calculate_profile_params(sim_info,
        #                          log10_min_mass=9,
        #                          log10_max_mass=11)

    plot_data(sim_info, output_name_list)
