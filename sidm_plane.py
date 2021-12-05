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
    "lines.markersize": 0.3,
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

        if name == 'DML025N752SigmaVelDep20Isotropic':

            filename = f"{sim_info.output_path}/Profile_params_" + name + ".hdf5"
            with h5py.File(filename, "r") as file:
                log10_M200c = file["Data/M200c"][:]
                Type = file["Data/StructureType"][:]
                r0 = file["Data/r0_isothermal_fit"][:]
                rho0 = file["Data/rho0_isothermal_fit"][:]
                v0 = file["Data/velocity_dispersion_0"][:]
                sigma0 = file["Data/cross_section_0"][:]

            sample = np.where((log10_M200c >= 9.8) & (log10_M200c <= 11) & (r0 > 0.01))[0]
            cen = np.where(Type[sample] == 10)[0]
            cen = sample[cen]

            r0 = r0[cen]
            rho0 = rho0[cen]
            v0 = v0[cen]
            sigma0 = sigma0[cen]

            mass_bins = np.arange(9, 11, 0.1)
            median_rho0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=rho0, statistic="median", bins=mass_bins, )
            median_r0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=r0, statistic="median", bins=mass_bins, )
            median_v0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=v0, statistic="median", bins=mass_bins, )
            median_s0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=sigma0, statistic="median", bins=mass_bins, )

            print(median_s0)
            print(median_r0)
            print(median_rho0)
            print(mass_bins)

            x = r0
            y = 10**rho0
            #plot(x, y,'o',color='tab:orange',label=name)
            plot(median_r0,10**median_rho0,'o',color='white',ms=3)
            plot(median_r0,10**median_rho0,'o',color='tab:orange',ms=2)

            # x = np.log10(x)
            # y = np.log10(y)
            # popt, pcov = curve_fit(fit_model, x, y, p0=[1, 10])
            # print(popt, 'x=sigma0 x r0, y=rho0 (Full-ISO)')
            # xrange = np.arange(0, 2, 0.1)
            # yrange = fit_model(xrange, *popt)
            # plt.plot(10 ** xrange, 10 ** yrange, ':',color='black')

        if name == 'DML025N752SigmaConstant01':

            filename = f"{sim_info.output_path}/Profile_params_" + name + ".hdf5"
            with h5py.File(filename, "r") as file:
                log10_M200c = file["Data/M200c"][:]
                Type = file["Data/StructureType"][:]
                r0 = file["Data/r0_isothermal_fit"][:]
                rho0 = file["Data/rho0_isothermal_fit"][:]
                v0 = file["Data/velocity_dispersion_0"][:]

            sample = np.where((log10_M200c >= 9.8) & (log10_M200c <= 11) & (r0 > 0.01))[0]
            cen = np.where(Type[sample] == 10)[0]
            cen = sample[cen]

            r0 = r0[cen]
            rho0 = rho0[cen]
            v0 = v0[cen]

            mass_bins = np.arange(9, 11, 0.1)
            median_rho0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=rho0, statistic="median", bins=mass_bins, )
            median_r0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=r0, statistic="median", bins=mass_bins, )
            median_v0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=v0, statistic="median", bins=mass_bins, )

            x = r0 * v0
            y = 10 ** rho0
            #plot(x, y, 'o',color='tab:blue', label=name)
            plot(median_r0, 10 ** median_rho0, 'o',color='white', ms=3)
            plot(median_r0, 10 ** median_rho0, 'o',color='tab:blue', ms=2)

        if name == 'DML025N752SigmaConstant10':

            filename = f"{sim_info.output_path}/Profile_params_" + name + ".hdf5"
            with h5py.File(filename, "r") as file:
                log10_M200c = file["Data/M200c"][:]
                Type = file["Data/StructureType"][:]
                r0 = file["Data/r0_pseudo_isothermal_fit"][:]
                rho0 = file["Data/rho0_pseudo_isothermal_fit"][:]
                v0 = file["Data/velocity_dispersion_0"][:]
                # r0 = file["Data/r0_isothermal_fit"][:]
                # rho0 = file["Data/rho0_isothermal_fit"][:]
                # v0 = file["Data/velocity_dispersion_0"][:]

            sample = np.where((log10_M200c >= 9.8) & (log10_M200c <= 11) & (r0 > 0.001))[0]
            cen = np.where(Type[sample] == 10)[0]
            cen = sample[cen]

            r0 = r0[cen]
            rho0 = rho0[cen]
            v0 = v0[cen]

            mass_bins = np.arange(9, 11, 0.1)
            median_rho0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=rho0, statistic="median", bins=mass_bins, )
            median_r0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=r0, statistic="median", bins=mass_bins, )
            median_v0, _, _ = stat.binned_statistic(x=log10_M200c[cen], values=v0, statistic="median", bins=mass_bins, )

            x = r0
            y = 10 ** rho0
            plot(x, y, 'o',color='tab:red', label=name,alpha=0.5)
            plot(median_r0, 10 ** median_rho0, 'o',color='white', ms=3)
            plot(median_r0, 10 ** median_rho0, 'o',color='tab:red', ms=2)

        if name == 'DML025N752SigmaConstant00':

            filename = f"{sim_info.output_path}/Profile_params_" + name + ".hdf5"
            with h5py.File(filename, "r") as file:
                log10_M200c = file["Data/M200c"][:]
                Type = file["Data/StructureType"][:]
                r0 = file["Data/r0_pseudo_isothermal_fit"][:]
                rho0 = file["Data/rho0_pseudo_isothermal_fit"][:]
                v0 = file["Data/velocity_dispersion_0"][:]
                # r0 = file["Data/r0_isothermal_fit"][:]
                # rho0 = file["Data/rho0_isothermal_fit"][:]
                # v0 = file["Data/velocity_dispersion_0"][:]

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

            x = r0 * v0
            y = 10 ** rho0
            #plot(x, y, 'o',color='tab:green', label=name)
            plot(median_r0, 10 ** median_rho0, 'o',color='white', ms=3)
            plot(median_r0, 10 ** median_rho0, 'o',color='tab:green', ms=2)


    xscale('log')
    yscale('log')
    xlabel(r'($\sigma/m$)$\times$ r$_{0}$ [(cm$^{2}$/g)$\times$ kpc]')
    ylabel(r'$\rho_{0}$ [M$_{\odot}$/kpc$^{3}$]')
    #axis([1e-1, 5, 1e7, 1e8])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="lower right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/SIDM_plane.png", dpi=200)
    plt.close()

def plot_relations(sim_info, output_name_list):


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
    "lines.markersize": 0.3,
    "lines.linewidth": 1.,
    "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    for name in output_name_list:

        filename = f"{sim_info.output_path}/Profile_params_" + name + ".hdf5"
        with h5py.File(filename, "r") as file:
            index = file["Data/ID"][:]
            Type = file["Data/StructureType"][:]
            r0 = file["Data/r0_isothermal_fit"][:]
            rho0 = file["Data/rho0_isothermal_fit"][:]
            r0pse = file["Data/r0_pseudo_isothermal_fit"][:]
            rho0pse = file["Data/rho0_pseudo_isothermal_fit"][:]
            log10_M200c = file["Data/M200c"][:]

        cen = np.where((Type == 10) & (r0>0.01))[0]
        index = index[cen]
        r0 = r0[cen]
        rho0 = rho0[cen]
        r0pse = r0pse[cen]
        rho0pse = rho0pse[cen]
        log10_M200c = log10_M200c[cen]

        rs = []
        for i in range(len(index)):
            select = sim_info.halo_data.halo_index == index[i]
            rs = np.append(rs, sim_info.halo_data.scale_radius[select] * 1e3)

        mass_bins = np.arange(9, 11, 0.1)
        median_rho0, _, _ = stat.binned_statistic(x=log10_M200c, values=rho0, statistic="median", bins=mass_bins, )
        median_r0, _, _ = stat.binned_statistic(x=log10_M200c, values=r0, statistic="median", bins=mass_bins, )
        median_rs, _, _ = stat.binned_statistic(x=log10_M200c, values=rs, statistic="median", bins=mass_bins, )
        median_M200c, _, _ = stat.binned_statistic(x=log10_M200c, values=log10_M200c, statistic="median", bins=mass_bins, )
        median_rho0pse, _, _ = stat.binned_statistic(x=log10_M200c, values=rho0pse, statistic="median", bins=mass_bins, )
        median_r0pse, _, _ = stat.binned_statistic(x=log10_M200c, values=r0pse, statistic="median", bins=mass_bins, )

        figure()
        ax = plt.subplot(1, 1, 1)
        grid(True)

        plot(rs, r0, 'o', color='tab:orange', alpha=0.5)
        plot(median_rs, median_r0, 'o',color='white', ms=3)
        plot(median_rs, median_r0, 'o',color='tab:orange', ms=2,label='Isothermal fit')

        plot(rs, r0pse, 'o', color='tab:blue', alpha=0.5)
        plot(median_rs, median_r0pse, 'o',color='white', ms=3)
        plot(median_rs, median_r0pse, 'o',color='tab:blue', ms=2,label='Pseudo-Isothermal fit')

        xscale('log')
        yscale('log')
        ylabel(r'$r_{0}$ [kpc]')
        xlabel(r'$r_{\mathrm{scale}}$ [kpc]')
        axis([1e-1, 1e2, 1e-1, 1e1])
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc="lower right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
        plt.savefig(f"{sim_info.output_path}/centrals_r0_rscale_relation_"+name+".png", dpi=200)
        plt.close()

        figure()
        ax = plt.subplot(1, 1, 1)
        grid(True)

        plot(10**log10_M200c, r0, 'o', color='tab:orange', alpha=0.5)
        plot(10**median_M200c, median_r0, 'o',color='white', ms=3)
        plot(10**median_M200c, median_r0, 'o',color='tab:orange', ms=2,label='Isothermal fit')

        plot(10**log10_M200c, r0pse, 'o', color='tab:blue', alpha=0.5)
        plot(10**median_M200c, median_r0pse, 'o',color='white', ms=3)
        plot(10**median_M200c, median_r0pse, 'o',color='tab:blue', ms=2,label='Pseudo-Isothermal fit')

        xscale('log')
        yscale('log')
        ylabel(r'$r_{0}$ [kpc]')
        xlabel(r'$M_{200c}$ [M$_{\odot}$]')
        axis([1e9,1e11, 1e-1, 1e1])
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc="lower right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
        plt.savefig(f"{sim_info.output_path}/centrals_r0_M200c_relation_"+name+".png", dpi=200)
        plt.close()

        figure()
        ax = plt.subplot(1, 1, 1)
        grid(True)

        plot(10**log10_M200c, 10**rho0, 'o', color='tab:orange', alpha=0.5)
        plot(10**median_M200c, 10**median_rho0, 'o',color='white', ms=3)
        plot(10**median_M200c, 10**median_rho0, 'o',color='tab:orange', ms=2,label='Isothermal fit')

        plot(10**log10_M200c, 10**rho0pse, 'o', color='tab:blue', alpha=0.5)
        plot(10**median_M200c, 10**median_rho0pse, 'o',color='white', ms=3)
        plot(10**median_M200c, 10**median_rho0pse, 'o',color='tab:blue', ms=2,label='Pseudo-Isothermal fit')

        xscale('log')
        yscale('log')
        ylabel(r'$\rho_{0}$ [M$_{\odot}$/kpc$^{3}$]')
        xlabel(r'$M_{200c}$ [M$_{\odot}$]')
        axis([1e9, 1e11, 1e5, 1e9])
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc="lower right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
        plt.savefig(f"{sim_info.output_path}/centrals_rho0_mass_relation_"+name+".png", dpi=200)
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
        #                          log10_max_mass=9.01)

    #plot_data(sim_info, output_name_list)
    plot_relations(sim_info, output_name_list)
