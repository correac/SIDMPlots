"""
Description here
"""
from argumentparser import ArgumentParser
from object import simulation_data
from time import time
from tqdm import tqdm
import h5py
import numpy as np
from merger_tree.build_tree import build_tree
from pylab import *
import matplotlib.pyplot as plt
import matplotlib as mpl

def read_rotation_curve_data(sim_info, min_halo_mass, max_halo_mass):

    filename = f"{sim_info.output_path}/Rotation_data_" + sim_info.simulation_name + ".hdf5"
    with h5py.File(filename, "r") as file:
        M200c = file["Data/M200c"][:]
        Type = file["Data/StructureType"][:]
        v_circ_2 = file["Data/Vcirc2kpc"][:]
        v_circ_nfw_2 = file["Data/Vcirc_nfw_2kpc"][:]
        halo_index = file["Data/ID"][:]

    ratio = (v_circ_2 - v_circ_nfw_2) / v_circ_nfw_2

    select = Type > 10
    halo_index = halo_index[select]
    ratio = ratio[select]
    M200c = M200c[select]

    select = np.logical_and(M200c >= min_halo_mass,M200c <= max_halo_mass)
    halo_index = halo_index[select]
    ratio = ratio[select]

    return halo_index, ratio

if __name__ == "__main__":

    config = ArgumentParser()
    #time_start = time()

    # Fetch relevant input parameters from lists
    directory = config.directory_list[0]
    snapshot = config.snapshot_list[0]
    catalogue = config.catalogue_list[0]
    sim_name = config.name_list[0]
    output = config.output_directory

    # Load all data and save it in SimInfo class
    sim_info = simulation_data.SimInfo(
        directory=directory,
        snapshot=snapshot,
        catalogue=catalogue,
        name=sim_name,
        output=output,
    )

    halo_index, ratio = read_rotation_curve_data(sim_info,
                                                 min_halo_mass=9.0,
                                                 max_halo_mass=9.25)

    print("Tracking %i haloes.."%len(halo_index))

    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 3.5),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.87,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "lines.markersize": 1,
        "lines.linewidth": 1,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    num_sample = len(halo_index)

    norm = mpl.colors.Normalize(vmin=-1, vmax=3)
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.seismic)
    cmap.set_array([])

    mass_loss = np.zeros(num_sample)
    z_accretion = np.ones(num_sample) * 4

    for i in tqdm(range(num_sample)):
        progenitor_index, merger_mass_ratio, mass, type, z = build_tree(sim_info, halo_index[i])

        if ratio[i] > 0.5: plt.plot(1+z, mass,'--', c=cmap.to_rgba(ratio[i]))

        mass_loss[i] = (np.max(mass)-mass[0])/mass[0]

        # Let's highlight when satellite is not satellite
        central = np.where(type==10)[0]
        if len(central)>2:
            if ratio[i] > 0.5: plt.plot(1+z[central.min():], mass[central.min():],'-', c=cmap.to_rgba(ratio[i]))
            if central.min() > 1:
                z_accretion[i] = z[central.min() - 1]

    colorbar(cmap, label=r"$\Delta V_{\mathrm{circ}}$(2~kpc)/$V_{\mathrm{circ,NFW}}$")
    plt.axis([1, 5, 1e8, 1e10])
    plt.xlabel("$z$")
    plt.ylabel("$\log_{10}$ M [M$_{\odot}$]")
    plt.yscale('log')
    plt.xscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.xticks(np.arange(1,6), ['0','1','2','3','4'])
    plt.savefig(f"{sim_info.output_path}/Satellites_mass_evolution_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    # Compute how much time it took to run the script
    #time_end = time()
    #script_runtime = time_end - time_start
    #m, s = divmod(script_runtime, 60)
    #h, m = divmod(m, 60)
    #print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

    ################
    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (5, 3.5),
        "figure.subplot.left": 0.15,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "lines.markersize": 1.5,
        "lines.linewidth": 1,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(ratio, mass_loss, 'o')

    plt.axis([-1, 2, 0.1, 100])
    plt.xlabel(r"$\Delta V_{\mathrm{circ}}$(2~kpc)/$V_{\mathrm{circ,NFW}}$")
    plt.ylabel("$\Delta M/M(z=0)$")
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{sim_info.output_path}/Satellites_mass_loss_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(ratio, z_accretion, 'o')

    plt.axis([-1, 2, 0, 4.1])
    plt.xlabel(r"$\Delta V_{\mathrm{circ}}$(2~kpc)/$V_{\mathrm{circ,NFW}}$")
    plt.ylabel("Accretion redshift")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{sim_info.output_path}/Satellites_accretion_time_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

