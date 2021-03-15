"""
Description
"""
import h5py
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scatter_rate import calc_density, sigma_1D
import scipy.stats as stat
from scipy.special import spence


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


def plot_relations(siminfo,output_path):

    # Load data:
    g = h5py.File(siminfo.halo_properties, "r")
    mass = g["Mass_200crit"][:] * 1e10  # convert to Msun
    c200 = g["cNFW_200crit"][:]
    Vmax = g["Vmax"][:] #km/s
    subtype = g["Structuretype"][:]
    centrals = subtype == 10

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

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    xrange = np.arange(7,15,0.2)
    xrange = 10**xrange
    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass[centrals], c200[centrals], xrange)

    plt.plot(mass[centrals], c200[centrals], 'o', color='tab:blue')
    plt.plot(xvalues, yvalues, '-', color='tab:blue')
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:blue')
    plt.axis([1e9,1e15,1,20])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("c$_{200}$")
    plt.xscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(output_path+"cM_relation.png", dpi=200)
    plt.close()

    ###########
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    xrange = np.arange(7, 15, 0.2)
    xrange = 10 ** xrange
    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass[centrals], Vmax[centrals], xrange)

    plt.plot(mass[centrals], Vmax[centrals], 'o', color='tab:blue')
    plt.plot(xvalues, yvalues, '-', color='tab:blue')
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:blue')
    plt.axis([1e9, 1e15, 0, 600])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("V$_{\mathrm{max}}$ [km/s]")
    plt.xscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(output_path + "VmaxM_relation.png", dpi=200)
    plt.close()


def bin_volumes(radial_bins):
    """Returns the volumes of the bins. """

    single_vol = lambda x: (4.0 / 3.0) * np.pi * x ** 3
    outer = single_vol(radial_bins[1:])
    inner = single_vol(radial_bins[:-1])
    return outer - inner

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)

def analyse_halo(mass, pos, vel):
    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    SumMasses, _, _ = stat.binned_statistic(x=r, values=np.ones(len(r)) * mass[0], statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3

    std_vel_x, _, _ = stat.binned_statistic(x=r, values=vel[:, 0], statistic="std", bins=radial_bins, )
    std_vel_y, _, _ = stat.binned_statistic(x=r, values=vel[:, 1], statistic="std", bins=radial_bins, )
    std_vel_z, _, _ = stat.binned_statistic(x=r, values=vel[:, 2], statistic="std", bins=radial_bins, )
    velocity = np.sqrt(std_vel_x ** 2 + std_vel_y ** 2 + std_vel_z ** 2) / np.sqrt(3.)
    velocity[np.where(np.isnan(velocity))[0]] = 0
    return density, velocity


def read_data(siminfo):

    with h5py.File(siminfo.latest_snapshot, "r") as hf:
        mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        pos = hf['PartType1/Coordinates'][:][:] * siminfo.a
        vel = hf['PartType1/Velocities'][:][:]
        # Read units
        unit_length_in_cgs = hf["/Units"].attrs["Unit length in cgs (U_L)"]
        unit_time_in_cgs = hf["/Units"].attrs["Unit time in cgs (U_t)"]
        vel *= unit_length_in_cgs / unit_time_in_cgs  # cm/s
        vel *= 1e-5  # km/s

    snapshot_file = h5py.File(siminfo.latest_snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")
    properties_file = h5py.File(siminfo.halo_properties, "r")

    c200c = properties_file["cNFW_200crit"][:]
    c200c[c200c == 0] = 1
    m200c = properties_file["Mass_200crit"][:] * 1e10
    R200c = properties_file["R_200crit"][:] * 1e3 * siminfo.a #kpc
    xCoP = properties_file["Xcminpot"][:] * siminfo.a
    yCoP = properties_file["Ycminpot"][:] * siminfo.a
    zCoP = properties_file["Zcminpot"][:] * siminfo.a
    m200c = np.log10(m200c)
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP

    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    density = np.zeros((len(centers), 3))
    sig_density = np.zeros((len(centers), 3))
    velocity = np.zeros((len(centers), 3))
    sig_velocity = np.zeros((len(centers), 3))
    rs = np.zeros(3)
    M200 = np.zeros(3)

    for i in range(0, 3):
        if i == 0: select_halos = np.where((m200c >= 9.9) & (m200c <= 10.1))[0]  # >10 star parts
        if i == 1: select_halos = np.where((m200c >= 10.9) & (m200c <= 11.1))[0]  # >10 star parts
        if i == 2: select_halos = np.where((m200c >= 11.9) & (m200c <= 12.1))[0]  # >10 star parts

        if len(select_halos) >= 20:
            select_random = np.random.random_integers(len(select_halos)-1, size=(20))
            select_halos = select_halos[select_random]

        rs[i] = np.median(R200c[select_halos] / c200c[select_halos]) # kpc
        M200[i] = np.median(10 ** m200c[select_halos])

        num_halos = len(select_halos)
        density_all = np.zeros((len(centers), num_halos))
        velocity_all = np.zeros((len(centers), num_halos))

        for halo in range(0, num_halos):

            halo_j = select_halos[halo]

            # Grab the start position in the particles file to read from
            halo_start_position = group_file["Offset"][halo_j]
            halo_end_position = group_file["Offset"][halo_j + 1]
            particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]
            particle_ids_from_snapshot = snapshot_file["PartType1/ParticleIDs"][...]

            _, indices_v, indices_p = np.intersect1d(particle_ids_in_halo,
                                                     particle_ids_from_snapshot,
                                                     assume_unique=True,
                                                     return_indices=True, )

            particles_mass = mass[indices_p].copy()
            particles_pos = pos[indices_p, :].copy()
            particles_pos -= CoP[halo_j, :]  # centering
            particles_pos *= 1e3  # kpc
            particles_vel = vel[indices_p, :].copy()

            density_halo, velocity_halo = analyse_halo(particles_mass, particles_pos, particles_vel)
            density_all[:, halo] = density_halo
            velocity_all[:, halo] = velocity_halo

        density[:, i] = np.mean(density_all[:, :], axis=1)
        sig_density[:, i] = np.std(density_all[:, :], axis=1)
        velocity[:, i] = np.mean(velocity_all[:, :], axis=1)
        sig_velocity[:, i] = np.std(velocity_all[:, :], axis=1)

    return rs, M200, density, sig_density, velocity, sig_velocity



def plot_halo_profiles(siminfo,output_path):

    rs, M200, density, sig_density, velocity, sig_velocity = read_data(siminfo)

    # Plot parameters
    params = {
        "font.size": 14,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (7, 8),
        "figure.subplot.left": 0.12,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.1,
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

    z = 0  # Redshift
    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    for i in range(3):

        j = i*2 + 1

        ax = plt.subplot(3, 2, j)
        plt.grid("True")

        if i==0:
            label = '$10^{11}$M$_{\odot}$ halo'
            color = 'tab:green'
        if i==1:
            label = '$10^{12}$M$_{\odot}$ halo'
            color = 'tab:orange'
        if i==2:
            label = '$10^{13}$M$_{\odot}$ halo'
            color = 'tab:blue'

        NFWsig1D = sigma_1D(centers, M200[i], rs[i], z, siminfo)  # km/s
        NFWrho = calc_density(centers, M200[i], rs[i], z, siminfo)  # Msun/kpc^3

        plt.plot(centers, NFWrho, lw=1, color='black', label="NFW profile")
        plt.plot(centers, density[:, i], lw=2, color=color)
        plt.fill_between(centers, density[:, i] - sig_density[:, i] / 2,
                         density[:, i] + sig_density[:, i] / 2, alpha=0.4,
                         color=color, label=label)

        xarray = np.array([2.3, 2.3])
        yarray = np.array([1e3, 1e9])
        plt.plot(xarray, yarray, '--', color='grey')
        plt.ylim(1e3, 1e9)
        plt.xlim(1, 1e3)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Radius [kpc]")
        plt.ylabel("Density [M$_{\odot}$/kpc$^{3}$]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

        ######################
        ax = plt.subplot(3, 2, j+1)
        plt.grid("True")

        plt.plot(centers, NFWsig1D, lw=1, color='black')
        plt.plot(centers, velocity[:, i], lw=2, color=color)
        plt.fill_between(centers, velocity[:, i] - sig_velocity[:, i] / 2,
                         velocity[:, i] + sig_velocity[:, i] / 2,
                         alpha=0.4, color=color)

        #plt.ylim(0, 80)
        plt.xlim(1, 1e3)
        plt.xscale("log")
        plt.xlabel("Radius [kpc]")
        plt.ylabel("Velocity dispersion [km/s]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    plt.savefig(output_path + "Density_profiles.png", dpi=200)
    plt.close()