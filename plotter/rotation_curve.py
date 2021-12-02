import numpy as np
import h5py
from scipy import interpolate
import scipy.stats as stat
from pylab import *
from tqdm import tqdm
from object import particle_data
import matplotlib.pyplot as plt

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


def calc_NFW(x, log10_M, c, siminfo):
    z = siminfo.z
    Om = siminfo.Omega_m
    Ol = siminfo.Omega_l
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    a = R200 / c

    # NFW profile #
    delta = 200. / 3.
    delta *= c ** 3 / (np.log(1. + c) - c / (1. + c))
    f = np.zeros(len(x))
    for i in range(0, len(x)): f[i] = rho_crit * delta * 1e-9 / ((x[i] / a) * (1. + x[i] / a) ** 2)
    return f

def calc_NFW_vcirc(r, log10_M, c, siminfo):

    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21
    G = 6.67408e-11 # m^3 kg^-1 s^-2
    G *= (1e2)**3/1e3 # cm^3 g^-1 s^-2
    G /= kpc_in_cgs**3
    G *= Msun_in_cgs # kpc^3 Msun^-1 s^-2


    z = siminfo.z
    Om = siminfo.Omega_m
    Ol = siminfo.Omega_l
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    x = r / R200
    gcx = np.log(1. + c * x) - c * x / (1. + c * x)
    gc = np.log(1. + c) - c / (1. + c)
    mass = 200 * 4 * np.pi * R200**3 * rho_crit * 1e-9 * gcx / (3. * gc) # Msun

    Vcirc = G * mass / r # kpc^2/s^2
    Vcirc *= (kpc_in_cgs/1e5)**2 # km^2/s^2
    Vcirc = np.sqrt(Vcirc) # km/s
    return Vcirc

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


def calculate_Vcirc(mass, pos, radial_bins):
    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21
    G = 6.67408e-11 #m^3 kg^-1 s^-2
    G *= (1e2)**3/1e3 #cm^3 g^-1 s^-2
    G /= kpc_in_cgs**3
    G *= Msun_in_cgs #kpc^3 Msun^-1 s^-2

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))
    centers = bin_centers(radial_bins)  # kpc

    SumMasses, _, _ = stat.binned_statistic(x=r, values=mass, statistic="sum", bins=radial_bins, )
    M_within_r = np.cumsum(SumMasses)
    V_circ_2 = G * M_within_r / centers #kpc^2 s^-2
    V_circ_2 *= (kpc_in_cgs/1e5)**2 #km^2/s^2
    V_circ = np.sqrt(V_circ_2)

    return V_circ

def calculate_density(mass, radius):

    radial_bins = np.arange(-1, 2, 0.1)
    radial_bins = 10 ** radial_bins

    SumMasses, _, _ = stat.binned_statistic(x=radius, values= mass, statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3
    return density


def make_rotation_curve_data(sim_info, log10_min_mass, log10_max_mass):

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0.2, 25, 0.25)
    centers = bin_centers(radial_bins) # kpc

    sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= log10_min_mass) &
        (sim_info.halo_data.log10_halo_mass <= log10_max_mass)
        )[0]

    M200c = sim_info.halo_data.log10_halo_mass[sample]
    c200c = sim_info.halo_data.concentration[sample]
    structure_type = sim_info.halo_data.structure_type[sample]
    index = sim_info.halo_data.halo_index[sample]

    num_halos = len(sample)
    v_circ_1 = np.zeros(num_halos)
    v_circ_2 = np.zeros(num_halos)
    v_circ_5 = np.zeros(num_halos)
    v_circ_nfw_1 = np.zeros(num_halos)
    v_circ_nfw_2 = np.zeros(num_halos)
    v_circ_nfw_5 = np.zeros(num_halos)
    Vmax = np.zeros(num_halos)

    for i in tqdm(range(num_halos)):

        halo_indx = sim_info.halo_data.halo_index[sample[i]]
        part_data = particle_data.load_particle_data(sim_info, halo_indx, sample[i])

        circular_velocity = calculate_Vcirc(part_data.masses.value,
                                       part_data.coordinates.value,
                                       radial_bins)

        f = interpolate.interp1d(centers, circular_velocity)
        v_circ_1[i] = f(1.0)
        v_circ_2[i] = f(2.0)
        v_circ_5[i] = f(5.0)

        v_nfw = calc_NFW_vcirc(centers, M200c[i], c200c[i], sim_info)
        f = interpolate.interp1d(centers, v_nfw)
        v_circ_nfw_1[i] = f(1.0)
        v_circ_nfw_2[i] = f(2.0)
        v_circ_nfw_5[i] = f(5.0)
        Vmax[i] = np.max(circular_velocity)

    # Output data
    filename = f"{sim_info.output_path}/Rotation_data_" + sim_info.simulation_name + ".hdf5"
    data_file = h5py.File(filename, 'w')
    f = data_file.create_group('Data')
    MH = f.create_dataset('ID', data=index)
    MH = f.create_dataset('StructureType', data=structure_type)
    MH = f.create_dataset('M200c', data=M200c)
    MH = f.create_dataset('c200c', data=c200c)
    MH = f.create_dataset('Vcirc1kpc', data=v_circ_1)
    MH = f.create_dataset('Vcirc2kpc', data=v_circ_2)
    MH = f.create_dataset('Vcirc5kpc', data=v_circ_5)
    MH = f.create_dataset('Vcirc_nfw_1kpc', data=v_circ_nfw_1)
    MH = f.create_dataset('Vcirc_nfw_2kpc', data=v_circ_nfw_2)
    MH = f.create_dataset('Vcirc_nfw_5kpc', data=v_circ_nfw_5)
    MH = f.create_dataset('Vmax', data=Vmax)
    data_file.close()
    return

def plot_rotation_curve(sim_info, log10_min_mass, log10_max_mass, structure_type):

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0.2, 25, 0.25)
    centers = bin_centers(radial_bins) # kpc

    select_sub_sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= log10_min_mass) &
        (sim_info.halo_data.log10_halo_mass <= log10_max_mass))[0]

    if structure_type == 10:
        select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] == structure_type)[0]
    else:
        select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] > 10)[0]

    sample = select_sub_sample[select_type]
    M200c = np.median(sim_info.halo_data.log10_halo_mass[sample])
    c200c = np.median(sim_info.halo_data.concentration[sample])

    if len(sample) >= 30:
        select_random = np.random.random_integers(len(sample) - 1, size=(30))
        sample = sample[select_random]

    num_halos = len(sample)
    v_circ = np.zeros((len(centers), num_halos))

    for i in tqdm(range(num_halos)):

        halo_indx = sim_info.halo_data.halo_index[sample[i]]
        part_data = particle_data.load_particle_data(sim_info, halo_indx, sample[i])

        circular_velocity = calculate_Vcirc(part_data.masses.value,
                                            part_data.coordinates.value,
                                            radial_bins)

        v_circ[:, i] = circular_velocity


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

    for i in range(len(sample)):
        plt.plot(centers, v_circ[:, i], '-')

    NFW_circ = calc_NFW_vcirc(centers, M200c, c200c, sim_info)
    plt.plot(centers, NFW_circ,'--',lw=2,color='black')

    plt.axis([0, 20, 0, 100])
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Circular velocity [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if structure_type == 10:
        plt.savefig(f"{sim_info.output_path}/Rotation_curve_"+sim_info.simulation_name+"_centrals.png", dpi=200)
    else:
        plt.savefig(f"{sim_info.output_path}/Rotation_curve_"+sim_info.simulation_name+"_satellites.png", dpi=200)
    plt.close()

def plot_rotation_curve_data(sim_info, output_name_list):

    for name in output_name_list:
        filename = f"{sim_info.output_path}/Rotation_data_" + name + ".hdf5"
        with h5py.File(filename, "r") as file:
            M200c = file["Data/M200c"][:]
            Type = file["Data/StructureType"][:]
            v_circ_1 = file["Data/Vcirc1kpc"][:]
            v_circ_2 = file["Data/Vcirc2kpc"][:]
            v_circ_5 = file["Data/Vcirc5kpc"][:]
            v_circ_nfw_1 = file["Data/Vcirc_nfw_1kpc"][:]
            v_circ_nfw_2 = file["Data/Vcirc_nfw_2kpc"][:]
            v_circ_nfw_5 = file["Data/Vcirc_nfw_5kpc"][:]
            Vmax = file["Data/Vmax"][:]

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
        "figure.subplot.wspace": 0.45,
        "figure.subplot.hspace": 0.35,
        "lines.markersize": 1,
        "lines.linewidth": 1.5,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    ratio = (v_circ_1 - v_circ_nfw_1) / v_circ_nfw_1
    cen = Type == 10
    sat = Type > 10
    plt.plot(M200c[cen], ratio[cen], 'o',color='tab:orange',label='Centrals')
    plt.plot(M200c[sat], ratio[sat], 'o',color='tab:blue',label='Satellites')

    x_range = np.arange(9, 11, 0.05)
    y_cen = stat.binned_statistic(x=M200c[cen], values=ratio[cen], statistic="median", bins=x_range, )[0]
    y_sat = stat.binned_statistic(x=M200c[sat], values=ratio[sat], statistic="median", bins=x_range, )[0]
    x_center = bin_centers(x_range)
    plt.plot(x_center, np.zeros(len(x_center)), '--', lw=1, color='black')
    plt.plot(x_center, y_cen, '-',lw=2, color='white')
    plt.plot(x_center, y_sat, '-',lw=2, color='white')
    plt.plot(x_center, y_cen, '-', color='tab:orange')
    plt.plot(x_center, y_sat, '-', color='tab:blue')

    plt.axis([9, 11, -1, 3])
    plt.xlabel("$\log_{10}$ M [M$_{\odot}$]")
    plt.ylabel("(V$_{\mathrm{c}}$(1kpc)-V$_{\mathrm{c,NFW}}$)/V$_{\mathrm{c,NFW}}$")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/M200_Vcirc_ratio_1_kpc_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    ratio = (v_circ_2 - v_circ_nfw_2) /v_circ_nfw_2
    cen = Type == 10
    sat = Type > 10
    plt.plot(M200c[cen], ratio[cen], 'o',color='tab:orange',label='Centrals')
    plt.plot(M200c[sat], ratio[sat], 'o',color='tab:blue',label='Satellites')

    x_range = np.arange(9, 11, 0.05)
    y_cen = stat.binned_statistic(x=M200c[cen], values=ratio[cen], statistic="median", bins=x_range, )[0]
    y_sat = stat.binned_statistic(x=M200c[sat], values=ratio[sat], statistic="median", bins=x_range, )[0]
    x_center = bin_centers(x_range)
    plt.plot(x_center, np.zeros(len(x_center)), '--', lw=1, color='black')
    plt.plot(x_center, y_cen, '-',lw=2, color='white')
    plt.plot(x_center, y_sat, '-',lw=2, color='white')
    plt.plot(x_center, y_cen, '-', color='tab:orange')
    plt.plot(x_center, y_sat, '-', color='tab:blue')

    plt.axis([9, 11, -1, 3])
    plt.xlabel("$\log_{10}$ M [M$_{\odot}$]")
    plt.ylabel("(V$_{\mathrm{c}}$(2kpc)-V$_{\mathrm{c,NFW}}$)/V$_{\mathrm{c,NFW}}$")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/M200_Vcirc_ratio_2_kpc_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    ratio = (v_circ_5 - v_circ_nfw_5) / v_circ_nfw_5
    cen = Type == 10
    sat = Type > 10
    plt.plot(M200c[cen], ratio[cen], 'o',color='tab:orange',label='Centrals')
    plt.plot(M200c[sat], ratio[sat], 'o',color='tab:blue',label='Satellites')

    x_range = np.arange(9, 11, 0.05)
    y_cen = stat.binned_statistic(x=M200c[cen], values=ratio[cen], statistic="median", bins=x_range, )[0]
    y_sat = stat.binned_statistic(x=M200c[sat], values=ratio[sat], statistic="median", bins=x_range, )[0]
    x_center = bin_centers(x_range)
    plt.plot(x_center, np.zeros(len(x_center)), '--', lw=1, color='black')
    plt.plot(x_center, y_cen, '-',lw=2, color='white')
    plt.plot(x_center, y_sat, '-',lw=2, color='white')
    plt.plot(x_center, y_cen, '-', color='tab:orange')
    plt.plot(x_center, y_sat, '-', color='tab:blue')

    plt.axis([9, 11, -1, 3])
    plt.xlabel("$\log_{10}$ M [M$_{\odot}$]")
    plt.ylabel("(V$_{\mathrm{c}}$(5kpc)-V$_{\mathrm{c,NFW}}$)/V$_{\mathrm{c,NFW}}$")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/M200_Vcirc_ratio_5_kpc_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    # figure()
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plt.plot(yrange, yrange, '--', lw=1, color='grey')
    # plt.plot(NFW_circ_max, NFW_circ_2kpc, '-', lw=2, color='black')
    #
    # v_circ_median = stat.binned_statistic(x=Vmax, values=v_circ_2, statistic="median", bins=Vmax_range, )[0]
    # v_circ_16 = stat.binned_statistic(x=Vmax, values=v_circ_2,
    #                                   statistic=lambda v_circ_2: percentile(v_circ_2, 16), bins=Vmax_range, )[0]
    # v_circ_84 = stat.binned_statistic(x=Vmax, values=v_circ_2,
    #                                   statistic=lambda v_circ_2: percentile(v_circ_2, 84), bins=Vmax_range, )[0]
    # plt.plot(Vmax_center, v_circ_median, '-')
    # plt.fill_between(Vmax_center, v_circ_16, v_circ_84, alpha=0.4)
    #
    # plt.plot(Vmax, v_circ_2, 'o')
    #
    # plt.axis([0, 100, 0, 100])
    # plt.xlabel("Maximum Circular Velocity [km/s]")
    # plt.ylabel("Circular velocity at 2kpc [km/s]")
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # if galaxy_type == "centrals":
    #     plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_" + siminfo.name + "_2kpc_centrals.png", dpi=200)
    # else:
    #     plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_" + siminfo.name + "_2kpc_satellites.png", dpi=200)
    # plt.close()
    #
    # figure()
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plt.plot(yrange, yrange, '--', lw=1, color='grey')
    # plt.plot(NFW_circ_max, NFW_circ_5kpc, '-', lw=2, color='black')
    #
    # v_circ_median = stat.binned_statistic(x=Vmax, values=v_circ_5, statistic="median", bins=Vmax_range, )[0]
    # v_circ_16 = stat.binned_statistic(x=Vmax, values=v_circ_5,
    #                                   statistic=lambda v_circ_5: percentile(v_circ_5, 16), bins=Vmax_range, )[0]
    # v_circ_84 = stat.binned_statistic(x=Vmax, values=v_circ_5,
    #                                   statistic=lambda v_circ_5: percentile(v_circ_5, 84), bins=Vmax_range, )[0]
    # plt.plot(Vmax_center, v_circ_median, '-')
    # plt.fill_between(Vmax_center, v_circ_16, v_circ_84, alpha=0.4)
    #
    # plt.plot(Vmax, v_circ_5, 'o')
    #
    # plt.axis([0, 100, 0, 100])
    # plt.xlabel("Maximum Circular Velocity [km/s]")
    # plt.ylabel("Circular velocity at 5kpc [km/s]")
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # if galaxy_type == "centrals":
    #     plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_" + siminfo.name + "_5kpc_centrals.png", dpi=200)
    # else:
    #     plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_" + siminfo.name + "_5kpc_satellites.png", dpi=200)
    # plt.close()
    #
    # ###################
    #
    # figure()
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plt.plot(M200c, v_circ_1, 'o')
    # v_median, _, _ = stat.binned_statistic(x=M200c, values=v_circ_1, statistic="median", bins=mass_range, )
    # mass = bin_centers(mass_range)
    # plt.plot(mass, v_median, '-', lw=2, color='tab:blue')
    # plt.plot(mass_range, NFW_circ_1kpc, '--', lw=2, color='black')
    #
    # plt.axis([8, 11, 0, 50])
    # # plt.xscale('log')
    # # plt.yscale('log')
    # plt.xlabel("$\log_{10}$ M$_{200}$ [M$_{\odot}$]")
    # plt.ylabel("Circular velocity at 1kpc [km/s]")
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # if galaxy_type == "centrals":
    #     plt.savefig(f"{siminfo.output_path}/M200_Vcirc_" + siminfo.name + "_1kpc_centrals.png", dpi=200)
    # else:
    #     plt.savefig(f"{siminfo.output_path}/M200_Vcirc_" + siminfo.name + "_1kpc_satellites.png", dpi=200)
    # plt.close()
    #
    # figure()
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plt.plot(M200c, v_circ_2, 'o')
    # v_median, _, _ = stat.binned_statistic(x=M200c, values=v_circ_2, statistic="median", bins=mass_range, )
    # plt.plot(mass, v_median, '-', lw=2, color='tab:blue')
    # plt.plot(mass_range, NFW_circ_2kpc, '--', lw=2, color='black')
    #
    # plt.axis([8, 11, 0, 50])
    # # plt.xscale('log')
    # # plt.yscale('log')
    # plt.xlabel("$\log_{10}$ M$_{200}$ [M$_{\odot}$]")
    # plt.ylabel("Circular velocity at 2kpc [km/s]")
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # if galaxy_type == "centrals":
    #     plt.savefig(f"{siminfo.output_path}/M200_Vcirc_" + siminfo.name + "_2kpc_centrals.png", dpi=200)
    # else:
    #     plt.savefig(f"{siminfo.output_path}/M200_Vcirc_" + siminfo.name + "_2kpc_satellites.png", dpi=200)
    # plt.close()
    #
    # figure()
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # plt.plot(M200c, v_circ_5, 'o')
    # v_median, _, _ = stat.binned_statistic(x=M200c, values=v_circ_5, statistic="median", bins=mass_range, )
    # plt.plot(mass, v_median, '-', lw=2, color='tab:blue')
    #
    # plt.plot(mass_range, NFW_circ_5kpc, '--', lw=2, color='black')
    #
    # plt.axis([8, 11, 10, 100])
    # # plt.xscale('log')
    # # plt.yscale('log')
    # plt.xlabel("$\log_{10}$ M$_{200}$ [M$_{\odot}$]")
    # plt.ylabel("Circular velocity at 5kpc [km/s]")
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # if galaxy_type == "centrals":
    #     plt.savefig(f"{siminfo.output_path}/M200_Vcirc_" + siminfo.name + "_5kpc_centrals.png", dpi=200)
    # else:
    #     plt.savefig(f"{siminfo.output_path}/M200_Vcirc_" + siminfo.name + "_5kpc_satellites.png", dpi=200)
    # plt.close()