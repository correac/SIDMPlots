import numpy as np
import h5py
from scipy import interpolate
import scipy.stats as stat
from scipy.optimize import curve_fit
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

def NFW_curve_fiducial(M200c, siminfo):

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0.0, 26, 0.25)
    r = bin_centers(radial_bins)  # kpc

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
    R200 = 10**M200c / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    c = c_M_relation(M200c)

    Vmax = np.zeros(len(R200))
    V_fid = np.zeros(len(R200))
    for i in range(len(R200)):

        x = r / R200[i]
        gcx = np.log(1. + c[i] * x) - c[i] * x / (1. + c[i] * x)
        gc = np.log(1. + c[i]) - c[i] / (1. + c[i])
        mass = 200 * 4 * np.pi * R200[i]**3 * rho_crit * 1e-9 * gcx / (3. * gc) # Msun

        Vcirc = G * mass / r # kpc^2/s^2
        Vcirc *= (kpc_in_cgs/1e5)**2 # km^2/s^2
        Vcirc = np.sqrt(Vcirc) # km/s
        Vmax[i] = np.max(Vcirc)

        r_fid = Vmax[i] / 10.
        if r_fid < 2.: r_fid = 2.
        if r_fid > 25.: r_fid = 25.

        f = interpolate.interp1d(r, Vcirc)
        V_fid[i] = f(r_fid)

    return Vmax, V_fid


def NFW_curve(M200c, siminfo):

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0.2, 25, 0.25)
    r = bin_centers(radial_bins)  # kpc

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
    R200 = 10**M200c / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    c = c_M_relation(M200c)

    Vmax = np.zeros(len(R200))
    Vcirc_2kpc = np.zeros(len(R200))
    for i in range(len(R200)):

        x = r / R200[i]
        gcx = np.log(1. + c[i] * x) - c[i] * x / (1. + c[i] * x)
        gc = np.log(1. + c[i]) - c[i] / (1. + c[i])
        mass = 200 * 4 * np.pi * R200[i]**3 * rho_crit * 1e-9 * gcx / (3. * gc) # Msun

        Vcirc = G * mass / r # kpc^2/s^2
        Vcirc *= (kpc_in_cgs/1e5)**2 # km^2/s^2
        Vcirc = np.sqrt(Vcirc) # km/s
        Vmax[i] = np.max(Vcirc)

        f = interpolate.interp1d(r, Vcirc)
        Vcirc_2kpc[i] = f(2.0)

    return Vmax, Vcirc_2kpc


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
    radial_bins = np.arange(0.2, 26, 0.25)
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
    v_fid = np.zeros(num_halos)
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

        Vmax[i] = np.max(circular_velocity)

        r_fid = Vmax[i] / 10.
        if r_fid < 2.: r_fid = 2.
        if r_fid > 25.: r_fid = 25.
        v_fid[i] = f(r_fid)

        v_nfw = calc_NFW_vcirc(centers, M200c[i], c200c[i], sim_info)
        f = interpolate.interp1d(centers, v_nfw)
        v_circ_nfw_1[i] = f(1.0)
        v_circ_nfw_2[i] = f(2.0)
        v_circ_nfw_5[i] = f(5.0)

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
    MH = f.create_dataset('V_fiducial', data=v_fid)
    MH = f.create_dataset('Vmax', data=Vmax)
    data_file.close()
    return

def plot_rotation_curve(sim_info, log10_min_mass, log10_max_mass, structure_type):

    mean_mass = 0.5 * (log10_max_mass + log10_min_mass)

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

    plt.xlim([0, 20])
    #plt.axis([0, 20, 0, 100])
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Circular velocity [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    output_file = f"{sim_info.output_path}/Rotation_curve_"+sim_info.simulation_name
    if structure_type == 10:
        plt.savefig(output_file+"_mass_%.1f_"%mean_mass+"centrals.png", dpi=200)
    else:
        plt.savefig(output_file+"_mass_%.1f_"%mean_mass+"satellites.png", dpi=200)
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
            v_fid = file["Data/V_fiducial"][:]
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

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    cen = Type == 10
    sat = Type > 10
    plt.plot(Vmax[cen], v_circ_2[cen], 'o', color='tab:orange', label='Centrals')
    plt.plot(Vmax[sat], v_circ_2[sat], 'o', color='tab:blue', label='Satellites')

    NFW_M200c = np.arange(8,12.1,0.1)
    NFW_Vmax, NFW_Vcirc_2 = NFW_curve(NFW_M200c, sim_info)
    plt.plot(NFW_Vmax, NFW_Vcirc_2,'-', lw=1, color='black',label='NFW')
    x_range = np.arange(0, 85, 5)
    plt.plot(x_range, x_range,'--', lw=1, color='grey')

    plt.axis([10, 80, 10, 80])
    plt.xlabel("V$_{\mathrm{max}}$ [kpc]")
    plt.ylabel("V$_{\mathrm{c}}$(2kpc) [kpc]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/Vmax_Vcirc_2_kpc_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    cen = Type == 10
    sat = Type > 10
    plt.plot(Vmax[cen], v_fid[cen], 'o', color='tab:orange', label='Centrals')
    plt.plot(Vmax[sat], v_fid[sat], 'o', color='tab:blue', label='Satellites')

    NFW_Vmax, NFW_fid = NFW_curve_fiducial(NFW_M200c, sim_info)
    plt.plot(NFW_Vmax, NFW_fid, '-', lw=1, color='black', label='NFW')
    x_range = np.arange(5, 85, 5)
    plt.plot(x_range, x_range, '--', lw=1, color='grey')

    plt.axis([0, 80, 0, 80])
    plt.xlabel("V$_{\mathrm{max}}$ [km/s]")
    plt.ylabel("V$_{\mathrm{fiducial}}$ [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/Vmax_V_fiducial_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    ratio = np.zeros(len(M200c))
    for i in range(len(M200c)):
        f = interpolate.interp1d(NFW_M200c, NFW_fid)
        ratio[i] = v_fid[i] / f(M200c[i])

    cen = Type == 10
    sat = Type > 10
    plt.plot(M200c[cen], ratio[cen], 'o', color='tab:orange', label='Centrals')
    plt.plot(M200c[sat], ratio[sat], 'o', color='tab:blue', label='Satellites')

    plt.axis([9, 11, 0.4, 2.5])
    plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
    plt.ylabel("$V_{\mathrm{fid}}/V_{\mathrm{fid-NFW}}$")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/M200c_Vfid_Vmax_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    cen = Type == 10
    sat = Type > 10
    plt.plot(M200c[cen], v_circ_2[cen], 'o', color='tab:orange', label='Centrals')
    plt.plot(M200c[sat], v_circ_2[sat], 'o', color='tab:blue', label='Satellites')

    plt.plot(NFW_M200c, NFW_Vcirc_2, '-', lw=1, color='black', label='NFW')

    plt.axis([9, 11, 0, 80])
    plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
    plt.ylabel("V$_{\mathrm{c}}$(2kpc) [kpc]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc='upper left',labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/M200c_Vcirc_2_kpc_" + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

def function(x,a,b,c,d):
    f = a + b*x + c*x**2 + d*x**3
    return f

def plot_rotation_relative_to_CDM(sim_info, output_name_list):

    name = 'DML025N752SigmaConstant00'
    filename = f"{sim_info.output_path}/Rotation_data_" + name + ".hdf5"
    with h5py.File(filename, "r") as file:
        CDM_M200c = file["Data/M200c"][:]
        CDM_Type = file["Data/StructureType"][:]
        CDM_v_fid = file["Data/V_fiducial"][:]

    cen = np.where(CDM_Type == 10)[0]
    sat = np.where(CDM_Type > 10)[0]

    x_range = np.arange(8.9, 12.1, 0.2)
    ydata_50 = stat.binned_statistic(x=CDM_M200c[cen], values=CDM_v_fid[cen], statistic="median", bins=x_range, )[0]
    ydata_16 = stat.binned_statistic(x=CDM_M200c[cen], values=CDM_v_fid[cen],
                                     statistic=lambda y: np.percentile(y, 16), bins=x_range, )[0]
    ydata_84 = stat.binned_statistic(x=CDM_M200c[cen], values=CDM_v_fid[cen],
                                     statistic=lambda y: np.percentile(y, 84), bins=x_range, )[0]
    ydata_50_cen = ydata_50.copy()
    xdata = bin_centers(x_range)

    no_nan = np.isnan(ydata_50) == False
    popt_50, pcov = curve_fit(function, xdata[no_nan]-9.0, ydata_50[no_nan], p0=[0.5, 1, 2, 1])
    popt_16, pcov = curve_fit(function, xdata[no_nan]-9.0, ydata_16[no_nan], p0=[0.5, 1, 2, 1])
    popt_84, pcov = curve_fit(function, xdata[no_nan]-9.0, ydata_84[no_nan], p0=[0.5, 1, 2, 1])

    ydata_50 = stat.binned_statistic(x=CDM_M200c[sat], values=CDM_v_fid[sat], statistic="median", bins=x_range, )[0]
    ydata_16 = stat.binned_statistic(x=CDM_M200c[sat], values=CDM_v_fid[sat],
                                     statistic=lambda y: np.percentile(y, 16), bins=x_range, )[0]
    ydata_84 = stat.binned_statistic(x=CDM_M200c[sat], values=CDM_v_fid[sat],
                                     statistic=lambda y: np.percentile(y, 84), bins=x_range, )[0]
    ydata_50_sat = ydata_50.copy()

    no_nan = np.isnan(ydata_50) == False
    popt_50_sat, pcov = curve_fit(function, xdata[no_nan]-9.0, ydata_50[no_nan], p0=[0.5, 1, 2, 1])
    popt_16_sat, pcov = curve_fit(function, xdata[no_nan]-9.0, ydata_16[no_nan], p0=[0.5, 1, 2, 1])
    popt_84_sat, pcov = curve_fit(function, xdata[no_nan]-9.0, ydata_84[no_nan], p0=[0.5, 1, 2, 1])


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
        "lines.markersize": 0.5,
        "lines.linewidth": 1.5,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(CDM_M200c[cen], CDM_v_fid[cen], 'o', color='tab:orange', label='Centrals')
    plt.plot(CDM_M200c[sat], CDM_v_fid[sat], 'o', color='tab:blue', label='Satellites')

    plt.fill_between(xdata, function(xdata-9.0,*popt_16), function(xdata-9.0,*popt_84), alpha=0.2, color='tab:orange')
    plt.plot(xdata, function(xdata-9.0,*popt_50),'-',lw=2,color='white')
    plt.plot(xdata, function(xdata-9.0,*popt_50),'-',lw=1,color='tab:orange')
    plt.plot(xdata, ydata_50_cen, '-', lw=1, color='tab:red')

    plt.fill_between(xdata, function(xdata-9.0,*popt_16_sat), function(xdata-9.0,*popt_84_sat), alpha=0.2, color='tab:blue')
    plt.plot(xdata, function(xdata-9.0,*popt_50_sat),'-',lw=2,color='white')
    plt.plot(xdata, function(xdata-9.0,*popt_50_sat),'-',lw=1,color='tab:blue')
    plt.plot(xdata, ydata_50_sat, '-', lw=1, color='tab:purple')

    plt.axis([9, 12, 1, 500])
    plt.yscale('log')
    plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
    plt.ylabel("$V_{\mathrm{fid}}$ [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{sim_info.output_path}/M200c_Vfid_test_CDM.png", dpi=200)
    plt.close()

    #######
    # Let's look at CDM scatter
    ratio = np.zeros(len(cen))
    for i in range(len(cen)):
        ratio[i] = CDM_v_fid[cen[i]] / function(CDM_M200c[cen[i]] - 9.0, *popt_50)

    ydata_1_cen = stat.binned_statistic(x=CDM_M200c[cen], values=ratio,
                                     statistic=lambda y: np.percentile(y, 1), bins=x_range, )[0]
    ydata_99_cen = stat.binned_statistic(x=CDM_M200c[cen], values=ratio,
                                     statistic=lambda y: np.percentile(y, 99), bins=x_range, )[0]
    ydata_10_cen = stat.binned_statistic(x=CDM_M200c[cen], values=ratio,
                                     statistic=lambda y: np.percentile(y, 5), bins=x_range, )[0]
    ydata_90_cen = stat.binned_statistic(x=CDM_M200c[cen], values=ratio,
                                     statistic=lambda y: np.percentile(y, 95), bins=x_range, )[0]

    ratio = np.zeros(len(sat))
    for i in range(len(CDM_M200c[sat])):
        ratio[i] = CDM_v_fid[sat[i]] / function(CDM_M200c[sat[i]] - 9.0, *popt_50_sat)

    ydata_1_sat = stat.binned_statistic(x=CDM_M200c[sat], values=ratio,
                                     statistic=lambda y: np.percentile(y, 1), bins=x_range, )[0]
    ydata_99_sat = stat.binned_statistic(x=CDM_M200c[sat], values=ratio,
                                     statistic=lambda y: np.percentile(y, 99), bins=x_range, )[0]
    ydata_10_sat = stat.binned_statistic(x=CDM_M200c[sat], values=ratio,
                                     statistic=lambda y: np.percentile(y, 5), bins=x_range, )[0]
    ydata_90_sat = stat.binned_statistic(x=CDM_M200c[sat], values=ratio,
                                     statistic=lambda y: np.percentile(y, 95), bins=x_range, )[0]

    #######
    for name in output_name_list:
        filename = f"{sim_info.output_path}/Rotation_data_" + name + ".hdf5"
        with h5py.File(filename, "r") as file:
            M200c = file["Data/M200c"][:]
            Type = file["Data/StructureType"][:]
            v_fid = file["Data/V_fiducial"][:]

        cen = np.where(Type == 10)[0]
        sat = np.where(Type > 10)[0]

        figure()
        ax = plt.subplot(1, 1, 1)
        plt.grid("True")

        plt.plot(M200c[cen], v_fid[cen], 'o', color='tab:orange', label='Centrals')
        plt.plot(M200c[sat], v_fid[sat], 'o', color='tab:blue', label='Satellites')

        plt.fill_between(xdata, function(xdata - 9.0, *popt_16), function(xdata - 9.0, *popt_84), alpha=0.2,
                         color='tab:orange')
        plt.plot(xdata, function(xdata - 9.0, *popt_50), '-', lw=2, color='white')
        plt.plot(xdata, function(xdata - 9.0, *popt_50), '-', lw=1, color='tab:orange')
        plt.plot(xdata, ydata_50_cen, '-', lw=1, color='tab:red')

        plt.fill_between(xdata, function(xdata - 9.0, *popt_16_sat), function(xdata - 9.0, *popt_84_sat), alpha=0.2,
                         color='tab:blue')
        plt.plot(xdata, function(xdata - 9.0, *popt_50_sat), '-', lw=2, color='white')
        plt.plot(xdata, function(xdata - 9.0, *popt_50_sat), '-', lw=1, color='tab:blue')
        plt.plot(xdata, ydata_50_sat, '-', lw=1, color='tab:purple')

        plt.axis([9, 12, 1, 500])
        plt.yscale('log')
        plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
        plt.ylabel("$V_{\mathrm{fid}}$ [km/s]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
        plt.savefig(f"{sim_info.output_path}/M200c_Vfid_test_" + sim_info.simulation_name + ".png", dpi=200)
        plt.close()

        # Plot parameters
        params = {
            "font.size": 12,
            "font.family": "Times",
            "text.usetex": True,
            "figure.figsize": (7, 3.5),
            "figure.subplot.left": 0.1,
            "figure.subplot.right": 0.95,
            "figure.subplot.bottom": 0.18,
            "figure.subplot.top": 0.95,
            "figure.subplot.wspace": 0.25,
            "figure.subplot.hspace": 0.25,
            "lines.markersize": 0.5,
            "lines.linewidth": 1.5,
        }
        rcParams.update(params)

        figure()
        ax = plt.subplot(1, 2, 1)
        plt.grid("True")

        ratio = np.zeros(len(cen))
        for i in range(len(cen)):
            ratio[i] = v_fid[cen[i]] / function(M200c[cen[i]]-9.0,*popt_50)

        plt.plot(M200c[cen], ratio, 'o', color='tab:orange', label='Centrals')#,alpha=0.5)
        ydata = stat.binned_statistic(x=M200c[cen], values=ratio, statistic="median", bins=x_range, )[0]
        plt.plot(np.array([9,12]), np.array([1,1]), '-', lw=1, color='darkblue')
        plt.plot(xdata, ydata, '--', lw=1, color='black')

        #plt.fill_between(xdata, ydata_1_cen, ydata_99_cen, alpha=0.3, color='lightgrey', zorder=2)
        plt.fill_between(xdata, ydata_10_cen, ydata_90_cen, alpha=0.3, color='grey', zorder=2,label='CDM 5-95 percentiles')

        plt.axis([9, 11, 0, 2.5])
        plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
        plt.ylabel("$V_{\mathrm{fid}}/V_{\mathrm{fid-CDM}}$")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

        ax = plt.subplot(1, 2, 2)
        plt.grid("True")

        ratio = np.zeros(len(sat))
        for i in range(len(M200c[sat])):
            ratio[i] = v_fid[sat[i]] / function(M200c[sat[i]] - 9.0, *popt_50_sat)

        plt.plot(M200c[sat], ratio, 'o', color='tab:blue', label='Satellites')#, alpha=0.5)
        ydata = stat.binned_statistic(x=M200c[sat], values=ratio, statistic="median", bins=x_range, )[0]
        plt.plot(np.array([9,12]), np.array([1,1]), '-', lw=1, color='darkblue')
        plt.plot(xdata, ydata, '--', lw=1, color='black')

        #plt.fill_between(xdata, ydata_1_sat, ydata_99_sat, alpha=0.3, color='lightgrey', zorder=2)
        plt.fill_between(xdata, ydata_10_sat, ydata_90_sat, alpha=0.3, color='grey', zorder=2,label='CDM 5-95 percentiles')

        plt.axis([9, 11, 0, 2.5])
        plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
        plt.ylabel("$V_{\mathrm{fid}}/V_{\mathrm{fid-CDM}}$")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
        plt.savefig(f"{sim_info.output_path}/M200c_Vfid_Vfid_CDM_" + sim_info.simulation_name + ".png", dpi=200)
        plt.close()

        #Let's plot ratios

        figure()
        ax = plt.subplot(1, 2, 1)
        plt.grid("True")

        ratio = np.zeros(len(cen))
        for i in range(len(cen)):
            ratio[i] = (v_fid[cen[i]]-function(M200c[cen[i]] - 9.0, *popt_50)) / function(M200c[cen[i]] - 9.0, *popt_50)

        plt.plot(M200c[cen], ratio, 'o', color='tab:orange', label='Centrals')  # ,alpha=0.5)
        ydata = stat.binned_statistic(x=M200c[cen], values=ratio, statistic="median", bins=x_range, )[0]
        plt.plot(np.array([9, 12]), np.array([0, 0]), '-', lw=1, color='darkblue')
        plt.plot(xdata, ydata, '--', lw=1, color='black')

        #plt.fill_between(xdata, ydata_10_cen, ydata_90_cen, alpha=0.3, color='grey', zorder=2,
        #                 label='CDM 5-95 percentiles')

        plt.axis([9, 11, -1, 2])
        plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
        plt.ylabel("$(V_{\mathrm{fid}}-V_{\mathrm{fid-CDM}})/V_{\mathrm{fid-CDM}}$")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

        ax = plt.subplot(1, 2, 2)
        plt.grid("True")

        ratio = np.zeros(len(sat))
        for i in range(len(M200c[sat])):
            ratio[i] = (v_fid[sat[i]]-function(M200c[sat[i]] - 9.0, *popt_50_sat)) / function(M200c[sat[i]] - 9.0, *popt_50_sat)

        plt.plot(M200c[sat], ratio, 'o', color='tab:blue', label='Satellites')  # , alpha=0.5)
        ydata = stat.binned_statistic(x=M200c[sat], values=ratio, statistic="median", bins=x_range, )[0]
        plt.plot(np.array([9, 12]), np.array([0, 0]), '-', lw=1, color='darkblue')
        plt.plot(xdata, ydata, '--', lw=1, color='black')

        # plt.fill_between(xdata, ydata_1_sat, ydata_99_sat, alpha=0.3, color='lightgrey', zorder=2)
        #plt.fill_between(xdata, ydata_10_sat, ydata_90_sat, alpha=0.3, color='grey', zorder=2,
        #                 label='CDM 5-95 percentiles')

        plt.axis([9, 11, -1, 2])
        plt.xlabel("$\log_{10}M_{200c}$ [M$_{\odot}$]")
        plt.ylabel("$(V_{\mathrm{fid}}-V_{\mathrm{fid-CDM}})/V_{\mathrm{fid-CDM}}$")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc='upper right', labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
        plt.savefig(f"{sim_info.output_path}/M200c_ratio_Vfid_Vfid_CDM_" + sim_info.simulation_name + ".png", dpi=200)
        plt.close()