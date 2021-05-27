"""
Description
"""
import h5py
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scatter_rate import calc_density, sigma_1D
import scipy.stats as stat


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

def output_cM_vMax_relations(siminfo):

    # Load data:
    g = h5py.File(siminfo.halo_properties, "r")
    mass = g["Mass_200crit"][:] * 1e10  # convert to Msun
    c200 = g["cNFW_200crit"][:]
    Vmax = g["Vmax"][:]  # km/s
    subtype = g["Structuretype"][:]
    centrals = subtype == 10

    mass = mass[centrals]
    c200 = c200[centrals]
    Vmax = Vmax[centrals]

    np.savetxt(f"{siminfo.output_path}/cMVmax_" + siminfo.name + ".txt",
               np.transpose([mass, c200, Vmax]))

def velocity_dependence(x,w0):
    f = 2. * w0**4 / x**4
    f *= (2*np.log(1+0.5*x**2/w0**2)-np.log(1+x**2/w0**2))
    return f

def sigma(x,mx,mphi,alpha):
    w0 = 300. * (mphi/10.) * (10./mx)
    sigma0 = 274.85 * (alpha/0.01)**2 * (mx/10.) * (10./mphi)**4
    sigmam = sigma0
    y = velocity_dependence(x,w0)
    y *= sigmam #cm^2/gr
    y *= 2
    return y

def output_particles_cross_section(siminfo):

    # Load data:
    g = h5py.File(siminfo.snapshot, "r")
    unit_length = g["Units"].attrs["Unit length in cgs (U_L)"]
    unit_mass = g["Units"].attrs["Unit mass in cgs (U_M)"]

    velocity_dispersion = g["PartType1/Velocity_dispersion"][:][:] #km/s
    cross_section = g["PartType1/Cross_section"][:] * unit_length**2 / unit_mass #cm^2/g
    probability = g["PartType1/SIDM_probability"][:] #internal units
    timestep = g["PartType1/Time_step_size"][:] #internal units
    probability *= timestep
    search_radius = g["PartType1/SIDM_search_radius"][:] * 1e3 #kpc

    velocity_dispersion *= 4 / np.sqrt(np.pi)  #Assume vel. dispersion relates to relative particle velocities

    xrange = np.arange(0, 3, 0.2)
    xrange = 10 ** xrange
    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(velocity_dispersion, cross_section, xrange)

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


    #plt.plot(mass, c200, 'o', color='tab:blue', alpha=0.2)
    plt.plot(xvalues, yvalues, '-', color='white', zorder=10)
    plt.plot(xvalues, yvalues, '-', lw=1.5, color='tab:blue', zorder=10,label=siminfo.name)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:blue')

    plt.plot(xvalues, sigma(xvalues, 3.0,0.3,6.74e-6), '-', lw=2, color='tab:green',
             label=r'$m_{x}{=}3.0$GeV, $m_{\phi}{=}0.3$MeV, $\alpha{=}6.74e{-}6$')


    plt.axis([1, 1e3, 0, 60])
    plt.xlabel("Velocity dispersion [km/s]")
    plt.ylabel("Cross section [cm$^{2}$/g]")
    plt.xscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{siminfo.output_path}/particle_data_1.png", dpi=200)
    plt.close()

    np.savetxt(f"{siminfo.output_path}/particles_relation_1_" + siminfo.name + ".txt",
               np.transpose([xvalues, yvalues, yvalues_err_down, yvalues_err_up]))

    # xrange = np.arange(0, 3, 0.2)
    # xrange = 10 ** xrange
    # xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(search_radius, probability, xrange)

    # figure()
    # ax = plt.subplot(1, 1, 1)
    # plt.grid("True")
    #
    # # plt.plot(mass, c200, 'o', color='tab:blue', alpha=0.2)
    # plt.plot(xvalues, yvalues, '-', color='white', zorder=10)
    # plt.plot(xvalues, yvalues, '-', lw=1.5, color='tab:blue', zorder=10)
    # plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:blue')
    #
    # plt.axis([1, 14, 1e-8, 1e-2])
    # plt.xlabel("Search Radius [kpc]")
    # plt.ylabel("Probability")
    # plt.yscale('log')
    # ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    # #plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    # plt.savefig(f"{siminfo.output_path}/particle_data_2.png", dpi=200)
    # plt.close()
    #
    # np.savetxt(f"{siminfo.output_path}/particles_relation_2_" + siminfo.name + ".txt",
    #            np.transpose([xvalues, yvalues, yvalues_err_down, yvalues_err_up]))



def plot_relations(siminfo, name_list):

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

    color = ['tab:blue', 'tab:orange']
    i = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/cMVmax_" + name + ".txt")

        mass = data[:, 0]
        c200 = data[:, 1]

        xrange = np.arange(6,15,0.2)
        xrange = 10**xrange
        xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass, c200, xrange)

        plt.plot(mass, c200, 'o', color=color[i],alpha=0.2)
        plt.plot(xvalues, yvalues, '-', color='white',zorder=10)
        plt.plot(xvalues, yvalues, '-', lw=1.5, color=color[i], label=name,zorder=10)
        plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color=color[i])
        i += 1

    plt.axis([1e6,1e15,1,50])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("c$_{200}$")
    plt.xscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{siminfo.output_path}/cM_relation.png", dpi=200)
    plt.close()

    ###########
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    color = ['tab:blue', 'tab:orange']
    i = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/cMVmax_" + name + ".txt")

        mass = data[:, 0]
        Vmax = data[:, 2]

        xrange = np.arange(6, 15, 0.2)
        xrange = 10 ** xrange
        xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(mass, Vmax, xrange)

        plt.plot(mass, Vmax, 'o', color=color[i], alpha=0.2)
        plt.plot(xvalues, yvalues, '-', color='white',zorder=10)
        plt.plot(xvalues, yvalues, '-', lw=1.5, color=color[i], label=name,zorder=10)
        plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color=color[i])
        i += 1

    plt.axis([1e6, 1e15, 1, 1e3])
    plt.xlabel("M${}_{200,\mathrm{crit}}$ ($M_\odot$)")
    plt.ylabel("V$_{\mathrm{max}}$ [km/s]")
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig(f"{siminfo.output_path}/VmaxM_relation.png", dpi=200)
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

def analyse_halo(mass, pos, vel, sigma):
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

    nozero = sigma > 0
    sigma, _, _ = stat.binned_statistic(x=r[nozero], values=sigma[nozero], statistic="median", bins=radial_bins, )
    sigma[np.where(np.isnan(sigma))[0]] = 0
    return density, velocity, sigma


def read_data(siminfo, option):

    with h5py.File(siminfo.snapshot, "r") as hf:
        mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        pos = hf['PartType1/Coordinates'][:][:] * siminfo.a
        vel = hf['PartType1/Velocities'][:][:]
        cross_section = hf["PartType1/Cross_section"][:]
        # Read units
        unit_length_in_cgs = hf["/Units"].attrs["Unit length in cgs (U_L)"]
        unit_time_in_cgs = hf["/Units"].attrs["Unit time in cgs (U_t)"]
        unit_mass_in_cgs = hf["Units"].attrs["Unit mass in cgs (U_M)"]
        vel *= unit_length_in_cgs / unit_time_in_cgs  # cm/s
        vel *= 1e-5  # km/s
        cross_section *= unit_length_in_cgs ** 2 / unit_mass_in_cgs  # cm^2/g

    snapshot_file = h5py.File(siminfo.snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")
    properties_file = h5py.File(siminfo.halo_properties, "r")

    c200c = properties_file["cNFW_200crit"][:]
    c200c[c200c == 0] = 1
    m200c = properties_file["Mass_200crit"][:] * 1e10
    m200c[m200c <= 0] = 1
    R200c = properties_file["R_200crit"][:] * 1e3 #kpc
    xCoP = properties_file["Xcminpot"][:]
    yCoP = properties_file["Ycminpot"][:]
    zCoP = properties_file["Zcminpot"][:]
    m200c = np.log10(m200c)
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP
    Type = properties_file["Structuretype"][:]

    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    density = np.zeros((len(centers), 3))
    sig_density = np.zeros((len(centers), 3))
    velocity = np.zeros((len(centers), 3))
    sig_velocity = np.zeros((len(centers), 3))
    sigma = np.zeros((len(centers), 3))
    sig_sigma = np.zeros((len(centers), 3))
    rs = np.zeros(3)
    M200 = np.zeros(3)

    for i in range(0, 3):
        if i == 0: select_halos = np.where((m200c >= 8.9) & (m200c <= 9.1))[0]  # >10 star parts
        if i == 1: select_halos = np.where((m200c >= 9.9) & (m200c <= 10.1))[0]  # >10 star parts
        if i == 2: select_halos = np.where((m200c >= 10.9) & (m200c <= 11.1))[0]  # >10 star parts

        if option == "centrals":
            centrals = np.where(Type[select_halos] == 10)[0]
            select_halos = select_halos[centrals]
        else :
            satellites = np.where(Type[select_halos] > 10)[0]
            select_halos = select_halos[satellites]


        if len(select_halos) >= 10:
            select_random = np.random.random_integers(len(select_halos)-1, size=(10))
            select_halos = select_halos[select_random]

        rs[i] = np.median(R200c[select_halos] / c200c[select_halos]) # kpc
        M200[i] = np.median(10 ** m200c[select_halos])

        num_halos = len(select_halos)
        density_all = np.zeros((len(centers), num_halos))
        velocity_all = np.zeros((len(centers), num_halos))
        sigma_all = np.zeros((len(centers), num_halos))

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
            particles_sigma = cross_section[indices_p].copy()
            if len(particles_mass) == 0: continue

            density_halo, velocity_halo, sigma_halo = analyse_halo(particles_mass, particles_pos, particles_vel, particles_sigma)
            density_all[:, halo] = density_halo
            velocity_all[:, halo] = velocity_halo
            sigma_all[:, halo] = sigma_halo

        density[:, i] = np.mean(density_all[:, :], axis=1)
        sig_density[:, i] = np.std(density_all[:, :], axis=1)
        velocity[:, i] = np.mean(velocity_all[:, :], axis=1)
        sig_velocity[:, i] = np.std(velocity_all[:, :], axis=1)
        sigma[:, i] = np.mean(sigma_all[:, :], axis=1)
        sig_sigma[:, i] = np.std(sigma_all[:, :], axis=1)

    return rs, M200, density, sig_density, velocity, sig_velocity, sigma, sig_sigma


def output_halo_profiles(siminfo):

    rs, M200, density, sig_density, velocity, sig_velocity, sigma, sig_sigma = read_data(siminfo, "centrals")

    np.savetxt(f"{siminfo.output_path}/Halo_details_" + siminfo.name + ".txt",
               np.transpose([rs, M200]))

    np.savetxt(f"{siminfo.output_path}/Density_profiles_" + siminfo.name + ".txt",
               np.transpose([density[:,0], density[:,1], density[:,2],
                             sig_density[:,0], sig_density[:,1], sig_density[:,2]]))

    np.savetxt(f"{siminfo.output_path}/Velocity_profiles_" + siminfo.name + ".txt",
               np.transpose([velocity[:,0], velocity[:,1], velocity[:,2],
                             sig_velocity[:,0], sig_velocity[:,1], sig_velocity[:,2]]))

    np.savetxt(f"{siminfo.output_path}/Cross_section_profiles_" + siminfo.name + ".txt",
               np.transpose([sigma[:,0], sigma[:,1], sigma[:,2],
                             sig_sigma[:,0], sig_sigma[:,1], sig_sigma[:,2]]))


    rs, M200, density, sig_density, velocity, sig_velocity, sigma, sig_sigma = read_data(siminfo, "satellites")

    np.savetxt(f"{siminfo.output_path}/Subhalo_details_" + siminfo.name + ".txt",
               np.transpose([rs, M200]))

    np.savetxt(f"{siminfo.output_path}/Density_profiles_sub_" + siminfo.name + ".txt",
               np.transpose([density[:,0], density[:,1], density[:,2],
                             sig_density[:,0], sig_density[:,1], sig_density[:,2]]))

    np.savetxt(f"{siminfo.output_path}/Velocity_profiles_sub_" + siminfo.name + ".txt",
               np.transpose([velocity[:,0], velocity[:,1], velocity[:,2],
                             sig_velocity[:,0], sig_velocity[:,1], sig_velocity[:,2]]))

    np.savetxt(f"{siminfo.output_path}/Cross_section_profiles_sub_" + siminfo.name + ".txt",
               np.transpose([sigma[:,0], sigma[:,1], sigma[:,2],
                             sig_sigma[:,0], sig_sigma[:,1], sig_sigma[:,2]]))

def plot_halo_profiles(siminfo, name_list, option):

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

    z = siminfo.z  # Redshift
    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    for i in range(3):

        j = i*2 + 1

        ax = plt.subplot(3, 2, j)
        plt.grid("True")

        color = ['tab:blue', 'tab:orange']
        k = 0
        for name in name_list:
            if option == "centrals":
                data = np.loadtxt(f"{siminfo.output_path}/Density_profiles_" + name + ".txt")
            else :
                data = np.loadtxt(f"{siminfo.output_path}/Density_profiles_sub_" + name + ".txt")

            density = data[:,0:3]
            sig_density = data[:,3:6]

            plt.plot(centers, density[:, i], lw=2, color=color[k], label=name)
            plt.fill_between(centers, density[:, i] - sig_density[:, i] / 2,
                             density[:, i] + sig_density[:, i] / 2, alpha=0.4,
                             color=color[k])

            k += 1

        if option == "centrals":
            data = np.loadtxt(f"{siminfo.output_path}/Halo_details_" + name + ".txt")
        else :
            data = np.loadtxt(f"{siminfo.output_path}/Subhalo_details_" + name + ".txt")
        M200 = data[:, 1]
        rs = data[:, 0]

        NFWrho = calc_density(centers, M200[i], rs[i], z, siminfo)  # Msun/kpc^3
        plt.plot(centers, NFWrho, lw=1, color='black', label="NFW profile")

        mass = np.log10(M200[i])
        plt.text(0.05,0.05,'$10^{%0.1f'%mass+'}$M$_{\odot}$ halo', transform=ax.transAxes)

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

        k = 0
        for name in name_list:
            if option == "centrals":
                data = np.loadtxt(f"{siminfo.output_path}/Velocity_profiles_" + name + ".txt")
            else :
                data = np.loadtxt(f"{siminfo.output_path}/Velocity_profiles_sub_" + name + ".txt")
            velocity = data[:, 0:3]
            sig_velocity = data[:, 3:6]

            plt.plot(centers, velocity[:, i], lw=2, color=color[k], label=name)
            plt.fill_between(centers, velocity[:, i] - sig_velocity[:, i] / 2,
                             velocity[:, i] + sig_velocity[:, i] / 2,
                             alpha=0.4, color=color[k])
            k += 1

        if option == "centrals":
            data = np.loadtxt(f"{siminfo.output_path}/Halo_details_" + name + ".txt")
        else:
            data = np.loadtxt(f"{siminfo.output_path}/Subhalo_details_" + name + ".txt")

        M200 = data[:, 1]
        rs = data[:, 0]
        NFWsig1D = sigma_1D(centers, M200[i], rs[i], z, siminfo)  # km/s
        plt.plot(centers, NFWsig1D, lw=1, color='black')

        #plt.ylim(0, 80)
        plt.xlim(1, 1e3)
        plt.xscale("log")
        plt.xlabel("Radius [kpc]")
        plt.ylabel("Velocity dispersion [km/s]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    if option == "centrals":
        plt.savefig(f"{siminfo.output_path}/Density_profiles_halos.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/Density_profiles_subhalos.png", dpi=200)
    plt.close()

def plot_cross_section_profiles(siminfo, name_list, option):

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

    z = siminfo.z  # Redshift
    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    for i in range(3):

        j = i*2 + 1

        ax = plt.subplot(3, 2, j)
        plt.grid("True")

        color = ['tab:blue', 'tab:orange']
        k = 0
        for name in name_list:
            if option == "centrals":
                data = np.loadtxt(f"{siminfo.output_path}/Cross_section_profiles_" + name + ".txt")
            else:
                data = np.loadtxt(f"{siminfo.output_path}/Cross_section_profiles_sub_" + name + ".txt")
            density = data[:,0:3]
            sig_density = data[:,3:6]

            plt.plot(centers, density[:, i], lw=2, color=color[k], label=name)
            plt.fill_between(centers, density[:, i] - sig_density[:, i] / 2,
                             density[:, i] + sig_density[:, i] / 2, alpha=0.4,
                             color=color[k])

            k += 1

        if option == "centrals":
            data = np.loadtxt(f"{siminfo.output_path}/Halo_details_" + name + ".txt")
        else:
            data = np.loadtxt(f"{siminfo.output_path}/Subhalo_details_" + name + ".txt")

        M200 = data[:, 1]
        mass = np.log10(M200[i])
        plt.text(0.05,0.05,'$10^{%0.1f'%mass+'}$M$_{\odot}$ halo', transform=ax.transAxes)

        xarray = np.array([2.3, 2.3])
        yarray = np.array([0,90])
        plt.plot(xarray, yarray, '--', color='grey')
        plt.ylim(0, 50)
        plt.xlim(1, 1e3)
        plt.xscale("log")
        plt.xlabel("Radius [kpc]")
        plt.ylabel("Cross section [cm$^{2}$/g]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

        ######################
        ax = plt.subplot(3, 2, j+1)
        plt.grid("True")

        k = 0
        for name in name_list:
            if option == "centrals":
                data = np.loadtxt(f"{siminfo.output_path}/Velocity_profiles_" + name + ".txt")
            else:
                data = np.loadtxt(f"{siminfo.output_path}/Velocity_profiles_sub_" + name + ".txt")
            velocity = data[:, 0:3]
            sig_velocity = data[:, 3:6]

            plt.plot(centers, velocity[:, i], lw=2, color=color[k], label=name)
            plt.fill_between(centers, velocity[:, i] - sig_velocity[:, i] / 2,
                             velocity[:, i] + sig_velocity[:, i] / 2,
                             alpha=0.4, color=color[k])
            k += 1

        if option == "centrals":
            data = np.loadtxt(f"{siminfo.output_path}/Halo_details_" + name + ".txt")
        else:
            data = np.loadtxt(f"{siminfo.output_path}/Subhalo_details_" + name + ".txt")
        M200 = data[:, 1]
        rs = data[:, 0]
        NFWsig1D = sigma_1D(centers, M200[i], rs[i], z, siminfo)  # km/s
        plt.plot(centers, NFWsig1D, lw=1, color='black')

        #plt.ylim(0, 80)
        plt.xlim(1, 1e3)
        plt.xscale("log")
        plt.xlabel("Radius [kpc]")
        plt.ylabel("Velocity dispersion [km/s]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    if option == "centrals":
        plt.savefig(f"{siminfo.output_path}/Cross_section_profiles_halos.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/Cross_section_profiles_subhalos.png", dpi=200)
    plt.close()


def read_individual_profiles(siminfo):

    with h5py.File(siminfo.snapshot, "r") as hf:
        mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        pos = hf['PartType1/Coordinates'][:][:] * siminfo.a
        vel = hf['PartType1/Velocities'][:][:]
        # Read units
        unit_length_in_cgs = hf["/Units"].attrs["Unit length in cgs (U_L)"]
        unit_time_in_cgs = hf["/Units"].attrs["Unit time in cgs (U_t)"]
        vel *= unit_length_in_cgs / unit_time_in_cgs  # cm/s
        vel *= 1e-5  # km/s

    snapshot_file = h5py.File(siminfo.snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")
    properties_file = h5py.File(siminfo.halo_properties, "r")

    c200c = properties_file["cNFW_200crit"][:]
    c200c[c200c == 0] = 1
    m200c = properties_file["Mass_200crit"][:] * 1e10
    m200c[m200c <= 0] = 1
    R200c = properties_file["R_200crit"][:] * 1e3 #kpc
    xCoP = properties_file["Xcminpot"][:]
    yCoP = properties_file["Ycminpot"][:]
    zCoP = properties_file["Zcminpot"][:]
    subtype = properties_file["Structuretype"]
    m200c = np.log10(m200c)
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP

    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    rs = np.zeros(3)
    M200 = np.zeros(3)

    for i in range(0, 3):
        if i == 0:
            select_halos = np.where((m200c >= 8.8) & (m200c <= 9.2))[0]  # >10 star parts
        if i == 1:
            select_halos = np.where((m200c >= 9.8) & (m200c <= 10.2))[0]  # >10 star parts
        if i == 2:
            select_halos = np.where((m200c >= 10.8) & (m200c <= 11.2))[0]  # >10 star parts

        select_main = np.where(subtype[select_halos] == 10)[0]
        select_sub = np.where(subtype[select_halos] > 10)[0]

        if len(select_main) >= 30:
            select_random = np.random.random_integers(len(select_main)-1, size=(30))
            select_main = select_main[select_random]

        if len(select_sub) >= 30:
            select_random = np.random.random_integers(len(select_sub)-1, size=(30))
            select_sub = select_sub[select_random]

        rs[i] = np.median(R200c[select_halos] / c200c[select_halos]) # kpc
        M200[i] = np.median(10 ** m200c[select_halos])

        num_halos = len(select_main)
        if i == 0: density_main_10 = np.zeros((len(centers), num_halos))
        if i == 1: density_main_11 = np.zeros((len(centers), num_halos))
        if i == 2: density_main_12 = np.zeros((len(centers), num_halos))

        for halo in range(0, num_halos-1):

            halo_j = select_halos[select_main[halo]]

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
            if len(particles_mass) == 0: continue

            density_halo, _ = analyse_halo(particles_mass, particles_pos, particles_vel)

            if i == 0: density_main_10[:, halo] = density_halo
            if i == 1: density_main_11[:, halo] = density_halo
            if i == 2: density_main_12[:, halo] = density_halo


        num_halos = len(select_sub)

        if i == 0: density_sub_10 = np.zeros((len(centers), num_halos))
        if i == 1: density_sub_11 = np.zeros((len(centers), num_halos))
        if i == 2: density_sub_12 = np.zeros((len(centers), num_halos))

        # Check
        if num_halos == 0:continue

        for halo in range(0, num_halos-1):

            halo_j = select_halos[select_sub[halo]]

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
            if len(particles_mass) == 0: continue

            density_halo, _ = analyse_halo(particles_mass, particles_pos, particles_vel)
            if i == 0: density_sub_10[:, halo] = density_halo
            if i == 1: density_sub_11[:, halo] = density_halo
            if i == 2: density_sub_12[:, halo] = density_halo

    return rs, M200, density_main_10, density_main_11, density_main_12, \
           density_sub_10, density_sub_11, density_sub_12

def plot_individual_profiles(siminfo,output_path):

    rs, M200, density_main_10, density_main_11, density_main_12, \
    density_sub_10, density_sub_11, density_sub_12 = read_individual_profiles(siminfo)

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
            label = '$10^{9}$M$_{\odot}$ halo'
            color = 'tab:green'
            density = density_main_10
            density_sub = density_sub_10
        if i==1:
            label = '$10^{10}$M$_{\odot}$ halo'
            color = 'tab:orange'
            density = density_main_11
            density_sub = density_sub_11
        if i==2:
            label = '$10^{11}$M$_{\odot}$ halo'
            color = 'tab:blue'
            density = density_main_12
            density_sub = density_sub_12

        NFWrho = calc_density(centers, M200[i], rs[i], z, siminfo)  # Msun/kpc^3
        plt.plot(centers, NFWrho, lw=1, color='black', label="NFW profile")

        for k in range(len(density[0,:])):
            #nozero = np.where(density[:,k]>0)[0]
            plt.plot(centers, density[:, k], lw=0.5, color=color)

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

        plt.plot(centers, NFWrho, lw=1, color='black', label="NFW profile")

        for k in range(len(density_sub[0,:])):
            plt.plot(centers, density_sub[:, k], lw=0.5, color=color)

        plt.ylim(1e3, 1e9)
        plt.xlim(1, 1e3)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Radius [kpc]")
        plt.ylabel("Density [M$_{\odot}$/kpc$^{3}$]")
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

    plt.savefig(output_path + "Density_profiles_all.png", dpi=200)
    plt.close()
