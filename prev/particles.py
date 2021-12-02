import h5py
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import scipy.stats as stat

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

def analyse_halo(pos, property):
    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))
    median_prop, _, _ = stat.binned_statistic(x=r, values=property, statistic="median", bins=radial_bins, )
    return median_prop

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


def read_data(siminfo):

    with h5py.File(siminfo.snapshot, "r") as hf:
        pos = hf['PartType1/Coordinates'][:][:] * siminfo.a
        cross_section = hf["PartType1/Cross_section"][:]
        #cross_section = np.ones(len(pos[:,0]))
        max_sidm_event_per_timestep = hf["PartType1/Max_SIDM_events_per_timestep"][:]
        ngb = hf["PartType1/Number_of_neighbours"][:]
        sidm_events = hf["PartType1/SIDM_events"][:]
        sidm_prob = hf["PartType1/SIDM_probability"][:]
        sidm_h = hf["PartType1/SIDM_search_radius"][:] * siminfo.a
        timestep = hf["PartType1/Time_step_size"][:]
        sidm_prob *= timestep

        # Read units
        unit_length_in_cgs = hf["/Units"].attrs["Unit length in cgs (U_L)"]
        unit_time_in_cgs = hf["/Units"].attrs["Unit time in cgs (U_t)"]
        unit_mass_in_cgs = hf["Units"].attrs["Unit mass in cgs (U_M)"]
        kpc_in_cgs = 3.086e21
        cross_section *= unit_length_in_cgs ** 2 / unit_mass_in_cgs  # cm^2/g
        sidm_h *= unit_length_in_cgs / kpc_in_cgs #kpc
        pos *= unit_length_in_cgs / kpc_in_cgs #kpc

    snapshot_file = h5py.File(siminfo.snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")
    properties_file = h5py.File(siminfo.halo_properties, "r")

    m200c = properties_file["Mass_200crit"][:] * 1e10
    m200c[m200c <= 0] = 1
    m200c = np.log10(m200c)
    xCoP = properties_file["Xcminpot"][:] * unit_length_in_cgs / kpc_in_cgs #kpc
    yCoP = properties_file["Ycminpot"][:] * unit_length_in_cgs / kpc_in_cgs #kpc
    zCoP = properties_file["Zcminpot"][:] * unit_length_in_cgs / kpc_in_cgs #kpc
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP
    Type = properties_file["Structuretype"][:]
    select_halos = np.where((m200c >= 9.9) & (m200c <= 10.1))[0]  # >10 star parts
    select = np.where(Type[select_halos] == 10)[0]
    select_halos = select_halos[select]
    select_random = np.random.random_integers(len(select_halos) - 1, size=(10))
    select_halos = select_halos[select_random]

    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    num_halos = len(select_halos)
    mean_properties = np.zeros((len(centers), 7))
    properties = np.zeros((len(centers), num_halos, 7))

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

        if len(indices_p) == 0: continue

        particles_pos = pos[indices_p, :].copy()
        particles_pos -= CoP[halo_j, :]  # centering

        # Let's do cross section
        particles_prop = cross_section[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 0] = halo_prop

        # Let's do SIDM events
        particles_prop = max_sidm_event_per_timestep[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 1] = halo_prop

        # Let's do ngb
        particles_prop = ngb[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 2] = halo_prop

        # Let's do sidm_events
        particles_prop = sidm_events[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 3] = halo_prop

        # Let's do sidm_prob
        particles_prop = sidm_prob[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 4] = halo_prop

        # Let's do sidm_h
        particles_prop = sidm_h[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 5] = halo_prop

        # Let's do timestep
        particles_prop = timestep[indices_p].copy()
        halo_prop = analyse_halo(particles_pos, particles_prop)
        properties[:, halo, 6] = halo_prop

    for i in range(0,7):
        mean_properties[:, i] = np.mean(properties[:, :, i], axis=1)

    return centers, mean_properties

def output_particles_properties(siminfo):

    centers, properties = read_data(siminfo)
    np.savetxt(f"{siminfo.output_path}/Particles_properties_" + siminfo.name + ".txt",
               np.transpose([centers, properties[:, 0], properties[:, 1], properties[:, 2], properties[:, 3],
                             properties[:, 4], properties[:, 5], properties[:, 6]]))


def plot_particles_properties(siminfo, name_list):

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
    color = ['tab:blue', 'tab:orange']

    ######################
    figure()

    ax = plt.subplot(3, 2, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Particles_properties_" + name + ".txt")
        radius = data[:, 0]
        property = data[:, 1:7]

        plt.plot(radius, property[:, 0], lw=2, color=color[k], label=name)
        k += 1

    #plt.ylim(1e3, 1e9)
    plt.xlim(1, 1e3)
    plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Cross section")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)

    ######################
    ax = plt.subplot(3, 2, 2)
    plt.grid("True")

    k = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Particles_properties_" + name + ".txt")
        radius = data[:, 0]
        property = data[:, 1:8]

        plt.plot(radius, property[:, 6], lw=2, color=color[k])
        k += 1

    #plt.ylim(1e3, 1e9)
    plt.xlim(1, 1e3)
    plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Timestep size")

    ######################
    ax = plt.subplot(3, 2, 3)
    plt.grid("True")

    k = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Particles_properties_" + name + ".txt")
        radius = data[:, 0]
        property = data[:, 1:7]

        plt.plot(radius, property[:, 2], lw=2, color=color[k])
        k += 1

    # plt.ylim(1e3, 1e9)
    plt.xlim(1, 1e3)
    plt.xscale("log")
    # plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Number of neighbours")

    ######################
    ax = plt.subplot(3, 2, 4)
    plt.grid("True")

    k = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Particles_properties_" + name + ".txt")
        radius = data[:, 0]
        property = data[:, 1:7]

        plt.plot(radius, property[:, 3], lw=2, color=color[k])
        k += 1

    # plt.ylim(1e3, 1e9)
    plt.xlim(1, 1e3)
    plt.xscale("log")
    # plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("SIDM events")

    ######################
    ax = plt.subplot(3, 2, 5)
    plt.grid("True")

    k = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Particles_properties_" + name + ".txt")
        radius = data[:, 0]
        property = data[:, 1:7]

        plt.plot(radius, property[:, 4], lw=2, color=color[k])
        k += 1

    # plt.ylim(1e3, 1e9)
    plt.xlim(1, 1e3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("SIDM probability")

    ######################
    ax = plt.subplot(3, 2, 6)
    plt.grid("True")

    k = 0
    for name in name_list:
        data = np.loadtxt(f"{siminfo.output_path}/Particles_properties_" + name + ".txt")
        radius = data[:, 0]
        property = data[:, 1:7]

        plt.plot(radius, property[:, 5], lw=2, color=color[k])
        k += 1

    # plt.ylim(1e3, 1e9)
    plt.xlim(1, 1e3)
    plt.xscale("log")
    #plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("SIDM search radius [kpc]")

    plt.savefig(f"{siminfo.output_path}/Particles_property.png", dpi=200)