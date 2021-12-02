import numpy as np
import h5py
from scipy import interpolate
import scipy.stats as stat
from pylab import *
import matplotlib.pyplot as plt

def c_M_relation(M0):
    """
    Concentration-mass relation from Correa et al. (2015).
    This relation is most suitable for Planck cosmology.
    """
    log_M0 = np.log10(M0)
    z = 0
    # Best-fit params:
    alpha = 1.7543 - 0.2766 * (1. + z) + 0.02039 * (1. + z) ** 2
    beta = 0.2753 + 0.0035 * (1. + z) - 0.3038 * (1. + z) ** 0.0269
    gamma = -0.01537 + 0.02102 * (1. + z) ** (-0.1475)

    log_10_c200 = alpha + beta * log_M0 * (1. + gamma * log_M0 ** 2)
    c200 = 10 ** log_10_c200
    return c200


def calc_NFW(x, M, c, z, siminfo):
    Om = siminfo.Om
    Ol = 1. - Om
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol) #Msun/Mpc^3
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    a = R200 / c

    # NFW profile #
    delta = 200. / 3.
    delta *= c ** 3 / (np.log(1. + c) - c / (1. + c))
    f = np.zeros(len(x))
    for i in range(0, len(x)): f[i] = rho_crit * delta * 1e-9 / ((x[i] / a) * (1. + x[i] / a) ** 2)
    return f

def calc_NFW_vcirc(r, M, c, z, siminfo):

    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21
    G = 6.67408e-11 #m^3 kg^-1 s^-2
    G *= (1e2)**3/1e3 #cm^3 g^-1 s^-2
    G /= kpc_in_cgs**3
    G *= Msun_in_cgs #kpc^3 Msun^-1 s^-2


    Om = siminfo.Om
    Ol = 1. - Om
    rho_crit = siminfo.rhocrit0 * (Om * (1. + z) ** 3 + Ol) #Msun/Mpc^3
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    x = r / R200
    gcx = np.log(1. + c * x) - c * x / (1. + c * x)
    gc = np.log(1. + c) - c / (1. + c)
    mass = 200 * 4 * np.pi * R200**3 * rho_crit * 1e-9 * gcx / (3. * gc) #Msun

    Vcirc = G * mass / r #kpc^2/s^2
    Vcirc *= (kpc_in_cgs/1e5)**2 #km^2/s^2
    Vcirc = np.sqrt(Vcirc) #km/s
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


def analyse_rotation(mass,radius):
    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21
    G = 6.67408e-11 #m^3 kg^-1 s^-2
    G *= (1e2)**3/1e3 #cm^3 g^-1 s^-2
    G /= kpc_in_cgs**3
    G *= Msun_in_cgs #kpc^3 Msun^-1 s^-2

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0.2, 25, 0.25)
    centers = bin_centers(radial_bins)  # kpc

    SumMasses, _, _ = stat.binned_statistic(x=radius, values=mass, statistic="sum", bins=radial_bins, )
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


def rotation_curve(siminfo, galaxy_type):

    with h5py.File(siminfo.snapshot, "r") as hf:
        DMpart_mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        DMpart_pos = hf['PartType1/Coordinates'][:][:] * siminfo.a # Mpc

        #if siminfo.sim_type == "Hydro":
        #    Gaspart_mass = hf['PartType0/Masses'][:] * 1e10  # Msun
        #    Gaspart_pos = hf['PartType0/Coordinates'][:][:] * siminfo.a # Mpc
        #    Starpart_mass = hf['PartType0/Masses'][:] * 1e10  # Msun
        #    Starpart_pos = hf['PartType0/Coordinates'][:][:] * siminfo.a # Mpc


    with h5py.File(siminfo.halo_properties, "r") as properties_file:
        m200c = properties_file["Mass_200crit"][:] * 1e10 # Msun
        c200c = properties_file["cNFW_200crit"][:]
        xCoP = properties_file["Xcminpot"][:] # Mpc
        yCoP = properties_file["Ycminpot"][:] # Mpc
        zCoP = properties_file["Zcminpot"][:] # Mpc
        Type = properties_file["Structuretype"][:]
        Vmax = properties_file["Vmax"][:] # km/s

    c200c[c200c <= 0] = 1
    m200c[m200c <= 0] = 1
    m200c = np.log10(m200c)
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP
    #sample = np.where((Vmax >= 70) & (Vmax <= 80))[0]
    sample = np.where((m200c >= 9) & (m200c < 9.2))[0]

    if galaxy_type == "centrals":
        centrals = np.where(Type[sample] == 10)[0]
    else:
        centrals = np.where(Type[sample] > 10)[0]

    print(len(sample),len(centrals))
    sample = sample[centrals]

    Mmedian = np.median(m200c[sample])
    cmedian = np.median(c200c[sample])

    if len(sample) >= 50:
       select_random = np.random.random_integers(len(sample) - 1, size=(50))
       sample = sample[select_random]

    # Reduce
    CoP = CoP[sample, :]
    m200c = m200c[sample]
    radial_bins = np.arange(0.2, 25, 0.25)
    centers = bin_centers(radial_bins) # kpc

    v_circ = np.zeros((len(centers), len(sample)))

    #if siminfo.sim_type == "Hydro":
    #    num_parts = len(DMpart_mass) + len(Gaspart_mass) + len(Starpart_mass)
    #    particles_pos = np.zeros((num_parts,3))
    #    particles_pos[0:len(DMpart_mass),:] = DMpart_pos
    #    particles_pos[len(DMpart_mass):len(DMpart_mass)+len(Gaspart_mass),:] = Gaspart_pos
    #    particles_pos[len(DMpart_mass)+len(Gaspart_mass):len(DMpart_mass)+len(Gaspart_mass)+len(Starpart_mass),:] =  Starpart_pos

    #    mass = DMpart_mass.copy()
    #    mass = np.append(mass, Gaspart_mass)
    #    mass = np.append(mass, Starpart_mass)
    #else:
    #    particles_pos = DMpart_pos.copy()
    #    mass = DMpart_mass.copy()

    snapshot_file = h5py.File(siminfo.snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")

    particles_pos = DMpart_pos.copy()
    mass = DMpart_mass.copy()

    for i in range(len(sample)):
        halo_i = sample[i]

        # Grab the start position in the particles file to read from
        halo_start_position = group_file["Offset"][halo_i]
        halo_end_position = group_file["Offset"][halo_i + 1]
        particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]
        particle_ids_from_snapshot = snapshot_file["PartType1/ParticleIDs"][...]

        _, indices_v, indices_p = np.intersect1d(particle_ids_in_halo,
                                                 particle_ids_from_snapshot,
                                                 assume_unique=True,
                                                 return_indices=True, )

        all_pos = particles_pos[indices_p,:].copy()
        all_pos -= CoP[i, :]  # centering
        all_pos *= 1e3  # kpc

        radius = np.linalg.norm(all_pos[:, :3], axis=1)  # computing distances to the CoP
        vrot = analyse_rotation(mass[indices_p], radius)  #
        v_circ[:, i] = vrot


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
        plt.plot(centers, v_circ[:, i], '-', label='M$_{200}$=%.3f [$\log_{10}$M$_{\odot}$]'%m200c[i])

    NFW_circ = calc_NFW_vcirc(centers, 10**Mmedian, cmedian, siminfo.z, siminfo)
    plt.plot(centers, NFW_circ,'--',lw=2,color='black')

    plt.axis([0, 20, 0, 100])
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Circular velocity [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/Rotation_curve_"+siminfo.name+"_centrals.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/Rotation_curve_"+siminfo.name+"_satellites.png", dpi=200)
    plt.close()



def percentile(y,z):
   return np.percentile(y,z)

def output_rotation_data(siminfo):

    with h5py.File(siminfo.snapshot, "r") as hf:
        DMpart_mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        DMpart_pos = hf['PartType1/Coordinates'][:][:] * siminfo.a # Mpc
        #if siminfo.sim_type == "Hydro":
        #    Gaspart_mass = hf['PartType0/Masses'][:] * 1e10  # Msun
        #    Gaspart_pos = hf['PartType0/Coordinates'][:][:] * siminfo.a # Mpc
        #    Starpart_mass = hf['PartType0/Masses'][:] * 1e10  # Msun
        #    Starpart_pos = hf['PartType0/Coordinates'][:][:] * siminfo.a # Mpc


    with h5py.File(siminfo.halo_properties, "r") as properties_file:
        m200c = properties_file["Mass_200crit"][:] * 1e10 # Msun
        c200c = properties_file["cNFW_200crit"][:]
        xCoP = properties_file["Xcminpot"][:] # Mpc
        yCoP = properties_file["Ycminpot"][:] # Mpc
        zCoP = properties_file["Zcminpot"][:] # Mpc
        Type = properties_file["Structuretype"][:]
        ID = properties_file["ID"][:]

    c200c[c200c <= 0] = 1
    m200c[m200c <= 0] = 1
    m200c = np.log10(m200c)
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP
    sample = np.where((m200c >= 9) & (m200c < 10))[0]
    #sample = np.where(m200c >= 10)[0]
    print(len(sample))

    # Reduce
    CoP = CoP[sample, :]
    m200c = m200c[sample]
    c200c = c200c[sample]
    ID = ID[sample]
    Type = Type[sample]

    v_circ_1 = np.zeros(len(sample))
    v_circ_2 = np.zeros(len(sample))
    v_circ_5 = np.zeros(len(sample))
    Vmax = np.zeros(len(sample))

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0.2, 25, 0.25)
    centers = bin_centers(radial_bins)  # kpc

    #if siminfo.sim_type == "Hydro":
    #    num_parts = len(DMpart_mass) + len(Gaspart_mass) + len(Starpart_mass)
    #    particles_pos = np.zeros((num_parts,3))
    #    particles_pos[0:len(DMpart_mass),:] = DMpart_pos
    #    particles_pos[len(DMpart_mass):len(DMpart_mass)+len(Gaspart_mass),:] = Gaspart_pos
    #    particles_pos[len(DMpart_mass)+len(Gaspart_mass):len(DMpart_mass)+len(Gaspart_mass)+len(Starpart_mass),:] =  Starpart_pos

    #    mass = DMpart_mass.copy()
    #    mass = np.append(mass, Gaspart_mass)
    #    mass = np.append(mass, Starpart_mass)
    #else:

    snapshot_file = h5py.File(siminfo.snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")

    particles_pos = DMpart_pos.copy()
    mass = DMpart_mass.copy()

    for i in range(len(sample)):
        halo_i = sample[i]
        print(i)

        # Grab the start position in the particles file to read from
        halo_start_position = group_file["Offset"][halo_i]
        halo_end_position = group_file["Offset"][halo_i + 1]
        particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]
        particle_ids_from_snapshot = snapshot_file["PartType1/ParticleIDs"][...]

        _, indices_v, indices_p = np.intersect1d(particle_ids_in_halo,
                                                 particle_ids_from_snapshot,
                                                 assume_unique=True,
                                                 return_indices=True, )

        all_pos = particles_pos[indices_p, :].copy()
        all_pos -= CoP[i, :]  # centering
        all_pos *= 1e3  # kpc

        if len(indices_p) < 10: continue

        radius = np.linalg.norm(all_pos[:, :3], axis=1)  # computing distances to the CoP
        vrot = analyse_rotation(mass[indices_p], radius)  #
        f = interpolate.interp1d(centers, vrot)

        v_circ_1[i] = f(1.0)
        v_circ_2[i] = f(2.0)
        v_circ_5[i] = f(5.0)
        Vmax[i] = np.max(vrot)

    # Output data
    filename = f"{siminfo.output_path}/Rotation_data_" + siminfo.name + ".hdf5"
    data_file = h5py.File(filename, 'w')
    f = data_file.create_group('Data')
    MH = f.create_dataset('ID', data=ID)
    MH = f.create_dataset('StructureType', data=Type)
    MH = f.create_dataset('M200c', data=m200c)
    MH = f.create_dataset('c200c', data=c200c)
    MH = f.create_dataset('Vcirc1kpc', data=v_circ_1)
    MH = f.create_dataset('Vcirc2kpc', data=v_circ_2)
    MH = f.create_dataset('Vcirc5kpc', data=v_circ_5)
    MH = f.create_dataset('Vmax', data=Vmax)
    data_file.close()
    return



def plot_rotation_data(siminfo, galaxy_type):
    filename = f"{siminfo.output_path}/Rotation_data_" + siminfo.name + ".hdf5"

    with h5py.File(filename, "r") as file:
        M200c = file["Data/M200c"][:]
        Type = file["Data/StructureType"][:]
        v_circ_1 = file["Data/Vcirc1kpc"][:]
        v_circ_2 = file["Data/Vcirc2kpc"][:]
        v_circ_5 = file["Data/Vcirc5kpc"][:]
        Vmax = file["Data/Vmax"][:]

    if galaxy_type == "centrals":
        select = Type == 10
    else :
        select = Type > 10

    v_circ_1 = v_circ_1[select]
    v_circ_2 = v_circ_2[select]
    v_circ_5 = v_circ_5[select]
    Vmax = Vmax[select]
    M200c = M200c[select]

    mass_range = np.arange(6, 16, 0.2)
    c200_range = c_M_relation(10 ** mass_range)
    NFW_circ_1kpc = np.zeros(len(mass_range))
    NFW_circ_2kpc = np.zeros(len(mass_range))
    NFW_circ_5kpc = np.zeros(len(mass_range))
    NFW_circ_max = np.zeros(len(mass_range))
    radius = np.arange(0.1, 20, 0.1)
    for i in range(len(mass_range)):
        NFW_circ = calc_NFW_vcirc(radius, 10 ** mass_range[i], c200_range[i], 0, siminfo)
        f = interpolate.interp1d(radius, NFW_circ)
        NFW_circ_1kpc[i] = f(1.0)
        NFW_circ_2kpc[i] = f(2.0)
        NFW_circ_5kpc[i] = f(5.0)
        NFW_circ_max[i] = np.max(NFW_circ)

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
        "lines.markersize": 1.5,
        "lines.linewidth": 1.5,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    yrange = np.array([10,300])
    plt.plot(yrange,yrange,'--',lw=1,color='grey')
    plt.plot(NFW_circ_max, NFW_circ_1kpc, '-', lw=2, color='black')

    Vmax_range = np.arange(0,100,5)
    v_circ_median = stat.binned_statistic(x=Vmax, values=v_circ_1, statistic="median", bins=Vmax_range,)[0]
    v_circ_16 = stat.binned_statistic(x=Vmax, values=v_circ_1,
                                      statistic=lambda v_circ_1: percentile(v_circ_1, 16), bins=Vmax_range,)[0]
    v_circ_84 = stat.binned_statistic(x=Vmax, values=v_circ_1,
                                      statistic=lambda v_circ_1: percentile(v_circ_1, 84), bins=Vmax_range,)[0]
    Vmax_center = bin_centers(Vmax_range)
    plt.plot(Vmax_center, v_circ_median, '-')
    plt.fill_between(Vmax_center, v_circ_16, v_circ_84, alpha=0.4)

    plt.plot(Vmax, v_circ_1, 'o')

    plt.axis([0, 100, 0, 100])
    plt.xlabel("Maximum Circular Velocity [km/s]")
    plt.ylabel("Circular velocity at 1kpc [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_"+siminfo.name+"_1kpc_centrals.png", dpi=200)
    else :
        plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_"+siminfo.name+"_1kpc_satellites.png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(yrange,yrange,'--',lw=1,color='grey')
    plt.plot(NFW_circ_max, NFW_circ_2kpc, '-', lw=2, color='black')

    v_circ_median = stat.binned_statistic(x=Vmax, values=v_circ_2, statistic="median", bins=Vmax_range,)[0]
    v_circ_16 = stat.binned_statistic(x=Vmax, values=v_circ_2,
                                      statistic=lambda v_circ_2: percentile(v_circ_2, 16), bins=Vmax_range,)[0]
    v_circ_84 = stat.binned_statistic(x=Vmax, values=v_circ_2,
                                      statistic=lambda v_circ_2: percentile(v_circ_2, 84), bins=Vmax_range,)[0]
    plt.plot(Vmax_center, v_circ_median, '-')
    plt.fill_between(Vmax_center, v_circ_16, v_circ_84, alpha=0.4)

    plt.plot(Vmax, v_circ_2, 'o')

    plt.axis([0, 100, 0, 100])
    plt.xlabel("Maximum Circular Velocity [km/s]")
    plt.ylabel("Circular velocity at 2kpc [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_"+siminfo.name+"_2kpc_centrals.png", dpi=200)
    else :
        plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_"+siminfo.name+"_2kpc_satellites.png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(yrange,yrange,'--',lw=1,color='grey')
    plt.plot(NFW_circ_max, NFW_circ_5kpc, '-', lw=2, color='black')

    v_circ_median = stat.binned_statistic(x=Vmax, values=v_circ_5, statistic="median", bins=Vmax_range,)[0]
    v_circ_16 = stat.binned_statistic(x=Vmax, values=v_circ_5,
                                      statistic=lambda v_circ_5: percentile(v_circ_5, 16), bins=Vmax_range,)[0]
    v_circ_84 = stat.binned_statistic(x=Vmax, values=v_circ_5,
                                      statistic=lambda v_circ_5: percentile(v_circ_5, 84), bins=Vmax_range,)[0]
    plt.plot(Vmax_center, v_circ_median, '-')
    plt.fill_between(Vmax_center, v_circ_16, v_circ_84, alpha=0.4)

    plt.plot(Vmax, v_circ_5, 'o')

    plt.axis([0, 100, 0, 100])
    plt.xlabel("Maximum Circular Velocity [km/s]")
    plt.ylabel("Circular velocity at 5kpc [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_"+siminfo.name+"_5kpc_centrals.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/Vmax_Vcirc_"+siminfo.name+"_5kpc_satellites.png", dpi=200)
    plt.close()

    ###################

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(M200c, v_circ_1, 'o')
    v_median, _, _ = stat.binned_statistic(x=M200c, values=v_circ_1, statistic="median", bins=mass_range,)
    mass = bin_centers(mass_range)
    plt.plot(mass, v_median, '-',lw=2,color='tab:blue')
    plt.plot(mass_range, NFW_circ_1kpc,'--',lw=2,color='black')

    plt.axis([8, 11, 0, 50])
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel("$\log_{10}$ M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Circular velocity at 1kpc [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/M200_Vcirc_"+siminfo.name+"_1kpc_centrals.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/M200_Vcirc_"+siminfo.name+"_1kpc_satellites.png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(M200c, v_circ_2, 'o')
    v_median, _, _ = stat.binned_statistic(x=M200c, values=v_circ_2, statistic="median", bins=mass_range,)
    plt.plot(mass, v_median, '-',lw=2,color='tab:blue')
    plt.plot(mass_range, NFW_circ_2kpc,'--',lw=2,color='black')

    plt.axis([8, 11, 0, 50])
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel("$\log_{10}$ M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Circular velocity at 2kpc [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/M200_Vcirc_"+siminfo.name+"_2kpc_centrals.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/M200_Vcirc_"+siminfo.name+"_2kpc_satellites.png", dpi=200)
    plt.close()

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(M200c, v_circ_5, 'o')
    v_median, _, _ = stat.binned_statistic(x=M200c, values=v_circ_5, statistic="median", bins=mass_range,)
    plt.plot(mass, v_median, '-',lw=2,color='tab:blue')

    plt.plot(mass_range, NFW_circ_5kpc,'--',lw=2,color='black')

    plt.axis([8, 11, 10, 100])
    #plt.xscale('log')
    #plt.yscale('log')
    plt.xlabel("$\log_{10}$ M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Circular velocity at 5kpc [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/M200_Vcirc_"+siminfo.name+"_5kpc_centrals.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/M200_Vcirc_"+siminfo.name+"_5kpc_satellites.png", dpi=200)
    plt.close()


def plot_rotation_of_sample(siminfo, galaxy_type):
    filename = f"{siminfo.output_path}/Rotation_data_" + siminfo.name + ".hdf5"

    with h5py.File(filename, "r") as file:
        Type = file["Data/StructureType"][:]
        Vmax = file["Data/Vmax"][:]
        ID = file["Data/ID"][:]

    if galaxy_type == "centrals":
        select = Type == 10
    else :
        select = Type > 10

    Vmax = Vmax[select]
    ID = ID[select]

    select = np.where((Vmax >= 80) & (Vmax <= 85))[0]
    ID_sample = ID[select]
    print(len(ID_sample))

    with h5py.File(siminfo.snapshot, "r") as hf:
        DMpart_mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        DMpart_pos = hf['PartType1/Coordinates'][:][:] * siminfo.a # Mpc

    with h5py.File(siminfo.halo_properties, "r") as properties_file:
        m200c = properties_file["Mass_200crit"][:] * 1e10 # Msun
        c200c = properties_file["cNFW_200crit"][:]
        xCoP = properties_file["Xcminpot"][:] # Mpc
        yCoP = properties_file["Ycminpot"][:] # Mpc
        zCoP = properties_file["Zcminpot"][:] # Mpc
        ID = properties_file["ID"][:]

    _, index_j, sample = np.intersect1d(ID_sample, ID, assume_unique=True, return_indices=True, )

    CoP = np.zeros((len(sample), 3))
    CoP[:, 0] = xCoP[sample]
    CoP[:, 1] = yCoP[sample]
    CoP[:, 2] = zCoP[sample]


    Mmedian = np.median(m200c[sample])
    cmedian = np.median(c200c[sample])

    print(len(sample), np.log10(Mmedian), cmedian)

    # Reduce
    radial_bins = np.arange(0.2, 25, 0.25)
    centers = bin_centers(radial_bins) # kpc
    v_circ = np.zeros((len(centers), len(sample)))

    snapshot_file = h5py.File(siminfo.snapshot, "r")
    group_file = h5py.File(siminfo.catalog_groups, "r")
    particles_file = h5py.File(siminfo.catalog_particles, "r")

    particles_pos = DMpart_pos.copy()
    mass = DMpart_mass.copy()

    for i in range(len(sample)):
        halo_i = sample[i]
        print(i)

        # Grab the start position in the particles file to read from
        halo_start_position = group_file["Offset"][halo_i]
        halo_end_position = group_file["Offset"][halo_i + 1]
        particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]
        particle_ids_from_snapshot = snapshot_file["PartType1/ParticleIDs"][...]

        _, indices_v, indices_p = np.intersect1d(particle_ids_in_halo,
                                                 particle_ids_from_snapshot,
                                                 assume_unique=True,
                                                 return_indices=True, )

        all_pos = particles_pos[indices_p,:].copy()
        all_pos -= CoP[i, :]  # centering
        all_pos *= 1e3  # kpc

        if len(indices_p) < 10: continue

        radius = np.linalg.norm(all_pos[:, :3], axis=1)  # computing distances to the CoP
        vrot = analyse_rotation(mass[indices_p], radius)  #
        v_circ[:, i] = vrot


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
        plt.plot(centers, v_circ[:, i], '-', label='M$_{200}$=%.3f [$\log_{10}$M$_{\odot}$]'%m200c[i])

    NFW_circ = calc_NFW_vcirc(centers, Mmedian, cmedian, siminfo.z, siminfo)
    plt.plot(centers, NFW_circ,'--',lw=2,color='black')

    plt.axis([0, 20, 0, 100])
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Circular velocity [km/s]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    if galaxy_type == "centrals":
        plt.savefig(f"{siminfo.output_path}/Rotation_curve_Vmax_"+siminfo.name+"_centrals.png", dpi=200)
    else:
        plt.savefig(f"{siminfo.output_path}/Rotation_curve_Vmax_"+siminfo.name+"_satellites.png", dpi=200)
    plt.close()
