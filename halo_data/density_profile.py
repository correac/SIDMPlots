import scipy.stats as stat
import numpy as np
from tqdm import tqdm
from object import particle_data
from .circular_velocity import calculate_Vcirc, calculate_Vcirc_all_part_type

class make_profile_data:

    def __init__(self, log10_M200, rs, centers,
                 density, sig_density,
                 velocity, sig_velocity,
                 sigma, sig_sigma,
                 structure_type):
        self.centers = centers
        self.log10_M200 = log10_M200
        self.rs = rs
        self.density = density
        self.sig_density = sig_density
        self.velocity = velocity
        self.sig_velocity = sig_velocity
        self.sigma = sigma
        self.sig_sigma = sig_sigma
        self.structure_type = structure_type

    def add_data(self, density, sig_density, velocity, sig_velocity, sigma, sig_sigma):
        self.density = np.append(self.density, density)
        self.sig_density = np.append(self.sig_density, sig_density)
        self.velocity = np.append(self.velocity, velocity)
        self.sig_velocity = np.append(self.sig_velocity, sig_velocity)
        self.sigma = np.append(self.sigma, sigma)
        self.sig_sigma = np.append(self.sig_sigma, sig_sigma)


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

def calculate_profiles(mass, pos, vel, sigma, radial_bins):

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    SumMasses, _, _ = stat.binned_statistic(x=r, values=mass, statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3

    std_vel_x, _, _ = stat.binned_statistic(x=r, values=vel[:, 0], statistic="std", bins=radial_bins, )
    std_vel_y, _, _ = stat.binned_statistic(x=r, values=vel[:, 1], statistic="std", bins=radial_bins, )
    std_vel_z, _, _ = stat.binned_statistic(x=r, values=vel[:, 2], statistic="std", bins=radial_bins, )
    velocity = np.sqrt(std_vel_x ** 2 + std_vel_y ** 2 + std_vel_z ** 2) / np.sqrt(3.)
    velocity[np.where(np.isnan(velocity))[0]] = 0

    sigma_profile, _, _ = stat.binned_statistic(x=r, values=sigma, statistic="median", bins=radial_bins, )

    return density, velocity, sigma_profile


def calculate_profiles_hydro(mass, pos, vel, radial_bins):

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    SumMasses, _, _ = stat.binned_statistic(x=r, values=mass, statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3

    std_vel_x, _, _ = stat.binned_statistic(x=r, values=vel[:, 0], statistic="std", bins=radial_bins, )
    std_vel_y, _, _ = stat.binned_statistic(x=r, values=vel[:, 1], statistic="std", bins=radial_bins, )
    std_vel_z, _, _ = stat.binned_statistic(x=r, values=vel[:, 2], statistic="std", bins=radial_bins, )
    velocity = np.sqrt(std_vel_x ** 2 + std_vel_y ** 2 + std_vel_z ** 2) / np.sqrt(3.)
    velocity[np.where(np.isnan(velocity))[0]] = 0

    return density, velocity

def calculate_profiles_all_part_type(
        star_mass, star_pos, gas_mass, gas_pos, dm_mass, dm_pos, density_bins, velocity_bins):

    # Radial coordinates [kpc units]
    r_stars = np.sqrt(np.sum(star_pos ** 2, axis=1))
    r_gas = np.sqrt(np.sum(gas_pos ** 2, axis=1))
    r_dm = np.sqrt(np.sum(dm_pos ** 2, axis=1))

    r = np.array([r_stars, r_gas, r_dm])
    mass = np.array([star_mass, gas_mass, dm_mass])

    SumMasses, _, _ = stat.binned_statistic(x=r, values=mass, statistic="sum", bins=density_bins, )
    density = (SumMasses / bin_volumes(density_bins))  # Msun/kpc^3

    velocity = calculate_Vcirc_all_part_type(mass, r, velocity_bins)
    return density, velocity

def compute_density_profiles(sim_info, log10_min_mass, log10_max_mass, structure_type, profile_data):

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(-1, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    select_sub_sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= log10_min_mass) &
        (sim_info.halo_data.log10_halo_mass <= log10_max_mass))[0]

    if structure_type == 10:
        select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] == structure_type)[0]
    else:
        select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] > 10)[0]

    sample = select_sub_sample[select_type]
    log10_M200 = np.median(sim_info.halo_data.log10_halo_mass[sample])
    rs = np.median(sim_info.halo_data.scale_radius[sample]) * 1e3  # kpc

    num_halos = len(sample)
    sigma_all = np.zeros((len(centers), num_halos))
    density_all = np.zeros((len(centers), num_halos))
    velocity_all = np.zeros((len(centers), num_halos))

    for i in tqdm(range(num_halos)):

        halo_indx = sim_info.halo_data.halo_index[sample[i]]
        part_data = particle_data.load_particle_data(sim_info, halo_indx, sample[i])

        num_part = len(part_data.masses.value[part_data.bound_particles_only])
        if num_part < 10: continue

        density, velocity, sigma = calculate_profiles(part_data.masses.value[part_data.bound_particles_only],
                                                      part_data.coordinates.value[part_data.bound_particles_only, :],
                                                      part_data.velocities.value[part_data.bound_particles_only],
                                                      part_data.cross_section[part_data.bound_particles_only],
                                                      radial_bins)

        sigma_all[:, i] = sigma
        density_all[:, i] = density
        velocity_all[:, i] = velocity

    sigma = np.mean(sigma_all[:, :], axis=1)
    sig_sigma = np.std(sigma_all[:, :], axis=1)
    density = np.mean(density_all[:, :], axis=1)
    sig_density = np.std(density_all[:, :], axis=1)
    velocity = np.mean(velocity_all[:, :], axis=1)
    sig_velocity = np.std(velocity_all[:, :], axis=1)

    if profile_data == None:
        profile_data = make_profile_data(log10_M200, rs, centers, density, sig_density,
                                         velocity, sig_velocity,
                                         sigma, sig_sigma,
                                         structure_type)
    else:
        profile_data.add_data(density, sig_density, velocity, sig_velocity, sigma, sig_sigma)

    return profile_data


def calculate_halo_data(sim_info, halo_index, density_radial_bins, velocity_radial_bins):

    num_haloes = len(halo_index)

    centered_radial_bins = bin_centers(density_radial_bins)  # kpc
    density = np.zeros((len(centered_radial_bins), num_haloes))
    veldisp = np.zeros((len(centered_radial_bins), num_haloes))
    sigmaprofile = np.zeros((len(centered_radial_bins), num_haloes))

    centered_velocity_radial_bins = bin_centers(velocity_radial_bins)  # kpc
    velocity = np.zeros((len(centered_velocity_radial_bins), num_haloes))

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.dark_matter.ids)

        if len(bound_particles_only) < 10: continue

        density[:,i], veldisp[:,i], sigmaprofile[:,i] = \
            calculate_profiles(part_data.dark_matter.masses.value[bound_particles_only],
                               part_data.dark_matter.coordinates.value[bound_particles_only, :],
                               part_data.dark_matter.velocities.value[bound_particles_only],
                               part_data.dark_matter.cross_section[bound_particles_only],
                               density_radial_bins)

        velocity[:,i] = \
            calculate_Vcirc(part_data.dark_matter.masses.value[bound_particles_only],
                            part_data.dark_matter.coordinates.value[bound_particles_only, :],
                            velocity_radial_bins)

    return density, velocity, veldisp, sigmaprofile


def calculate_halo_data_hydro(sim_info, halo_index, density_radial_bins, velocity_radial_bins):

    num_haloes = len(halo_index)

    centered_radial_bins = bin_centers(density_radial_bins)  # kpc
    stars_density = np.zeros((len(centered_radial_bins), num_haloes))
    stars_veldisp = np.zeros((len(centered_radial_bins), num_haloes))
    gas_density = np.zeros((len(centered_radial_bins), num_haloes))
    gas_veldisp = np.zeros((len(centered_radial_bins), num_haloes))
    density = np.zeros((len(centered_radial_bins), num_haloes))

    centered_velocity_radial_bins = bin_centers(velocity_radial_bins)  # kpc
    stars_velocity = np.zeros((len(centered_velocity_radial_bins), num_haloes))
    gas_velocity = np.zeros((len(centered_velocity_radial_bins), num_haloes))
    velocity = np.zeros((len(centered_velocity_radial_bins), num_haloes))

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        stars_bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.stars.ids)
        gas_bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.gas.ids)
        dm_bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.dark_matter.ids)

        if (len(stars_bound_particles_only) < 10) or (len(gas_bound_particles_only) < 10): continue

        stars_density[:,i], stars_veldisp[:,i] = \
            calculate_profiles_hydro(part_data.stars.masses.value[stars_bound_particles_only],
                                     part_data.stars.coordinates.value[stars_bound_particles_only, :],
                                     part_data.stars.velocities.value[stars_bound_particles_only],
                                     density_radial_bins)

        gas_density[:, i], gas_veldisp[:, i] = \
            calculate_profiles_hydro(part_data.gas.masses.value[gas_bound_particles_only],
                                     part_data.gas.coordinates.value[gas_bound_particles_only, :],
                                     part_data.gas.velocities.value[gas_bound_particles_only],
                                     density_radial_bins)

        stars_velocity[:,i] = \
            calculate_Vcirc(part_data.stars.masses.value[stars_bound_particles_only],
                            part_data.stars.coordinates.value[stars_bound_particles_only, :],
                            velocity_radial_bins)

        gas_velocity[:,i] = \
            calculate_Vcirc(part_data.gas.masses.value[gas_bound_particles_only],
                            part_data.gas.coordinates.value[gas_bound_particles_only, :],
                            velocity_radial_bins)

        density[:,i], velocity[:,i] = \
            calculate_profiles_all_part_type(
                part_data.stars.masses.value[stars_bound_particles_only],
                part_data.stars.coordinates.value[stars_bound_particles_only, :],
                part_data.gas.masses.value[gas_bound_particles_only],
                part_data.gas.coordinates.value[gas_bound_particles_only, :],
                part_data.dark_matter.masses.value[dm_bound_particles_only],
                part_data.dark_matter.coordinates.value[dm_bound_particles_only, :],
                density_radial_bins, velocity_radial_bins)


    Hydro_data = {'stars_density': stars_density,
                  'stars_velocity': stars_velocity,
                  'stars_veldisp': stars_veldisp,
                  'gas_density': gas_density,
                  'gas_velocity': gas_velocity,
                  'gas_veldisp': gas_veldisp,
                  'density': density,
                  'velocity': velocity
                  }

    return Hydro_data


def calculate_velocity_dispersion(sim_info, halo_index, density_radial_bins):

    num_haloes = len(halo_index)

    centered_radial_bins = bin_centers(density_radial_bins)  # kpc
    veldisp = np.zeros((len(centered_radial_bins), num_haloes))
    sigmaprofile = np.zeros((len(centered_radial_bins), num_haloes))

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        if len(part_data.bound_particles_only) < 10: continue

        _, veldisp[:,i], sigmaprofile[:,i] = calculate_profiles(part_data.masses.value[part_data.bound_particles_only],
                                            part_data.coordinates.value[part_data.bound_particles_only, :],
                                            part_data.velocities.value[part_data.bound_particles_only],
                                            part_data.cross_section[part_data.bound_particles_only],
                                            density_radial_bins)

    return veldisp, sigmaprofile
