import numpy as np
import h5py
import scipy.stats as stat
from output_data.analyse_tree import calculate_ratio_with_M200c
from argumentparser import ArgumentParser
from object import simulation_data
from merger_tree.snapshot_data import snapshot_info
import time
from tqdm import tqdm
from object import particle_data
from halo_data.circular_velocity import calculate_Vcirc
from halo_data.density_profile import calculate_profiles

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)

def look_for_high_ratio_haloes(Mmin, Mmax, Type, sim_name):

    ratio, _, halo_indx_rotation_data, _ = calculate_ratio_with_M200c(Mmin, Mmax, Type, sim_name)
    select_high = np.where(ratio >= 0.3)[0]
    return halo_indx_rotation_data[select_high]


def clean_sample(mask, name):

    filename = "./output_data/Tree_data_" + name + "_Snap36.hdf5"
    with h5py.File(filename, "r") as file:
        mass = file["Assembly_history/Mass"][:][:]

    select = []
    for i in range(len(mask)):
        mass_i = mass[mask[i], :]
        delta_m = (mass_i[:-1] - mass_i[1:]) / mass_i[1:]
        if np.max(np.abs(delta_m)) < 10: select = np.append(select, i)

    mask = mask[select.astype('int')]
    return mask

def load_progenitor_indx_list(halo_indx, name):

    filename = "./output_data/Tree_data_" + name + "_Snap36.hdf5"
    with h5py.File(filename, "r") as file:
        indx_hist = file["Assembly_history/ID"][:]
        progenitors = file["Assembly_history/Progenitor_index"][:][:]

    _, mask_halo, mask_progenitor = np.intersect1d(halo_indx, indx_hist, assume_unique=True, return_indices=True,)

    mask_progenitor = clean_sample(mask_progenitor, name)
    progenitors = progenitors[mask_progenitor,:]
    return progenitors


def calculate_number_of_collisions(num_sidm_events, pos, radial_bins):

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    sum_sidm_events, _, _ = stat.binned_statistic(x=r, values=num_sidm_events, statistic="sum", bins=radial_bins, )
    return sum_sidm_events

def calculate_halo_data(sim_info, halo_index, density_radial_bins, velocity_radial_bins):

    num_haloes = len(halo_index)

    centered_radial_bins = bin_centers(density_radial_bins)  # kpc
    density = np.zeros((len(centered_radial_bins), num_haloes))
    veldisp = np.zeros((len(centered_radial_bins), num_haloes))

    centered_velocity_radial_bins = bin_centers(velocity_radial_bins)  # kpc
    velocity = np.zeros((len(centered_velocity_radial_bins), num_haloes))
    numcoll = np.zeros((len(centered_velocity_radial_bins), num_haloes))

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        if len(part_data.bound_particles_only) < 10: continue

        density[:,i], veldisp[:,i], _ = calculate_profiles(part_data.masses.value[part_data.bound_particles_only],
                                            part_data.coordinates.value[part_data.bound_particles_only, :],
                                            part_data.velocities.value[part_data.bound_particles_only],
                                            part_data.cross_section[part_data.bound_particles_only],
                                            density_radial_bins)

        velocity[:,i] = calculate_Vcirc(part_data.masses.value[part_data.bound_particles_only],
                                   part_data.coordinates.value[part_data.bound_particles_only, :],
                                   velocity_radial_bins)

        numcoll[:,i] = calculate_number_of_collisions(part_data.num_sidm_events[part_data.bound_particles_only],
                                   part_data.coordinates.value[part_data.bound_particles_only, :],
                                   velocity_radial_bins)

    return density, velocity, veldisp, numcoll

def load_sample_evolution(progenitor_index, sim_info):

    output_file = './output_data/data/High_ratio_sample_Tree_'+ sim_info.simulation_name + '.hdf5'

    data_file = h5py.File(output_file, 'w')
    f = data_file.create_group('Assembly_history')
    f.create_dataset('Progenitor_index_list', data=progenitor_index)

    initial_snap = sim_info.initial_snap
    final_snap = 20

    # Related to internal structure evolution
    radial_bins = np.arange(-1, 3, 0.1)
    radial_bins = 10**radial_bins
    centered_radial_bins = bin_centers(radial_bins)  # kpc

    vel_radial_bins = np.arange(0.2, 25, 0.25)
    centered_velocity_radial_bins = bin_centers(vel_radial_bins)  # kpc

    density, velocity, velocity_dispersion, num_collisions = calculate_halo_data(sim_info, progenitor_index[:,0], radial_bins, vel_radial_bins)

    data_file = h5py.File(output_file, 'a')
    f = data_file.create_group('Profile_evolution')
    f.create_dataset('Density_snapshot_%04i'%initial_snap, data=density)
    f.create_dataset('Velocity_snapshot_%04i'%initial_snap, data=velocity)
    f.create_dataset('VelDisp_snapshot_%04i' % initial_snap, data=velocity_dispersion)
    f.create_dataset('NumCollisions_snapshot_%04i' % initial_snap, data=num_collisions)

    for snap in range(initial_snap-1,final_snap,-1):

        print('Reading snapshot data', snap)
        start_time = time.time()

        i = initial_snap - snap
        snapshot_data = snapshot_info(sim_info, snap)
        density, velocity, velocity_dispersion, num_collisions = calculate_halo_data(snapshot_data, progenitor_index[:,i], radial_bins, vel_radial_bins)

        print("--- %s seconds ---" % (time.time() - start_time))

        f.create_dataset('Density_snapshot_%04i' % snap, data=density)
        f.create_dataset('Velocity_snapshot_%04i' % snap, data=velocity)
        f.create_dataset('VelDisp_snapshot_%04i' % snap, data=velocity_dispersion)
        f.create_dataset('NumCollisions_snapshot_%04i' % snap, data=num_collisions)

    f.create_dataset('Density_radial_bins', data=centered_radial_bins)
    f.create_dataset('Velocity_radial_bins', data=centered_velocity_radial_bins)
    data_file.close()

    return


def main(config: ArgumentParser):

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

    Mmin = 9
    Mmax = 10
    Type = 15

    halo_indx = look_for_high_ratio_haloes(Mmin, Mmax, Type, sim_info.simulation_name)

    progenitor_list = load_progenitor_indx_list(halo_indx, sim_info.simulation_name)

    load_sample_evolution(progenitor_list, sim_info)


if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)


