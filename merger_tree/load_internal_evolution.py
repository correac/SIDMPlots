import numpy as np
from halo_data import density_profile, dynamical_relaxation
from .snapshot_data import snapshot_info
import h5py
import time

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)


def load_internal_evolution(sim_info, progenitor_index, output_file):

    initial_snap = sim_info.initial_snap
    final_snap = 10

    # Related to internal structure evolution
    radial_bins = np.arange(-1, 3, 0.1)
    radial_bins = 10**radial_bins
    centered_radial_bins = bin_centers(radial_bins)  # kpc

    vel_radial_bins = np.arange(0.1, 25, 0.25)
    centered_velocity_radial_bins = bin_centers(vel_radial_bins)  # kpc

    density, velocity = \
        density_profile.calculate_halo_data(sim_info, progenitor_index[:,0], radial_bins, vel_radial_bins)

    data_file = h5py.File(output_file, 'a')
    f = data_file.create_group('Profile_evolution')
    f.create_dataset('Density_snapshot_%04i'%initial_snap, data=density)
    f.create_dataset('Velocity_snapshot_%04i'%initial_snap, data=velocity)

    for snap in range(initial_snap-1,final_snap,-1):

        print('Reading snapshot data', snap)
        start_time = time.time()

        i = initial_snap - snap
        snapshot_data = snapshot_info(sim_info, snap)
        density, velocity = \
            density_profile.calculate_halo_data(snapshot_data, progenitor_index[:,i], radial_bins, vel_radial_bins)

        print("--- %s seconds ---" % (time.time() - start_time))

        f.create_dataset('Density_snapshot_%04i' % snap, data=density)
        f.create_dataset('Velocity_snapshot_%04i' % snap, data=velocity)

    f.create_dataset('Density_radial_bins', data=centered_radial_bins)
    f.create_dataset('Velocity_radial_bins', data=centered_velocity_radial_bins)
    data_file.close()

    return


def load_velocity_dispersion(sim_info, progenitor_index, output_file):

    initial_snap = sim_info.initial_snap
    final_snap = 10

    # Related to internal structure evolution
    radial_bins = np.arange(-1, 3, 0.1)
    radial_bins = 10**radial_bins
    centered_radial_bins = bin_centers(radial_bins)  # kpc

    veldisp, sigmaprofile = density_profile.calculate_velocity_dispersion(sim_info, progenitor_index[:,0], radial_bins)

    snapshot_data = snapshot_info(sim_info, initial_snap)
    relaxation = dynamical_relaxation.calculate_relaxation_criteria(snapshot_data, progenitor_index[:,0])

    data_file = h5py.File(output_file, 'a')
    f = data_file.create_group('Additional_data')
    f.create_dataset('Velocity_Dispersion_snapshot_%04i'%initial_snap, data=veldisp)
    f.create_dataset('Sigma_profile_snapshot_%04i'%initial_snap, data=sigmaprofile)
    f.create_dataset('Dynamical_relaxation_snapshot_%04i'%initial_snap, data=relaxation)

    for snap in range(initial_snap-1,final_snap,-1):

        print('Reading snapshot data', snap)
        start_time = time.time()

        i = initial_snap - snap
        snapshot_data = snapshot_info(sim_info, snap)

        veldisp, sigmaprofile = \
            density_profile.calculate_velocity_dispersion(snapshot_data, progenitor_index[:, i], radial_bins)

        relaxation = dynamical_relaxation.calculate_relaxation_criteria(snapshot_data, progenitor_index[:, i])

        print("--- %s seconds ---" % (time.time() - start_time))

        f.create_dataset('Velocity_Dispersion_snapshot_%04i' % snap, data=veldisp)
        f.create_dataset('Sigma_profile_snapshot_%04i' % snap, data=sigmaprofile)
        f.create_dataset('Dynamical_relaxation_snapshot_%04i' % snap, data=relaxation)

    f.create_dataset('VelDisp_radial_bins', data=centered_radial_bins)
    data_file.close()

    return


