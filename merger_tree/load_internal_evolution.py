import numpy as np
from halo_data import density_profile
from .snapshot_data import snapshot_info
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
    density = np.zeros((len(centered_radial_bins),initial_snap - final_snap))

    vel_radial_bins = np.arange(0.2, 25, 0.25)
    centered_velocity_radial_bins = bin_centers(vel_radial_bins)  # kpc
    velocity = np.zeros((len(centered_velocity_radial_bins),initial_snap - final_snap))

    density[:,0], velocity[:,0] = density_profile.calculate_halo_data(sim_info, progenitor_index[0],
                                                                      radial_bins, vel_radial_bins)

    for snap in range(initial_snap-1,final_snap,-1):

        i = initial_snap - snap
        print('Reading snapshot data', snap)
        start_time = time.time()

        snapshot_data = snapshot_info(sim_info, snap)
        density[:, i], velocity[:, i] = \
            density_profile.calculate_halo_data(snapshot_data, progenitor_index[i], radial_bins, vel_radial_bins)

        print("--- %s seconds ---" % (time.time() - start_time))

        if np.sum(density[:, initial_snap-snap]) == 0.:break

    evolution_data = {'density': density,'velocity': velocity,
                      'density_radial_bins': centered_radial_bins,
                      'velocity_radial_bins': centered_velocity_radial_bins}

    return evolution_data