import numpy as np
import h5py
from tqdm import tqdm

from .build_tree import build_tree
from .load_internal_evolution import load_internal_evolution

def make_tree_data(sim_info):

    select_sub_sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= 9) &
        (sim_info.halo_data.log10_halo_mass <= 9.01))[0]

    select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] > 10)[0]

    sample = select_sub_sample[select_type]
    num_halos = len(sample)
    halo_index = sim_info.halo_data.halo_index[sample]

    # Output data
    filename = f"{sim_info.output_path}/Tree_data_" + sim_info.simulation_name + ".hdf5"
    data_file = h5py.File(filename, 'w')
    f = data_file.create_group('Data')
    f.create_dataset('ID', data=halo_index)

    for i in tqdm(range(num_halos)):

        tree_data = build_tree(sim_info, halo_index[i])

        # Write data to file while it is being calculated..
        f.create_dataset('Mass_%04i'%i, data=tree_data['M200crit'])
        f.create_dataset('Redshift_%04i'%i, data=tree_data['redshift'])
        f.create_dataset('Structure_Type_%04i'%i, data=tree_data['structure_type'])
        f.create_dataset('Merger_mass_ratio_%04i'%i, data=tree_data['merger_mass_ratio'])
        f.create_dataset('Progenitor_index_%04i'%i, data=tree_data['progenitor_index'])

        # evolution_data = load_internal_evolution(sim_info, tree_data['progenitor_index'])
        # f.create_dataset('Density_%04i'%i, data=evolution_data['density'])
        # f.create_dataset('Velocity_%04i'%i, data=evolution_data['velocity'])


    # f.create_dataset('Density_radial_bins', data=evolution_data['density_radial_bins'])
    # f.create_dataset('Velocity_radial_bins', data=evolution_data['velocity_radial_bins'])
    data_file.close()
