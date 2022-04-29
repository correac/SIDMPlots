import numpy as np
import h5py
from tqdm import tqdm

from .build_tree import build_tree
from .load_internal_evolution import load_internal_evolution

def make_tree_data(sim_info):

    select_sub_sample = np.where(
        (sim_info.halo_data.log10_halo_mass > 9.5) &
        (sim_info.halo_data.log10_halo_mass <= 10.0))[0]

    select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] > 10)[0]

    sample = select_sub_sample[select_type]
    halo_index = sim_info.halo_data.halo_index[sample]

    # Output data
    output_file = f"{sim_info.output_path}/Tree_data_" + sim_info.simulation_name + "_satellites_95_10.hdf5"
    progenitor_index = build_tree(sim_info, halo_index, output_file)

    load_internal_evolution(sim_info, progenitor_index, output_file)
