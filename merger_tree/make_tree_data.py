import numpy as np
import h5py
from .build_tree import build_tree
from .load_internal_evolution import load_internal_evolution, load_velocity_dispersion

def make_tree_data(sim_info):

    select_sub_sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= 9.3) &
        (sim_info.halo_data.log10_halo_mass < 9.7))[0]


    #select_sub_sample = np.where(sim_info.halo_data.log10_halo_mass >= 9.0)[0]

    select_type = np.where(sim_info.halo_data.structure_type[select_sub_sample] == 10)[0]

    sample = select_sub_sample[select_type]
    halo_index = sim_info.halo_data.halo_index[sample]

    # Output data
    output_file = f"{sim_info.output_path}/Tree_data_Centrals_" + sim_info.simulation_name + "_93_97.hdf5"
    progenitor_index = build_tree(sim_info, halo_index, output_file)

    load_internal_evolution(sim_info, progenitor_index, output_file)


def add_tree_data(sim_info, input_file):

    # Output data

    output_file = sim_info.output_path+"/"+input_file+".hdf5"
    with h5py.File(output_file, "r") as file:
        progenitor_index = file["Assembly_history/Progenitor_index"][:][:]
        M200c = file["Assembly_history/Mass"][:][:]
        M200c = np.log10(M200c[:,0])

    select_sub_sample = np.where(M200c <= 9.5)[0]

    output_file = sim_info.output_path+"/"+input_file+"_90_95.hdf5" 
    load_velocity_dispersion(sim_info, progenitor_index[select_sub_sample,:], output_file)
