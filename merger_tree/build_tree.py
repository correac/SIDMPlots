import numpy as np
import velociraptor
import h5py
from .snapshot_data import snapshot_info
from .tree_dataset import TreeCatalogue
import time

def build_tree(sim_info, halo_index):

    initial_snap = sim_info.initial_snap
    final_snap = 10

    # Let's collect some data from the halo that we are following,
    progenitor_index = np.zeros(initial_snap - final_snap)
    progenitor_index[0] = halo_index
    z = np.zeros(initial_snap - final_snap)

    # Related to mass assembly history
    merger_mass_ratio = np.zeros(initial_snap - final_snap)
    mass = np.zeros(initial_snap - final_snap)
    type = np.zeros(initial_snap - final_snap)
    host_distance = np.zeros(initial_snap - final_snap)

    catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_base_name}" + "_%04i.properties" % initial_snap
    catalogue = velociraptor.load(catalogue_file)
    m200c = catalogue.masses.mass_200crit.to("Msun").value
    mass_descendant = m200c[halo_index]
    z[0] = catalogue.z
    mass[0] = m200c[halo_index]
    type[0] = catalogue.structure_type.structuretype[halo_index]

    #print("read z=0 tree data")
    start_time = time.time()
    tree_file = f"{sim_info.directory}/merger_tree/MergerTree.snapshot_0%i.VELOCIraptor.tree" % initial_snap
    tree_data = TreeCatalogue(tree_file)

    halo = tree_data.catalogue.ProgenOffsets.value[halo_index]
    num_progenitors = tree_data.catalogue.NumProgen.value[halo_index]
    progenitor_list = np.arange(num_progenitors)+halo
    proID = tree_data.catalogue.Progenitors.value[progenitor_list]
    #print("--- %s seconds ---" % (time.time() - start_time))

    for snap in range(initial_snap-1,final_snap,-1):

        i = initial_snap - snap
        #print('snapshot', snap)
        snapshot_data = snapshot_info(sim_info, snap)
        path_to_catalogue_file = f"{snapshot_data.directory}/{snapshot_data.catalogue_name}"
        catalogue = velociraptor.load(path_to_catalogue_file)
        m200c = catalogue.masses.mass_200crit.to("Msun").value
        z[i] = catalogue.z

        #print("read z=0 tree data")
        start_time = time.time()
        tree_file = f"{sim_info.directory}/merger_tree/MergerTree.snapshot_0%i.VELOCIraptor.tree" % snap
        tree_data = TreeCatalogue(tree_file)

        #print("--- %s seconds ---" % (time.time() - start_time))

        #print("Look for progenitor")
        start_time = time.time()
        if num_progenitors > 10: proID = proID[0:10]

        _, indx_ID, indx_proID = np.intersect1d(tree_data.catalogue.ID.value,
                                                proID, assume_unique=True, return_indices=True, )
        if len(indx_ID) == 0: break

        largest_mass_progenitor = np.where(m200c[indx_ID] == np.max(m200c[indx_ID]))[0]
        other_progenitors = np.where(m200c[indx_ID] != np.max(m200c[indx_ID]))[0]

        if len(largest_mass_progenitor) > 1:
            other_progenitors = largest_mass_progenitor[1:]
            largest_mass_progenitor = largest_mass_progenitor[0]

        connect = tree_data.catalogue.ProgenOffsets.value[indx_ID[largest_mass_progenitor]]
        num_progenitors = tree_data.catalogue.NumProgen.value[indx_ID[largest_mass_progenitor]]
        progenitor_list = np.arange(num_progenitors) + connect
        proID = tree_data.catalogue.Progenitors.value[progenitor_list.astype('int')]
        progenitor_index[i] = indx_ID[largest_mass_progenitor][0]

        if len(indx_ID) > 1:
            merger_mass_ratio[i] = np.max(m200c[indx_ID[other_progenitors]] / mass_descendant)

        mass_descendant = m200c[progenitor_index[i].astype('int')]
        mass[i] = m200c[progenitor_index[i].astype('int')]
        type[i] = catalogue.structure_type.structuretype[progenitor_index[i].astype('int')]
        #print("--- %s seconds ---" % (time.time() - start_time))

        if num_progenitors == 0: break

    tree_data = {'progenitor_index': progenitor_index,
                 'merger_mass_ratio': merger_mass_ratio,
                 'M200crit': mass,
                 'structure_type': type,
                 'redshift': z}

    return tree_data
