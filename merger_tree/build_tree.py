import numpy as np
import velociraptor
import h5py
from .snapshot_data import snapshot_info
from .tree_dataset import TreeCatalogue
import time
from tqdm import tqdm

def look_for_progenitor_index(progenitor_offset, num_progenitors, progenitors_ID, ID, m200c):

    num_haloes = len(num_progenitors)
    pro_indx = np.ones(num_haloes) * (-1)
    merger_ratio = np.zeros(num_haloes)

    for i in range(num_haloes):

        if num_progenitors[i] == 0:
            continue

        num = num_progenitors[i]
        if num > 10 : num = 10
        progenitor_list = np.arange(num) + progenitor_offset[i]
        proID = progenitors_ID[progenitor_list]

        _, indx_ID, _ = np.intersect1d(ID, proID, assume_unique=True, return_indices=True, )

        if len(indx_ID) == 0:
            continue

        largest_mass_progenitor = np.argmax(m200c[indx_ID])
        pro_indx[i] = indx_ID[largest_mass_progenitor]

        if len(indx_ID) > 1:
            other_progenitors = np.where(m200c[indx_ID] != m200c[pro_indx[i].astype('int')])[0]

            if len(other_progenitors) > 1:
                mass_progenitors = m200c[indx_ID[other_progenitors].astype('int')]
                ratio = mass_progenitors / m200c[pro_indx[i].astype('int')]
                merger_ratio[i] = np.max(ratio)


    return pro_indx, merger_ratio

def build_tree(sim_info, halo_index, output_file):

    data_file = h5py.File(output_file, 'w')
    f = data_file.create_group('Assembly_history')
    f.create_dataset('ID', data=halo_index)

    initial_snap = sim_info.initial_snap
    final_snap = 10

    num_haloes = len(halo_index)

    # Let's collect some data from the halo that we are following,
    progenitor_index = np.zeros((num_haloes,initial_snap - final_snap))
    progenitor_index[:, 0] = halo_index
    z = np.zeros(initial_snap - final_snap)

    # Related to mass assembly history
    merger_mass_ratio = np.zeros((num_haloes,initial_snap - final_snap))
    mass = np.zeros((num_haloes,initial_snap - final_snap))
    type = np.zeros((num_haloes,initial_snap - final_snap))
    #host_distance = np.zeros(initial_snap - final_snap)

    catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_base_name}" + "_%04i.properties" % initial_snap
    catalogue = velociraptor.load(catalogue_file)
    m200c = catalogue.masses.mass_200crit.to("Msun").value
    z[0] = catalogue.z
    mass[:, 0] = m200c[halo_index]
    type[:, 0] = catalogue.structure_type.structuretype[halo_index]

    #print("read z=0 tree data")
    #start_time = time.time()
    tree_file = f"{sim_info.directory}/merger_tree/MergerTree.snapshot_0%i.VELOCIraptor.tree" % initial_snap
    tree_data = TreeCatalogue(tree_file)

    halo_offset = tree_data.catalogue.ProgenOffsets.value[halo_index]
    num_progenitors = tree_data.catalogue.NumProgen.value[halo_index]
    progenitors_ID = tree_data.catalogue.Progenitors.value
    #print("--- %s seconds ---" % (time.time() - start_time))

    for snap in tqdm(range(initial_snap-1,final_snap,-1)):

        i = initial_snap - snap
        #print('snapshot', snap)
        snapshot_data = snapshot_info(sim_info, snap)
        path_to_catalogue_file = f"{snapshot_data.directory}/{snapshot_data.catalogue_name}"
        catalogue = velociraptor.load(path_to_catalogue_file)
        m200c = catalogue.masses.mass_200crit.to("Msun").value
        z[i] = catalogue.z

        #print("read z=0 tree data")
        #start_time = time.time()
        tree_file = f"{sim_info.directory}/merger_tree/MergerTree.snapshot_0%i.VELOCIraptor.tree" % snap
        tree_data = TreeCatalogue(tree_file)

        #print("--- %s seconds ---" % (time.time() - start_time))

        #print("Look for progenitor")
        #start_time = time.time()

        pro_indx, merger_ratio = look_for_progenitor_index(
            halo_offset,
            num_progenitors,
            progenitors_ID,
            tree_data.catalogue.ID.value,
            m200c,
        )

        halo_offset = tree_data.catalogue.ProgenOffsets.value[pro_indx.astype('int')]
        num_progenitors = tree_data.catalogue.NumProgen.value[pro_indx.astype('int')]
        progenitors_ID = tree_data.catalogue.Progenitors.value

        progenitor_index[:, i] = pro_indx
        merger_mass_ratio[:, i] = merger_ratio
        mass[:, i] = m200c[pro_indx.astype('int')]
        type[:, i] = catalogue.structure_type.structuretype[pro_indx.astype('int')]
        #print("--- %s seconds ---" % (time.time() - start_time))

    # Write data to file while it is being calculated..
    f.create_dataset('Mass', data=mass)
    f.create_dataset('Redshift', data=z)
    f.create_dataset('Structure_Type', data=type)
    f.create_dataset('Merger_mass_ratio', data=merger_mass_ratio)
    f.create_dataset('Progenitor_index', data=progenitor_index)
    data_file.close()

    # tree_data = {'progenitor_index': progenitor_index,
    #              'merger_mass_ratio': merger_mass_ratio,
    #              'M200crit': mass,
    #              'structure_type': type,
    #              'redshift': z}
    # return tree_data
    return progenitor_index
