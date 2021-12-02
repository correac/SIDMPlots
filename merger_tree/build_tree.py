import numpy as np
import h5py
from velociraptor import load

def build_tree(sim_info, halo_index):

    initial_snap = sim_info.n_snapshots - 1
    final_snap = 10

    progenitor_index = np.zeros(initial_snap - final_snap)
    progenitor_index[0] = halo_index
    merger_mass_ratio = np.zeros(initial_snap - final_snap)
    z = np.zeros(initial_snap - final_snap)
    mass = np.zeros(initial_snap - final_snap)
    type = np.zeros(initial_snap - final_snap)

    catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_base_name}" + "_%04i.properties" % initial_snap
    catalogue = load(catalogue_file)
    m200c = catalogue.masses.mass_200crit.to("Msun").value
    mass_descendant = m200c[halo_index]
    z[0] = catalogue.z
    mass[0] = m200c[halo_index]
    type[0] = catalogue.structure_type.structuretype[halo_index]

    tree_file = f"{sim_info.directory}/merger_tree/MergerTree.snapshot_0%i.VELOCIraptor.tree" % initial_snap
    with h5py.File(tree_file, "r") as hf:
        Progenitors = hf["Progenitors"][:]
        NumProgenitors = hf["NumProgen"][:]
        ProgenOffset = hf["ProgenOffsets"][:]

    halo = ProgenOffset[halo_index]
    num_progenitors = NumProgenitors[halo_index]
    progenitor_list = np.arange(num_progenitors)+halo
    proID = Progenitors[progenitor_list]

    for snap in range(initial_snap-1,final_snap,-1):

        catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_base_name}" + "_%04i.properties" % snap
        catalogue = load(catalogue_file)
        m200c = catalogue.masses.mass_200crit.to("Msun").value
        z[initial_snap-snap] = catalogue.z

        tree_file = f"{sim_info.directory}/merger_tree/MergerTree.snapshot_0%i.VELOCIraptor.tree" % snap
        with h5py.File(tree_file, "r") as hf:
            Progenitors = hf["Progenitors"][:]
            NumProgenitors = hf["NumProgen"][:]
            ProgenOffset = hf["ProgenOffsets"][:]
            ID = hf["ID"][:]

        if num_progenitors > 10: proID = proID[0:10]

        _, indx_ID, indx_proID = np.intersect1d(ID, proID, assume_unique=True, return_indices=True, )
        if len(indx_ID) == 0: break

        largest_mass_progenitor = np.where(m200c[indx_ID] == np.max(m200c[indx_ID]))[0]
        other_progenitors = np.where(m200c[indx_ID] != np.max(m200c[indx_ID]))[0]

        if len(largest_mass_progenitor) > 1:
            other_progenitors = largest_mass_progenitor[1:]
            largest_mass_progenitor = largest_mass_progenitor[0]

        connect = ProgenOffset[indx_ID[largest_mass_progenitor]]
        num_progenitors = NumProgenitors[indx_ID[largest_mass_progenitor]]
        progenitor_list = np.arange(num_progenitors) + connect
        proID = Progenitors[progenitor_list.astype('int')]
        progenitor_index[initial_snap - snap] = indx_ID[largest_mass_progenitor]

        if len(indx_ID) > 1:
            merger_mass_ratio[initial_snap - snap] = np.max(m200c[indx_ID[other_progenitors]] / mass_descendant)

        mass_descendant = m200c[indx_ID[largest_mass_progenitor]]
        mass[initial_snap-snap] = m200c[indx_ID[largest_mass_progenitor]]
        type[initial_snap-snap] = catalogue.structure_type.structuretype[indx_ID[largest_mass_progenitor]]
        #print(num_progenitors, proID, np.log10(mass[initial_snap-snap]), snap)

        if num_progenitors == 0: break

    #print(merger_mass_ratio)
    return progenitor_index, merger_mass_ratio, mass, type, z
