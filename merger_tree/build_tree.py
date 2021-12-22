import numpy as np
import velociraptor
import h5py
from ..halo_data.density_profile import calculate_halo_data
from .snapshot_data import snapshot_info

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)

def build_tree(sim_info, halo_index):

    initial_snap = sim_info.n_snapshots - 1
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

    # Related to internal structure evolution
    radial_bins = np.arange(-1, 3, 0.1)
    radial_bins = 10**radial_bins
    density = np.zeros((len(radial_bins),initial_snap - final_snap))
    vel_radial_bins = np.arange(0.2, 25, 0.25)
    velocity = np.zeros((len(radial_bins),initial_snap - final_snap))

    catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_base_name}" + "_%04i.properties" % initial_snap
    catalogue = velociraptor.load(catalogue_file)
    m200c = catalogue.masses.mass_200crit.to("Msun").value
    mass_descendant = m200c[halo_index]
    z[0] = catalogue.z
    mass[0] = m200c[halo_index]
    type[0] = catalogue.structure_type.structuretype[halo_index]
    density[:,0], velocity[:,0] = calculate_halo_data(sim_info, halo_index, radial_bins, vel_radial_bins)

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

        snapshot_data = snapshot_info(sim_info, snap)
        path_to_catalogue_file = f"{snapshot_data.directory}/{snapshot_data.catalogue_name}"
        catalogue = velociraptor.load(path_to_catalogue_file)
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
        density[:, initial_snap-snap], velocity[:, initial_snap-snap] = \
            calculate_halo_data(snapshot_data, indx_ID[largest_mass_progenitor], radial_bins, vel_radial_bins)

        if num_progenitors == 0: break

    density_radial_bins = bin_centers(radial_bins)  # kpc
    velocity_radial_bins = bin_centers(vel_radial_bins)  # kpc

    tree_data = {'progenitor_index': progenitor_index,
                 'merger_mass_ratio': merger_mass_ratio,
                 'M200crit': mass,
                 'structure_type': type,
                 'redshift': z,
                 'density': density,
                 'velocity': velocity,
                 'density_radial_bins': density_radial_bins,
                 'velocity_radial_bins': velocity_radial_bins}

    return tree_data
