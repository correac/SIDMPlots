import numpy as np
from tqdm import tqdm
from object import particle_data
import velociraptor
import unyt
import h5py


class Catalogue:
    """
    General class containing halo properties
    """

    def __init__(self, path_to_catalogue: str):
        """
        Parameters
        ----------
        path_to_catalogue: str
        Path to the catalogue with halo properties
        dm_particle_mass: unyt.array.unyt_quantity
        Minimum dark matter particle mass in units of Msun. Haloes that contain less than
        1000 dark mattter particles are disregarded
        """

        self.path_to_catalogue = path_to_catalogue

        # Load catalogue using velociraptor python library
        catalogue = velociraptor.load(self.path_to_catalogue)

        # Selecting haloes that contain at less 10 DM particles
        dm_particle_mass = 9.7e6 #Msun

        mask = np.where(
            catalogue.masses.mass_200crit.to("Msun").value >= unyt.unyt_quantity(10 * dm_particle_mass, "Msun")
        )[0]

        # Compute the number of haloes following the selection mask
        self.number_of_haloes = len(mask)

        # Structure type
        self.structure_type = catalogue.structure_type.structuretype[mask]

        # Log10 halo mass in units of Msun
        self.log10_halo_mass = np.log10(
            catalogue.masses.mass_200crit.to("Msun").value[mask]
        )

        # Ids of haloes satisfying the selection criterion
        self.halo_index = mask.copy()

        self.xminpot = catalogue.positions.xcminpot.to("Mpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("Mpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("Mpc").value[mask]

        self.id_mbp = catalogue.ids.id_mbp.value[mask]



def select_haloes(sim_info, cdm_ids):
    """
    Select particles that are gravitationally bound to halo
    Parameters
    ----------
    halo_id: int
    Halo id from the catalogue
    Returns
    -------
    Output: Tuple[np.ndarray, np.ndarray]
    A tuple containing ids of the stellar particles and gas particles
    """
    particles_file = h5py.File(f"{sim_info.directory}/{sim_info.catalogue_particles}", "r")
    group_file = h5py.File(f"{sim_info.directory}/{sim_info.catalogue_groups}", "r")

    halo_position = group_file["Offset"][:]
    particle_ids_sim = particles_file["Particle_IDs"][:]

    match_sim = np.zeros(len(cdm_ids))
    for i in range(len(cdm_ids)):
        select = np.where(particle_ids_sim == cdm_ids[i])[0]
        if len(select) == 1:
            match_sim[i] = select

    # Remove zero case
    match_sim = match_sim[match_sim>0]
    match_sim = match_sim.astype('int')

    min_match = np.min(match_sim)
    max_match = np.max(match_sim)
    select_min_range = np.where(min_match >= halo_position)[0]
    select_max_range = np.where(max_match <= halo_position)[0]

    if select_max_range[0]-1 > select_min_range[-1]:
        halo = np.arange(select_min_range[-1],select_max_range[0]-1,1)
    else:
        halo = np.array([select_min_range[-1]])

    return halo

def look_for_particle_ids(sim_info, halo_index):

    part_data = particle_data.load_particle_data(sim_info, halo_index)
    dm_bound_particles_only = part_data.select_bound_particles(sim_info, halo_index, part_data.dark_matter.ids)

    pos = part_data.dark_matter.coordinates.value[dm_bound_particles_only, :]

    r = np.sqrt(np.sum(pos ** 2, axis=1))
    within_2_kpc = r <= 2.
    dm_part_ids = part_data.dark_matter.ids[dm_bound_particles_only]
    dm_part_ids = dm_part_ids[within_2_kpc]

    if len(dm_part_ids) > 10:
        dm_part_ids = dm_part_ids[0:11]  # Select random 10 dm parts within 2kpc region

    return dm_part_ids

def look_for_min_distance_and_mass(cdm_info, sample, sidm_info, sidm_haloes):

    cdm_halo_mass = cdm_info.halo_data.log10_halo_mass[sample]
    sidm_haloes_all = sidm_info.halo_data.halo_index

    _, sample_sidm, _ = np.intersect1d(
        sidm_haloes_all, sidm_haloes,
        assume_unique=True,
        return_indices=True,
    )

    sidm_halo_mass = sidm_info.halo_data.log10_halo_mass[sample_sidm]
    sidm_mass_diff = sidm_halo_mass - cdm_halo_mass

    select = np.where(np.abs(sidm_mass_diff) == np.min(np.abs(sidm_mass_diff)))[0]

    if len(select) == 1:
        match = sidm_haloes[select]

    else:
        cdm_center = np.array([
            cdm_info.halo_data.xminpot[sample],
            cdm_info.halo_data.yminpot[sample],
            cdm_info.halo_data.zminpot[sample],
        ])

        sidm_center = np.zeros((3,len(select)))
        sidm_center[0,:] = sidm_info.halo_data.xminpot[sample_sidm[select]]
        sidm_center[1,:] = sidm_info.halo_data.yminpot[sample_sidm[select]]
        sidm_center[2,:] = sidm_info.halo_data.zminpot[sample_sidm[select]]

        # sidm_center = np.array([
        #     sidm_info.halo_data.xminpot[sample_sidm[select]],
        #     sidm_info.halo_data.yminpot[sample_sidm[select]],
        #     sidm_info.halo_data.zminpot[sample_sidm[select]],
        # ])

        r = sidm_center - cdm_center
        r = np.sqrt(np.sum(r ** 2, axis=1))
        select_pos = np.where(r == np.min(r))[0]
        match = sidm_haloes[select[select_pos[0]]]

    return match # returning SIDM halo index


def match_simulations(cdm_info, sidm_info):
    """
    :param cdm_info: class for CDM simulation
    :param sidm_info: class data for SIDM simulation
    :return: output file
    """

    #sample = np.where(cdm_info.halo_data.log10_halo_mass >= 10)[0]
    sample = np.where(cdm_info.halo_data.log10_stellar_mass >= 8.5)[0]

    halo_index = cdm_info.halo_data.halo_index[sample]
    num_haloes = len(sample)

    matched_halo_sidm = np.zeros(num_haloes)

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        cdm_part_ids = look_for_particle_ids(cdm_info, halo_index[i])

        sidm_haloes = select_haloes(sidm_info, cdm_part_ids)

        if len(sidm_haloes) == 1:
            matched_halo_sidm[i] = sidm_haloes

        else:
            matched_halo_sidm[i] = look_for_min_distance_and_mass(
                cdm_info, sample[i], sidm_info, sidm_haloes)

        # print('sidm haloes', matched_halo_sidm[i], halo_index[i])

    # Output data
    output_file = f"{cdm_info.output_path}/Halo_match_" + sidm_info.simulation_name + ".hdf5"
    data_file = h5py.File(output_file, 'w')
    f = data_file.create_group('Halo_matched')
    f.create_dataset('CDM_HaloIndex', data=halo_index)
    f.create_dataset('SIDM_HaloIndex', data=matched_halo_sidm)
    data_file.close()


