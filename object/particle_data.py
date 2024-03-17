import swiftsimio as sw
import numpy
import unyt
import h5py
from typing import Tuple
import velociraptor
from astropy.cosmology import Planck13 as cosmo

class load_dark_matter_particle_data:
    """
    Class containing particles properties
    """

    def __init__(
            self, sim_info, CoP, CoV, size
    ):
        """
        Parameters
        ----------
        """
        mask = sw.mask(f"{sim_info.directory}/{sim_info.snapshot_name}")

        # region is a 3x2 list [[left, right], [bottom, top], [front, back]]
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, CoP)]

        # Constrain the mask
        mask.constrain_spatial(region)

        # Now load the snapshot with this mask
        data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=mask)

        self.ids = data.dark_matter.particle_ids.value

        self.coordinates = data.dark_matter.coordinates.to("Mpc") * sim_info.a
        self.coordinates -= CoP * sim_info.a  # centering
        self.coordinates *= 1e3  # to kpc

        self.velocities = data.dark_matter.velocities.to("km/s") - CoV
        self.masses = data.dark_matter.masses.to("Msun")

        if (hasattr(data.dark_matter, 'cross_section')):
            unit_length_in_cgs = data.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
            unit_mass_in_cgs = data.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
            self.cross_section = data.dark_matter.cross_section.value
            self.cross_section *= unit_length_in_cgs ** 2 / unit_mass_in_cgs  # cm^2/g

            check = numpy.isinf(self.cross_section) == True
            self.cross_section[check] = 0

            self.num_sidm_events = data.dark_matter.sidm_events.value

        else:
            self.cross_section = numpy.zeros(sim_info.num_dm_particles)
            self.num_sidm_events = numpy.zeros(sim_info.num_dm_particles)


class load_star_particle_data:
    """
    Class containing particles properties
    """

    def __init__(
            self, sim_info, CoP, CoV, size
    ):
        """
        Parameters
        ----------
        """
        mask = sw.mask(f"{sim_info.directory}/{sim_info.snapshot_name}")

        # region is a 3x2 list [[left, right], [bottom, top], [front, back]]
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, CoP)]

        # Constrain the mask
        mask.constrain_spatial(region)

        # Now load the snapshot with this mask
        data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=mask)

        self.ids = data.stars.particle_ids.value
        num_stars = len(self.ids)

        self.coordinates = data.stars.coordinates.to("Mpc") * sim_info.a
        self.coordinates -= CoP * sim_info.a  # centering
        self.coordinates *= 1e3  # to kpc

        self.velocities = data.stars.velocities.to("km/s") - CoV
        self.masses = data.stars.masses.to("Msun")

        self.luminosities_r_band = data.stars.luminosities.GAMA_r.value
        self.luminosities_u_band = data.stars.luminosities.GAMA_u.value
        self.luminosities_i_band = data.stars.luminosities.GAMA_i.value
        self.luminosities_g_band = data.stars.luminosities.GAMA_g.value
        self.luminosities_K_band = data.stars.luminosities.GAMA_K.value

        self.metal_mass_fractions = data.stars.metal_mass_fractions.value
        self.smoothed_metal_mass_fractions = data.stars.smoothed_metal_mass_fractions.value

        stars_birthz = (
                1.0 / data.stars.birth_scale_factors.value - 1.0
        )

        cosmic_age_z0 = cosmo.age(0.0).value
        cosmic_age = cosmo.age(stars_birthz).value
        self.age = numpy.ones(num_stars) * cosmic_age_z0 - cosmic_age



class load_gas_particle_data:
    """
    Class containing particles properties
    """

    def __init__(
            self, sim_info, CoP, CoV, size
    ):
        """
        Parameters
        ----------
        """
        mask = sw.mask(f"{sim_info.directory}/{sim_info.snapshot_name}")

        # region is a 3x2 list [[left, right], [bottom, top], [front, back]]
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, CoP)]

        # Constrain the mask
        mask.constrain_spatial(region)

        # Now load the snapshot with this mask
        data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=mask)

        self.ids = data.gas.particle_ids.value

        self.coordinates = data.gas.coordinates.to("Mpc") * sim_info.a
        self.coordinates -= CoP * sim_info.a  # centering
        self.coordinates *= 1e3  # to kpc

        self.velocities = data.gas.velocities.to("km/s") - CoV
        self.masses = data.gas.masses.to("Msun")



class load_particle_data:
    """
    Class containing particles properties
    """
    def __init__(
        self, sim_info, halo_index, index = None,
    ):
        """
        Parameters
        ----------
        """

        # The full metadata object is available from within the mask
        size = unyt.unyt_array([1, 1, 1], 'Mpc')

        if index == None:

            catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_name}"
            catalogue = velociraptor.load(catalogue_file)
            x = catalogue.positions.xcminpot.to("Mpc").value[halo_index.astype('int')]
            y = catalogue.positions.ycminpot.to("Mpc").value[halo_index.astype('int')]
            z = catalogue.positions.zcminpot.to("Mpc").value[halo_index.astype('int')]

        else :
            x = sim_info.halo_data.xminpot[index]
            y = sim_info.halo_data.yminpot[index]
            z = sim_info.halo_data.zminpot[index]
        
        CoP = unyt.unyt_array([x, y, z], 'Mpc') / sim_info.a #to comoving

        if index == None:

            catalogue_file = f"{sim_info.directory}/{sim_info.catalogue_name}"
            catalogue = velociraptor.load(catalogue_file)
            vx = catalogue.velocities.vxcminpot.to("km/s").value[halo_index.astype('int')]
            vy = catalogue.velocities.vycminpot.to("km/s").value[halo_index.astype('int')]
            vz = catalogue.velocities.vzcminpot.to("km/s").value[halo_index.astype('int')]

        else:
            vx = sim_info.halo_data.vxminpot[index]
            vy = sim_info.halo_data.vyminpot[index]
            vz = sim_info.halo_data.vzminpot[index]

        CoV = unyt.unyt_array([vx, vy, vz], "km/s")

        self.dark_matter = load_dark_matter_particle_data(sim_info, CoP, CoV, size)

        if sim_info.simulation_type == 'Hydro':

            self.stars = load_star_particle_data(sim_info, CoP, CoV, size)

            self.gas = load_gas_particle_data(sim_info, CoP, CoV, size)


    def select_bound_particles(self, sim_info, halo_index, ids):
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

        halo_start_position = group_file["Offset"][halo_index.astype('int')]
        try:
            halo_end_position = group_file["Offset"][halo_index.astype('int') + 1]
        except IndexError:
            return numpy.array([-1])

        particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]

        _, _, mask = numpy.intersect1d(
            particle_ids_in_halo,
            ids,
            assume_unique=True,
            return_indices=True,
        )

        # Ensure that there are no negative indices
        mask = mask[mask > 0]

        return mask
