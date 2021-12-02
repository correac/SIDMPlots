import swiftsimio as sw
import numpy
import unyt


class load_particle_data:
    """
    Class containing particles properties
    """
    def __init__(
        self, sim_info, halo_index, index,
    ):
        """
        Parameters
        ----------
        """

        mask = sw.mask(f"{sim_info.directory}/{sim_info.snapshot_name}")

        # The full metadata object is available from within the mask
        size = unyt.unyt_array([0.5, 0.5, 0.5], 'Mpc')

        x = sim_info.halo_data.xminpot[index]
        y = sim_info.halo_data.yminpot[index]
        z = sim_info.halo_data.zminpot[index]
        origin = unyt.unyt_array([x, y, z], 'Mpc')

        # region is a 3x2 list [[left, right], [bottom, top], [front, back]]
        region = [[-0.5 * b + o, 0.5 * b + o] for b, o in zip(size, origin)]

        # Constrain the mask
        mask.constrain_spatial(region)

        # Now load the snapshot with this mask
        data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=mask)

        self.coordinates = data.dark_matter.coordinates.to("Mpc") - origin
        self.coordinates *= sim_info.a
        self.coordinates *= 1e3 #to kpc

        vx = sim_info.halo_data.vxminpot[index]
        vy = sim_info.halo_data.vyminpot[index]
        vz = sim_info.halo_data.vzminpot[index]
        origin = unyt.unyt_array([vx, vy, vz], "km/s")

        self.velocities = data.dark_matter.velocities.to("km/s") - origin
        self.masses = data.dark_matter.masses.to("Msun")

        if (hasattr(data.dark_matter,'cross_section')):
            unit_length_in_cgs = data.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
            unit_mass_in_cgs = data.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
            self.cross_section = data.dark_matter.cross_section.value
            self.cross_section *= unit_length_in_cgs ** 2 / unit_mass_in_cgs  # cm^2/g

            check = numpy.isinf(self.cross_section) == True
            self.cross_section[check] = 0
        else:
            self.cross_section = np.zeros(len(self.masses))