import numpy as np
import unyt
from velociraptor import load


class HaloCatalogue:
    """
    General class containing halo properties
    """

    def __init__(
        self, path_to_catalogue: str, dm_particle_mass: float,
    ):
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
        catalogue = load(self.path_to_catalogue)

        # Selecting haloes that contain at less 1000 DM particles
        mask = catalogue.masses.mass_200crit.to("Msun").value >= unyt.unyt_quantity(1e3 * dm_particle_mass, "Msun")

        # Compute the number of haloes following the selection mask
        self.number_of_haloes = mask.sum()

        # Structure type
        self.structure_type = catalogue.structure_type.structuretype[mask]

        # Log10 halo mass in units of Msun
        self.log10_halo_mass = np.log10(
            catalogue.masses.mass_200crit.to("Msun").value[mask]
        )

        self.concentration = catalogue.concentration.cnfw.value[mask]
        self.virial_radius = catalogue.radii.r_200crit.to("Mpc").value[mask]

        self.scale_radius = self.virial_radius / self.concentration

        # Ids of haloes satisfying the selection criterion
        self.halo_index = np.array([i for i in range(len(mask)) if mask[i] == True])

        self.xminpot = catalogue.positions.xcminpot.to("Mpc").value[mask]
        self.yminpot = catalogue.positions.ycminpot.to("Mpc").value[mask]
        self.zminpot = catalogue.positions.zcminpot.to("Mpc").value[mask]

        self.vxminpot = catalogue.velocities.vxcminpot.to("km/s").value[mask]
        self.vyminpot = catalogue.velocities.vycminpot.to("km/s").value[mask]
        self.vzminpot = catalogue.velocities.vzcminpot.to("km/s").value[mask]

        self.vmax = catalogue.velocities.vmax.to("km/s").value[mask]