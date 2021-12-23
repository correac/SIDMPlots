from typing import List, Union, Tuple, Dict
import unyt
import glob
import os
import swiftsimio


class snapshot_info:

    def __init__(self, sim_info, snap):

        self.directory = sim_info.directory
        self.catalogue_name = sim_info.catalogue_base_name + "_%04i.properties" % snap
        self.snapshot_name = sim_info.snapshot_base_name + "_%04i.hdf5" % snap
        self.output_path = sim_info.output_path

        # Find the group and particle catalogue files
        #self.__find_groups_and_particles_catalogues()

        # Load snapshot via swiftsimio
        #self.snapshot = swiftsimio.load(f"{self.directory}/{self.snapshot_name}")

        # Load conversion units
        #self.__load_conversion()

        # Load cosmological params
        #self.__load_cosmology()
        return

    def __load_conversion(self) -> None:

        # Conversion from internal units to kpc
        kpc = 3.08567758e21
        Msun = 1.9891e33
        yr = 3.1556926e7
        Myr = yr * 1e6

        self.to_kpc_units = (
                self.snapshot.metadata.internal_code_units["Unit length in cgs (U_L)"][0]
                / kpc
        )

        # Conversion from internal units to Msun
        self.to_Msun_units = (
                self.snapshot.metadata.internal_code_units["Unit mass in cgs (U_M)"][0]
                / Msun
        )

        # Conversion from internal units to Myr
        self.to_Myr_units = (
                self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
                / Myr
        )

        # Conversion from internal units to yr
        self.to_yr_units = (
                self.snapshot.metadata.internal_code_units["Unit time in cgs (U_t)"][0]
                / yr
        )

    def __load_cosmology(self) -> None:

        # Box size of the simulation in kpc
        self.boxSize = self.snapshot.metadata.boxsize.to("kpc").value[0]

        # Cosmic scale factor
        self.a = self.snapshot.metadata.scale_factor

        self.hubble_time_Gyr = self.snapshot.metadata.time.to("Gyr").value

        self.Omega_m = self.snapshot.metadata.cosmology.Om0

        self.Omega_l = 1 - self.Omega_m

        self.h = self.snapshot.metadata.cosmology.H0.value / 100.

        self.z = self.snapshot.metadata.redshift

        self.rhocrit0 = 2.7754e11 * self.h ** 2  # Msun / Mpc^3

        # Maximum softening
        self.softening = (
                self.snapshot.metadata.gravity_scheme[
                    "Maximal physical DM softening length (Plummer equivalent) [internal units]"
                ][0] * self.to_kpc_units
        )

        # Maximum softening for baryons
        self.baryon_max_soft = (
                self.snapshot.metadata.gravity_scheme[
                    "Maximal physical baryon softening length  [internal units]"
                ][0]
                * self.to_kpc_units
        )

        self.dm_particle_mass = 1e6 # Need to fix, table mass is empty

        self.num_dm_particles = self.snapshot.metadata.n_dark_matter


    def __find_groups_and_particles_catalogues(self) -> None:
        """
        Finds paths to the fields with particles catalogue and groups catalogue
        """

        catalogue_num = "".join([s for s in self.catalogue_name if s.isdigit()])
        catalogue_groups_paths: List[str] = glob.glob(
            f"{self.directory}/*{catalogue_num}.catalog_groups*"
        )
        catalogue_particles_paths: List[str] = glob.glob(
            f"{self.directory}/*{catalogue_num}.catalog_particles*"
        )

        # We expect one file for particle groups
        if len(catalogue_groups_paths) == 1:
            self.catalogue_groups = catalogue_groups_paths[0].split("/")[-1]
        else:
            raise IOError("Couldn't find catalogue_groups file")

        # We expect two files: one for bound and the other for unbound particles
        if len(catalogue_particles_paths) == 2:
            for path in catalogue_particles_paths:
                if path.find("unbound") == -1:
                    self.catalogue_particles = path.split("/")[-1]
        else:
            raise IOError("Couldn't find catalogue_particles file")

        return
