"""
Description here
"""

import os
import h5py
import glob
#from scatter_rate import plot_cosmic_scatter_rate
from HMF import make_HMF, HMF_output_files
from plotter import output_cM_vMax_relations, \
    plot_relations, output_halo_profiles, plot_halo_profiles

class SimInfo:
    def __init__(self,folder,snap,output_path,name):
        self.name = name
        self.folder = folder
        self.output_path = output_path

        properties = os.path.join(folder, "subhalo_%04i.properties" % snap)
        if os.path.exists(properties):
            self.halo_properties = os.path.join(folder, "subhalo_%04i.properties" % snap)
        else :
            self.halo_properties = os.path.join(folder, "halo_%04i.properties" % snap)

        catalog = os.path.join(folder,"subhalo_%04i.catalog_groups" % snap)
        if os.path.exists(catalog):
            self.catalog_groups = os.path.join(folder,"subhalo_%04i.catalog_groups"%snap)
        else :
            self.catalog_groups = os.path.join(folder, "halo_%04i.catalog_groups" % snap)

        particles = os.path.join(folder, "subhalo_%04i.catalog_particles" % snap)
        if os.path.exists(particles):
            self.catalog_particles = os.path.join(folder, "subhalo_%04i.catalog_particles" % snap)
        else :
            self.catalog_particles = os.path.join(folder, "halo_%04i.catalog_particles" % snap)

        self.latest_snapshot = os.path.join(self.folder,"snapshot_%04i.hdf5"%snap)
        self.n_snapshots = int(snap)

        snapshot_file = h5py.File(self.latest_snapshot, "r")
        self.boxSize = snapshot_file["/Header"].attrs["BoxSize"][0] # Mpc
        self.a = snapshot_file["/Cosmology"].attrs["Scale-factor"]
        self.z = 1. / self.a - 1.
        self.h = snapshot_file["/Cosmology"].attrs["h"]
        self.Om = snapshot_file["/Cosmology"].attrs["Omega_m"]
        self.Ol = snapshot_file["/Cosmology"].attrs["Omega_lambda"]
        self.rhocrit0 = 2.7754e11 * self.h**2 # Msun / Mpc^3
        self.sigma = snapshot_file["/SIDMScheme"].attrs["SIDM cross section [cgs units]"]
        self.dmpart_mass =  snapshot_file["PartType1/Masses"][0] * 1e10 # Msun

    def snapshot(self,snap):
        return os.path.join(self.folder,"snapshot_%04i.hdf5"%snap)


def sidmplots(siminfo):

    #plot_cosmic_scatter_rate(siminfo,output_path)
    #HMF_output_files(siminfo)
    #output_cM_vMax_relations(siminfo)
    output_halo_profiles(siminfo)
    #plot_individual_profiles(siminfo)


if __name__ == '__main__':
    from utils import *

    # Load SIDMplots production details
    output_path = args.output
    number_of_inputs = len(args.snapshot)
    directory_list = args.directory
    snapshot_list = args.snapshot

    name_list = (
        args.run_names
        if args.run_names is not None
        else [None] * number_of_inputs
    )

    # Loop over simulation list
    for sims in range(number_of_inputs):
        directory = directory_list[sims]
        snap_number = int(snapshot_list[sims])
        sim_name = name_list[sims]
        siminfo = SimInfo(directory, snap_number, output_path, sim_name)

        # Run SIDMplots
        sidmplots(siminfo)

    #make_HMF(siminfo, name_list)
    #plot_relations(siminfo, name_list)
    plot_halo_profiles(siminfo, name_list)

