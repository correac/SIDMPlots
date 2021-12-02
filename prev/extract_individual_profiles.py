import numpy as np
import h5py
import scipy.stats as stat

def bin_volumes(radial_bins):
    """Returns the volumes of the bins. """

    single_vol = lambda x: (4.0 / 3.0) * np.pi * x ** 3
    outer = single_vol(radial_bins[1:])
    inner = single_vol(radial_bins[:-1])
    return outer - inner


def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)


def analyse_halo(mass, pos, vel):
    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    SumMasses, _, _ = stat.binned_statistic(x=r, values=np.ones(len(r)) * mass[0], statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3

    std_vel_x, _, _ = stat.binned_statistic(x=r, values=vel[:, 0], statistic="std", bins=radial_bins, )
    std_vel_y, _, _ = stat.binned_statistic(x=r, values=vel[:, 1], statistic="std", bins=radial_bins, )
    std_vel_z, _, _ = stat.binned_statistic(x=r, values=vel[:, 2], statistic="std", bins=radial_bins, )
    velocity = np.sqrt(std_vel_x ** 2 + std_vel_y ** 2 + std_vel_z ** 2) / np.sqrt(3.)
    velocity[np.where(np.isnan(velocity))[0]] = 0
    return density, velocity

def read_data(which_halos,snap,folder,output_path,name,mass_select):

    radial_bins = np.arange(0, 3, 0.1)
    radial_bins = 10**radial_bins
    centers = bin_centers(radial_bins) #kpc

    with h5py.File(folder+"/snapshot_00%02i.hdf5"%snap,"r") as hf:
        a = hf["/Header"].attrs["Scale-factor"]
        mass = hf['PartType1/Masses'][:] * 1e10 #Msun
        pos = hf['PartType1/Coordinates'][:][:] * a
        vel = hf['PartType1/Velocities'][:][:]
        # Read units
        unit_length_in_cgs = hf["/Units"].attrs["Unit length in cgs (U_L)"]
        unit_time_in_cgs = hf["/Units"].attrs["Unit time in cgs (U_t)"]
        vel *= unit_length_in_cgs / unit_time_in_cgs  # cm/s
        vel *= 1e-5  # km/s

    snapshot_file = h5py.File(folder+"/snapshot_00%02i.hdf5"%snap, "r")
    group_file = h5py.File(folder+"/halo_00%02i.catalog_groups"%snap, "r")
    particles_file = h5py.File(folder+"/halo_00%02i.catalog_particles"%snap, "r")
    properties_file = h5py.File(folder+"/halo_00%02i.properties"%snap, "r")

    c200c = properties_file["cNFW_200crit"][:]
    m200c = properties_file["Mass_200crit"][:] * 1e10
    m200c[m200c <= 0] = 1
    m200c = np.log10(m200c)
    CoP = np.zeros((len(m200c), 3))
    CoP[:, 0] = properties_file["Xcminpot"][:]
    CoP[:, 1] = properties_file["Ycminpot"][:]
    CoP[:, 2] = properties_file["Zcminpot"][:]
    subtype = properties_file["Structuretype"][:]

    min_mass = mass_select-0.2
    max_mass = mass_select+0.2

    select_halos = np.where((m200c >= min_mass) & (m200c <= max_mass))[0]  # >10 star parts

    # Checking sample
    if which_halos == 'subhalos':
        select = np.where(subtype[select_halos] > 10)[0]
        select_halos = select_halos[select]
    else:
        select = np.where(subtype[select_halos] == 10)[0]
        select_halos = select_halos[select]

    if len(select_halos) >= 5:
        select_random = np.random.random_integers(len(select_halos) - 1, size=(5))
        select_halos = select_halos[select_random]


    M200 = m200c[select_halos]
    c200 = c200c[select_halos]
    num_halos = len(select_halos)

    if which_halos == 'subhalos':
        file = output_path + "Profile_subhalos"
        file += "_" + name + ".txt"
        output = np.zeros((num_halos, 2))
        output[:, 0] = M200
        output[:, 1] = c200
        np.savetxt(file, output, fmt="%s")

    else:
        file = output_path + "Profile_halos"
        file += "_" + name + ".txt"
        output = np.zeros((num_halos, 2))
        output[:, 0] = M200
        output[:, 1] = c200
        np.savetxt(file, output, fmt="%s")

    for halo in range(0, num_halos):
        halo_j = select_halos[halo]

        # Grab the start position in the particles file to read from
        halo_start_position = group_file["Offset"][halo_j]
        halo_end_position = group_file["Offset"][halo_j + 1]
        particle_ids_in_halo = particles_file["Particle_IDs"][halo_start_position:halo_end_position]
        particle_ids_from_snapshot = snapshot_file["PartType1/ParticleIDs"][...]

        _, indices_v, indices_p = np.intersect1d(particle_ids_in_halo,
                                                 particle_ids_from_snapshot,
                                                 assume_unique=True,
                                                 return_indices=True, )

        particles_mass = mass[indices_p].copy()
        particles_pos = pos[indices_p, :].copy()
        particles_pos -= CoP[halo_j, :]  # centering
        particles_pos *= 1e3  # kpc
        particles_vel = vel[indices_p, :].copy()
        if len(particles_mass) == 0 :continue

        density_halo, velocity_halo = analyse_halo(particles_mass, particles_pos, particles_vel)

        output = np.zeros((len(centers), 3))
        output[:, 0] = centers
        output[:, 1] = density_halo
        output[:, 2] = velocity_halo

        if which_halos == 'subhalos':
            file = output_path + "Profile_subhalos"
            file += "_"+name+"_%02i"%halo +".txt"
            np.savetxt(file, output, fmt="%s")

        else:
            file = output_path + "Profile_halos"
            file += "_"+name+"_%02i"%halo +".txt"
            np.savetxt(file, output, fmt="%s")






if __name__ == '__main__':

    from utils import *

    output_path = args.output
    folder = args.directory
    snapshot = args.snapshot
    name = args.name

    mass = 9.5
    read_data("halos",snapshot,folder,output_path,name,mass)
    read_data("subhalos",snapshot,folder,output_path,name,mass)

    mass = 10
    read_data("halos",snapshot,folder,output_path,name,mass)
    read_data("subhalos",snapshot,folder,output_path,name,mass)

    #mass = 11
    #read_data("halos",snapshot,folder,output_path,name,mass)
    #read_data("subhalos",snapshot,folder,output_path,name,mass)
