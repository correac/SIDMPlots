import numpy as np
import h5py
from .density_profile import calculate_halo_data, calculate_halo_data_hydro

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)

def calculate_relaxation(sim_info, sample):

    MoP = np.array([sim_info.halo_data.xminpot[sample],
                    sim_info.halo_data.yminpot[sample],
                    sim_info.halo_data.zminpot[sample]])
    CoM = np.array([sim_info.halo_data.xcom[sample],
                    sim_info.halo_data.ycom[sample],
                    sim_info.halo_data.zcom[sample]])
    pos = MoP - CoM
    r = np.sqrt( pos[0,:]**2 + pos[1,:]**2 + pos[2,:]**2 )

    relaxation = r / (0.07 * sim_info.halo_data.virial_radius[sample])

    return relaxation

def load_profiles(sim_info, halo_index, output_file):

    # Related to internal structure evolution
    radial_bins = np.arange(-1, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centered_radial_bins = bin_centers(radial_bins)  # kpc

    vel_radial_bins = np.arange(0.25, 100, 0.25)
    centered_velocity_radial_bins = bin_centers(vel_radial_bins)  # kpc

    density, velocity, veldisp, sigmaprofile = \
        calculate_halo_data(sim_info, halo_index, radial_bins, vel_radial_bins)

    data_file = h5py.File(output_file, 'a')
    f = data_file.create_group('Profile_data')
    f.create_dataset('Dark_matter_Density_profile', data=density)
    f.create_dataset('Dark_matter_Circular_Velocity', data=velocity)
    f.create_dataset('Dark_matter_Velocity_dispersion', data=veldisp)
    f.create_dataset('Dark_matter_Sigma_profile', data=sigmaprofile)
    f.create_dataset('Density_radial_bins', data=centered_radial_bins)
    f.create_dataset('Velocity_radial_bins', data=centered_velocity_radial_bins)

    if sim_info.simulation_type == 'Hydro':

        hydro_data = calculate_halo_data_hydro(sim_info, halo_index, radial_bins, vel_radial_bins)
        f.create_dataset('Stars_Density_profile', data=hydro_data['stars_density'])
        f.create_dataset('Stars_Circular_Velocity', data=hydro_data['stars_velocity'])
        f.create_dataset('Stars_Velocity_dispersion', data=hydro_data['stars_veldisp'])

        f.create_dataset('Gas_Density_profile', data=hydro_data['gas_density'])
        f.create_dataset('Gas_Circular_Velocity', data=hydro_data['gas_velocity'])
        f.create_dataset('Gas_Velocity_dispersion', data=hydro_data['gas_veldisp'])

        f.create_dataset('Density_profile', data=hydro_data['density'])
        f.create_dataset('Circular_Velocity', data=hydro_data['velocity'])

    data_file.close()
    return

def make_halo_data(sim_info):

    sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= 9) &
        (sim_info.halo_data.log10_halo_mass < 12))[0]

    halo_index = sim_info.halo_data.halo_index[sample]
    M200c = sim_info.halo_data.log10_halo_mass[sample]
    c200c = sim_info.halo_data.concentration[sample]
    R200c = sim_info.halo_data.virial_radius[sample]
    structure_type = sim_info.halo_data.structure_type[sample]
    relaxation = calculate_relaxation(sim_info, sample)
    Vmax = sim_info.halo_data.vmax[sample]

    # Output data
    output_file = f"{sim_info.output_path}/Halo_data_" + sim_info.simulation_name + ".hdf5"
    data_file = h5py.File(output_file, 'w')
    f = data_file.create_group('Halo_data')
    f.create_dataset('ID', data=halo_index)
    f.create_dataset('StructureType', data=structure_type)
    f.create_dataset('M200c', data=M200c)
    f.create_dataset('c200c', data=c200c)
    f.create_dataset('R200c', data=R200c)
    f.create_dataset('Vmax', data=Vmax)
    f.create_dataset('Dynamical_relaxation', data=relaxation)

    if sim_info.simulation_type == 'Hydro':

        Mstar = sim_info.halo_data.log10_stellar_mass[sample]
        Mgas = sim_info.halo_data.log10_gas_mass[sample]
        GalaxySize = sim_info.halo_data.galaxy_size[sample]
        SFR = sim_info.halo_data.sfr[sample]
        Metallicity = sim_info.halo_data.metallicity_stars[sample]

        f.create_dataset('Mstar', data=Mstar)
        f.create_dataset('Mgas', data=Mgas)
        f.create_dataset('GalaxySize', data=GalaxySize)
        f.create_dataset('SFR', data=SFR)
        f.create_dataset('Metallicity', data=Metallicity)


    data_file.close()

    load_profiles(sim_info, halo_index, output_file)
