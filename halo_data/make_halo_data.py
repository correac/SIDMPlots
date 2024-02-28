import numpy as np
import h5py
from .density_profile import calculate_halo_data, calculate_halo_data_hydro
from .morphology import calculate_morphology
from .galaxy_image import make_galaxy_images
from object import particle_data
import swiftsimio as sw
from astropy.cosmology import Planck13

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
    radial_bins = np.arange(-1, 3.1, 0.1)
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

    if sim_info.simulation_type == 'Hydro':
        sample = np.where((sim_info.halo_data.log10_halo_mass >= 10) & (sim_info.halo_data.log10_stellar_mass >= 8.5))[0]
    else:
        sample = np.where(sim_info.halo_data.log10_halo_mass >= 10)[0]
    
    centrals = np.where(sim_info.halo_data.structure_type[sample] == 10)[0]
    sample = sample[centrals]

    halo_index = sim_info.halo_data.halo_index[sample]
    M200c = sim_info.halo_data.log10_halo_mass[sample]
    c200c = sim_info.halo_data.concentration[sample]
    R200c = sim_info.halo_data.virial_radius[sample]
    structure_type = sim_info.halo_data.structure_type[sample]
    relaxation = calculate_relaxation(sim_info, sample)
    Vmax = sim_info.halo_data.vmax[sample]

    # # Morphology/shape estimations
    data = calculate_morphology(sim_info, sample)

    kappa = data['kappa']
    smomentum = data['smomentum']
    Lmomentum = data['Lmomentum']
    GalaxyHalfLightRadius = data['GalaxyHalfLightRadius']
    GalaxyProjectedHalfLightRadius = data['GalaxyProjectedHalfLightRadius']
    GalaxyLuminosity = data['GalaxyLuminosity']


    # make_galaxy_images(sim_info, halo_index, sample, smomentum, kappa)
    # make_galaxy_images(sim_info, halo_index[0:25], sample[0:25], smomentum[:,0:25], kappa[0:25])


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
    f.create_dataset('DynamicalRelaxation', data=relaxation)
    f.create_dataset('CrossSection', data=data['cross_section'])
    f.create_dataset('AxisRadius', data=data['radial_bins'])
    f.create_dataset('DMMajorAxis_a', data=data['DM_a_axis'])
    f.create_dataset('DMMinorAxis_b', data=data['DM_b_axis'])
    f.create_dataset('DMMinorAxis_c', data=data['DM_c_axis'])
    f.create_dataset('DMNparticlesWithinAxisRadius', data=data['DMNparts'])

    if sim_info.simulation_type == 'Hydro':

        # f.create_dataset('StarsMajorAxis_a', data=data['Stars_a_axis'])
        # f.create_dataset('StarsMinorAxis_b', data=data['Stars_b_axis'])
        # f.create_dataset('StarsMinorAxis_c', data=data['Stars_c_axis'])
        # f.create_dataset('StarsNparticlesWithinAxisRadius', data=data['StarsNparts'])
        #
        # f.create_dataset('GasMajorAxis_a', data=data['Gas_a_axis'])
        # f.create_dataset('GasMinorAxis_b', data=data['Gas_b_axis'])
        # f.create_dataset('GasMinorAxis_c', data=data['Gas_c_axis'])
        # f.create_dataset('GasNparticlesWithinAxisRadius', data=data['GasNparts'])

        Mstar = sim_info.halo_data.log10_stellar_mass[sample]
        Mgas = sim_info.halo_data.log10_gas_mass[sample]
        GalaxyHalfMassRadius = sim_info.halo_data.half_mass_radius_star[sample]
        GalaxyProjectedHalfMassRadius = sim_info.halo_data.half_mass_projected_radius_star[sample]
        SFR = sim_info.halo_data.sfr[sample]
        Metallicity = sim_info.halo_data.metallicity_stars[sample]

        f.create_dataset('Mstar', data=Mstar)
        f.create_dataset('Mgas', data=Mgas)
        f.create_dataset('GalaxyHalfMassRadius', data=GalaxyHalfMassRadius)
        f.create_dataset('GalaxyProjectedHalfMassRadius', data=GalaxyProjectedHalfMassRadius)
        f.create_dataset('GalaxyHalfLightRadius', data=GalaxyHalfLightRadius)
        f.create_dataset('GalaxyProjectedHalfLightRadius', data=GalaxyProjectedHalfLightRadius)
        f.create_dataset('GalaxyLuminosity', data=GalaxyLuminosity)
        f.create_dataset('SFR', data=SFR)
        f.create_dataset('Metallicity', data=Metallicity)
        f.create_dataset('kappa', data=kappa)
        f.create_dataset('SpecificAngularMomentum', data=Lmomentum)

    data_file.close()

    load_profiles(sim_info, halo_index, output_file)


def load_particle_data(sim_info, output_file):

    data = sw.load(f"{sim_info.directory}/{sim_info.snapshot_name}")

    data_file = h5py.File(output_file, 'a')

    if sim_info.simulation_type == 'Hydro':

        star_coordinates = data.stars.coordinates.to("Mpc").value * sim_info.a * 1e3 # to kpc
        star_masses = data.stars.masses.to("Msun").value / 1e10 # in units / 1e10 Msun
        star_metallicity = data.stars.metal_mass_fractions.value / 0.0134 # in units / Zun
        star_age_scale_factor = data.stars.birth_scale_factors.value
        star_age_redshift = 1. / star_age_scale_factor - 1.0
        tH = Planck13.age(0).value
        star_age_cosmic_time = Planck13.age(star_age_redshift).value
        star_age_lookback_time = tH - star_age_cosmic_time

        f = data_file.create_group('StarParticleData')
        data_part = f.create_dataset('Coordinates', data=star_coordinates)
        data_part.attrs["Description"] = np.string_("Particle coordinates in units of [kpc]")

        data_part = f.create_dataset('Masses', data=star_masses)
        data_part.attrs["Description"] = np.string_("Particle masses in units of [1e10 Msun]")

        data_part = f.create_dataset('Metallicity', data=star_metallicity)
        data_part.attrs["Description"] = np.string_("Particle metallicity in units of [Zsun], with Zsun=0.0134 (Asplund et al. 2009)")

        data_part = f.create_dataset('Age', data=star_age_lookback_time)
        data_part.attrs["Description"] = np.string_("Particle age indicated as lookback time in units of [Gyr]")

    DM_coordinates = data.dark_matter.coordinates.to("Mpc").value * sim_info.a * 1e3 # to kpc
    DM_masses = data.dark_matter.masses.to("Msun").value / 1e10 # in units / 1e10 Msun

    f = data_file.create_group('DMParticleData')
    data_part = f.create_dataset('Coordinates', data=DM_coordinates)
    data_part.attrs["Description"] = np.string_("Particle coordinates in units of [kpc]")

    data_part = f.create_dataset('Masses', data=DM_masses)
    data_part.attrs["Description"] = np.string_("Particle masses in units of [1e10 Msun]")

    data_file.close()
    return


def make_processed_data(sim_info):

    if sim_info.simulation_type == 'Hydro':
        sample = np.where((sim_info.halo_data.log10_halo_mass >= 10) & (sim_info.halo_data.log10_stellar_mass >= 6))[0]
    else:
        sample = np.where(sim_info.halo_data.log10_halo_mass >= 10)[0]

    centrals = np.where(sim_info.halo_data.structure_type[sample] == 10)[0]
    sample = sample[centrals]

    halo_index = sim_info.halo_data.halo_index[sample]
    M200c = sim_info.halo_data.log10_halo_mass[sample]
    c200c = sim_info.halo_data.concentration[sample]
    R200c = sim_info.halo_data.virial_radius[sample]
    Vmax = sim_info.halo_data.vmax[sample]
    xminpot = sim_info.halo_data.xminpot[sample]
    yminpot = sim_info.halo_data.yminpot[sample]
    zminpot = sim_info.halo_data.zminpot[sample]
    xcom = sim_info.halo_data.xcom[sample]
    ycom = sim_info.halo_data.ycom[sample]
    zcom = sim_info.halo_data.zcom[sample]

    # Morphology/shape estimations
    # data = calculate_morphology(sim_info, sample)

    # kappa = data['kappa']
    # Lmomentum = data['Lmomentum']
    # GalaxyHalfLightRadius = data['GalaxyHalfLightRadius']
    # GalaxyProjectedHalfLightRadius = data['GalaxyProjectedHalfLightRadius']
    # GalaxyLuminosity = data['GalaxyLuminosity']

    # Output data
    output_file = f"{sim_info.output_path}/TangoSIDM_" + sim_info.simulation_name + ".hdf5"
    data_file = h5py.File(output_file, 'w')

    description = "TangoSIDM is a simulation project with the goal of exploring the impact of Self-Interacting Dark Matter "
    description += "(SIDM) on galaxy formation. This file is a dataset of the TangoSIDM simulations and contains the Halo "
    description += "Catalogue, Star particle data and Dark Matter particle data of the simulations at redshift 0."
    data_file.create_dataset('Description', data=np.string_(description))

    contact = "Dataset generated by Camila Correa (CEA Paris-Saclay). Email: camila.correa@cea.fr,"
    contact += " website: www.camilacorrea.com. TangoSIDM: www.tangosidm.com."
    data_file.create_dataset('Contact', data=np.string_(contact))

    Reference = np.string_(
        ['Correa, C. A., et al., (2022), MNRAS, Vol 517, Issue 2, pp. 3045. The publication of this '
         'specific dataset is currently in preparation.'])
    data_file.create_dataset('Reference', data=Reference)

    f = data_file.create_group('HaloCatalogue')
    f.create_dataset('HaloID', data=halo_index)

    data_halo = f.create_dataset('M200c', data=M200c)
    data_halo.attrs["Description"] = np.string_("Halo mass in units of [log10 Msun]")

    f.create_dataset('c200c', data=c200c)

    data_halo = f.create_dataset('R200c', data=R200c)
    data_halo.attrs["Description"] = np.string_("Halo virial radius in units of [Mpc]")

    f.create_dataset('Vmax', data=Vmax)

    data_halo = f.create_dataset('XminimumPotential', data=xminpot)
    data_halo.attrs["Description"] = np.string_("Minimum of the potential, x coordinate in units of [Mpc]")

    data_halo = f.create_dataset('YminimumPotential', data=yminpot)
    data_halo.attrs["Description"] = np.string_("Minimum of the potential, y coordinate in units of [Mpc]")

    data_halo = f.create_dataset('ZminimumPotential', data=zminpot)
    data_halo.attrs["Description"] = np.string_("Minimum of the potential, z coordinate in units of [Mpc]")

    data_halo = f.create_dataset('XCenterOfMass', data=xcom)
    data_halo.attrs["Description"] = np.string_("Center of mass, x coordinate in units of [Mpc]")

    data_halo = f.create_dataset('YCenterOfMass', data=ycom)
    data_halo.attrs["Description"] = np.string_("Center of mass, y coordinate in units of [Mpc]")

    data_halo = f.create_dataset('ZCenterOfMass', data=zcom)
    data_halo.attrs["Description"] = np.string_("Center of mass, z coordinate in units of [Mpc]")

    if sim_info.simulation_type == 'Hydro':

        Mstar = sim_info.halo_data.log10_stellar_mass[sample]
        Mgas = sim_info.halo_data.log10_gas_mass[sample]
        GalaxyHalfMassRadius = sim_info.halo_data.half_mass_radius_star[sample]
        GalaxyProjectedHalfMassRadius = sim_info.halo_data.half_mass_projected_radius_star[sample]
        SFR = sim_info.halo_data.sfr[sample]
        Metallicity = sim_info.halo_data.metallicity_stars[sample]

        data_halo = f.create_dataset('GalaxyStellarMass', data=Mstar)
        data_halo.attrs["Description"] = np.string_("Galaxy stellar mass under an aperture of 50 kpc. Units [log10 Msun]")

        data_halo = f.create_dataset('GalaxyGasMass', data=Mgas)
        data_halo.attrs["Description"] = np.string_("Galaxy gas mass under an aperture of 50 kpc. Units [log10 Msun]")

        data_halo = f.create_dataset('GalaxyHalfMassRadius', data=GalaxyHalfMassRadius)
        data_halo.attrs["Description"] = np.string_("Galaxy half mass radius in units [kpc]")

        data_halo = f.create_dataset('GalaxyProjectedHalfMassRadius', data=GalaxyProjectedHalfMassRadius)
        data_halo.attrs["Description"] = np.string_("Galaxy projected half mass radius in units [kpc]")

        # f.create_dataset('GalaxyHalfLightRadius', data=GalaxyHalfLightRadius)
        # f.create_dataset('GalaxyProjectedHalfLightRadius', data=GalaxyProjectedHalfLightRadius)
        # f.create_dataset('GalaxyLuminosity', data=GalaxyLuminosity)
        f.create_dataset('GalaxySFR', data=SFR)
        f.create_dataset('GalaxyMetallicity', data=Metallicity)
        # f.create_dataset('GalaxyKappa', data=kappa)
        # f.create_dataset('GalaxySpecificAngularMomentum', data=Lmomentum)

    data_file.close()

    load_particle_data(sim_info, output_file)
