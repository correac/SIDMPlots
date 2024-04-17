import numpy as np
import h5py
from object import particle_data
from tqdm import tqdm

def calculate_magnitude(sim_info, sample_satellites):

    num_haloes = len(sample_satellites)
    GalaxyLuminosity =  np.zeros(num_haloes)

    halo_index = sim_info.halo_data.halo_index[sample_satellites]

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.stars.ids)

        if len(bound_particles_only) < 1: continue
        GalaxyLuminosity[i] = np.sum(part_data.stars.luminosities_r_band[bound_particles_only])

    Lum_Sun = 3.828e26
    GalaxyMr = -2.5 * np.log10(GalaxyLuminosity) + 2.5 * np.log10(Lum_Sun) + 4.76
    GalaxyMv = GalaxyMr + 0.2
    return GalaxyMv, GalaxyMr, GalaxyLuminosity


def look_for_satellites(sim_info, sample_host):

    host_index = []
    satellites_index = []
    satellites_sample = []
    num_sample = len(sample_host)

    for i in range(num_sample):

        x = sim_info.halo_data.xminpot - sim_info.halo_data.xminpot[sample_host[i]]
        y = sim_info.halo_data.yminpot - sim_info.halo_data.yminpot[sample_host[i]]
        z = sim_info.halo_data.zminpot - sim_info.halo_data.zminpot[sample_host[i]]

        distance = np.sqrt(x**2 + y**2 + z**2)
        Rvir = sim_info.halo_data.virial_radius[sample_host[i]]

        select = np.where((distance>0) & (distance<= Rvir))[0]
        flag = sim_info.halo_data.satellite_flag[select]
        sub_select = np.where(flag == True)[0]
        select = select[sub_select]

        sub_select = np.where(sim_info.halo_data.log10_halo_mass[select] > 7)[0]
        select = select[sub_select]
        num_satellites = len(select)

        satellites_sample = np.append(satellites_sample, select)
        satellites_index = np.append(satellites_index, sim_info.halo_data.halo_index[select])
        host_index = np.append(host_index, np.ones(num_satellites) * sim_info.halo_data.halo_index[sample_host[i]] )

    satellites_sample = satellites_sample.astype('int')
    return satellites_sample, satellites_index, host_index

def make_satellites_catalog(sim_info):

    sample_host = np.where((sim_info.halo_data.log10_halo_mass >= 11.9) &
                           (sim_info.halo_data.log10_halo_mass <= 12.1))[0]
    centrals = np.where(sim_info.halo_data.structure_type[sample_host] == 10)[0]
    sample_host = sample_host[centrals]
    print('num host',len(sample_host))

    sample_host = sample_host[0:1]

    sample_satellites, index_satellites, index_host = look_for_satellites(sim_info, sample_host)

    print('num satellites',len(sample_satellites))

    M200c = sim_info.halo_data.log10_halo_mass[sample_satellites]
    M200c_host = sim_info.halo_data.log10_halo_mass[sample_host]
    Vmax = sim_info.halo_data.vmax[sample_satellites]

    # # Morphology/shape estimations
    Mv_satellite, Mr_satellite, Luminosity_satellite = calculate_magnitude(sim_info, sample_satellites)


    # Output data
    output_file = f"{sim_info.output_path}/Satellites_data_" + sim_info.simulation_name + ".hdf5"
    data_file = h5py.File(output_file, 'w')
    f = data_file.create_group('Halo_data')
    f.create_dataset('IDs_Satellite', data=index_satellites)
    f.create_dataset('IDs_Host', data=index_host)
    f.create_dataset('M200c_Host', data=M200c_host)
    f.create_dataset('M200c', data=M200c)
    f.create_dataset('Vmax', data=Vmax)
    f.create_dataset('Mv', data=Mv_satellite)
    f.create_dataset('Mr', data=Mr_satellite)
    f.create_dataset('Luminosity_r_band', data=Luminosity_satellite)

    Mstar = sim_info.halo_data.log10_stellar_mass[sample_satellites]
    Mgas = sim_info.halo_data.log10_gas_mass[sample_satellites]
    GalaxyHalfMassRadius = sim_info.halo_data.half_mass_radius_star[sample_satellites]
    GalaxyProjectedHalfMassRadius = sim_info.halo_data.half_mass_projected_radius_star[sample_satellites]
    SFR = sim_info.halo_data.sfr[sample_satellites]
    Metallicity = sim_info.halo_data.metallicity_stars[sample_satellites]

    f.create_dataset('Mstar', data=Mstar)
    f.create_dataset('Mgas', data=Mgas)
    f.create_dataset('GalaxyHalfMassRadius', data=GalaxyHalfMassRadius)
    f.create_dataset('GalaxyProjectedHalfMassRadius', data=GalaxyProjectedHalfMassRadius)
    f.create_dataset('SFR', data=SFR)
    f.create_dataset('Metallicity', data=Metallicity)

    data_file.close()


