"""
Description here
"""
from argumentparser import ArgumentWithInputFiles
import h5py
from fit_profile.fit_density_and_shape import fit_density_and_shape
from fit_profile.fit_profile import fit_density
from fit_profile.ratio_analysis import calculate_ratio


if __name__ == "__main__":

    config = ArgumentWithInputFiles()

    # Fetch relevant input parameters from lists
    sim_name = config.name_list[0]
    output = config.output_directory
    input_file = config.input_file_list
    # input_file_CDM = config.input_file_list_CDM
    # sim_type = config.sim_type_list[0]
    print(sim_name)
    print(input_file)

    # Gather data
    # halo_mask, halo_ratio, halo_mass, halo_type = calculate_ratio(output, input_file, input_file_CDM, sim_type)
    # n0_pse, r0_pse, rho0_pse, r0_iso, rho0_iso = fit_density(halo_mask, output, input_file)

    rs, rhos, gamma, q, halo_mask, halo_mass, halo_type = fit_density_and_shape(input_file)

    output_file = "DensityFitParams_" + sim_name + ".hdf5"
    data_file = h5py.File(output_file, 'w')
    f = data_file.create_group('HaloData')
    f.create_dataset('HaloID', data=halo_mask)
    # f.create_dataset('VfidRatio', data=halo_ratio)
    f.create_dataset('M200c', data=halo_mass)
    f.create_dataset('StructureType', data=halo_type)

    f = data_file.create_group('BestFitProfile')
    f.create_dataset('rs', data=rs)
    f.create_dataset('rhos', data=rhos)
    f.create_dataset('gamma', data=gamma)
    f.create_dataset('q', data=q)


    # f = data_file.create_group('IsothermalProfile')
    # f.create_dataset('r0', data=r0_iso)
    # f.create_dataset('rho0', data=rho0_iso)

    # f = data_file.create_group('PseudoIsothermalProfile')
    # f.create_dataset('n0', data=n0_pse)
    # f.create_dataset('r0', data=r0_pse)
    # f.create_dataset('rho0', data=rho0_pse)

    data_file.close()
