"""
Description here
"""
from argumentparser import ArgumentWithInputFiles
import h5py
from fit_profile.fit_profile import fit_density
from fit_profile.ratio_analysis import calculate_ratio


if __name__ == "__main__":

    config = ArgumentWithInputFiles()

    # Fetch relevant input parameters from lists
    sim_name = config.name_list[0]
    output = config.output_directory
    input_file = config.input_file_list
    input_file_CDM = config.input_file_list_CDM

    # Gather data
    halo_mask, halo_ratio, halo_mass, halo_type = calculate_ratio(output, input_file, input_file_CDM)

    n0_pse, r0_pse, rho0_pse, r0_iso, rho0_iso = fit_density(halo_mask, output, input_file)

    output_file = "DensityFitParams_" + sim_name + ".hdf5"
    data_file = h5py.File(output_file, 'a')
    f = data_file.create_group('HaloData')
    f.create_dataset('HaloID', data=halo_mask)
    f.create_dataset('VfidRatio', data=halo_ratio)
    f.create_dataset('M200c', data=halo_mass)
    f.create_dataset('StructureType', data=halo_type)

    f = data_file.create_group('IsothermalProfile')
    f.create_dataset('r0', data=r0_iso)
    f.create_dataset('rho0', data=rho0_iso)

    f = data_file.create_group('PseudoIsothermalProfile')
    f.create_dataset('n0', data=n0_pse)
    f.create_dataset('r0', data=r0_pse)
    f.create_dataset('rho0', data=rho0_pse)

    data_file.close()
