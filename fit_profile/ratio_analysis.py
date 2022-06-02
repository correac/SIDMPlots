import h5py
import numpy as np
from scipy import interpolate

def calculate_r_fid(M200c):
    r_fid = 2.0 * (10**M200c / 1e9) ** 0.23
    return r_fid


def calculate_vcirc_at_fiducial_radius(radius, velocity, M200c):

    r_fid = calculate_r_fid(M200c)
    f = interpolate.interp1d(radius, velocity)
    v_fid = f(r_fid)
    return v_fid


def calculate_ratio(output_path, input_file, input_file_CDM):

    # Let's read some data :
    filename = output_path+"/"+input_file+".hdf5"
    with h5py.File(filename, "r") as file:
        M200c = file["Halo_data/M200c"][:]
        M200c = np.log10(M200c)
        Structure_type = file["Halo_data/StructureType"][:]
        Velocity = file["Profile_data/Circular_Velocity"][:][:]
        radial_bins = file["Profile_data/Velocity_radial_bins"][:]

    filename = output_path+"/"+input_file_CDM+".hdf5"
    with h5py.File(filename, "r") as file:
        CDM_M200c = file["Halo_data/M200c"][:]
        CDM_M200c = np.log10(CDM_M200c)
        CDM_Structure_type = file["Halo_data/StructureType"][:]
        CDM_Velocity = file["Profile_data/Circular_Velocity"][:][:]

    mass_bins = np.arange(9, 12, 0.1)
    num_mass_bins = len(mass_bins)-1

    halo_mass = []
    halo_ratio = []
    halo_mask = []
    halo_type = []

    for i in range(num_mass_bins):

        sample = np.where((CDM_M200c >= mass_bins[i]) & (CDM_M200c <= mass_bins[i+1]))[0]
        if len(sample) > 5:

            subsample = np.where(CDM_Structure_type[sample] == 10)[0]

            if len(subsample) > 5:
                median_CDM_Velocity_Centrals = np.median(CDM_Velocity[:, sample[subsample]], axis=1)
                median_CDM_M200c_Centrals = np.median(CDM_M200c[sample[subsample]])
                CDM_v_fid_centrals = calculate_vcirc_at_fiducial_radius(radial_bins,
                                                                        median_CDM_Velocity_Centrals,
                                                                        median_CDM_M200c_Centrals)
            else :
                CDM_v_fid_centrals = -1

            subsample = np.where(CDM_Structure_type[sample] > 10)[0]

            if len(subsample) > 5:
                median_CDM_Velocity_Satellites = np.median(CDM_Velocity[:, sample[subsample]], axis=1)
                median_CDM_M200c_Satellites = np.median(CDM_M200c[sample[subsample]])
                CDM_v_fid_satellites = calculate_vcirc_at_fiducial_radius(radial_bins,
                                                                          median_CDM_Velocity_Satellites,
                                                                          median_CDM_M200c_Satellites)

            else :
                CDM_v_fid_satellites = -1

        else :
            CDM_v_fid_centrals = -1
            CDM_v_fid_satellites = -1

        # Generate another sample array
        sample = np.where((M200c >= mass_bins[i]) & (M200c <= mass_bins[i+1]))[0]
        if len(sample) == 0: continue

        ratio_sub = np.zeros(len(sample))
        for j in range(len(sample)):

            v_fid = calculate_vcirc_at_fiducial_radius(radial_bins, Velocity[:, sample[j]], M200c[sample[j]])

            if Structure_type[sample[j]] == 10:
                ratio_sub[j] = v_fid - CDM_v_fid_centrals
                ratio_sub[j] /= CDM_v_fid_centrals
            else:
                ratio_sub[j] = v_fid - CDM_v_fid_satellites
                ratio_sub[j] /= CDM_v_fid_satellites


        halo_mass = np.append(halo_mass, M200c[sample])
        halo_ratio = np.append(halo_ratio, ratio_sub)
        halo_mask = np.append(halo_mask, sample)
        halo_type = np.append(halo_type, Structure_type[sample])


    return halo_mask, halo_ratio, halo_mass, halo_type
