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
    filename = f"{output_path}/f{input_file}.hdf5"
    with h5py.File(filename, "r") as file:
        M200c = file["Assembly_history/Mass"][:][:]
        M200c = np.log10(M200c[:,0])
        Velocity = file["Profile_evolution/Velocity_snapshot_0036"][:][:]
        radial_bins = file["Profile_evolution/Velocity_radial_bins"][:]

    filename = f"{output_path}/f{input_file_CDM}.hdf5"
    with h5py.File(filename, "r") as file:
        CDM_M200c = file["Assembly_history/Mass"][:][:]
        CDM_M200c = np.log10(CDM_M200c[:,0])
        CDM_Velocity = file["Profile_evolution/Velocity_snapshot_0036"][:][:]

    mass_bins = np.arange(9, 11.2, 0.1)
    num_mass_bins = len(mass_bins)-1

    halo_mass = []
    halo_ratio = []
    halo_mask = []

    for i in range(num_mass_bins):

        sample = np.where((CDM_M200c >= mass_bins[i]) & (CDM_M200c <= mass_bins[i+1]))[0]
        median_CDM_Velocity = np.median(CDM_Velocity[:, sample], axis=1)
        median_CDM_M200c = np.median(CDM_M200c[sample])
        CDM_v_fid = calculate_vcirc_at_fiducial_radius(radial_bins, median_CDM_Velocity, median_CDM_M200c)

        # Generate another sample array
        sample = np.where((M200c >= mass_bins[i]) & (M200c <= mass_bins[i+1]))[0]
        ratio_sub = np.zeros(len(sample))
        for j in range(len(sample)):

            v_fid = calculate_vcirc_at_fiducial_radius(radial_bins, Velocity[:,sample[j]], M200c[sample[j]])
            ratio_sub[j] = v_fid - CDM_v_fid
            ratio_sub[j] /= CDM_v_fid

        halo_mass = np.append(halo_mass, M200c[sample])
        halo_ratio = np.append(halo_ratio, ratio_sub)
        halo_mask = np.append(halo_mask, sample)

    return halo_mask, halo_ratio, halo_mass