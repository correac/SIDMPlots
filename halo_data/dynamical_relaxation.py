import numpy as np

def calculate_relaxation_criteria(sim_info, halo_index):

    halo_data = sim_info.load_halo_catalogue(halo_index)

    MoP = np.array([halo_data.xminpot, halo_data.yminpot, halo_data.zminpot])
    CoM = np.array([halo_data.xcom, halo_data.ycom, halo_data.zcom])
    pos = MoP - CoM
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    relaxation = r / (0.07 * halo_data.virial_radius)

    return relaxation