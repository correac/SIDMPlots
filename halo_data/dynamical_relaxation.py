import numpy as np

def calculate_relaxation_criteria(sim_info, halo_index):

    sim_info.load_halo_catalogue(halo_index)

    MoP = np.array([sim_info.xminpot, sim_info.yminpot, sim_info.zminpot])
    CoM = np.array([sim_info.xcom, sim_info.ycom, sim_info.zcom])
    pos = MoP - CoM
    r = np.sqrt( pos[0,:]**2 + pos[1,:]**2 + pos[2,:]**2 )

    relaxation = r / (0.07 * sim_info.virial_radius)

    return relaxation
