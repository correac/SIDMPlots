import numpy as np
import scipy.stats as stat

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)

def calculate_Vcirc(mass, pos, radial_bins):
    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21
    G = 6.67408e-11 #m^3 kg^-1 s^-2
    G *= (1e2)**3/1e3 #cm^3 g^-1 s^-2
    G /= kpc_in_cgs**3
    G *= Msun_in_cgs #kpc^3 Msun^-1 s^-2

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))
    centers = bin_centers(radial_bins)  # kpc

    SumMasses, _, _ = stat.binned_statistic(x=r, values=mass, statistic="sum", bins=radial_bins, )
    M_within_r = np.cumsum(SumMasses)
    V_circ_2 = G * M_within_r / centers #kpc^2 s^-2
    V_circ_2 *= (kpc_in_cgs/1e5)**2 #km^2/s^2
    V_circ = np.sqrt(V_circ_2)

    return V_circ