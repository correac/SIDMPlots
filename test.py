import numpy as np
import h5py
from scipy import interpolate
import scipy.stats as stat
from pylab import *
from tqdm import tqdm
from object import particle_data
import matplotlib.pyplot as plt

def c_M_relation(log10_M0):
    """
    Concentration-mass relation from Correa et al. (2015).
    This relation is most suitable for Planck cosmology.
    """
    z = 0

    # Best-fit params:
    alpha = 1.7543 - 0.2766 * (1. + z) + 0.02039 * (1. + z) ** 2
    beta = 0.2753 + 0.0035 * (1. + z) - 0.3038 * (1. + z) ** 0.0269
    gamma = -0.01537 + 0.02102 * (1. + z) ** (-0.1475)

    log_10_c200 = alpha + beta * log10_M0 * (1. + gamma * log10_M0 ** 2)
    c200 = 10 ** log_10_c200
    return c200

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)


def calc_NFW_vcirc(r, log10_M):

    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21
    G = 6.67408e-11 # m^3 kg^-1 s^-2
    G *= (1e2)**3/1e3 # cm^3 g^-1 s^-2
    G /= kpc_in_cgs**3
    G *= Msun_in_cgs # kpc^3 Msun^-1 s^-2


    z = 0
    Om = 0.306
    Ol = 1-Om
    h = 0.6777
    rho_crit = 2.7754e11 * h ** 2 * (Om * (1. + z) ** 3 + Ol)
    R200 = 10**log10_M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    c = c_M_relation(log10_M)
    x = r / R200
    gcx = np.log(1. + c * x) - c * x / (1. + c * x)
    gc = np.log(1. + c) - c / (1. + c)
    mass = 200 * 4 * np.pi * R200**3 * rho_crit * 1e-9 * gcx / (3. * gc) # Msun

    Vcirc = G * mass / r # kpc^2/s^2
    Vcirc *= (kpc_in_cgs/1e5)**2 # km^2/s^2
    Vcirc = np.sqrt(Vcirc) # km/s
    return Vcirc


# Plot parameters
params = {
    "font.size": 12,
    "font.family": "Times",
    "text.usetex": True,
    "figure.figsize": (5, 4),
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.45,
    "figure.subplot.hspace": 0.35,
    "lines.markersize": 4,
    "lines.linewidth": 2,
}
rcParams.update(params)

figure()
ax = plt.subplot(1, 1, 1)
plt.grid("True")

M200c = np.arange(9,11,0.5)
radial_bins = np.arange(0.2, 25, 0.25)
centers = bin_centers(radial_bins) # kpc

for i in range(len(M200c)):
    v_circ = calc_NFW_vcirc(centers, M200c[i])
    Vmax = np.max(v_circ)
    r_fid = Vmax / 10
    print(r_fid,M200c[i],Vmax)
    f = interpolate.interp1d(centers, v_circ)
    v_fid = f(r_fid)
    plt.plot(centers, v_circ, '-')
    plt.plot([r_fid], [v_fid], 'o')
    plt.plot([r_fid], [Vmax], '+')

plt.plot(np.ones(2)*2., np.array([0,60]), '--',lw=1,color='grey')

plt.axis([0, 25, 0, 60])
plt.xlabel("Radius [kpc]")
plt.ylabel("Circular velocity [km/s]")
ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
plt.savefig("test.png", dpi=200)
plt.close()