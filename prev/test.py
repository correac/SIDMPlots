import numpy as np
from scipy.special import spence
from colossus.cosmology import cosmology
cosmology.setCosmology('planck13');
from colossus.lss import mass_function


def calc_density(x, M, a, z):
    h = 0.6777
    Om = 0.307
    rho_crit0 = 2.7754e11 * h**2 # Msun / Mpc^3
    Ol = 1. - Om
    rho_crit = rho_crit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    # NFW profile #
    c = R200 / a
    delta = 200. / 3.
    delta *= c ** 3 / (np.log(1. + c) - c / (1. + c))
    f = np.zeros(len(x))
    for i in range(0, len(x)): f[i] = rho_crit * delta * 1e-9 / ((x[i] / a) * (1. + x[i] / a) ** 2)
    return f


def sigma_1D(x, M, a, z):
    h = 0.6777
    Om = 0.307
    rho_crit0 = 2.7754e11 * h**2 # Msun / Mpc^3
    Ol = 1. - Om
    rho_crit = rho_crit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    G = 4.3e-6  # kpc km^2 Msun^-1 s^-2
    ff = np.zeros(len(x))
    for i in range(0, len(x)):
        s = x[i] / R200
        c = R200 / a
        u = 1. + c * s
        gc = 1. / (np.log(1. + c) - c / (1. + c))
        f = 0.5 * c ** 2 * gc * s * u ** 2 * G * M / R200
        u = 1. + c * s
        Li = spence(u)
        f *= (np.pi ** 2 - np.log(c * s) - 1. / (c * s) - 1. / u ** 2 - 6. / u + (
                    1. + 1. / (c ** 2 * s ** 2) - 4. / (c * s) - 2. / u) * np.log(u) + 3 * (np.log(u)) ** 2 + 6. * Li)
        ff[i] = f
    return np.sqrt(ff)


def gamma(M, z):

    sigma = 10
    h = 0.6777
    Om = 0.307
    rho_crit0 = 2.7754e11 * h**2 # Msun / Mpc^3
    Ol = 1. - Om
    rho_crit = rho_crit0 * (Om * (1. + z) ** 3 + Ol)
    R200 = M / (4. * np.pi * 200 * rho_crit / 3.)
    R200 = R200 ** (1. / 3.)  # Mpc
    R200 *= 1e3  # kpc

    c200 = 5.72 * (M / 1e14) ** (-0.081) * (1. + z) ** (-0.71)

    a = R200 / c200  # r200/ r200/rs [kpc]

    radial_bins = np.arange(-5, np.log10(R200), 0.1)
    x = 10 ** radial_bins
    deltax = x[1:] - x[:1]
    deltax = np.append(deltax, deltax[-1])
    Msun_in_cgs = 1.98848e33
    kpc_in_cgs = 3.08567758e21

    sig1D = sigma_1D(x, M, a, z)  # km/s
    sig1D[np.isnan(sig1D)] = 0
    sig1D *= 1e5  # cm/s

    rho = calc_density(x, M, a, z)  # Msun/kpc^3
    rho *= Msun_in_cgs / kpc_in_cgs ** 3  # g/cm^2

    f = sig1D * rho * (4. / np.sqrt(np.pi)) * sigma  # 1/s
    Gyr = 1e9 * 365.25 * 24 * 3600

    f *= Gyr  # Gyr^-1 particle^-1

    rho = calc_density(x, M, a, z)  # Msun/kpc^3
    vol = 4 * np.pi * x ** 2 * deltax * rho  # Msun
    Total = np.sum(f * vol) / M
    return Total


def analytic_scatter(z, Mmin, Mmax):
    h = 0.6777
    Om = 0.307
    rho_crit0 = 2.7754e11 * h**2 # Msun / Mpc^3
    rho_crit = rho_crit0 #* ( Om * (1.+ z)**3 + Ol)

    M = np.arange(Mmin, Mmax, 0.25)
    M = 10 ** M
    mfunc = mass_function.massFunction(M / h, z, mdef='200c', model='tinker08', q_out='dndlnM')
    mfunc /= h**3
    mfunc *= M
    mfunc /= rho_crit
    f = mfunc

    G = np.zeros(len(M))
    for i in range(0, len(M)): G[i] = gamma(M[i], z)  # [Gyr^-1 particle^-1]

    n = f * G
    n = np.sum(n * 0.25 / np.log(10.))
    return n  # [Gyr^-1 particle^-1]