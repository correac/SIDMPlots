import h5py
import numpy as np
import emcee
from multiprocessing import Pool
import time
from scipy.optimize import minimize
from pylab import *
import matplotlib.pyplot as plt

def make_plot(xdata, ydata, popt, halo):
    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.16,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.4,
        "figure.subplot.hspace": 0.4,
        "lines.markersize": 6,
        "lines.linewidth": 1.0,
        "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    ######################
    figure()

    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(xdata, ydata, lw=2, color='tab:blue')

    model = fit_density_model(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
    plt.plot(xdata, 10**model, lw=1, color='black')

    plt.ylim(1e4, 1e9)
    plt.xlim(0.1, 1e3)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Radius [kpc]")
    plt.ylabel("Density [M$_{\odot}$/kpc$^{3}$]")
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig('output_data_%i.png'%halo, dpi=300)

def fit_density_model(xdata, gamma, beta, rs, rhos, q):

    R = xdata * np.sqrt((1. + q**2)/2.)
    #R = xdata * np.sqrt((1. + 1./q**2)/2.)

    yrange = rhos - gamma * np.log10(R / rs)
    # yrange -= (3. - gamma) * np.log10(1. + R / rs)
    yrange -= beta * np.log10(1. + R / rs)


    return yrange

def log_prior(theta):
    """
    The natural logarithm of the prior probability.
    It sets prior to 1 (log prior to 0) if params are in range, and zero (-inf) otherwise.
    Args: theta (tuple): a sample containing individual parameter values
    """
    gamma, beta, rs, rhos, q = theta
    log_prior = -np.inf

    if 0 < gamma < 10 and 0 < beta < 10 and 0.1 < rs < 50 and 5 < rhos < 12 and 0 < q < 1 : log_prior = 0.0

    return log_prior

def log_posterior(theta, x, y, yerr):
    """
    The natural logarithm of the joint posterior.

    Args:
        theta (tuple): a sample containing individual parameter values
        x (array): values over which the data/model is defined
        y (array): the set of data
        yerr (array): the standard deviation of the data points
    """
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

def log_likelihood(theta, x, y, yerr):
    """
    The natural logarithm of the joint likelihood.

    Args:
        theta (tuple): a sample containing individual parameter values
        x (array): values over which the data/model is defined
        y (array): the set of data
        yerr (array): the standard deviation of the data points
    """
    gamma, beta, rs, rhos, q = theta
    model = fit_density_model(x, gamma, beta, rs, rhos, q)
    sigma2 = yerr**2
    log_l = -0.5 * np.sum((y - model) ** 2 / sigma2)
    return log_l

def run_mcmc(x, y, yerr, soln):

    pos = soln.x + 1e-4 * np.random.randn(32, 5)
    nwalkers, ndim = pos.shape

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x, y, yerr), pool=pool)
        start = time.time()
        sampler.run_mcmc(pos, 500, progress=True)
        end = time.time()
        multi_time = end - start
        print("Multiprocessing took {0:.1f} minutes".format(multi_time / 60))

    samples = sampler.get_chain(discard=100, thin=15, flat=True)
    gamma = np.median(samples[:, 0])
    beta = np.median(samples[:, 1])
    rs = np.median(samples[:, 2])
    rhos = np.median(samples[:, 3])
    q = np.median(samples[:, 4])

    print("Mean autocorrelation time: {0:.3f} steps".format(np.mean(sampler.get_autocorr_time(quiet=True))))
    return gamma, rs, rhos, q

def fit_profile(x, y, R200):

    nozero = y > 0
    x = x[nozero]
    y = y[nozero]

    # Let's define inner region as everything within 25 kpc
    inner = np.where((x <= R200) & (x >= 2.0))[0]

    x = x[inner]
    y = np.log10(y[inner])
    yerr = np.ones(len(y)) * 0.1

    # First fit profile based on Isothermal model
    np.random.seed(42)
    nll = lambda *args: -log_likelihood(*args)
    initial = np.array([1.0, 3.0, 2, 8, 0.5])
    soln = minimize(nll, initial, args=(x, y, yerr))
    gamma, beta, rs, rhos, q = soln.x
    print('=======')
    print('Gamma, Beta, Rs, Rhos, q')
    print(gamma, beta, rs, rhos, q)

    # gamma, rs, rhos, q = run_mcmc(x, y, yerr, soln)
    # print(gamma, rs, rhos, q)

    sol = np.array([gamma, beta, rs, rhos, q])

    return sol


def fit_density_and_shape(filename):

    # Let's read some data :
    with h5py.File(filename, "r") as file:
        M200c = file["Halo_data/M200c"][:]
        R200c = file["Halo_data/R200c"][:] * 1e3
        Structure_type = file["Halo_data/StructureType"][:]
        radial_bins = file["Profile_data/Density_radial_bins"][:]
        Density = file["Profile_data/Dark_matter_Density_profile"][:][:]

    select = np.where(M200c >= 10)[0]
    select_type = np.where(Structure_type[select] ==10)[0]
    halo_mask = select[select_type[0:10]]

    num_sample = len(halo_mask)
    gamma = np.zeros(num_sample)
    beta = np.zeros(num_sample)
    rs = np.zeros(num_sample)
    rhos = np.zeros(num_sample)
    q = np.zeros(num_sample)

    for j in range(num_sample):

        print('halo ',j,'out of ',num_sample)
        Density_halo_j = Density[:,halo_mask[j].astype('int')]
        popt = fit_profile(radial_bins, Density_halo_j, R200c[halo_mask[j]])

        gamma[j] = popt[0]
        beta[j] =  popt[1]
        rs[j] = popt[2]
        rhos[j] = popt[3]
        q[j] = popt[4]

        make_plot(radial_bins, Density_halo_j, popt, j)

    return rs, rhos, gamma, q, halo_mask, M200c[halo_mask], Structure_type[halo_mask]