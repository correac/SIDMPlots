import h5py
import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
import emcee
from multiprocessing import Pool
import time
from scipy.optimize import minimize

def log_prior_pse(theta):
    """
    The natural logarithm of the prior probability.
    It sets prior to 1 (log prior to 0) if params are in range, and zero (-inf) otherwise.
    Args: theta (tuple): a sample containing individual parameter values
    """
    r0, rho0, n0 = theta
    log_prior = -np.inf
    if 0.001 < r0 < 5 and 5 < rho0 < 9 and 0 < n0 < 2 : log_prior = 0.0

    return log_prior

def log_prior(theta):
    """
    The natural logarithm of the prior probability.
    It sets prior to 1 (log prior to 0) if params are in range, and zero (-inf) otherwise.
    Args: theta (tuple): a sample containing individual parameter values
    """
    r0, rho0 = theta
    log_prior = -np.inf
    if 0.001 < r0 < 5 and 5 < rho0 < 9 : log_prior = 0.0

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

def log_posterior_pse(theta, x, y, yerr):
    """
    The natural logarithm of the joint posterior.

    Args:
        theta (tuple): a sample containing individual parameter values
        x (array): values over which the data/model is defined
        y (array): the set of data
        yerr (array): the standard deviation of the data points
    """
    lp = log_prior_pse(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_pse(theta, x, y, yerr)

def log_likelihood(theta, x, y, yerr):
    """
    The natural logarithm of the joint likelihood.

    Args:
        theta (tuple): a sample containing individual parameter values
        x (array): values over which the data/model is defined
        y (array): the set of data
        yerr (array): the standard deviation of the data points
    """
    r0, rho0 = theta
    model = fit_isothermal_model(x, r0, rho0)
    sigma2 = yerr**2
    log_l = -0.5 * np.sum((y - model) ** 2 / sigma2)
    return log_l

def log_likelihood_pse(theta, x, y, yerr):
    """
    The natural logarithm of the joint likelihood.

    Args:
        theta (tuple): a sample containing individual parameter values
        x (array): values over which the data/model is defined
        y (array): the set of data
        yerr (array): the standard deviation of the data points
    """
    r0, rho0, n0 = theta
    model = fit_pseudo_isothermal_model(x, r0, rho0, n0)
    sigma2 = yerr**2
    log_l = -0.5 * np.sum((y - model) ** 2 / sigma2)
    return log_l

def diff_isothermal_equation(f,x,n):
    """
    Differential equation that describes the isothermal profile
    """
    y, z = f
    dfdx = [z,-(n+2)*(1./x)*z-n*(n+1)*(1./x**2)-(1./x**n)*np.exp(y)]
    return dfdx

def fit_pseudo_isothermal_model(xdata, a, b, n):
    """
    For this isothermal form I let the slope to vary.
    """
    if n<0: return 0
    xrange = np.arange(-5, 5, 0.01)
    xrange = 10**xrange
    xrange = xrange / a
    y0 = [0, 0]

    sol = odeint(diff_isothermal_equation, y0, xrange, args=(n,))
    yrange = np.exp(sol[:, 0])
    yrange = np.log10(yrange)
    finterpolate = interp1d(xrange, yrange)
    x = xdata / a

    max_x = 1e6
    if len(x) > 1: max_x = np.max(x)
    if a<=0 or max_x>1e5 : return 0

    ydata = finterpolate(x)
    f = b + ydata
    return f

def fit_isothermal_model(xdata, a, b):
    xrange = np.arange(-5, 5, 0.01)
    xrange = 10**xrange
    xrange = xrange / a
    y0 = [0, 0]
    n = 0

    sol = odeint(diff_isothermal_equation, y0, xrange, args=(n,))
    yrange = np.exp(sol[:, 0])
    yrange = np.log10(yrange)
    finterpolate = interp1d(xrange, yrange)
    x = xdata / a

    max_x = 1e6
    if len(x) > 1: max_x = np.max(x)
    if a<=0 or max_x>1e5 : return 0

    ydata = finterpolate(x)
    f = b + ydata
    return f

def run_mcmc(x, y, yerr, soln):

    pos = soln.x + 1e-4 * np.random.randn(32, 2)
    nwalkers, ndim = pos.shape

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x, y, yerr), pool=pool)
        start = time.time()
        sampler.run_mcmc(pos, 500, progress=True)
        end = time.time()
        multi_time = end - start
        print("Multiprocessing took {0:.1f} minutes".format(multi_time / 60))

    samples = sampler.get_chain(discard=100, thin=15, flat=True)
    r0 = np.median(samples[:, 0])
    rho0 = np.median(samples[:, 1])

    print("Mean autocorrelation time: {0:.3f} steps".format(np.mean(sampler.get_autocorr_time(quiet=True))))
    return r0, rho0

def run_mcmc_pse(x, y, yerr, soln):

    pos = soln.x + 1e-4 * np.random.randn(32, 3)
    nwalkers, ndim = pos.shape

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_pse, args=(x, y, yerr), pool=pool)
        start = time.time()
        sampler.run_mcmc(pos, 500, progress=True)
        end = time.time()
        multi_time = end - start
        print("Multiprocessing took {0:.1f} minutes".format(multi_time / 60))

    samples = sampler.get_chain(discard=100, thin=15, flat=True)
    r0 = np.median(samples[:, 0])
    rho0 = np.median(samples[:, 1])
    n0 = np.median(samples[:, 2])

    print("Mean autocorrelation time: {0:.3f} steps".format(np.mean(sampler.get_autocorr_time(quiet=True))))
    return r0, rho0, n0


def fit_profile(x, y):

    nozero = y > 0
    x = x[nozero]
    y = y[nozero]

    # Let's define inner region as everything within 25 kpc
    inner = np.where((x <= 5) & (x >= 0.5))[0]

    x = x[inner]
    y = np.log10(y[inner])
    yerr = np.ones(len(y)) * 0.1

    # First fit profile based on Isothermal model
    np.random.seed(42)
    nll = lambda *args: -log_likelihood(*args)
    initial = np.array([1, 7])
    soln = minimize(nll, initial, args=(x, y, yerr))
    r0, rho0 = soln.x
    print('=======')
    print(r0, rho0)

    r0, rho0 = run_mcmc(x, y, yerr, soln)
    print(r0, rho0)
    sol_iso = np.array([r0, rho0])

    # Let's also consider pseudo-Isothermal model
    np.random.seed(42)
    nll = lambda *args: -log_likelihood_pse(*args)
    initial = np.array([1, 7, 0])
    soln = minimize(nll, initial, args=(x, y, yerr))
    r0, rho0, n0 = soln.x
    print('=======')
    print(r0, rho0, n0)

    r0, rho0, n0 = run_mcmc_pse(x, y, yerr, soln)
    print(r0, rho0, n0)
    sol_pse = np.array([r0, rho0, n0])
    print('=======')

    return sol_iso, sol_pse


def fit_density(halo_mask, output_path, input_file):

    num_sample = len(halo_mask)
    n0_pse = np.zeros(num_sample)
    r0_pse = np.zeros(num_sample)
    rho0_pse = np.zeros(num_sample)
    r0_iso = np.zeros(num_sample)
    rho0_iso = np.zeros(num_sample)

    filename = output_path+"/"+input_file+".hdf5"
    with h5py.File(filename, "r") as file:
        radial_bins = file["Profile_data/Density_radial_bins"][:]
        Density = file["Profile_data/Dark_matter_Density_profile"][:][:]

    for j in range(num_sample):
        print('halo ',j,'out of ',num_sample)
        Density_halo_j = Density[:,halo_mask[j].astype('int')]
        popt_iso, popt_pse = fit_profile(radial_bins, Density_halo_j)

        n0_pse[j] = popt_pse[2]
        r0_pse[j] = popt_pse[0]
        rho0_pse[j] = popt_pse[1]
        r0_iso[j] = popt_iso[0]
        rho0_iso[j] = popt_iso[1]

    return n0_pse, r0_pse, rho0_pse, r0_iso, rho0_iso
