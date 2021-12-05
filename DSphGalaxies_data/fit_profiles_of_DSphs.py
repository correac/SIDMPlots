import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.optimize import minimize
import emcee
import corner

from multiprocessing import Pool
import time


from pylab import *
import matplotlib.pyplot as plt


def read_data(file):
    data = np.loadtxt(file)
    # Remove potential zeros before we begin
    nozero = data[:,2] > 0
    data = data[nozero,:]
    min_rho = np.log10(data[:,1]) >= 5.
    data = data[min_rho, :]
    min_r = data[:,0] >= 1e-1
    data = data[min_r, :]

    # Next:
    x = data[:,0]
    y = np.log10(data[:,1])
    yerr = np.sqrt( (y-np.log10(data[:,2]))**2 + (np.log10(data[:,3])-y)**2)
    errorbar = np.log10( np.sqrt( (data[:,1]-data[:,2])**2 + (data[:,3]-data[:,1])**2) )

    return x, y, yerr, errorbar

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
    xrange = np.arange(-5, 5, 0.01)
    xrange = 10**xrange
    xrange = xrange / a
    y0 = [0, 0]

    sol = odeint(diff_isothermal_equation, y0, xrange, args=(n,))
    yrange = np.exp(sol[:, 0])
    yrange = np.log10(yrange)
    finterpolate = interp1d(xrange, yrange)
    x = xdata / a
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
    ydata = finterpolate(x)
    f = b + ydata
    return f

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

def log_likelihood_pseudo(theta, x, y, yerr):
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

def log_prior(theta):
    """
    The natural logarithm of the prior probability.

    It sets prior to 1 (log prior to 0) if params are in range, and zero (-inf) otherwise.

    Args:
        theta (tuple): a sample containing individual parameter values
    """
    r0, rho0 = theta

    log_prior = -np.inf

    if 0 < r0 < 5 and 6. < rho0 < 10:
        log_prior = 0

    return log_prior

def log_prior_pseudo(theta):
    """
    The natural logarithm of the prior probability.

    It sets prior to 1 (log prior to 0) if params are in range, and zero (-inf) otherwise.

    Args:
        theta (tuple): a sample containing individual parameter values
    """
    r0, rho0, n0 = theta

    log_prior = -np.inf

    if 0 < r0 < 5 and 6. < rho0 < 10 and -1 < n0 < 1:
        log_prior = 0

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

def log_posterior_pseudo(theta, x, y, yerr):
    """
    The natural logarithm of the joint posterior.

    Args:
        theta (tuple): a sample containing individual parameter values
        x (array): values over which the data/model is defined
        y (array): the set of data
        yerr (array): the standard deviation of the data points
    """
    lp = log_prior_pseudo(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood_pseudo(theta, x, y, yerr)


def fit_profile(x, y, yerr, dwarf, check):

    # Let's define inner region as everything within 25 kpc
    # inner = x <= 100

    # np.random.seed(42)
    # nll = lambda *args: -log_likelihood(*args)
    # initial = np.array([1, 10, 0])
    # soln = minimize(nll, initial, args=(x[inner], y[inner], yerr[inner]))
    # r0_full_iso, rho0_full_iso, n0 = soln.x

    # First fit profile based on Isothermal model
    # popt, pcov = curve_fit(fit_isothermal_model, x[inner], y[inner], p0=[2, 8.5])
    # r0 = popt[0]
    # rho0 = popt[1]

    # popt, pcov = curve_fit(fit_pseudo_isothermal_model, x[inner], y[inner], p0=[0.5, 8.5, 0])
    # r0_pse_iso = popt[0]
    # rho0_pse_iso = popt[1]
    # ns0_pse_iso = popt[2]


    pos = np.array([0.1, 9]) + 1e-3 * np.random.randn(64, 2)
    nwalkers, ndim = pos.shape

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=(x, y, yerr), pool=pool)
        start = time.time()
        sampler.run_mcmc(pos, 1000, progress=True)
        end = time.time()
        multi_time = end - start
        print("Multiprocessing took {0:.1f} minutes".format(multi_time / 60))

    labels = ["$r_{0}$", r"log$_{10}\rho_{0}$"]
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

    print(dwarf+" ==================")
    print("Fits for isothermal model")
    r0_iso = np.median(flat_samples[:, 0])
    low_r0_iso = np.percentile(flat_samples[:,0],16)
    high_r0_iso = np.percentile(flat_samples[:,0],84)
    print('r0',r0_iso, low_r0_iso, high_r0_iso)
    rho0_iso = np.median(flat_samples[:, 1])
    low_rho0_iso = np.percentile(flat_samples[:,1],16)
    high_rho0_iso = np.percentile(flat_samples[:,1],84)
    print('rho0',rho0_iso, low_rho0_iso, high_rho0_iso)


    # Plot parameters
    params = {
        "font.size": 12,
        "font.family": "Times",
        "text.usetex": True,
        "figure.figsize": (4, 3),
        "figure.subplot.left": 0.18,
        "figure.subplot.right": 0.95,
        "figure.subplot.bottom": 0.18,
        "figure.subplot.top": 0.95,
        "figure.subplot.wspace": 0.25,
        "figure.subplot.hspace": 0.25,
        "lines.markersize": 6,
        "lines.linewidth": 1.5,
        "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    # Make the base corner plot
    fig_iso = corner.corner(
        flat_samples, labels=labels, quantiles=[0.16, 0.84],
        truths=[r0_iso, rho0_iso], show_titles=True,
        title_kwargs={"fontsize": 16}
    )

    plt.savefig('corner_plot_isothermal_model_'+dwarf+'.png', dpi=200)
    plt.close()

    pos = np.array([0.1, 9, 0]) + 1e-3 * np.random.randn(64, 3)
    nwalkers, ndim = pos.shape

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_pseudo, args=(x, y, yerr), pool=pool)
        start = time.time()
        sampler.run_mcmc(pos, 1000, progress=True)
        end = time.time()
        multi_time = end - start
        print("Multiprocessing took {0:.1f} minutes".format(multi_time / 60))

    labels = ["$r_{0}$", r"log$_{10}\rho_{0}$", "$n_{0}$"]
    flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)

    print("Fits for pseudo-isothermal model")
    r0_pse = np.median(flat_samples[:, 0])
    low_r0_pse = np.percentile(flat_samples[:,0],16)
    high_r0_pse = np.percentile(flat_samples[:,0],84)
    print('r0',r0_pse, low_r0_pse, high_r0_pse)
    rho0_pse = np.median(flat_samples[:, 1])
    low_rho0_pse = np.percentile(flat_samples[:,1],16)
    high_rho0_pse = np.percentile(flat_samples[:,1],84)
    print('rho0',rho0_pse, low_rho0_pse, high_rho0_pse)
    n0_pse = np.median(flat_samples[:, 2])
    low_n0_pse = np.percentile(flat_samples[:,2],16)
    high_n0_pse = np.percentile(flat_samples[:,2],84)
    print('n0',n0_pse, low_n0_pse, high_n0_pse)

    # Make the base corner plot
    fig_pse = corner.corner(
        flat_samples, labels=labels, quantiles=[0.16, 0.84],
        truths=[r0_pse, rho0_pse, n0_pse], show_titles=True,
        title_kwargs={"fontsize": 16}
    )

    plt.savefig('corner_plot_pseudo_isothermal_model_'+dwarf+'.png', dpi=200)
    plt.close()



    if check == 'yes':

        # Plot parameters
        params = {
            "font.size": 12,
            "font.family": "Times",
            "text.usetex": True,
            "figure.figsize": (4, 3),
            "figure.subplot.left": 0.18,
            "figure.subplot.right": 0.95,
            "figure.subplot.bottom": 0.18,
            "figure.subplot.top": 0.95,
            "figure.subplot.wspace": 0.25,
            "figure.subplot.hspace": 0.25,
            "lines.markersize": 4,
            "lines.linewidth": 1.5,
            "figure.max_open_warning": 0,
        }
        rcParams.update(params)

        #######################
        # Plot the density profile
        figure()
        ax = plt.subplot(1, 1, 1)
        grid(True)

        plot(x, y, '-', color='tab:orange',label=dwarf)
        plt.fill_between(x, y - yerr / 2, y + yerr / 2, alpha=0.4, color='grey')

        xrange = np.arange(-2,1,0.2)
        xrange = 10**xrange
        plot(xrange, fit_isothermal_model(xrange, r0_iso, rho0_iso),'--',
            color='tab:red',label='Isothermal fit')

        plot(xrange, fit_pseudo_isothermal_model(xrange, r0_pse, rho0_pse, n0_pse),'--',
             color='tab:blue',label='Pseudo-Isothermal fit')

        xscale('log')
        xlabel(r'r [kpc]')
        ylabel(r'$\log_{10}\rho$ [M$_{\odot}$/kpc$^{3}$]')
        axis([1e-2, 1e1, 5, 9])
        ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
        plt.legend(loc='lower left')
        plt.savefig("DSh_"+dwarf+"_fit.png", dpi=200)
        plt.close()

    output_file = "DSh_"+dwarf+"_fit.txt"
    output = np.zeros((2, 9))
    output[0, :] = np.array([r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso, 0, 0, 0])
    output[1, :] = np.array([r0_pse, high_r0_pse, low_r0_pse,
                             rho0_pse, high_rho0_pse, low_rho0_pse,
                             n0_pse, high_n0_pse, low_n0_pse])
    np.savetxt(output_file, output, fmt="%s")
    return

def read_fits(file):
    data = np.loadtxt(file)
    r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso = data[0,0:6]
    r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse = data[1,0:6]
    return r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso,\
           r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse

def plot_data(DSphs):


    # Plot parameters
    params = {
    "font.size": 10,
    "font.family": "Times",
    "text.usetex": True,
    "figure.figsize": (4, 3),
    "figure.subplot.left": 0.18,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.18,
    "figure.subplot.top": 0.95,
    "figure.subplot.wspace": 0.25,
    "figure.subplot.hspace": 0.25,
    "lines.markersize": 2,
    "lines.linewidth": 1.,
    "figure.max_open_warning": 0,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    grid(True)

    i = 0
    color=['tab:blue', 'tab:red', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'black']
    for dwarf in DSphs:

        file = "DSh_"+dwarf+"_fit.txt"
        r0_iso, high_r0_iso, low_r0_iso, rho0_iso, high_rho0_iso, low_rho0_iso,\
            r0_pse, high_r0_pse, low_r0_pse, rho0_pse, high_rho0_pse, low_rho0_pse= read_fits(file)

        yerr = np.array([(10**rho0_iso-10**low_rho0_iso, 10**high_rho0_iso-10**rho0_iso)]).T
        xerr = np.array([(r0_iso-low_r0_iso, high_r0_iso-r0_iso)]).T
        plt.errorbar([r0_iso], [10**rho0_iso], xerr=xerr, yerr=yerr, fmt='',color=color[i])
        plt.plot([r0_iso], [10**rho0_iso],'o',ms=4,color='white')
        plt.plot([r0_iso], [10**rho0_iso],'o',ms=3,color=color[i],label=dwarf)

        yerr = np.array([(10**rho0_pse-10**low_rho0_pse, 10**high_rho0_pse-10**rho0_pse)]).T
        xerr = np.array([(r0_pse-low_r0_pse, high_r0_pse-r0_pse)]).T
        plt.errorbar([r0_pse], [10**rho0_pse], xerr=xerr, yerr=yerr, fmt='',color=color[i])
        plt.plot([r0_pse], [10**rho0_pse],'v',ms=4,color='white')
        plt.plot([r0_pse], [10**rho0_pse],'v',ms=3,color=color[i])

        i+=1

    yscale('log')
    xscale('log')
    xlabel(r'r$_{0}$ [kpc]')
    ylabel(r'$\rho_{0}$ [M$_{\odot}$/kpc$^{3}$]')
    axis([1e-2, 1, 1e7, 5e9])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc="upper right", labelspacing=0.2, handlelength=1.5, handletextpad=0.4, frameon=False)
    plt.savefig("SIDM_plane.png", dpi=200)
    plt.close()

if __name__ == "__main__":

    #DSphs = ['Carina','Draco','Fornax','LeoI','Sculptor','Sextans','UMi']
    DSphs = ['UMi']

    # for dwarf in DSphs:
    #
    #     file = './data/DSph_'+dwarf+'.txt'
    #     x, y, yerr, errorbar = read_data(file)
    #     fit_profile(x, y, yerr, dwarf, "yes")


    DSphs = ['Carina','Draco','Fornax','LeoI','Sculptor','Sextans','UMi']
    plot_data(DSphs)