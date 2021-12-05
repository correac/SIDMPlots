import numpy as np
import h5py
import scipy.stats as stat
import os
from pylab import *
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.optimize import root
from tqdm import tqdm
from object import particle_data


class make_profile_params:

    def __init__(self,mass):
        self.log10_M200c = mass
        self.rho0_isothermal = []
        self.r0_isothermal = []
        self.n0_isothermal= []
        self.rho0_pse_isothermal = []
        self.r0_pse_isothermal = []
        self.n0_pse_isothermal= []
        self.v0 = []
        self.sigma0 = []

    def append_data(self, rho0_isothermal, r0_isothermal, n0_isothermal,
                    rho0_pse_isothermal, r0_pse_isothermal, n0_pse_isothermal,
                    v0, sigma0):
        self.rho0_isothermal = np.append(self.rho0_isothermal, rho0_isothermal)
        self.r0_isothermal = np.append(self.r0_isothermal, r0_isothermal)
        self.n0_isothermal = np.append(self.n0_isothermal, n0_isothermal)

        self.rho0_pse_isothermal = np.append(self.rho0_pse_isothermal, rho0_pse_isothermal)
        self.r0_pse_isothermal = np.append(self.r0_pse_isothermal, r0_pse_isothermal)
        self.n0_pse_isothermal = np.append(self.n0_pse_isothermal, n0_pse_isothermal)
        self.v0 = np.append(self.v0, v0)
        self.sigma0 = np.append(self.sigma0, sigma0)


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
    if a<=0 or np.max(x)>1e5 : return 0
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
    if a<=0 or np.max(x)>1e5 : return 0

    ydata = finterpolate(x)
    f = b + ydata
    return f


def fit_profile(x, y, velocity, sigma, halo_index, sim_info):

    nozero = y > 0
    x = x[nozero]
    y = y[nozero]
    velocity = velocity[nozero]
    sigma = sigma[nozero]

    # Let's define inner region as everything within 25 kpc
    inner = x <= 25

    # First fit profile based on Isothermal model
    popt, pcov = curve_fit(fit_isothermal_model, x[inner], np.log10(y[inner]), p0=[1, 10])
    r0_full_iso = popt[0]
    rho0_full_iso = popt[1]
    ns0_full_iso = 0

    # Let's also consider pseudo-Isothermal model
    popt, pcov = curve_fit(fit_pseudo_isothermal_model, x[inner], np.log10(y[inner]), p0=[1, 10, 0])
    r0_pse_iso = popt[0]
    rho0_pse_iso = popt[1]
    ns0_pse_iso = popt[2]

    # Median velocity dispersion value in core
    v0 = np.median(velocity[x <= 10])

    # Median velocity dispersion value in core
    sigma0 = np.median(sigma[x <= 10])

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

    plot(x[inner], np.log10(y[inner]), '-', color='tab:orange')

    xrange = np.arange(-1,np.log10(25),0.2)
    xrange = 10**xrange
    plot(xrange, fit_isothermal_model(xrange, r0_full_iso, rho0_full_iso),'--',
         color='tab:blue',label='Isothermal fit')
    plot(xrange, fit_pseudo_isothermal_model(xrange, r0_pse_iso, rho0_pse_iso, ns0_pse_iso),'--',
         color='tab:red',label='Pseudo-Isothermal fit')

    plt.plot([r0_full_iso], [rho0_full_iso], 'o', color='tab:blue')
    plt.plot([r0_pse_iso], [rho0_pse_iso], 'v', color='tab:red')

    plt.plot([r0_full_iso], [fit_isothermal_model(r0_full_iso, r0_full_iso, rho0_full_iso)], '*', color='tab:blue')
    plt.plot([r0_pse_iso], [fit_pseudo_isothermal_model(r0_pse_iso, r0_pse_iso, rho0_pse_iso, ns0_pse_iso)], '>', color='tab:red')

    xscale('log')
    xlabel(r'r [kpc]')
    ylabel(r'$\log_{10}\rho$ [M$_{\odot}$/kpc$^{3}$]')
    # axis([0, 50, 5, 10])
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.legend(loc='lower left')
    plt.savefig(f"{sim_info.output_path}/Test_fit/halo_%i_"%halo_index + sim_info.simulation_name + ".png", dpi=200)
    plt.close()

    return rho0_full_iso, r0_full_iso, ns0_full_iso, rho0_pse_iso, r0_pse_iso, ns0_pse_iso, v0, sigma0


def bin_volumes(radial_bins):
    """Returns the volumes of the bins. """

    single_vol = lambda x: (4.0 / 3.0) * np.pi * x ** 3
    outer = single_vol(radial_bins[1:])
    inner = single_vol(radial_bins[:-1])
    return outer - inner


def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)


def calculate_profiles(mass, pos, vel, sigma, radial_bins):

    # Radial coordinates [kpc units]
    r = np.sqrt(np.sum(pos ** 2, axis=1))

    SumMasses, _, _ = stat.binned_statistic(x=r, values=mass, statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3

    std_vel_x, _, _ = stat.binned_statistic(x=r, values=vel[:, 0], statistic="std", bins=radial_bins, )
    std_vel_y, _, _ = stat.binned_statistic(x=r, values=vel[:, 1], statistic="std", bins=radial_bins, )
    std_vel_z, _, _ = stat.binned_statistic(x=r, values=vel[:, 2], statistic="std", bins=radial_bins, )
    velocity = np.sqrt(std_vel_x ** 2 + std_vel_y ** 2 + std_vel_z ** 2) / np.sqrt(3.)
    velocity[np.where(np.isnan(velocity))[0]] = 0

    sigma_profile, _, _ = stat.binned_statistic(x=r, values=sigma, statistic="median", bins=radial_bins, )

    return density, velocity, sigma_profile


def calculate_profile_params(sim_info, log10_min_mass, log10_max_mass):

    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(-0.3, 3, 0.1)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    sample = np.where(
        (sim_info.halo_data.log10_halo_mass >= log10_min_mass) &
        (sim_info.halo_data.log10_halo_mass <= log10_max_mass))[0]

    log10_M200 = sim_info.halo_data.log10_halo_mass[sample]
    type = sim_info.halo_data.structure_type[sample]

    num_halos = len(sample)
    profile_params = make_profile_params(log10_M200)

    for i in tqdm(range(num_halos)):

        halo_indx = sim_info.halo_data.halo_index[sample[i]]
        part_data = particle_data.load_particle_data(sim_info, halo_indx, sample[i])

        density, velocity, sigma = calculate_profiles(part_data.masses.value,
                                               part_data.coordinates.value,
                                               part_data.velocities.value,
                                               part_data.cross_section,
                                               radial_bins)


        rho0_fi, r0_fi, n0_fi, \
        rho0_ps, r0_ps, n0_ps, v0, sigma0 = fit_profile(centers, density, velocity, sigma, halo_indx, sim_info)

        profile_params.append_data(rho0_fi, r0_fi, n0_fi, rho0_ps, r0_ps, n0_ps, v0, sigma0)

    # Output data
    filename = f"{sim_info.output_path}/Test_fit/Profile_params_" + sim_info.simulation_name + ".hdf5"
    data_file = h5py.File(filename, 'w')
    f = data_file.create_group('Data')
    MH = f.create_dataset('ID', data=sim_info.halo_data.halo_index[sample])
    MH = f.create_dataset('StructureType', data=type)
    MH = f.create_dataset('M200c', data=log10_M200)
    MH = f.create_dataset('r0_isothermal_fit', data=profile_params.r0_isothermal)
    MH = f.create_dataset('rho0_isothermal_fit', data=profile_params.rho0_isothermal)
    MH = f.create_dataset('r0_pseudo_isothermal_fit', data=profile_params.r0_pse_isothermal)
    MH = f.create_dataset('rho0_pseudo_isothermal_fit', data=profile_params.rho0_pse_isothermal)
    MH = f.create_dataset('n0_pseudo_isothermal_fit', data=profile_params.n0_pse_isothermal)
    MH = f.create_dataset('velocity_dispersion_0', data=profile_params.v0)
    MH = f.create_dataset('cross_section_0', data=profile_params.sigma0)
    data_file.close()

    return