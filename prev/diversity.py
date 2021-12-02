import numpy as np
import h5py
from scipy import interpolate
import scipy.stats as stat
from pylab import *
import matplotlib.pyplot as plt


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

def analyse_halo(mass, radius):
    # Define radial bins [log scale, kpc units]
    radial_bins = np.arange(-0.5, 2, 0.2)
    radial_bins = 10 ** radial_bins
    centers = bin_centers(radial_bins)  # kpc

    SumMasses, _, _ = stat.binned_statistic(x=radius, values=mass, statistic="sum", bins=radial_bins, )
    density = (SumMasses / bin_volumes(radial_bins))  # Msun/kpc^3
    nonzero = density == 0     # clean..
    density[nonzero] = 1
    return centers, np.log10(density)

def density_slope(x,y):
    dx = x[1:]-x[:-1]
    dy = 10**y[1:]-10**y[:-1]
    g = dy * x[:-1] / (dx * 10**y[:-1])
    mask = x[:-1] <= 5
    g5 = np.median(g[mask])
    mask = x[:-1] <= 10
    g10 = np.median(g[mask])
    mask = x[:-1] <= 30
    g30 = np.median(g[mask])
    return g5, g10, g30

def make_diversity_data(siminfo):

    with h5py.File(siminfo.snapshot, "r") as hf:
        mass = hf['PartType1/Masses'][:] * 1e10  # Msun
        pos = hf['PartType1/Coordinates'][:][:] * siminfo.a # Mpc

    with h5py.File(siminfo.halo_properties, "r") as properties_file:
        c200c = properties_file["cNFW_200crit"][:]
        m200c = properties_file["Mass_200crit"][:] * 1e10 # Msun
        R200c = properties_file["R_200crit"][:] * 1e3 # kpc
        xCoP = properties_file["Xcminpot"][:] # Mpc
        yCoP = properties_file["Ycminpot"][:] # Mpc
        zCoP = properties_file["Zcminpot"][:] # Mpc
        Type = properties_file["Structuretype"][:]
        ID = properties_file["ID"][:]

    c200c[c200c == 0] = 1
    m200c[m200c <= 0] = 1
    m200c = np.log10(m200c)
    CoP = np.zeros((len(xCoP), 3))
    CoP[:, 0] = xCoP
    CoP[:, 1] = yCoP
    CoP[:, 2] = zCoP
    sample = np.where((m200c >= 9) & (m200c < 11))[0]

    # Reduce
    CoP = CoP[sample,:]
    m200c = m200c[sample]
    c200c = c200c[sample]
    R200c = R200c[sample]
    ID = ID[sample]
    Type = Type[sample]

    Density1kpc = np.zeros(len(sample))
    Density1p5kpc = np.zeros(len(sample))
    Density2kpc = np.zeros(len(sample))
    Density2p5kpc = np.zeros(len(sample))
    Density3kpc = np.zeros(len(sample))
    Density3p5kpc = np.zeros(len(sample))
    Density5kpc = np.zeros(len(sample))
    Density10kpc = np.zeros(len(sample))
    Gamma5kpc = np.zeros(len(sample))
    Gamma10kpc = np.zeros(len(sample))
    Gamma30kpc = np.zeros(len(sample))

    for i in range(len(sample)):
        particles_pos = pos.copy()
        particles_pos -= CoP[i, :]  # centering
        particles_pos *= 1e3  # kpc

        radius = np.linalg.norm(particles_pos[:, :3], axis=1)  # computing distances to the CoP
        x, y = analyse_halo(mass, radius)  # radius [kpc], log10 density [Msun/kpc^3]
        f = interpolate.interp1d(x, y)

        Density1kpc[i] = f(1.0)
        Density1p5kpc[i] = f(1.5)
        Density2kpc[i] = f(2.0)
        Density2p5kpc[i] = f(2.5)
        Density3kpc[i] = f(3.0)
        Density3p5kpc[i] = f(3.5)
        Density5kpc[i] = f(5.0)
        Density10kpc[i] = f(10.0)
        Gamma5kpc[i], Gamma10kpc[i], Gamma30kpc[i] = density_slope(x, y)

    filename = f"{siminfo.output_path}/Diversity_data_" + siminfo.name + ".hdf5"
    data_file = h5py.File(filename, 'w')
    f = data_file.create_group('Data')
    MH = f.create_dataset('ID', data=ID)
    MH = f.create_dataset('StructureType', data=Type)
    MH = f.create_dataset('M200c', data=m200c)
    MH = f.create_dataset('R200c', data=R200c)
    MH = f.create_dataset('c200c', data=c200c)
    MH = f.create_dataset('DMDensity1kpc', data=Density1kpc)
    MH = f.create_dataset('DMDensity1p5kpc', data=Density1p5kpc)
    MH = f.create_dataset('DMDensity2kpc', data=Density2kpc)
    MH = f.create_dataset('DMDensity2p5kpc', data=Density2p5kpc)
    MH = f.create_dataset('DMDensity3kpc', data=Density3kpc)
    MH = f.create_dataset('DMDensity3p5kpc', data=Density3p5kpc)
    MH = f.create_dataset('DMDensity5kpc', data=Density5kpc)
    MH = f.create_dataset('DMDensity10kpc', data=Density10kpc)
    MH = f.create_dataset('DMGamma5kpc', data=Gamma5kpc)
    MH = f.create_dataset('DMGamma10kpc', data=Gamma10kpc)
    MH = f.create_dataset('DMGamma30kpc', data=Gamma30kpc)
    data_file.close()

def median_relations(xdata, ydata, xrange):

    nonan = np.where(np.isnan(ydata)==False)[0]
    x = xdata[nonan]
    y = ydata[nonan]

    xvalues = np.ones(len(xrange) - 1) * (-10)
    yvalues = np.zeros(len(xrange) - 1)
    yvalues_err_down = np.zeros(len(xrange) - 1)
    yvalues_err_up = np.zeros(len(xrange) - 1)

    perc = [16, 84]

    for i in range(0, len(xrange) - 2):
        mask = (x > xrange[i]) & (x < xrange[i + 1])
        if len(x[mask]) > 4:
            xvalues[i] = np.median(x[mask])
            yvalues[i] = np.median(y[mask])
            yvalues_err_down[i], yvalues_err_up[i] = np.transpose(np.percentile(y[mask], perc))

    mask = xvalues>-10
    xvalues = xvalues[mask]
    yvalues = yvalues[mask]
    yvalues_err_down = yvalues_err_down[mask]
    yvalues_err_up = yvalues_err_up[mask]

    return xvalues, yvalues, yvalues_err_down, yvalues_err_up


def plot_central_density(x, y1, y2, y3, siminfo):

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
        "lines.markersize": 2,
        "lines.linewidth": 2,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    xrange = np.arange(7,15,0.2)
    xrange = 10**xrange

    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(x, y1, xrange)
    #plt.plot(x, y1, 'o', color='tab:blue',mfc='none',alpha=0.5)
    plt.plot(xvalues, yvalues, '-', color='tab:blue',label='Density at 1 kpc',zorder=10)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:blue',zorder=10)

    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(x, y2, xrange)
    plt.plot(x, y2, 'o', color='tab:orange',mfc='none',alpha=0.2)
    plt.plot(xvalues, yvalues, '-', color='tab:orange',label='Density at 2 kpc',zorder=10)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:orange',zorder=10)

    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(x, y3, xrange)
    #plt.plot(x, y3, 'o', color='tab:green',mfc='none',alpha=0.5)
    plt.plot(xvalues, yvalues, '-', color='tab:green',label='Density at 3 kpc',zorder=10)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:green',zorder=10)

    plt.axis([1e9,1e11,5,8])
    plt.xlabel("M$_{200}$ [$M_\odot$]")
    plt.ylabel("Density [log$_{10}$ M$_{\odot}$/kpc$^{3}$]")
    plt.xscale('log')
    plt.legend()
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{siminfo.output_path}/Density_"+siminfo.name+".png", dpi=200)
    plt.close()

def plot_gamma(x, y5, y10, y30, siminfo):

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
        "lines.markersize": 2,
        "lines.linewidth": 2,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    xrange = np.arange(7,15,0.2)
    xrange = 10**xrange

    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(x, y5, xrange)
    plt.plot(x, y5, 'o', color='tab:blue',mfc='none',alpha=0.2)
    plt.plot(xvalues, yvalues, '-', color='tab:blue',label='Gamma within 5 kpc',zorder=10)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:blue',zorder=10)

    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(x, y10, xrange)
    #plt.plot(x, y10, 'o', color='tab:orange',mfc='none',alpha=0.5)
    plt.plot(xvalues, yvalues, '-', color='tab:orange',label='Gamma within 10 kpc',zorder=10)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:orange',zorder=10)

    xvalues, yvalues, yvalues_err_down, yvalues_err_up = median_relations(x, y30, xrange)
    #plt.plot(x, y30, 'o', color='tab:green',mfc='none',alpha=0.5)
    plt.plot(xvalues, yvalues, '-', color='tab:green',label='Gamma within 30 kpc',zorder=10)
    plt.fill_between(xvalues, yvalues_err_down, yvalues_err_up, alpha=0.2, color='tab:green',zorder=10)
    plt.axis([1e9,1e11,-2,0.5])

    plt.xlabel("M$_{200}$ [$M_\odot$]")
    plt.ylabel("$\Gamma$ (Log. density slope within 5kpc)")
    plt.xscale('log')
    plt.legend()
    ax.tick_params(direction='in', axis='both', which='both', pad=4.5)
    plt.savefig(f"{siminfo.output_path}/DensitySlope_"+siminfo.name+".png", dpi=200)
    plt.close()

def plot_diversity(siminfo):

    filename = f"{siminfo.output_path}/Diversity_data_" + siminfo.name + ".hdf5"
    with h5py.File(filename, 'r') as f:
        mass = f['Data/M200c'][:]
        density1kpc = f['Data/DMDensity1kpc'][:]
        density2kpc = f['Data/DMDensity2kpc'][:]
        density3kpc = f['Data/DMDensity3kpc'][:]
        gamma5kpc = f['Data/DMGamma5kpc'][:]
        gamma10kpc = f['Data/DMGamma10kpc'][:]
        gamma30kpc = f['Data/DMGamma30kpc'][:]


    plot_gamma(10**mass, gamma5kpc, gamma10kpc, gamma30kpc, siminfo)
    plot_central_density(10**mass, density1kpc, density2kpc, density3kpc, siminfo)

def plot_diversity_comparison(siminfo, name_list):

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
        "lines.markersize": 2,
        "lines.linewidth": 2,
    }
    rcParams.update(params)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    color = np.array(['tab:orange','tab:blue'])
    xrange = np.arange(9, 11, 0.2)

    k = 0
    for name in name_list:

        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            density3kpc = f['Data/DMDensity3kpc'][:]
            type = f['Data/StructureType'][:]

        centrals = type == 10
        plt.plot(mass[centrals],density3kpc[centrals],'o', color=color[k],label=name,alpha=0.1)
        x, y, _, _ = median_relations(mass[centrals],density3kpc[centrals], xrange)
        plt.plot(x,y,'-',color=color[k],zorder=10)

        k += 1

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("DM Density at 3kpc")
    plt.legend()
    plt.axis([9,11,5,8])
    plt.savefig(f"{siminfo.output_path}/M200_DMDensity3kpc_centrals_comparison.png",dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            density3kpc = f['Data/DMDensity3kpc'][:]
            type = f['Data/StructureType'][:]

        satellites = type > 10
        plt.plot(mass[satellites], density3kpc[satellites],'o',  color=color[k], label=name,alpha=0.1)
        x, y, _, _ = median_relations(mass[satellites], density3kpc[satellites], xrange)
        plt.plot(x,y,'-',color=color[k],zorder=10)

        k += 1

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("DM Density at 3kpc")
    plt.legend()
    plt.axis([9,11,5,8])
    plt.savefig(f"{siminfo.output_path}/M200_DMDensity3kpc_satellites_comparison.png", dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    color = np.array(['tab:orange','tab:blue'])
    xrange = np.arange(9, 11, 0.2)

    k = 0
    for name in name_list:

        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            density1kpc = f['Data/DMDensity1kpc'][:]
            type = f['Data/StructureType'][:]

        centrals = type == 10
        plt.plot(mass[centrals],density1kpc[centrals],'o', color=color[k],label=name,alpha=0.1)
        x, y, _, _ = median_relations(mass[centrals],density1kpc[centrals], xrange)
        plt.plot(x,y,'-',color=color[k],zorder=10)

        k += 1

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("DM Density at 1kpc")
    plt.legend()
    plt.axis([9,11,5,8])
    plt.savefig(f"{siminfo.output_path}/M200_DMDensity1kpc_centrals_comparison.png",dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            density1kpc = f['Data/DMDensity1kpc'][:]
            type = f['Data/StructureType'][:]

        satellites = type > 10
        plt.plot(mass[satellites], density1kpc[satellites],'o',  color=color[k], label=name,alpha=0.1)
        x, y, _, _ = median_relations(mass[satellites], density1kpc[satellites], xrange)
        plt.plot(x,y,'-',color=color[k],zorder=10)

        k += 1

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("DM Density at 1kpc")
    plt.legend()
    plt.axis([9,11,5,8])
    plt.savefig(f"{siminfo.output_path}/M200_DMDensity1kpc_satellites_comparison.png", dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            gamma5kpc = f['Data/DMGamma5kpc'][:]
            type = f['Data/StructureType'][:]

        satellites = type > 10
        plt.plot(mass[satellites], gamma5kpc[satellites], 'o', color=color[k], label=name,alpha=0.1)
        x, y, _, _ = median_relations(mass[satellites], gamma5kpc[satellites], xrange)
        plt.plot(x,y,'-',color=color[k],zorder=10)

        k += 1
    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Density slope within 5kpc")
    plt.axis([9,11,-1.5,1])
    plt.legend()
    plt.savefig(f"{siminfo.output_path}/M200_gamma5kpc_satellites_comparison.png", dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            gamma5kpc = f['Data/DMGamma5kpc'][:]
            type = f['Data/StructureType'][:]

        centrals = type == 10
        plt.plot(mass[centrals], gamma5kpc[centrals], 'o', color=color[k], label=name,alpha=0.1)
        x, y, _, _ = median_relations(mass[centrals], gamma5kpc[centrals], xrange)
        plt.plot(x,y,'-',color=color[k],zorder=10)
        k += 1

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Density slope within 5kpc")
    plt.axis([9,11,-1.5,1])
    plt.legend()
    plt.savefig(f"{siminfo.output_path}/M200_gamma5kpc_centrals_comparison.png", dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            gamma5kpc = f['Data/DMGamma30kpc'][:]
            type = f['Data/StructureType'][:]

        satellites = type > 10
        plt.plot(mass[satellites], gamma5kpc[satellites], 'o', color=color[k], label=name, alpha=0.1)
        x, y, _, _ = median_relations(mass[satellites], gamma5kpc[satellites], xrange)
        plt.plot(x, y, '-', color=color[k], zorder=10)

        k += 1
    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Density slope within 30kpc")
    plt.axis([9, 11, -1.5, 1])
    plt.legend()
    plt.savefig(f"{siminfo.output_path}/M200_gamma30kpc_satellites_comparison.png", dpi=200)

    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    k = 0
    for name in name_list:
        filename = f"{siminfo.output_path}/Diversity_data_" + name + ".hdf5"
        with h5py.File(filename, 'r') as f:
            mass = f['Data/M200c'][:]
            gamma5kpc = f['Data/DMGamma30kpc'][:]
            type = f['Data/StructureType'][:]

        centrals = type == 10
        plt.plot(mass[centrals], gamma5kpc[centrals], 'o', color=color[k], label=name, alpha=0.1)
        x, y, _, _ = median_relations(mass[centrals], gamma5kpc[centrals], xrange)
        plt.plot(x, y, '-', color=color[k], zorder=10)
        k += 1

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Density slope within 30kpc")
    plt.axis([9, 11, -1.5, 1])
    plt.legend()
    plt.savefig(f"{siminfo.output_path}/M200_gamma30kpc_centrals_comparison.png", dpi=200)


def plot_diversity_with_stellar_mass(siminfo):

    filename = f"{siminfo.output_path}/Diversity_data_" + siminfo.name + ".hdf5"
    with h5py.File(filename, 'r') as f:
        gamma5kpc = f['Data/DMGamma5kpc'][:]
        gamma30kpc = f['Data/DMGamma30kpc'][:]
        density1kpc = f['Data/DMDensity1kpc'][:]
        density3kpc = f['Data/DMDensity3kpc'][:]
        type = f['Data/StructureType'][:]
        ID = f['Data/ID'][:]

    properties_file = h5py.File(siminfo.halo_properties, "r")
    prop_mass = properties_file["M_star_30kpc"][:] * 1e10
    prop_ID = properties_file["ID"][:]

    _, indices_v, indices_p = np.intersect1d(ID, prop_ID, assume_unique=True, return_indices=True, )
    mass = prop_mass[indices_p]
    type = type[indices_v]
    gamma5kpc = gamma5kpc[indices_v]
    gamma30kpc = gamma30kpc[indices_v]
    density1kpc = density1kpc[indices_v]
    density3kpc = density3kpc[indices_v]

    select = mass > 1e7
    mass = np.log10(mass[select])
    type = type[select]
    gamma5kpc = gamma5kpc[select]
    gamma30kpc = gamma30kpc[select]
    density1kpc = density1kpc[select]
    density3kpc = density3kpc[select]

    centrals = type == 10
    satellites = type > 10
    xrange = np.arange(7, 10, 0.2)

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
        "lines.markersize": 2,
        "lines.linewidth": 2,
    }
    rcParams.update(params)

    ################
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(mass[centrals],density3kpc[centrals],'o', color='tab:red',label="centrals",alpha=0.1)
    plt.plot(mass[satellites],density3kpc[satellites],'o', color='tab:blue',label="satellites",alpha=0.1)
    x, y, _, _ = median_relations(mass[centrals],density3kpc[centrals], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:red')
    x, y, _, _ = median_relations(mass[satellites],density3kpc[satellites], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:blue')

    plt.xlabel("log$_{10}$M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("DM Density at 3kpc")
    plt.legend()
    plt.axis([7,10,5,8])
    plt.savefig(f"{siminfo.output_path}/Mstellar_DMDensity3kpc.png",dpi=200)

    ################
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(mass[centrals],density1kpc[centrals],'o', color='tab:red',label="centrals",alpha=0.1)
    plt.plot(mass[satellites],density1kpc[satellites],'o', color='tab:blue',label="satellites",alpha=0.1)
    x, y, _, _ = median_relations(mass[centrals],density1kpc[centrals], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:red')
    x, y, _, _ = median_relations(mass[satellites],density1kpc[satellites], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:blue')

    plt.xlabel("log$_{10}$M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("DM Density at 1kpc")
    plt.legend()
    plt.axis([7,10,5,8])
    plt.savefig(f"{siminfo.output_path}/Mstellar_DMDensity1kpc.png",dpi=200)


    ################
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(mass[centrals],gamma5kpc[centrals],'o', color='tab:red',label="centrals",alpha=0.1)
    plt.plot(mass[satellites],gamma5kpc[satellites],'o', color='tab:blue',label="satellites",alpha=0.1)
    x, y, _, _ = median_relations(mass[centrals],gamma5kpc[centrals], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:red')
    x, y, _, _ = median_relations(mass[satellites],gamma5kpc[satellites], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:blue')

    plt.xlabel("log$_{10}$M$_{200}$ [M$_{\odot}$]")
    plt.ylabel("Density slope within 5kpc")
    plt.axis([7,10,-1.5,1])
    plt.legend()
    plt.savefig(f"{siminfo.output_path}/Mstellar_gamma5kpc.png", dpi=200)


    ################
    figure()
    ax = plt.subplot(1, 1, 1)
    plt.grid("True")

    plt.plot(mass[centrals],gamma30kpc[centrals],'o', color='tab:red',label="centrals",alpha=0.1)
    plt.plot(mass[satellites],gamma30kpc[satellites],'o', color='tab:blue',label="satellites",alpha=0.1)
    x, y, _, _ = median_relations(mass[centrals],gamma30kpc[centrals], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:red')
    x, y, _, _ = median_relations(mass[satellites],gamma30kpc[satellites], xrange)
    plt.plot(x,y,'-',zorder=10,color='tab:blue')

    plt.xlabel("log$_{10}$M$_{*}$ [M$_{\odot}$]")
    plt.ylabel("Density slope within 30kpc")
    plt.axis([7, 10, -1.5, 1])
    plt.legend()
    plt.savefig(f"{siminfo.output_path}/Mstellar_gamma30kpc.png", dpi=200)