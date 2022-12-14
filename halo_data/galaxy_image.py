from velociraptor import load as load_catalogue
import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib
from matplotlib import gridspec
import unyt

from swiftsimio import load, mask
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.projection import project_gas
from swiftsimio.visualisation.projection import project_pixel_grid
from swiftsimio.visualisation.slice import kernel_gamma
from swiftsimio.visualisation.smoothing_length_generation import (
    generate_smoothing_lengths,
)
from swiftsimio import swift_cosmology_to_astropy

from unyt import pc, kpc, msun, cm, yr
from unyt import proton_mass_cgs as mH
from unyt import boltzmann_constant_cgs as kB
from unyt import gravitational_constant_cgs as G
from unyt import unyt_array

from astropy.visualization import make_lupton_rgb

from scipy.optimize import curve_fit

# Define function for string formatting of scientific notation
def sci_notation(num, decimal_digits=1, precision=None, exponent=None):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext, with specified number of significant
    decimal digits and precision (number of decimal digits
    to show). The exponent to be used can also be specified
    explicitly.
    """
    try:
        if exponent is None:
            exponent = int(np.floor(np.log10(abs(num))))
        coeff = round(num / float(10 ** exponent), decimal_digits)
        if precision is None:
            precision = decimal_digits

        return r"${0:.{2}f}\times10^{{{1:d}}}$".format(coeff, exponent, precision)

    except:
        return "0"

def exponential(x, Sigma0, H, offset):
    return Sigma0 * np.exp(-np.abs(x + offset) / H)

def calculate_scaleheight_fit(mass_map, r_img_kpc):
    r_abs_max_kpc = 7.5
    xx = np.linspace(
        -r_img_kpc.value, r_img_kpc.value, len(mass_map[:, 0]), endpoint=True
    )
    z = (np.tile(xx, (len(xx), 1))).T
    z_1D = np.ravel(z[:, (np.abs(xx) < r_abs_max_kpc)])
    S_1D = np.ravel(mass_map[:, (np.abs(xx) < r_abs_max_kpc)])

    popt, pcov = curve_fit(
        exponential, z_1D[np.isfinite(S_1D)], S_1D[np.isfinite(S_1D)]
    )

    return popt[1]


def get_angular_momentum_vector(rparticles, vparticles, rgalaxy, vgalaxy, mparticles):
    # make sure everything is in the same unit system
    rparticles = rparticles.to_physical()
    vparticles = vparticles.to_physical()

    r = rparticles - rgalaxy
    v = vparticles - vgalaxy
    m = mparticles

    d = np.linalg.norm(r, axis=1)

    # mask out the innermost x% of star particles
    dmin = np.quantile(d, 0.3, axis=0)
    dmax = np.quantile(d, 0.5, axis=0)
    m[d < dmin] = 0.0
    m[d > dmax] = 0.0

    L = np.cross(r, v)
    L[:, 0] *= m
    L[:, 1] *= m
    L[:, 2] *= m

    Ltotal = np.sum(L, axis=0)

    Ltotal = Ltotal / np.linalg.norm(Ltotal)

    face_on_rotation_matrix = rotation_matrix_from_vector(Ltotal)
    edge_on_rotation_matrix = rotation_matrix_from_vector(Ltotal, axis="y")

    return face_on_rotation_matrix, edge_on_rotation_matrix

def make_rotation_matrix(ang_momentum):

    face_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum)
    edge_on_rotation_matrix = rotation_matrix_from_vector(ang_momentum, axis="y")
    return face_on_rotation_matrix, edge_on_rotation_matrix

def get_stars_surface_brightness_map(
    sim_info, halo_id, momentum, size, npix, r_img_kpc
):
    catalogue = load_catalogue(sim_info.halo_data.path_to_catalogue)

    # center of the halo
    x = catalogue.positions.xcmbp[halo_id]
    y = catalogue.positions.ycmbp[halo_id]
    z = catalogue.positions.zcmbp[halo_id]

    # angular momentum of the stars (for projection)
    # lx = catalogue.angular_momentum.lx_star[halo_id]
    # ly = catalogue.angular_momentum.ly_star[halo_id]
    # lz = catalogue.angular_momentum.lz_star[halo_id]
    #
    # angular_momentum_vector = np.array([lx.value, ly.value, lz.value])
    # angular_momentum_vector /= np.linalg.norm(angular_momentum_vector)
    angular_momentum_vector = momentum
    angular_momentum_vector /= np.linalg.norm(angular_momentum_vector)

    # needs to be in comoving coordinates for the mask
    region = [
        [
            x / catalogue.scale_factor - size / catalogue.scale_factor,
            x / catalogue.scale_factor + size / catalogue.scale_factor,
        ],
        [
            y / catalogue.scale_factor - size / catalogue.scale_factor,
            y / catalogue.scale_factor + size / catalogue.scale_factor,
        ],
        [
            z / catalogue.scale_factor - size / catalogue.scale_factor,
            z / catalogue.scale_factor + size / catalogue.scale_factor,
        ],
    ]

    visualise_region = [x - 0.5 * size, x + 0.5 * size, y - 0.5 * size, y + 0.5 * size]

    data_mask = mask(f"{sim_info.directory}/{sim_info.snapshot_name}")
    data_mask.constrain_spatial(region)
    data = load(f"{sim_info.directory}/{sim_info.snapshot_name}", mask=data_mask)
    data.stars.coordinates = data.stars.coordinates.to_physical()

    data.stars.smoothing_lengths = generate_smoothing_lengths(
        coordinates=data.stars.coordinates,
        boxsize=data.metadata.boxsize,
        kernel_gamma=kernel_gamma,
        neighbours=11,
        speedup_fac=1,
        dimension=3,
    )

    # face_on_rotation_matrix, edge_on_rotation_matrix = get_angular_momentum_vector(
    #     data.stars.coordinates,
    #     data.stars.velocities,
    #     [x, y, z],
    #     [
    #         catalogue.velocities.vxcmbp[halo_id],
    #         catalogue.velocities.vycmbp[halo_id],
    #         catalogue.velocities.vzcmbp[halo_id],
    #     ],
    #     data.stars.masses,
    # )

    face_on_rotation_matrix, edge_on_rotation_matrix = make_rotation_matrix(
        angular_momentum_vector
    )

    luminosities = [
        data.stars.luminosities.GAMA_i,
        data.stars.luminosities.GAMA_r,
        data.stars.luminosities.GAMA_g,
    ]
    rgb_image_face = np.zeros((npix, npix, len(luminosities)))

    for ilum in range(len(luminosities)):

        # Face on projection
        data.stars.usermass = luminosities[ilum]
        pixel_grid = project_pixel_grid(
            data.stars,
            resolution=int(npix),
            project="usermass",
            parallel=True,
            region=visualise_region,
            rotation_center=unyt.unyt_array([x, y, z]),
            rotation_matrix=face_on_rotation_matrix,
            boxsize=data.metadata.boxsize,
            backend="subsampled",
        )

        x_range = visualise_region[1] - visualise_region[0]
        y_range = visualise_region[3] - visualise_region[2]
        units = 1.0 / (x_range * y_range)
        # Unfortunately this is required to prevent us from {over,under}flowing
        # the units...
        units.convert_to_units(1.0 / (x_range.units * y_range.units))

        mass_map_face = unyt_array(pixel_grid, units=units)
        mass_map_face.convert_to_units(1.0 / pc ** 2)
        try:
            mass_map_face[mass_map_face == 0.0] = mass_map_face[
                mass_map_face > 0.0
            ].min()
        except:
            mass_map_face[mass_map_face == 0.0] = 1.0e-10

        rgb_image_face[:, :, ilum] = mass_map_face.T

    image_face = make_lupton_rgb(
        rgb_image_face[:, :, 0],
        rgb_image_face[:, :, 1],
        rgb_image_face[:, :, 2],
        Q=10,
        stretch=0.5,
    )
    # mask with circle
    lx, ly = mass_map_face.shape
    X, Y = np.ogrid[0:lx, 0:ly]
    mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
    #image_face[mask_circle, :] = 255

    H_kpc_gri = np.zeros(len(luminosities))
    rgb_image_edge = np.zeros((npix, npix, len(luminosities)))

    for ilum in range(len(luminosities)):
        # Face on projection
        data.stars.usermass = luminosities[ilum]
        pixel_grid = project_pixel_grid(
            data.stars,
            resolution=int(npix),
            project="usermass",
            parallel=True,
            region=visualise_region,
            rotation_center=unyt.unyt_array([x, y, z]),
            rotation_matrix=edge_on_rotation_matrix,
            boxsize=data.metadata.boxsize,
            backend="subsampled",
        )

        x_range = visualise_region[1] - visualise_region[0]
        y_range = visualise_region[3] - visualise_region[2]
        units = 1.0 / (x_range * y_range)
        # Unfortunately this is required to prevent us from {over,under}flowing
        # the units...
        units.convert_to_units(1.0 / (x_range.units * y_range.units))

        mass_map_edge = unyt_array(pixel_grid, units=units)
        mass_map_edge.convert_to_units(1.0 / pc ** 2)
        try:
            mass_map_edge[mass_map_edge == 0.0] = mass_map_edge[
                mass_map_edge > 0.0
            ].min()
        except:
            mass_map_edge[mass_map_edge == 0.0] = 1.0e-10

        try:
            H_kpc_gri[ilum] = calculate_scaleheight_fit(mass_map_edge.T, r_img_kpc)
        except:
            H_kpc_gri[ilum] = -1.0
        rgb_image_edge[:, :, ilum] = mass_map_edge.T

    print("H (gri): ", H_kpc_gri)

    image_edge = make_lupton_rgb(
        rgb_image_edge[:, :, 0],
        rgb_image_edge[:, :, 1],
        rgb_image_edge[:, :, 2],
        Q=10,
        stretch=0.5,
    )
    # mask with circle
    lx, ly = mass_map_edge.shape
    X, Y = np.ogrid[0:lx, 0:ly]
    mask_circle = (X - lx / 2) ** 2 + (Y - ly / 2) ** 2 > lx * ly / 4
    # image_edge[mask_circle, :] = 255

    return image_face, image_edge, visualise_region, x, y, -1.0, H_kpc_gri

def make_galaxy_images(sim_info, halo_index, halo_sample, momentum, kappa):

    #################################
    # Miscellaneous
    #################################

    nr_total_plots = 3

    npix = int(512 / 2)

    r_img_kpc = 30.0 * kpc
    lbar_kpc = 15.0 * kpc
    ypos_bar = 20.0 * kpc

    size = 2.0 * r_img_kpc

    ######################

    for ihalo, halo_id in enumerate(halo_index):

        ####
        # General information
        text = sim_info.simulation_name + ", z = %.1f" % (sim_info.z) + "\n\n"
        text += r"${\bf" + "VR\ halo\ id:\ \ \ %3.3i" % (halo_id) + r"}$" + "\n"
        text += (
                r"M$_{\mathrm{200,crit}}$ = "
                + sci_notation(10**sim_info.halo_data.log10_halo_mass[halo_id])
                + r" M$_{\odot}$"
                + "\n"
        )
        text += (
                r"M$_{\mathrm{*,30kpc}}$ = "
                + sci_notation(10**sim_info.halo_data.log10_stellar_mass[halo_id])
                + r" M$_{\odot}$"
                + "\n"
        )
        text += (
                r"M$_{\mathrm{gas,30kpc}}$ = "
                + sci_notation(10**sim_info.halo_data.log10_gas_mass[halo_id])
                + r" M$_{\odot}$"
                + "\n"
        )
        text += (
                r"$\kappa_{\mathrm{co}}$ = %.3f " % kappa[ihalo]
                + "\n"
        )
        print(sci_notation(10**sim_info.halo_data.log10_halo_mass[halo_id]))

        fig = plt.figure(figsize=(6.0, 3.5))
        fig.subplots_adjust(left=0.01, right=0.95, top=0.85, bottom=0.12)
        gs = gridspec.GridSpec(1, 3, wspace=0.0, hspace=0.0)
        #     3, 1, wspace=0.0, hspace=0.15, height_ratios=[0.05, 1.0]
        # )

        ax = plt.subplot(gs[0])
        ax.set_aspect("equal")
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(
            0.05, 0.95, text, ha="left", va="top", transform=ax.transAxes, fontsize=8
        )
        ax.text(
            0.05,
            1.20,
            "%i Galaxy" % (ihalo + 1),
            ha="left",
            va="bottom",
            transform=ax.transAxes,
            fontsize=14,
        )

        # Stars gri face-on
        ax = plt.subplot(gs[1])
        ax.set_title("Stars (gri) - face")
        (
            mass_map_face,
            mass_map_edge,
            visualise_region,
            x,
            y,
            totalmass,
            H_kpc_gri,
        ) = get_stars_surface_brightness_map(sim_info,
            halo_sample[ihalo], momentum[ihalo,:], size, npix, r_img_kpc
        )
        mass_map_face_plot = mass_map_face
        mass_map_edge_plot = mass_map_edge
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        im = ax.imshow(mass_map_face_plot, extent=visualise_region)
        # circle = plt.Circle(
        #     (x, y),
        #     (0.99 * r_img_kpc.value) / 1000.0,
        #     color="black",
        #     fill=False,
        #     linewidth=2,
        # )

        #ax.add_artist(circle)
        # ax.plot(
        #     [x - lbar_kpc / 2.0, x + lbar_kpc / 2.0],
        #     [y + ypos_bar, y + ypos_bar],
        #     color="white",
        #     linewidth=2,
        #     linestyle="solid",
        # )
        # ax.text(
        #     x,
        #     y + ypos_bar,
        #     "%i kpc" % (int(lbar_kpc.value)),
        #     color="white",
        #     verticalalignment="bottom",
        #     horizontalalignment="center",
        # )

        # Stars gri edge-on
        ax = plt.subplot(gs[2])
        ax.set_title("Stars (gri) - edge")
        ax.tick_params(labelleft=False, labelbottom=False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_axis_off()
        im = ax.imshow(mass_map_edge_plot, extent=visualise_region)
        # circle = plt.Circle(
        #     (x, y),
        #     (0.99 * r_img_kpc.value) / 1000.0,
        #     color="black",
        #     fill=False,
        #     linewidth=2,
        # )
        # ax.add_artist(circle)

        #ax.text(
        #    0.5,
        #    0.2,
        #    r"H$_{r}$ = %.2f kpc" % (H_kpc_gri[1]),
        #    ha="center",
        #    va="top",
        #    color="white",
        #    transform=ax.transAxes,
        #    fontsize=8,
        #)

        fig.savefig(
            f"{sim_info.output_path}" + "/surface_overview_halo%3.3i_" % (ihalo) + sim_info.simulation_name + ".png",
            dpi=300,
        )
        plt.close()