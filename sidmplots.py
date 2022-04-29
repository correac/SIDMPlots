"""
Description here
"""
from argumentparser import ArgumentParser
from plotter import html
from object import simulation_data
from time import time
from halo_data.density_profile import compute_density_profiles
from plotter.plot_profiles import plot_profiles
from plotter.scatter_rate import plot_cosmic_scatter_rate, compute_scatter_rate
from plotter.rotation_curve import plot_rotation_curve, make_rotation_curve_data
from plotter.rotation_curve import plot_rotation_curve_data
from plotter.rotation_curve import plot_rotation_relative_to_CDM
from plotter.HMF import make_HMF, plot_HMF
from plotter.plot_halo_data import store_halo_data, plot_relations
from plotter.core_expansion import plot_t0
from merger_tree.make_tree_data import make_tree_data

# from html import make_web, add_web_section, render_web, add_metadata_to_web
# from loadplots import loadPlots

def main(config: ArgumentParser):

    time_start = time()
    web = None
    profile_data_centrals_9 = None
    profile_data_satellites_9 = None
    profile_data_centrals_10 = None
    profile_data_satellites_10 = None
    profile_data_centrals_11 = None
    profile_data_satellites_11 = None
    cosmic_scatter_rate = None
    HMF_data = None
    halo_data = None
    output_name_list = []

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        # Fetch relevant input parameters from lists
        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]
        output = config.output_directory

        # Load all data and save it in SimInfo class
        sim_info = simulation_data.SimInfo(
            directory=directory,
            snapshot=snapshot,
            catalogue=catalogue,
            name=sim_name,
            output=output,
        )

        output_name_list.append(sim_info.simulation_name)

        # Make initial part of the webpage
        if sim == 0:
            web = html.make_web(sim_info.snapshot)
        elif web is not None:
            html.add_metadata_to_web(web, sim_info.snapshot)

        make_tree_data(sim_info)


        # profile_data_centrals_9 = compute_density_profiles(
        #     sim_info=sim_info,
        #     log10_min_mass=8.8,
        #     log10_max_mass=9.2,
        #     structure_type=10,
        #     profile_data=profile_data_centrals_9
        # )
        #
        # profile_data_satellites_9 = compute_density_profiles(
        #     sim_info=sim_info,
        #     log10_min_mass=8.8,
        #     log10_max_mass=9.2,
        #     structure_type=15,
        #     profile_data=profile_data_satellites_9
        # )
        #
        # profile_data_centrals_10 = compute_density_profiles(
        #     sim_info=sim_info,
        #     log10_min_mass=9.8,
        #     log10_max_mass=10.2,
        #     structure_type=10,
        #     profile_data=profile_data_centrals_10
        # )
        #
        # profile_data_satellites_10 = compute_density_profiles(
        #     sim_info=sim_info,
        #     log10_min_mass=9.8,
        #     log10_max_mass=10.2,
        #     structure_type=15,
        #     profile_data=profile_data_satellites_10
        # )
        #
        # profile_data_centrals_11 = compute_density_profiles(
        #     sim_info=sim_info,
        #     log10_min_mass=10.8,
        #     log10_max_mass=11.2,
        #     structure_type=10,
        #     profile_data=profile_data_centrals_11
        # )
        #
        # profile_data_satellites_11 = compute_density_profiles(
        #     sim_info=sim_info,
        #     log10_min_mass=10.8,
        #     log10_max_mass=11.2,
        #     structure_type=15,
        #     profile_data=profile_data_satellites_11
        # )
        #
        #cosmic_scatter_rate = compute_scatter_rate(sim_info, cosmic_scatter_rate)
        #
        # plot_rotation_curve(
        #     sim_info=sim_info,
        #     log10_min_mass=8.95,
        #     log10_max_mass=9.05,
        #     structure_type=10
        # )
        # plot_rotation_curve(
        #     sim_info=sim_info,
        #     log10_min_mass=8.95,
        #     log10_max_mass=9.05,
        #     structure_type=15
        # )
        #
        # plot_rotation_curve(
        #     sim_info=sim_info,
        #     log10_min_mass=9.95,
        #     log10_max_mass=10.05,
        #     structure_type=10
        # )
        #
        # plot_rotation_curve(
        #     sim_info=sim_info,
        #     log10_min_mass=9.95,
        #     log10_max_mass=10.05,
        #     structure_type=15
        # )
        #
        # HMF_data = make_HMF(sim_info=sim_info, HMF_data=HMF_data)
        #
        # halo_data = store_halo_data(sim_info=sim_info, halo_data=halo_data)

        #make_rotation_curve_data(sim_info=sim_info, log10_min_mass=9.0, log10_max_mass=12.0)
    
    #
    # plot_rotation_curve_data(sim_info, output_name_list=output_name_list)
    # plot_rotation_relative_to_CDM(sim_info, output_name_list=output_name_list)

    # plot_HMF(HMF_data, sim_info, output_name_list)
    #
    # plot_relations(halo_data, sim_info, output_name_list)

    # plot_t0(sim_info)

    #plot_cosmic_scatter_rate(cosmic_scatter_rate, sim_info, output_name_list)
    #
    # plot_profiles(profile_data=profile_data_centrals_9, sim_info=sim_info, output_name_list=output_name_list)
    # plot_profiles(profile_data=profile_data_satellites_9, sim_info=sim_info, output_name_list=output_name_list)
    #
    # plot_profiles(profile_data=profile_data_centrals_10, sim_info=sim_info, output_name_list=output_name_list)
    # plot_profiles(profile_data=profile_data_satellites_10, sim_info=sim_info, output_name_list=output_name_list)
    #
    # plot_profiles(profile_data=profile_data_centrals_11, sim_info=sim_info, output_name_list=output_name_list)
    # plot_profiles(profile_data=profile_data_satellites_11, sim_info=sim_info, output_name_list=output_name_list)

    # # After making individual plots finish up the website
    # # Load galaxy plots
    # loadPlots(web, siminfo)
    #
    # # Finish and output html file
    # render_web(web, siminfo.output_path)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)

