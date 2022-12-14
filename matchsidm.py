"""
Description here
"""
from argumentparser import ArgumentParser
from object import simulation_data
from time import time
from halo_data.matching_halos import match_simulations

def main(config: ArgumentParser):

    time_start = time()

    # Loop over simulation list
    for sim in range(config.number_of_inputs):

        # Fetch relevant input parameters from lists
        directory = config.directory_list[sim]
        snapshot = config.snapshot_list[sim]
        catalogue = config.catalogue_list[sim]
        sim_name = config.name_list[sim]
        output = config.output_directory
        sim_type = config.sim_type[sim]

        # Load all data and save it in SimInfo class
        if sim == 0:
            cdm_sim_info = simulation_data.SimInfo(
                directory=directory,
                snapshot=snapshot,
                catalogue=catalogue,
                name=sim_name,
                output=output,
                simtype=sim_type,
            )
        if sim == 1:
            sidm1_sim_info = simulation_data.SimInfo(
                directory=directory,
                snapshot=snapshot,
                catalogue=catalogue,
                name=sim_name,
                output=output,
                simtype=sim_type,
            )
        # if sim == 2:
        #     sidm2_sim_info = simulation_data.SimInfo(
        #         directory=directory,
        #         snapshot=snapshot,
        #         catalogue=catalogue,
        #         name=sim_name,
        #         output=output,
        #         simtype=sim_type,
        #     )

    match_simulations(cdm_sim_info, sidm1_sim_info)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)