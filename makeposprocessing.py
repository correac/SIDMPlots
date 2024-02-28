"""
Description here
"""
from argumentparser import ArgumentParser
from object import simulation_data
from time import time
from halo_data.make_halo_data import make_halo_data

def main(config: ArgumentParser):

    time_start = time()

    # Fetch relevant input parameters from lists
    directory = config.directory_list[0]
    snapshot = config.snapshot_list[0]
    catalogue = config.catalogue_list[0]
    sim_name = config.name_list[0]
    sim_type = config.sim_type[0]
    output = config.output_directory

    # Load all data and save it in SimInfo class
    sim_info = simulation_data.SimInfo(
        directory=directory,
        snapshot=snapshot,
        catalogue=catalogue,
        name=sim_name,
        output=output,
        simtype=sim_type,
    )

    # make_processed_data(sim_info)
    make_halo_data(sim_info)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)