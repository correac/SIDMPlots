"""
Description here
"""
from argumentparser import ArgumentWithInputFiles
from merger_tree.make_tree_data import add_tree_data
from object import simulation_data
from time import time

if __name__ == "__main__":

    time_start = time()

    config = ArgumentWithInputFiles()

    # Fetch relevant input parameters from lists
    directory = config.directory_list[0]
    snapshot = config.snapshot_list[0]
    catalogue = config.catalogue_list[0]
    sim_name = config.name_list[0]
    output = config.output_directory
    input_file = config.input_file_list

    # Load all data and save it in SimInfo class
    sim_info = simulation_data.SimInfo(
        directory=directory,
        snapshot=snapshot,
        catalogue=catalogue,
        name=sim_name,
        output=output,
    )

    add_tree_data(sim_info, input_file)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")
