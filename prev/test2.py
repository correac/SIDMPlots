import numpy as np
import h5py
import matplotlib.pyplot as plt
from argumentparser import ArgumentParser
from time import time
from halo_data.morphology import calculate_morphology, calculate_kappa_co
from object import simulation_data
from object import particle_data

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

    if sim_info.simulation_type == 'Hydro':
        sample = np.where((sim_info.halo_data.log10_halo_mass >= 12) & (sim_info.halo_data.log10_stellar_mass >= 10.5))[
            0]
    else:
        sample = np.where(sim_info.halo_data.log10_halo_mass >= 10)[0]

    centrals = np.where(sim_info.halo_data.structure_type[sample] == 10)[0]
    sample = sample[centrals]

    halo_index = sim_info.halo_data.halo_index[sample]
    num_sample = len(sample)

    # # Morphology/shape estimations
    data = calculate_morphology(sim_info, sample)

    kappa = data['kappa']
    print(kappa)

    for i in range(num_sample):

        if kappa[i] > 0.4:
            part_data = particle_data.load_particle_data(sim_info, halo_index[i])

            bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.stars.ids)

            if len(bound_particles_only) < 10: continue

            pos = part_data.stars.coordinates.value[bound_particles_only, :]
            mass = part_data.stars.masses.value[bound_particles_only]

            vel = part_data.stars.velocities.value[bound_particles_only]
            vel[:, 0] -= sim_info.halo_data.vxminpot[sample[i]]  # centering
            vel[:, 1] -= sim_info.halo_data.vyminpot[sample[i]]
            vel[:, 2] -= sim_info.halo_data.vzminpot[sample[i]]

            _,_,_, projected_pos = calculate_kappa_co(pos, vel, mass)

            fig = plt.figure(figsize=(6.0, 3.0))
            fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
            ax = plt.subplot(1, 2, 1)
            ax.tick_params(labelleft=False, labelbottom=False)
            plt.plot(projected_pos[:,0], projected_pos[:,1], '.', ms=0.05)
            ax = plt.subplot(1, 2, 2)
            ax.tick_params(labelleft=False, labelbottom=False)
            plt.plot(projected_pos[:,0], projected_pos[:,2], '.', ms=0.05)
            fig.savefig(
                f"{sim_info.output_path}" + "/face_on_halo%3.3i_" % (i) + sim_info.simulation_name + ".png",
                dpi=300,
            )
            plt.close()

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)