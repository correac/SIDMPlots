"""
Description here
"""
from argumentparser import ArgumentParser
from object import simulation_data
from time import time
import numpy as np
from tqdm import tqdm
from object import particle_data
from pNbody import *
from astropy import units as u


def make_pnbody_data(sim_info):

    # Define units
    u_Length = 1 * u.kpc
    u_Mass = 10 ** 10 * u.M_sun
    u_Velocity = 1 * u.km / u.s
    u_Time = u_Length / u_Velocity
    toMsol = u_Mass.to(u.M_sun).value

    # Select galaxy sample
    sample = np.where((sim_info.halo_data.log10_halo_mass >= 12) & (sim_info.halo_data.log10_stellar_mass >= 8.5))[0]
    centrals = np.where(sim_info.halo_data.structure_type[sample] == 10)[0]
    sample = sample[centrals]
    num_sample = len(sample)

    halo_index = sim_info.halo_data.halo_index[sample]

    for i in tqdm(range(num_sample)):

        output_file = f"{sim_info.output_path}/"+sim_info.simulation_name+"_galaxy_%i.hdf5" % halo_index[i]

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        mass = part_data.stars.masses.value[:] # Msun
        mass /= 1e10 # mass in 10^10 solar mass
        pos = part_data.stars.coordinates.value[:, :] # kpc
        vel = part_data.stars.velocities.value[:, :]
        vel[:, 0] -= sim_info.halo_data.vxminpot[sample[i]]  # centering
        vel[:, 1] -= sim_info.halo_data.vyminpot[sample[i]]  # in km/s
        vel[:, 2] -= sim_info.halo_data.vzminpot[sample[i]]
        age = part_data.stars.age[:]
        metallicity = part_data.stars.smoothed_metal_mass_fractions[:]

        Z_floor = 0.02 * 10**(-20.)
        zero_correction = np.where(metallicity <= Z_floor)[0]
        metallicity[zero_correction] = Z_floor
        metallicity = np.log10(metallicity / 0.02) # Floor of -20


        # create the pNbody object
        nb = Nbody(status='new',pos=pos,vel=vel,mass=mass,ftype='swift')

        # set all particles to stellar particles
        nb.set_tpe(4)

        nb.age                =   age            # <<< age of the stellar particles in Gyr
        nb.mh                 =   metallicity    # <<< metallicity  [M/H]

        # add units
        nb.UnitLength_in_cm         = u_Length.to(u.cm).value
        nb.UnitMass_in_g            = u_Mass.to(u.g).value
        nb.UnitVelocity_in_cm_per_s = u_Velocity.to(u.cm/u.s).value
        nb.Unit_time_in_cgs         = u_Time.to(u.s).value

        nb.hubblefactorcorrection      = False
        nb.comovingtoproperconversion  = False
        nb.atime                       = 1

        nb.rename(output_file)
        nb.write()



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

    make_pnbody_data(sim_info)

    # Compute how much time it took to run the script
    time_end = time()
    script_runtime = time_end - time_start
    m, s = divmod(script_runtime, 60)
    h, m = divmod(m, 60)
    print(f"The script was run in {h:.0f} hours, {m:.0f} minutes, and {s:.0f} seconds")

if __name__ == "__main__":

    config_parameters = ArgumentParser()
    main(config_parameters)