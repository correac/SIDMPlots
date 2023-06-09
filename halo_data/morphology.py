import numpy as np
from tqdm import tqdm
from object import particle_data
import scipy.stats as stat


def calculate_kappa_co(pos, vel, mass):
    # subhalo contain subhalo data and is strutured as follow
    # [ (0:3)CentreOfPotential[kpc]: (0)X | (1)Y | (2)Z  | (3:6)Velocity[km/s]: (3)Vx | (4)Vy | (5)Vz  | (6)R200c[kpc]]
    # partsDATA contains particles data and is structured as follow
    # [ (:3)Position[kpc]: (0)X | (1)Y | (2)Z  | (3)Mass[Msun] | (4:7)Velocity[km/s]: (4)Vx | (5)Vy | (6)Vz | (7)hsml]

    #particlesDATA = np.array(partsDATA).copy()  # isolating a copy

    # Centering onto subhalo CoP
    #particlesDATA[:, 0] -= halo_data.xminpot[halo_index]
    #particlesDATA[:, 1] -= halo_data.yminpot[halo_index]
    #particlesDATA[:, 2] -= halo_data.zminpot[halo_index]
    #particlesDATA[:, :3] += box_size / 2
    #particlesDATA[:, :3] %= box_size
    #particlesDATA[:, :3] -= box_size / 2  # end the unwrap

    # Center velocities on the subhalo CoM velocity
    #particlesDATA[:, 4] -= halo_data.vxminpot[halo_index]
    #particlesDATA[:, 5] -= halo_data.vyminpot[halo_index]
    #particlesDATA[:, 6] -= halo_data.vzminpot[halo_index]

    # Compute distances
    #distancesDATA = np.linalg.norm(particlesDATA[:, :3], axis=1)

    # Compute distances
    distancesDATA = np.sqrt(np.sum(pos ** 2, axis=1))

    # Restrict particles
    extract = distancesDATA < 30.0
    pos = pos[extract,:]
    vel = vel[extract,:]
    mass = mass[extract]

    #particlesDATA = particlesDATA[extract, :]
    distancesDATA = distancesDATA[extract]

    #Mstar = np.sum(particlesDATA[:, 3])  # compute total in-aperture stellar mass
    Mstar = np.sum(mass)

    # Compute 30kpc CoM to Sub CoM velocty offset & recenter
    #dvVmass = (
    #    np.sum(particlesDATA[:, 3][:, np.newaxis] * particlesDATA[:, 4:7], axis=0)
    #    / Mstar
    #)
    #particlesDATA[:, 4:7] -= dvVmass

    dvVmass = np.sum(mass[:, np.newaxis] * vel, axis=0) / Mstar
    vel -= dvVmass

    #particlesDATA[:, 4:7] -= dvVmass


    # Compute momentum
    smomentums = np.cross(pos, vel)
    momentum = np.sum(mass[:, np.newaxis] * smomentums, axis=0)

    #smomentums = np.cross(particlesDATA[:, :3], particlesDATA[:, 4:7])
    #momentum = np.sum(particlesDATA[:, 3][:, np.newaxis] * smomentums, axis=0)

    #extract = distancesDATA < 5.0
    #smomentum_inner_5kpc = np.cross(pos[extract,:] * vel[extract, :])
    #momentum_inner_5kpc = np.sum(mass[extract] * smomentum_inner_5kpc, axis=0)

    extract = distancesDATA < 5.0
    smomentum_inner_5kpc = np.cross(pos[extract,:], vel[extract, :])
    momentum_inner_5kpc = np.sum(mass[extract][:, np.newaxis] * smomentum_inner_5kpc, axis=0)

    #smomentum_inner_5kpc = np.cross(
    #    particlesDATA[extract, :3], particlesDATA[extract, 4:7]
    #)
    #momentum_inner_5kpc = np.sum(
    #    particlesDATA[extract, 3][:, np.newaxis] * smomentum_inner_5kpc, axis=0
    #)

    # Compute specific angular momentum
    sa_momentum = momentum / Mstar
    sa_momentum = np.linalg.norm(sa_momentum)

    # Compute rotational velocities
    smomentumz = np.sum(momentum * smomentums / np.linalg.norm(momentum), axis=1)
    cyldistances = ( distancesDATA ** 2- np.sum(momentum * pos / np.linalg.norm(momentum), axis=1)** 2 )

    # cyldistances = (
    #     distancesDATA ** 2
    #     - np.sum(momentum * particlesDATA[:, :3] / np.linalg.norm(momentum), axis=1)
    #     ** 2
    # )
    cyldistances = np.sqrt(np.abs(cyldistances))

    if len(cyldistances[cyldistances > 0]) > 0:
        cylmin = np.min(cyldistances[cyldistances > 0])
        cyldistances[cyldistances == 0] = cylmin
        vrots = smomentumz / cyldistances
    else:
        vrots = smomentumz

    # Compute kappa_co
    #Mvrot2 = np.sum((particlesDATA[:, 3] * vrots ** 2)[vrots > 0])
    #kappa_co = Mvrot2 / np.sum(
    #    particlesDATA[:, 3] * (np.linalg.norm(particlesDATA[:, 4:7], axis=1)) ** 2
    #)

    Mvrot2 = np.sum((mass * vrots ** 2)[vrots > 0])
    kappa_co = Mvrot2 / np.sum(mass * (np.linalg.norm(vel, axis=1)) ** 2 )

    # Apply rotation so that momentum vector corresponds to z-axis
    #momentum /= np.linalg.norm(momentum)
    #momentum_inner_5kpc /= np.linalg.norm(momentum_inner_5kpc)

    # Return
    return kappa_co, sa_momentum, momentum_inner_5kpc

def bin_centers(radial_bins):
    """Returns the centers of the bins. """

    outer = radial_bins[1:]
    inner = radial_bins[:-1]
    return 0.5 * (outer + inner)

def AxialRatios(pos, mass):
    """
    rs - CoM subtracted positions of *selected* particles in galactic units
    ms - *selected* particle masses in galactic units
    """

    radius = np.linalg.norm(pos[:, :3], axis=1)
    rs = pos[radius > 0, :].copy()
    ms = mass[radius > 0].copy()
    rs2 = radius[radius > 0] ** 2
    mst = np.sum(ms)

    # construct MoI tensor
    I_xx = ms * rs[:, 0] * rs[:, 0] / rs2
    I_xx = I_xx[np.isnan(I_xx) == False]  # remove nans
    I_xx = I_xx.sum()
    I_yy = ms * rs[:, 1] * rs[:, 1] / rs2
    I_yy = I_yy[np.isnan(I_yy) == False]
    I_yy = I_yy.sum()
    I_zz = ms * rs[:, 2] * rs[:, 2] / rs2
    I_zz = I_zz[np.isnan(I_zz) == False]
    I_zz = I_zz.sum()

    I_xy = ms * rs[:, 0] * rs[:, 1] / rs2
    I_xy = I_xy[np.isnan(I_xy) == False]
    I_xy = I_xy.sum()

    I_xz = ms * rs[:, 0] * rs[:, 2] / rs2
    I_xz = I_xz[np.isnan(I_xz) == False]
    I_xz = I_xz.sum()
    I_yz = ms * rs[:, 1] * rs[:, 2] / rs2
    I_yz = I_yz[np.isnan(I_yz) == False]
    I_yz = I_yz.sum()

    # I_xx = (
    #     rs2[:, [1, 2]].sum(axis=-1) / abs((rs2[:, [1, 2]].sum(axis=-1)) ** 0.5)
    # ) * ms
    # I_xx = I_xx[np.isnan(I_xx) == False]  # remove nans
    # I_xx = I_xx.sum()
    # I_yy = (
    #     rs2[:, [0, 2]].sum(axis=-1) / abs((rs2[:, [0, 2]].sum(axis=-1)) ** 0.5)
    # ) * ms
    # I_yy = I_yy[np.isnan(I_yy) == False]
    # I_yy = I_yy.sum()
    # I_zz = (
    #     rs2[:, [0, 1]].sum(axis=-1) / abs((rs2[:, [0, 1]].sum(axis=-1)) ** 0.5)
    # ) * ms
    # I_zz = I_zz[np.isnan(I_zz) == False]
    # I_zz = I_zz.sum()
    # I_xy = -((rs[:, 0] * rs[:, 1] / abs(rs[:, 0] * rs[:, 1]) ** 0.5) * ms)
    # I_xy = I_xy[np.isnan(I_xy) == False]
    # I_xy = I_xy.sum()
    # I_xz = -((rs[:, 0] * rs[:, 2] / abs(rs[:, 0] * rs[:, 2]) ** 0.5) * ms)
    # I_xz = I_xz[np.isnan(I_xz) == False]
    # I_xz = I_xz.sum()
    # I_yz = -((rs[:, 1] * rs[:, 2] / abs(rs[:, 1] * rs[:, 2]) ** 0.5) * ms)
    # I_yz = I_yz[np.isnan(I_yz) == False]
    # I_yz = I_yz.sum()
    I = np.array([[I_xx, I_xy, I_xz], [I_xy, I_yy, I_yz], [I_xz, I_yz, I_zz]])
    I /= mst

    # Get and order eigenvalues
    W, V = np.linalg.eig(I)
    W1, W2, W3 = np.sort(W)[::-1]

    # compute axes (unnormalised as we don't need absolute values)
    a = np.sqrt(W1)
    b = np.sqrt(W2)
    c = np.sqrt(W3)

    # a = np.sqrt(np.abs(W1 + W2 - W3))
    # b = np.sqrt(np.abs(W1 + W3 - W2))
    # c = np.sqrt(np.abs(W2 + W3 - W1))

    return a, b, c

def calculate_morphology(sim_info, sample):

    num_haloes = len(sample)

    radial_bins = np.arange(0, 3.1, 0.1)
    radial_bins = 10**radial_bins
    num_bins = len(radial_bins)

    DM_a_axis = np.zeros((num_bins, num_haloes))
    DM_b_axis = np.zeros((num_bins, num_haloes))
    DM_c_axis = np.zeros((num_bins, num_haloes))
    DMNparts = np.zeros((num_bins, num_haloes))
    cross_section = np.zeros((num_bins-1, num_haloes))

    kappa = np.zeros(num_haloes)
    ang_momentum = np.zeros((num_haloes, 3))
    smomentum = np.zeros((num_haloes, 3))

    Stars_a_axis = np.zeros((num_bins, num_haloes))
    Stars_b_axis = np.zeros((num_bins, num_haloes))
    Stars_c_axis = np.zeros((num_bins, num_haloes))
    StarsNparts = np.zeros((num_bins, num_haloes))

    Gas_a_axis = np.zeros((num_bins, num_haloes))
    Gas_b_axis = np.zeros((num_bins, num_haloes))
    Gas_c_axis = np.zeros((num_bins, num_haloes))
    GasNparts = np.zeros((num_bins, num_haloes))


    halo_index = sim_info.halo_data.halo_index[sample]

    for i in tqdm(range(num_haloes)):

        if halo_index[i] == -1: continue # no progenitor-found case

        part_data = particle_data.load_particle_data(sim_info, halo_index[i])

        bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.dark_matter.ids)

        if len(bound_particles_only) < 10: continue

        pos = part_data.dark_matter.coordinates.value[bound_particles_only, :]
        mass = part_data.dark_matter.masses.value[bound_particles_only]
        r = np.sqrt(np.sum(pos ** 2, axis=1))

        sigma = part_data.dark_matter.cross_section[bound_particles_only]
        cross_section[:, i], _, _ = stat.binned_statistic(x=r, values=sigma, statistic="median", bins=radial_bins, )

        for j in range(num_bins):

            select = np.where(r <= radial_bins[j])[0]

            DMNparts[j,i] = len(select)

            if len(select) < 10: continue

            DM_a_axis[j,i], DM_b_axis[j,i], DM_c_axis[j,i] = AxialRatios(pos[select, :], mass[select])

        if sim_info.simulation_type == 'Hydro':

            stars_bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.stars.ids)

            if (len(stars_bound_particles_only) < 10) : continue

            mass = part_data.stars.masses.value[stars_bound_particles_only]
            pos = part_data.stars.coordinates.value[stars_bound_particles_only, :]
            vel = part_data.stars.velocities.value[stars_bound_particles_only]
            vel[:, 0] -= sim_info.halo_data.vxminpot[sample[i]] # centering
            vel[:, 1] -= sim_info.halo_data.vyminpot[sample[i]]
            vel[:, 2] -= sim_info.halo_data.vzminpot[sample[i]]

            kappa[i], ang_momentum[i,:], smomentum[i,:] = calculate_kappa_co(pos, vel, mass)

            r = np.sqrt(np.sum(pos ** 2, axis=1))

            for j in range(num_bins):

                select = np.where(r <= radial_bins[j])[0]

                StarsNparts[j, i] = len(select)

                if len(select) < 10: continue

                Stars_a_axis[j, i], Stars_b_axis[j, i], Stars_c_axis[j, i] = AxialRatios(pos[select, :], mass[select])


            gas_bound_particles_only = part_data.select_bound_particles(sim_info, halo_index[i], part_data.gas.ids)

            if (len(gas_bound_particles_only) < 10) : continue

            mass = part_data.gas.masses.value[stars_bound_particles_only]
            pos = part_data.gas.coordinates.value[stars_bound_particles_only, :]
            r = np.sqrt(np.sum(pos ** 2, axis=1))

            for j in range(num_bins):

                select = np.where(r <= radial_bins[j])[0]

                GasNparts[j, i] = len(select)

                if len(select) < 10: continue

                Gas_a_axis[j, i], Gas_b_axis[j, i], Gas_c_axis[j, i] = AxialRatios(pos[select, :], mass[select])


    return {'DM_a_axis':DM_a_axis, 'DM_b_axis':DM_b_axis, 'DM_c_axis':DM_c_axis, 'DMNparts':DMNparts,
            'cross_section':cross_section, 'radial_bins':radial_bins, 'kappa':kappa,  'Lmomentum':ang_momentum, 'smomentum':smomentum,
            'Gas_a_axis': Gas_a_axis, 'Gas_b_axis': Gas_b_axis, 'Gas_c_axis': Gas_c_axis, 'GasNparts':GasNparts,
            'Stars_a_axis': Stars_a_axis, 'Stars_b_axis': Stars_b_axis, 'Stars_c_axis': Stars_c_axis, 'StarsNparts': StarsNparts}
