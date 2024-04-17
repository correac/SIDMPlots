#!/bin/bash

#python3 matchsidm.py \
#-d /projects/0/einf180/Tango_sims/L025N376/Hydro_Model_2/SigmaConstant00/ \
#   /projects/0/einf180/Tango_sims/L025N376/Hydro_Model_2/SigmaVelDep30Anisotropic/ \
#-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
#-c subhalo_0036.properties subhalo_0036.properties \
#-n RefL025N376SigmaConstant00 RefL025N376SigmaVel30 \
#-t Hydro Hydro \
#-o /home/ccorrea/SIDMPlots/output_data/

#python3 matchsidm.py \
# -d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant00/ \
#    /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaVelDep30Anisotropic/ \
# -s snapshot_0036.hdf5 snapshot_0036.hdf5 \
# -c subhalo_0036.properties subhalo_0036.properties \
# -n ReferenceL025N376SigmaConstant00 ReferenceL025N376SigmaVel30 \
# -t Hydro Hydro \
# -o /Users/cc276407/Simulation_data/home/TangoSIDM_Images/

python3 matchsidm.py \
 -d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant00/ \
    /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant10/ \
 -s snapshot_0036.hdf5 snapshot_0036.hdf5 \
 -c subhalo_0036.properties subhalo_0036.properties \
 -n ReferenceL025N376SigmaConstant00 ReferenceL025N376SigmaConstant10 \
 -t Hydro Hydro \
 -o /Users/cc276407/Simulation_data/home/TangoSIDM_Images/
