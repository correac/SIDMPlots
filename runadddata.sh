#!/bin/bash

# Runs SIDMPlots using the following variables:
python3 addsidmdata.py \
-d /projects/0/einf180/Tango_sims/L025N752/DMONLY/SigmaVelDep60Anisotropic  \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML025N752SigmaVelDep60Anisotropic \
-o ~/SIDMPlots/Simulation_data/ \
-f Tree_data_Satellites_DML025N752SigmaVelDep60Anisotropic

