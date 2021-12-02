#!/bin/bash

# Runs SIDMPlots using the following variables:
python3 -i sidm_plane.py \
-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant01 \
~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant10 \
-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
-c subhalo_0036.properties subhalo_0036.properties \
-n DML025N752SigmaConstant01 DML025N752SigmaConstant10 \
-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

