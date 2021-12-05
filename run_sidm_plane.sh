#!/bin/bash

# Runs SIDMPlots using the following variables:
python3 sidm_plane.py \
-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaVelDep20Isotropic \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML025N752SigmaVelDep20Isotropic \
-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

python3 sidm_plane.py \
-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant01 \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML025N752SigmaConstant01 \
-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

python3 sidm_plane.py \
-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant10 \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML025N752SigmaConstant10 \
-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

#python3 sidm_plane.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant00 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n DML025N752SigmaConstant00 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

#python3 sidm_plane.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant01 \
# ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant10 \
# ~/SimulationData/mahti/L025N752/DMONLY/SigmaVelDep20Isotropic \
# ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant00 \
#-s snapshot_0036.hdf5 snapshot_0036.hdf5 snapshot_0036.hdf5 snapshot_0036.hdf5 \
#-c subhalo_0036.properties subhalo_0036.properties subhalo_0036.properties subhalo_0036.properties \
#-n DML025N752SigmaConstant01 DML025N752SigmaConstant10 DML025N752SigmaVelDep20Isotropic DML025N752SigmaConstant00 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

