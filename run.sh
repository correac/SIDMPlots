#!/bin/bash

# Runs SIDMPlots using the following variables:
#python3 -i track_evolution.py \

python3 sidmplots.py \
-d ~/SimulationData/mahti/L006N188/DMONLY/SigmaConstant00 \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML006N188SigmaConstant00 \
-o ~/SimulationData/mahti/L006N188/DMONLY/plots/comparison/

python3 sidmplots.py \
-d ~/SimulationData/mahti/L006N188/DMONLY/SigmaVelDep100Anisotropic \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML006N188SigmaVelDep100Anisotropic \
-o ~/SimulationData/mahti/L006N188/DMONLY/plots/comparison/

python3 sidmplots.py \
-d ~/SimulationData/mahti/L006N188/DMONLY/SigmaVelDep100Isotropic \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n DML006N188SigmaVelDep100Isotropic \
-o ~/SimulationData/mahti/L006N188/DMONLY/plots/comparison/
#
#
#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant00 ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant01 \
#-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
#-c subhalo_0036.properties subhalo_0036.properties \
#-n DML025N752SigmaConstant00 DML025N752SigmaConstant01 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/
#
#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant00 ~/SimulationData/mahti/L025N752/DMONLY/SigmaVelDep20Isotropic \
#-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
#-c subhalo_0036.properties subhalo_0036.properties \
#-n DML025N752SigmaConstant00 DML025N752SigmaVelDep20Isotropic \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/
#
#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant00 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n DML025N752SigmaConstant00 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant01 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n DML025N752SigmaConstant01 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/
#
#
#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant10 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n DML025N752SigmaConstant10 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/
#
#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaVelDep20Isotropic \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n DML025N752SigmaVelDep20Isotropic \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant00 ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant10 \
#-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
#-c subhalo_0036.properties subhalo_0036.properties \
#-n DML025N752SigmaConstant00 DML025N752SigmaConstant10 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/

#python3 sidmplots.py \
#-d ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant01 ~/SimulationData/mahti/L025N752/DMONLY/SigmaConstant10 \
#-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
#-c subhalo_0036.properties subhalo_0036.properties \
#-n DML025N752SigmaConstant01 DML025N752SigmaConstant10 \
#-o ~/SimulationData/mahti/L025N752/DMONLY/plots/
