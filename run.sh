#!/bin/bash

#python3 makeposprocessing.py \
#-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant00/ \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n HaloCatalogueL025N376ReferenceSigmaConstant00 \
#-t Hydro \
#-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant00/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaConstant00/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376WeakStellarFBSigmaConstant00 \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaConstant00/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaVelDep30Anisotropic/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376WeakStellarFBSigmaVelDep30Anisotropic \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaVelDep30Anisotropic/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaVelDep60Anisotropic/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376WeakStellarFBSigmaVelDep60Anisotropic \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaVelDep60Anisotropic/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaVelDep30Anisotropic/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376ReferenceSigmaVelDep30Anisotropic \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaVelDep30Anisotropic/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaVelDep60Anisotropic/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376ReferenceSigmaVelDep60Anisotropic \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaVelDep60Anisotropic/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant10/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376ReferenceSigmaConstant10 \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/Reference/SigmaConstant10/

python3 makeposprocessing.py \
-d /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaConstant10/ \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n HaloCatalogueL025N376WeakStellarFBSigmaConstant10 \
-t Hydro \
-o /Users/cc276407/Simulation_data/snellius/TangoSIDM/L025N376/WeakStellarFB/SigmaConstant10/
