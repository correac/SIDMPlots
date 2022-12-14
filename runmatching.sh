#!/bin/bash

python3 matchsidm.py \
-d /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant00 /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant00 \
-s snapshot_0036.hdf5 snapshot_0036.hdf5 \
-c subhalo_0036.properties subhalo_0036.properties \
-n RefModel2L025N376SigmaConstant00 Test \
-t Hydro Hydro \
-o /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant00/output_data

#python3 $HOME/SIDMPlots/sidmdata.py \
#-d /projects/0/einf180/Tango_sims/L025N376/Hydro_Model_2/SigmaConstant00 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n RefModel2L025N376SigmaConstant00 \
#-t Hydro \
#-o $HOME/SIDMPlots/output_data/


#python3 analyse_high_ratio_haloes.py \
#-d /projects/0/einf180/Tango_sims/L025N752/DMONLY/SigmaVelDep20Isotropic  \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n DML025N752SigmaVelDep20 \
#-o ./