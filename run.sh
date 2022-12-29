#!/bin/bash

python3 sidmdata.py \
-d /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant10 \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n RefModel2SigmaConstant10 \
-t Hydro \
-o /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant10/output_data

#python3 sidmdata.py \
#-d /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant10 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n RefModel2SigmaConstant10 \
#-t Hydro \
#-o /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant10/output_data

#python3 sidmdata.py \
#-d /projects/0/einf180/Tango_sims/L025N376/Hydro_Model_2/SigmaVelDep60Anisotropic/ \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n RefModel2SigmaVel60 \
#-t Hydro \
#-o /home/ccorrea/SIDMPlots/output_data/GalaxyImages/L025N376Model2SigmaVel60

#python3 sidmdata.py \
#-d /projects/0/einf180/Tango_sims/L025N376/Hydro_Model_2/SigmaVelDep30Anisotropic/ \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n RefModel2SigmaVel30 \
#-t Hydro \
#-o /home/ccorrea/SIDMPlots/output_data/GalaxyImages/L025N376Model2SigmaVel30
