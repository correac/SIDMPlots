#!/bin/bash

#python3 sidmdata.py \
#-d /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant00 \
#-s snapshot_0036.hdf5 \
#-c subhalo_0036.properties \
#-n RefModel2SigmaConstant00 \
#-t Hydro \
#-o /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant00/output_data

python3 sidmdata.py \
-d /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant10 \
-s snapshot_0036.hdf5 \
-c subhalo_0036.properties \
-n RefModel2SigmaConstant10 \
-t Hydro \
-o /Users/camila/SimulationData/mahti/L025N376/Hydro/SigmaConstant10/output_data
