#!/bin/bash

# Runs SIDMPlots using the following variables:

#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant01/DMONLY /Users/Camila/SimulationData/cartesius/L025N376/SigmaVel45/DMONLY \
#-s 54 54 \
#-n DML025N376Sigma01 DML025N376SigmaVel45 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/Sigma_01_vel45_DMONLY_comparison/sidm_z0.5_plots

#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant01/DMONLY /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant10/DMONLY \
#-s 60 60 \
#-n DML025N376Sigma01 DML025N376Sigma10 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/Sigma_01_10_comparison
#
#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant00/DMONLY /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant10/DMONLY \
#-s 60 60 \
#-n DML025N376Sigma00 DML025N376Sigma10 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/Sigma_00_10_comparison
#
#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant00/DMONLY /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant01/DMONLY \
#-s 60 60 \
#-n DML025N376Sigma00 DML025N376Sigma01 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/Sigma_00_01_comparison
#
#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaVel45/DMONLY1node /Users/Camila/SimulationData/cartesius/L025N376/SigmaVel45/DMONLY \
#-s 55 55 \
#-n DML025N376SigmaVel451node DML025N376SigmaVel45 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/SigmaVel45_DMONLY_nodes_comparison
#

#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant01/Hydro /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant01/DMONLY \
#-s 54 54 \
#-n RefL025N376Sigma01 DML025N376Sigma01 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/Sigma_01_DMONLY_Hydro_comparison/sidm_z0.5_plots

#python sidmplots.py \
#-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaConstant01/Hydro /Users/Camila/SimulationData/cartesius/L025N376/SigmaVel45/Hydro \
#-s 54 54 \
#-n RefL025N376Sigma01 RefL025N376SigmaVel45 \
#-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/Sigma_01_vel45_Hydro_comparison/sidm_z0.5_plots

python sidmplots.py \
-d /Users/Camila/SimulationData/cartesius/L025N376/SigmaVel45/Hydro /Users/Camila/SimulationData/cartesius/L025N376/SigmaVel45/DMONLY \
-s 54 54 \
-n RefL025N376SigmaVel45 DML025N376SigmaVel45 \
-o /Users/camila/SimulationData/cartesius/L025N376/comparisons/SigmaVel45_DMONLY_Hydro_comparison/sidm_z0.5_plots

