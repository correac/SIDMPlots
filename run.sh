#!/bin/bash

# Runs SIDMPlots using the following variables:

#python sidmplots.py \
#-d ~/SimulationData/cartesius/L025N376/SigmaConstant01 ~/SimulationData/cartesius/L025N376/SigmaVel45 \
#-s 55 55 \
#-n Sigma01 SigmaVel45 \
#-o ~/SimulationData/cartesius/L025N376/comparisons/Sigma_01_vel45_comparison

python sidmplots.py \
-d ~/SimulationData/cartesius/L025N376/SigmaVel45 \
-s 60 \
-n SigmaVel45 \
-o ~/SimulationData/cartesius/L025N376/SigmaVel45/plots

#folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma15_w090"
#output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma15/"
#name="L025N256"
#snap="60"

#python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
#python extract_individual_profiles.py -d=$folder -s=$snap -n=$name -o=$output
#python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output


