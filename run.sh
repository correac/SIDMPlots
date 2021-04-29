#!/bin/bash

# Runs SIDMPlots using the following variables:

#folder="/Users/camila/SimulationData/cartesius/L025N256/SIDM/Sigma45VelDep/"
#output="/Users/camila/SimulationData/cartesius/L025N256/SIDM/Sigma45VelDep/plots/"
#snap="60"
#name="L025N256"

python sidmplots.py \
-d ~/SimulationData/cartesius/L025N376/SigmaConstant01 ~/SimulationData/cartesius/L025N376/SigmaConstant10  \
-s 60 60 \
-n Sigma01 Sigma10 \
-o ~/SimulationData/cartesius/L025N376/Sigma_01_10_comparison

#folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma15_w090"
#output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma15/"
#name="L025N256"
#snap="60"

#python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
#python extract_individual_profiles.py -d=$folder -s=$snap -n=$name -o=$output
#python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output


