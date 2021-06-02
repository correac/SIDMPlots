#!/bin/bash

# Runs SIDMPlots using the following variables:

python sidmplots.py \
-d /projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N376/L025N376_sigma_1/Hydro /projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N376/L025N376_sigma_1/DM-only/data \
-s 51 51 \
-n Sigma01Hydro Sigma01DMONLY \
-o /projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N376/L025N376_sigma_1/Comparison

#python sidmplots.py \
#-d ~/SimulationData/cartesius/L025N376/SigmaVel45 \
#-s 60 \
#-n SigmaVel45 \
#-o ~/SimulationData/cartesius/L025N376/SigmaVel45/plots

#folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma15_w090"
#output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma15/"
#name="L025N256"
#snap="60"

#python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
#python extract_individual_profiles.py -d=$folder -s=$snap -n=$name -o=$output
#python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output


