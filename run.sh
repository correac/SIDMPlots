#!/bin/bash

# Runs SIDMPlots using the following variables:

folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma10_w030/test_with_timestep_limiter/"
output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma10_test_with_limiter/"
name="L025N256"
snap="55"

folder="/Users/camila/TangoSIDM/SIDMPlots_data/Plots/L025N256_sigma10/data/"
output="/Users/camila/TangoSIDM/SIDMPlots_data/Plots/L025N256_sigma10/data/plots/"
snap="55"
name="L025N256"

#folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma15_w090"
#output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma15/"
#name="L025N256"
#snap="60"

#python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output


