#!/bin/bash

# Runs SIDMPlots using the following variables:

#folder="/Users/camila/TangoSIDM/SIDMPlots_data/L025N064"
#output="/Users/camila/TangoSIDM/SIDMPlots_data/L025N064/plots/"
#snap="60"
#name="L025N064"

folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma15_w090"
output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma15/"
name="L025N256"
snap="60"

python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
#python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output


