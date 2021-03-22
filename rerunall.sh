#!/bin/bash

# Runs SIDMPlots using the following variables:

#folder="/Users/camila/TangoSIDM/SIDMPlots_data/L025N064"
#output="/Users/camila/TangoSIDM/SIDMPlots_data/L025N064/plots/"
#snap="60"
#name="L025N064"

folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma10_w030"
output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma10/"
name="L025N256"
snap="60"

python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output

folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_vel/sigma20_w015"
output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma20/"
name="L025N256"
snap="60"

python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output

folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_10/"
output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma_cte10/"
name="L025N256"
snap="60"

python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output

folder="/projects/0/einf180/Cosmo_Volumes/L025boxSIDMTests/L025N256/L025N256_sigma_1/"
output="/home/ccorrea/SIDMPlots/plots/L025N256_sigma_cte01/"
name="L025N256"
snap="60"

python extract_profile_cosmo_box.py -d=$folder -s=$snap -n=$name -o=$output
python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output

