#!/bin/bash

# Runs SIDMPlots using the following variables:

#folder="/Users/camila/TangoSIDM/SIDMPlots_data/IsolatedHalo"
#output="/Users/camila/TangoSIDM/SIDMPlots_data/IsolatedHalo/plots/"
#snap="45"
#name="IsolatedHalo"

folder="/projects/0/einf180/Hernquist_halo/varible_hSIDM/L128/sigma_2423_notcheck/"
output="/home/ccorrea/SIDMPlots/plots/IsolatedHalo/"
snap="45"
name="IsolatedHalo"

python IsolatedHalo.py -d=$folder -s=$snap -n=$name -o=$output


