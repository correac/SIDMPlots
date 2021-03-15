#!/bin/bash

# Runs SIDMPlots using the following variables:

folder="/Users/camila/TangoSIDM/SIDMPlots_data"
output="/Users/camila/TangoSIDM/SIDMPlots_data/plots/"
snap="60"
name="L025N064"
 
python sidmplots.py -d=$folder -s=$snap -n=$name -o=$output


