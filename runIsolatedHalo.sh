#!/bin/bash

# Runs SIDMPlots using the following variables:

folder="/Users/camila/TangoSIDM/SIDMPlots_data/IsolatedHalo"
output="/Users/camila/TangoSIDM/SIDMPlots_data/IsolatedHalo/plots/"
snap="40"
name="IsolatedHalo"
 
python IsolatedHalo.py -d=$folder -s=$snap -n=$name -o=$output


