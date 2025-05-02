#!/bin/bash 
# This script converts a file in PLINK format (.ped and .map) to .bed format.

#############################
#### Plink ####

# To convert a file in PLINK format (.ped and .map) to .bed

ruta_file=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats/qmacd_ref_gen_rob.plk
output_name=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats/qmacd_ref_gen_rob
ruta_plink=~/bioinfo/Popgen_Qmacdougallii/bin/software/plink-1.07-x86_64/plink

# The --noweb option is used to run PLINK without attempting to check for updates online
$ruta_plink --file $ruta_file --noweb --recodeAD --out $output_name

# Convert to .bed format
$ruta_plink --file $ruta_file --noweb --make-bed --out $output_name