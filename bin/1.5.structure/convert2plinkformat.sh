#!/bin/bash 
#############################
#### Plink ####

#Para convertir archivo en formato PLINK (.ped y .map) a .bed

ruta_file=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats/qmacd_ref_gen_rob.plk
output_name=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats/qmacd_ref_gen_rob
ruta_plink=~/bioinfo/Popgen_Qmacdougallii/bin/software/plink-1.07-x86_64/plink

$ruta_plink --file $ruta_file --noweb --recodeAD --out $output_name


$ruta_plink --file $ruta_file --noweb --make-bed --out $output_name 
