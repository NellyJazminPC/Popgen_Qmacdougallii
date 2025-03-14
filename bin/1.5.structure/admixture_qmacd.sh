#!/bin/bash

# This script runs ADMIXTURE with cross-validation for K=2 to K=10 and saves the results in output files.

# Define paths
ruta_admixture=~/bioinfo/Popgen_Qmacdougallii/bin/software/admixture_linux-1.3.0
ruta_bed=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats
ruta_output=~/bioinfo/Popgen_Qmacdougallii/data/1.5.structure
ruta_final_output=~/bioinfo/Popgen_Qmacdougallii/data/1.5.structure/admixture_output

# Number of processors and random seed
num_procesadores=8
semilla=12345

# Change to the ADMIXTURE directory
cd $ruta_admixture

# Run ADMIXTURE with cross-validation for K=2 to K=10
for K in 1 2 3 4 5 6 7 8 9 10; 
do 
    ./admixture --cv -j$num_procesadores -s$semilla $ruta_bed/qmacd_ref_gen_rob.bed $K | tee $ruta_output/log${K}.out
done

# Extract cross-validation errors
grep CV $ruta_output/log*.out > $ruta_output/chooseK.txt

# Create the final output directory if it does not exist
mkdir -p $ruta_final_output

# Move the generated .P and .Q files to the final output directory
mv $ruta_admixture/*.P $ruta_final_output/
mv $ruta_admixture/*.Q $ruta_final_output/
mv $ruta_output/log*.out $ruta_final_output/
mv $ruta_output/chooseK.txt $ruta_final_output/

# Display the resulting file
cat $ruta_final_output/chooseK.txt