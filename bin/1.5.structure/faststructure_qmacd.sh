#!/bin/bash
# This script runs fastStructure with two different priors (simple and logistic) for a range of K values (1 to 10).
# It then selects the best K value using chooseK.py and moves the output files to a final output directory.

# Define paths
ruta_faststructure=/fastStructure-1.0
ruta_rel_bed=/workspace/data/structure_formats
ruta_output=/workspace/data/1.5.structure
ruta_final_output=/workspace/data/1.5.structure/faststructure_output

# Number of processors
num_procesadores=10

# Change to the fastStructure directory
cd $ruta_faststructure

##### S I M P L E  #######

for i in {1..10}; 
do 
    python structure.py -K $i --input=$ruta_rel_bed/qmacd_ref_gen_rob --output=$ruta_output/qmacd_ref_gen_rob.simple --full --seed=20 --prior=simple --format=bed
done

python chooseK.py --input=$ruta_output/qmacd_ref_gen_rob.simple > $ruta_output/chooseK_qmacd_ref_gen_rob.simple.txt

cat $ruta_output/chooseK_qmacd_ref_gen_rob.simple.txt

###### L O G I S T I C #######################

for i in {1..10}; 
do 
    python structure.py -K $i --input=$ruta_rel_bed/qmacd_ref_gen_rob --output=$ruta_output/qmacd_ref_gen_rob.logistic --full --seed=20 --prior=logistic --format=bed
done

python chooseK.py --input=$ruta_output/qmacd_ref_gen_rob.logistic > $ruta_output/chooseK_qmacd_ref_gen_rob.logistic.txt

cat $ruta_output/chooseK_qmacd_ref_gen_rob.logistic.txt

# Create the final output directory if it does not exist
mkdir -p $ruta_final_output

# Move all output files to the final output directory
mv $ruta_output/*.simple* $ruta_final_output/
mv $ruta_output/*.logistic* $ruta_final_output/