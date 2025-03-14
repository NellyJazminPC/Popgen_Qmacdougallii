#!/bin/bash

# Define paths
ruta_bayescan=/home/n311pc/bioinfo/Popgen_Qmacdougallii/bin/software/BayeScan2.1/binaries
ruta_output=/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.6.snps_outliers/bayescan_output
ruta_input=/home/n311pc/bioinfo/Popgen_Qmacdougallii/data/1.6.snps_outliers

# Number of available processors
num_procesadores=10

# Create the output directory if it does not exist
mkdir -p $ruta_output

# Change to the BayeScan directory
cd $ruta_bayescan

# Check if the BayeScan executable exists
if [[ ! -f ./BayeScan2.1_linux64bits ]]; then
    echo "Error: BayeScan2.1_linux64bits not found in $ruta_bayescan"
    exit 1
fi

# Run BayeScan for the first input file
output_dir_2pop=$ruta_output/bayescan_qmacd_ref_gen_qrob_2pop
mkdir -p $output_dir_2pop
./BayeScan2.1_linux64bits $ruta_input/bayescan_qmacd_ref_gen_qrob_2pop -od $output_dir_2pop -threads $num_procesadores

# Run BayeScan for the second input file
output_dir_9pop=$ruta_output/bayescan_qmacd_ref_gen_qrob_9pop
mkdir -p $output_dir_9pop
./BayeScan2.1_linux64bits $ruta_input/bayescan_qmacd_ref_gen_qrob_9pop -od $output_dir_9pop -threads $num_procesadores

# List the resulting files
echo "Results for 2pop:"
ls -lh $output_dir_2pop
echo "Results for 9pop:"
ls -lh $output_dir_9pop