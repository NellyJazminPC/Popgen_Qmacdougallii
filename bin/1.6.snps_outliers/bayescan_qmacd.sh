#!/bin/bash

# Define paths
ruta_bayescan=~/bioinfo/Popgen_Qmacdougallii/bin/software/BayeScan2.1/binaries
ruta_output=~/bioinfo/Popgen_Qmacdougallii/data/1.6.snps_outliers/bayescan_output
ruta_input=~/bioinfo/Popgen_Qmacdougallii/data/1.6.snps_outliers

# Number of available processors
num_procesadores=10

# Create the output directory if it does not exist
mkdir -p $ruta_output

# Change to the BayeScan directory
cd $ruta_bayescan

# Run BayeScan for the first input file
./Bayescan2.1_linux64bits $ruta_input/bayescan_qmacd_ref_gen_qrob_2pop -od $ruta_output -threads $num_procesadores

# Run BayeScan for the second input file
./Bayescan2.1_linux64bits $ruta_input/bayescan_qmacd_ref_gen_qrob_9pop -od $ruta_output -threads $num_procesadores

# List the resulting files
ls -lh $ruta_output