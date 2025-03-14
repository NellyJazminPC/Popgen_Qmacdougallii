# Workflow to Run ADMIXTURE and fastStructure

## 1. Conversion of PLINK Files
First, convert the files in PLINK format (.ped and .map) to .bed, .bim, and .fam using the script `convert2plinkformat.sh`.

```sh
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
```

## 2. Running ADMIXTURE
Run ADMIXTURE with cross-validation for K=2 to K=10 using the script `admixture_qmacd.sh`.

```sh
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
```

## 3. Running fastStructure
Run fastStructure with two different priors (simple and logistic) for a range of K values using the script `faststructure_qmacd.sh`.

```sh
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
```

## 3.1 Using fastStructure with Docker

To use fastStructure with Docker, follow these steps:

1. Pull the Docker image for fastStructure:

    ```bash
    docker pull fischuu/faststructure
    ```

2. Run the Docker container:

    ```bash
    docker run -v /path/to/your/data:/data fischuu/faststructure 
    ```

    Replace `/path/to/your/data` with the actual path to your data directory.

3. To run fastStructure for multiple K values, you can execute the `faststructure_qmacd.sh` script within the Docker container:

    ```sh
    docker run -v /path/to/your/data:/data fischuu/faststructure /data/faststructure_qmacd.sh
    ```

## References

- PLINK
Purcell S, Neale B, Todd-Brown K, Thomas L, Ferreira MA, Bender D, Maller J, Sklar P, de Bakker PI, Daly MJ, Sham PC. PLINK: a tool set for whole-genome association and population-based linkage analyses. Am J Hum Genet. 2007 Sep;81(3):559-75. doi: 10.1086/519795. Epub 2007 Jul 25. PMID: 17701901; PMCID: PMC1950838.

- ADMIXTURE
Alexander DH, Novembre J, Lange K. Fast model-based estimation of ancestry in unrelated individuals. Genome Res. 2009 Sep;19(9):1655-64. doi: 10.1101/gr.094052.109. Epub 2009 Aug 4. PMID: 19648217; PMCID: PMC2752134.

- fastStructure
Raj A, Stephens M, Pritchard JK. fastSTRUCTURE: Variational Inference of Population Structure in Large SNP Data Sets. Genetics. 2014 Nov;197(2):573-89. doi: 10.1534/genetics.114.164350. Epub 2014 Sep 17. PMID: 25143593; PMCID: PMC4231593.

- For more information, consult the [official ADMIXTURE documentation](https://dalexander.github.io/admixture/) and the [official fastStructure documentation](https://rajanil.github.io/fastStructure/).