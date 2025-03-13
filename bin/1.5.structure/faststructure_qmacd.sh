#!/bin/bash
# Para faststructure, structure prior= simple

# Definir rutas
ruta_faststructure=/fastStructure-1.0
ruta_rel_bed=/workspace/data/structure_formats
ruta_output=/workspace/data/1.5.structure
ruta_final_output=/workspace/data/1.5.structure/faststructure_output

# Número de procesadores
num_procesadores=10

# Cambiar al directorio de fastStructure
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

# Crear la carpeta de salida final si no existe
mkdir -p $ruta_final_output

# Mover todos los archivos de salida a la carpeta de salida final
mv $ruta_output/*.simple* $ruta_final_output/
mv $ruta_output/*.logistic* $ruta_final_output/