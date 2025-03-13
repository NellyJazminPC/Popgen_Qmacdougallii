#!/bin/bash

# Este script ejecuta ADMIXTURE con validación cruzada para K=2 a K=10 y guarda los resultados en archivos de salida.

# Definir rutas
ruta_admixture=~/bioinfo/Popgen_Qmacdougallii/bin/software/admixture_linux-1.3.0
ruta_bed=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats
ruta_output=~/bioinfo/Popgen_Qmacdougallii/data/1.5.structure
ruta_final_output=~/bioinfo/Popgen_Qmacdougallii/data/1.5.structure/admixture_output

# Número de procesadores y semilla aleatoria
num_procesadores=8
semilla=12345

# Cambiar al directorio de ADMIXTURE
cd $ruta_admixture

# Ejecutar ADMIXTURE con validación cruzada para K=2 a K=10
for K in 1 2 3 4 5 6 7 8 9 10; 
do 
    ./admixture --cv -j$num_procesadores -s$semilla $ruta_bed/qmacd_ref_gen_rob.bed $K | tee $ruta_output/log${K}.out
done

# Extraer los errores de validación cruzada
grep CV $ruta_output/log*.out > $ruta_output/chooseK.txt

# Crear la carpeta de salida final si no existe
mkdir -p $ruta_final_output

# Mover los archivos .P y .Q generados a la ruta de salida final
mv $ruta_admixture/*.P $ruta_final_output/
mv $ruta_admixture/*.Q $ruta_final_output/
mv $ruta_output/log*.out $ruta_final_output/
mv $ruta_output/chooseK.txt $ruta_final_output/

# Mostrar el archivo resultante
cat $ruta_final_output/chooseK.txt