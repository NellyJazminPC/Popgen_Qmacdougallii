# Ejecución del programa ADMIXTURE 1.3.0 en Ubuntu 24.04

## Requisitos previos

1. Instalación de [ADMIXTURE 1.3.0.](https://dalexander.github.io/admixture/download.html) 
2. Conversión de los datos en formato PLINK (.bed, .bim, .fam) 

    2.1 Se convirtip a `.ped`y `.map` con TASSEL

    2.2 Con el programa PLINK se convirtio a `.bed`, `.bim` y `.fam`


    ## Conversión de archivos `.ped` y `.map` a `.bed`, `.bim` y `.fam`

    Para convertir archivos `.ped` y `.map` a `.bed`, `.bim` y `.fam` utilizando [PLINK](https://zzz.bwh.harvard.edu/plink/download.shtml), puedes usar el siguiente comando:

    ```bash
    plink --file <input_file> --make-bed --out <output_file>
    ```

    - `<input_file>`: Nombre del archivo de entrada sin la extensión.
    - `<output_file>`: Nombre del archivo de salida sin la extensión.

    ### [Script para ejecutar plink](../1.5.structure/convert2plinkformat.sh)


### EXTRA: script para hacer los archivos por sitio y por zona.

Aunque para ello primero hay que separar el vcf en TASSEL.

```
#!/bin/bash 
#############################
#### Plink ####

#Para convertir archivo en formato PLINK (.ped y .map) a .bed


for i in CR.10.plk LS_04.plk MT_10.plk PZ.15.plk CY_10.plk MB_10.plk north.50.plk south.29.ind.plk CZ_10.plk MC_10.plk PZ_12.plk TZ_03.plk ; do


~/programs_bioinf/plink-1.07-x86_64/plink --file $i --recodeAD --out $i


~/programs_bioinf/plink-1.07-x86_64/plink --file $i --make-bed --out $i ;

done
```

## ADMIXTURE

## Comando básico

```bash
admixture <input_file>.bed <K>
```

- `<input_file>`: Nombre del archivo de entrada sin la extensión.
- `<K>`: Número de grupos ancestrales.

## Ejemplo de uso

```bash
admixture mydata.bed 3
```

Este comando ejecutará ADMIXTURE en el archivo `mydata.bed` asumiendo 3 grupos ancestrales.

## Opciones adicionales

- `-B`: Realiza bootstrapping.
- `--cv`: Realiza validación cruzada.

### Ejemplo con validación cruzada

```bash
admixture --cv mydata.bed 3
```

Este comando ejecutará ADMIXTURE con validación cruzada para 3 grupos ancestrales.

## Referencias

Para más información, consulta la [documentación oficial de ADMIXTURE](https://dalexander.github.io/admixture/).

### [Script para ejecutar ADMIXTURE](../1.5.structure/admixture_qmacd.sh)


```sh
#!/bin/bash

# Este script ejecuta ADMIXTURE con validación cruzada para K=2 a K=10 y guarda los resultados en archivos de salida.

# Definir rutas
ruta_admixture=~/bioinfo/Popgen_Qmacdougallii/bin/software/admixture_linux-1.3.0
ruta_bed=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats
ruta_output=~/bioinfo/Popgen_Qmacdougallii/data/1.5.structure

# Número de procesadores y semilla aleatoria
num_procesadores=8
semilla=12345

# Cambiar al directorio de ADMIXTURE
cd $ruta_admixture

# Ejecutar ADMIXTURE con validación cruzada para K=2 a K=10
for K in 2 3 4 5 6 7 8 9 10; 
do ./admixture --cv -j$num_procesadores -s$semilla $ruta_bed/qmacd_ref_gen_rob.bed $K | tee $ruta_output/log${K}.out; done

# Extraer los errores de validación cruzada
grep CV $ruta_output/log*.out > $ruta_output/chooseK.txt

# Mover los archivos .P y .Q generados a la ruta de salida
mv $ruta_admixture/*.P $ruta_output/
mv $ruta_admixture/*.Q $ruta_output/

# Mostrar el archivo resultante
cat $ruta_output/chooseK.txt
```

## fastStructure

Instalación:

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


Script
```sh
#!/bin/bash
# Para faststructure, structure prior= simple

# Definir rutas
ruta_faststructure=~/bioinfo/Popgen_Qmacdougallii/bin/software/fastStructure
ruta_rel_bed=~/bioinfo/Popgen_Qmacdougallii/data/structure_formats
ruta_output=~/bioinfo/Popgen_Qmacdougallii/data/1.5.structure

# Número de procesadores
num_procesadores=10

# Activar el entorno virtual
source $ruta_faststructure/venv/bin/activate

# Cambiar al directorio de fastStructure
cd $ruta_faststructure

##### S I M P L E  #######

for i in {1..8}; 
do 
    python3 structure.py -K $i --input=$ruta_rel_bed/qmacd_ref_gen_rob --output=$ruta_output/qmacd_ref_gen_rob.simple --full --seed=20 --prior=simple --threads=$num_procesadores
done

python3 chooseK.py --input=$ruta_output/qmacd_ref_gen_rob.simple > $ruta_output/chooseK_qmacd_ref_gen_rob.simple.txt

cat $ruta_output/chooseK_qmacd_ref_gen_rob.simple.txt

###### L O G I S T I C #######################

for i in {1..8}; 
do 
    python3 structure.py -K $i --input=$ruta_rel_bed/qmacd_ref_gen_rob --output=$ruta_output/qmacd_ref_gen_rob.logistic --full --seed=20 --prior=logistic --threads=$num_procesadores
done

python3 chooseK.py --input=$ruta_output/qmacd_ref_gen_rob.logistic > $ruta_output/chooseK_qmacd_ref_gen_rob.logistic.txt

cat $ruta_output/chooseK_qmacd_ref_gen_rob.logistic.txt

```

