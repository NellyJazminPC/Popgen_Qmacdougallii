# Ejecución del programa ADMIXTURE 1.3.0

## Requisitos previos

1. Instalación de [ADMIXTURE 1.3.0.](https://dalexander.github.io/admixture/download.html) 
2. Conversión de los datos en formato PLINK (.bed, .bim, .fam) 

    2.1 Se convirtip a `.ped`y `.map` con TASSEL

    2.2 Con el programa PLINK se convirtio a `.bed`, `.bim` y `.fam`


    ## Conversión de archivos `.ped` y `.map` a `.bed`, `.bim` y `.fam`

    Para convertir archivos `.ped` y `.map` a `.bed`, `.bim` y `.fam` utilizando PLINK, puedes usar el siguiente comando:

    ```bash
    plink --file <input_file> --make-bed --out <output_file>
    ```

    - `<input_file>`: Nombre del archivo de entrada sin la extensión.
    - `<output_file>`: Nombre del archivo de salida sin la extensión.

    ### Ejemplo de uso

    ```bash
    plink --file mydata --make-bed --out mydata_converted
    ```

    Este comando convertirá los archivos `mydata.ped` y `mydata.map` en `mydata_converted.bed`, `mydata_converted.bim` y `mydata_converted.fam`.

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

cv_error.sh

```sh
#!/bin/bash



for K in 1 2 3 4 5 6 7 8 9 10; 
do ./admixture --cv /media/nell/n311_pc/Quercus/Campos_Project_002_analisis/campos_002_BamHI-NsiI_4000000/output/Resources/admixture/var_filtro_4mill.bed  $K | tee log${K}.out; done
```

CV_error_per_sites_zones.sh

```sh
#!/bin/bash

### Script para hacer los análisis de admixture y faststructure simple y logistic


####################################################
######################## ADMIXTURE  ################

# -j número de procesadores
# -s random seed


ruta_bed=/home/nell/Bioinformatic/Qmacdougallii_genomics_and_environment/data/per_sites_and_zones

# PZ.15

for K in 1 2 3 4 5 6 7 8; 
do ./admixture --cv /home/nell/Bioinformatic/Qmacdougallii_genomics_and_environment/data/per_sites_and_zones/PZ.15.plk.bed  $K | tee log${K}.out; done

grep CV log*.out > PZ.15.chooseK.txt

cat PZ.15.chooseK.txt

# CR.10

for K in 1 2 3 4 5 6 7 8; 
do ./admixture --cv /home/nell/Bioinformatic/Qmacdougallii_genomics_and_environment/data/per_sites_and_zones/CR.10.plk.bed  $K | tee log${K}.out; done

grep CV log*.out > CR.10.chooseK.txt

cat CR.10.chooseK.txt
```

admixture_cv_error.sh

```sh
#!/bin/bash

ruta_abs_bed=/home/nell/Bioinformatic/Qmacdougallii_genomics_and_environment/data/var.79.inds.sorted.bed

# -j número de procesadores
# -s random seed
cd /home/nell/Bioinformatic/Qmacdougallii_genomics_and_environment/bin/admixture_linux-1.3.0

for K in 1 2 3 4 5 6 7 8 9 10; 
do ./admixture --cv $ruta_abs_bed  $K | tee log${K}.out; done

grep CV log*.out >chooseK.txt

cat chooseK.txt

```