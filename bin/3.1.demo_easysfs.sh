#!/bin/bash

# Crear la carpeta de salida
mkdir -p data/1.8.demography

# Explorar proyecciones posibles para cada popmap

# Para 9 poblaciones
python easySFS/easySFS.py -i data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf \
  -p data/structure_formats/popmap_9pop.txt -a

# Para 2 poblaciones (50 y 29)
python easySFS/easySFS.py -i data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf \
  -p data/structure_formats/popmap_2pop.txt -a

# Para 2 poblaciones (67 y 12)
python easySFS/easySFS.py -i data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf \
  -p data/structure_formats/popmap_2pop_PZ.txt -a

# Una vez seleccionadas las proyecciones óptimas, ejecuta los siguientes comandos reemplazando los valores de --proj

# Ejemplo para 9 poblaciones (reemplaza X1,X2,...,X9 por los valores elegidos)
# python easySFS/easySFS.py -i data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf \
#   -p data/structure_formats/popmap_9pop.txt \
#   -o data/1.8.demography/sfs_9pop -f --proj X1,X2,X3,X4,X5,X6,X7,X8,X9

# Ejemplo para 2 poblaciones (50 y 29)
# python easySFS/easySFS.py -i data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf \
#   -p data/structure_formats/popmap_2pop.txt \
#   -o data/1.8.demography/sfs_2pop -f --proj Y1,Y2

# Ejemplo para 2 poblaciones (67 y 12)
# python easySFS/easySFS.py -i data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf \
#   -p data/structure_formats/popmap_2pop_PZ.txt \
#   -o data/1.8.demography/sfs_2pop_PZ -f --proj Z1,Z2