Los comandos para ejecutar los análisis poblacionales con Stacks fueron los siguientes:

**Para el ensamble _de novo_**:

```sh
populations -V ../../data/1.3.assembly_variant_calling/denovo_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/denovo_9pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst_correction p_value --p_value_cutoff 0.05 -k --genepop --structure --plink
```

**Para el ensamble _reference genome - Quercus lobata_ (ref_gen_qlob)**:

```sh
populations -V ../../data/1.3.assembly_variant_calling/ref_gen_qlob_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/ref_gen_qlob_9pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst_correction p_value --p_value_cutoff 0.05 -k --genepop --structure --plink
```

**Para el ensamble _reference genome - Quercus robur_ (ref_gen_qrob)**:

```sh
populations -V ../../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/ref_gen_qrob_9pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst_correction p_value --p_value_cutoff 0.05 -k --genepop --structure --plink
```

Los resultados fueron similares pero seguimos con el archivo de **ref_gen_qrob**

