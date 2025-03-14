# Workflow - Stacks populations

The commands to run population analyses with Stacks were as follows:

### For the _de novo_ assembly:

```sh
populations -V ../../data/1.3.assembly_variant_calling/denovo_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/denovo_9pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst_correction p_value --p_value_cutoff 0.05 -k --genepop --structure --plink
```

### For the reference genome - Quercus lobata (ref_gen_qlob) assembly:

```sh
populations -V ../../data/1.3.assembly_variant_calling/ref_gen_qlob_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/ref_gen_qlob_9pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst_correction p_value --p_value_cutoff 0.05 -k --genepop --structure --plink
```

### For the reference genome - Quercus robur (ref_gen_qrob) assembly:

```sh
populations -V ../../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/ref_gen_qrob_9pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst_correction p_value --p_value_cutoff 0.05 -k --genepop --structure --plink
```

The results were similar but we continued with the **ref_gen_qrob** file.

## UPDATE - Corrections in the code line - February 2025

changed `--fst_correction p_value --p_value_cutoff` to `--fst-correction --p-value-cutoff`.

```sh
populations -V ../../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf -O ../../data/1.4.population_stacks/ref_gen_qrob_2pop_2zones -M popmap.txt -t 8 --min-maf 0.02 --hwe --fstats --fst-correction --p-value-cutoff 0.05 -k --genepop --structure --plink
```


The analysis was performed for the three situations:

- Assuming all individuals are from a single population.

- Assuming there are two populations, one from the northern zone and one from the southern zone.

- Assuming each site is a population (9 populations) and they are in two zones (north and south)

## References

- Catchen, J., Hohenlohe, P. A., Bassham, S., Amores, A., & Cresko, W. A. (2013). Stacks: an analysis tool set for population genomics. Molecular Ecology, 22(11), 3124-3140. doi:10.1111/mec.12354

- Rochette, N. C., Rivera-Colón, A. G., & Catchen, J. M. (2019). Stacks 2: Analytical methods for paired-end sequencing improve RADseq-based population genomics. Molecular Ecology, 28(21), 4737-4754. doi:10.1111/mec.15253