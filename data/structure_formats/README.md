# Structure input formats

This directory contains the population maps used in the demographic
analyses and serves as the local working directory for genotype formats
used in population-structure and SNP-outlier analyses.

## Population maps

The following population maps are retained:

- `popmap_2pop.txt`
- `popmap_2pop_PZ.txt`
- `popmap_one_pop.txt`

## Genotype formats

PED and MAP files were exported from the filtered dataset using TASSEL.
They were converted to PLINK formats with:

`bin/2.2.population_structure/2.2.4_prepare_plink.sh`

The downstream scripts expect the following locally generated files:

- `qmacd_ref_gen_rob.plk.ped`
- `qmacd_ref_gen_rob.plk.map`
- `qmacd_ref_gen_rob.bed`
- `qmacd_ref_gen_rob.bim`
- `qmacd_ref_gen_rob.fam`
- `qmacd_ref_gen_rob.raw`

These derived files are retained in the local working directory but are
not distributed through GitHub.

The filtered VCF used in this study is available from Zenodo:

`https://doi.org/10.5281/zenodo.20548053`
