# Data

This directory contains input documentation, selected intermediate
files, and lightweight data products used in the population genomic
analyses of *Quercus macdougallii*.

The filtered VCF used in this study is available from Zenodo:

`https://doi.org/10.5281/zenodo.20548053`

Large sequencing files, reference genomes, complete assembly outputs,
and locally generated program-specific formats are not distributed
through GitHub.

## Directory overview

### `raw/`

Documents the original genotyping-by-sequencing reads for the 79
individuals included in the study.

### `1.0.quality_analysis/`

Documents the initial quality assessment of the raw sequencing reads
with FastQC.

### `1.1.filter/`

Documents the read-trimming step and the three trimmed-read datasets
evaluated during assembly.

### `1.2.post-filter_quality_analysis/`

Documents the FastQC assessment performed after read trimming.

### `1.3.assembly_variant_calling/`

Contains documentation and selected summaries from the ipyrad assembly
and variant-calling stage.

The assembly based on the *Quercus robur* reference genome and the
`trim01` read dataset was selected for downstream analyses.

### `1.4.population_structure/`

Contains retained ADMIXTURE and fastStructure outputs used to evaluate
population structure.

### `1.6.snp_outlier_detection/`

Contains BayeScan input and output files used in candidate SNP
detection.

### `1.7.sequence_annotation/`

Contains representative sequences and BLAST query sequences for the
candidate loci.

### `1.5.demography/`

Contains site-frequency-spectrum files and Stairway Plot outputs used
in the demographic-history analyses.

### `1.8.outlier_allele_frequencies/`

Contains the curated table of coding candidate SNPs used to generate
the allele dosage and frequency plots.

The table reports the reference (`REF`) and alternate (`ALT`) alleles
recorded in the VCF.

### `reference_genomes/`

Documents the *Q. robur* and *Quercus lobata* reference genomes used
during assembly and sequence annotation.

The genome files and indexes are retained locally because of their
size.

### `structure_formats/`

Contains the population maps used in demographic analyses and
documents the locally generated genotype formats used by downstream
scripts.

## Additional directories

Sample metadata are stored under `metadata/`.

Generated tables, figures, logs, and intermediate outputs are written
under `results/`.

Specific files and preparation steps are described in the README within
each analysis directory.
