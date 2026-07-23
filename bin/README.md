# Analysis pipeline

This directory contains the scripts, Jupyter notebooks, configuration files, and workflow documentation used for the population genomic analyses of *Quercus macdougallii*.

The workflow is organized into numbered stages that broadly follow the order in which the analyses were performed. Some numbering gaps remain because exploratory or intermediate stages were removed from the public repository or consolidated into later workflow sections.

## Running the analyses

The refactored scripts in `2.1.genetic_diversity/` and `2.2.population_structure/` use repository-relative paths and automatically locate the repository root. They can therefore be launched from the repository root or from one of its subdirectories.

Examples:

```bash
Rscript bin/2.2.population_structure/2.2.1_pca.R
bash bin/2.2.population_structure/2.2.5_run_admixture.sh
```

Other scripts retained from earlier stages may still use relative paths based on their original workflow. Consult the README associated with each stage before running those analyses.

The local `.Rproj` file used during development is not included because it may contain machine-specific settings. Generated outputs are written to `results/`, which is intentionally excluded from version control unless a lightweight derived result is explicitly retained for reproducibility.

## 1.0 Initial quality assessment

### `1.0.quality_analysis.sh`

Performs the initial quality assessment of the 79 raw single-end GBS read files using FastQC.

## 1.1 Read trimming

### `1.1.filter_trimmomatic.sh`

Processes the raw reads with Trimmomatic using three alternative trimming strategies:

- `trim01`
- `trim02`
- `trim03`

These alternative datasets were evaluated before assembly and variant calling.

## 1.2 Post-filter quality assessment

### `1.2.post-filter_quality_analysis.sh`

Runs FastQC independently on the three trimmed-read datasets to compare their quality profiles and assess their suitability for downstream analyses.

## 1.3 Assembly and variant calling

### `1.3.assembly_variant_calling_ipyrad/`

Contains the Jupyter notebooks used for assembly and variant calling with ipyrad.

Three trimming datasets were evaluated using three assembly strategies:

- de novo assembly;
- reference-based assembly using the *Quercus lobata* genome;
- reference-based assembly using the *Quercus robur* genome.

The nine analyses were executed independently and are retained as separate notebooks. The reference-based assembly using the *Q. robur* genome and the `trim01` dataset was selected for downstream analyses.

## 2.1 Genetic diversity and differentiation

### `2.1.genetic_diversity/`

Contains the scripts used to calculate genetic diversity and
differentiation statistics.

#### `2.1.1_heterozygosity_fstatistics.R`

Calculates observed and expected heterozygosity, inbreeding
coefficients, and genetic differentiation statistics.

#### `2.1.2_snp_diversity_tajimasD.R`

Calculates SNP-based nucleotide diversity, Watterson's theta, and
Tajima's D after retaining strictly biallelic SNPs.

Summary tables are written under `results/`.

## 2.2 Population structure

### `2.2.population_structure/`

Contains the scripts used for principal component analysis,
discriminant analysis of principal components, minimum spanning
networks, ADMIXTURE, and fastStructure.

#### `2.2.1_pca.R`

Performs principal component analysis and exports individual scores,
explained variance, and manuscript figures.

#### `2.2.2_dapc.R`

Performs discriminant analysis of principal components, including
cross-validation and export of scores, assignments, and membership
probabilities.

#### `2.2.3_msn.R`

Calculates Nei's genetic distances and generates a minimum spanning
network.

#### `2.2.4_prepare_plink.sh`

Converts the PED and MAP files exported with TASSEL into the PLINK
formats required by downstream analyses.

#### `2.2.5_run_admixture.sh`

Runs ADMIXTURE across the evaluated values of K and records
cross-validation results.

#### `2.2.6_run_faststructure.sh`

Runs fastStructure using the simple and logistic prior models.

#### `2.2.7_plot_admixture_faststructure.R`

Generates model-selection summaries and ancestry-proportion plots.

Retained ADMIXTURE and fastStructure outputs are stored under:

`data/1.4.population_structure/`

Detailed requirements and execution instructions are provided in:

`bin/2.2.population_structure/README.md`

## 2.3 Demographic history

### `2.3.demographic_history/`

Contains the documentation and plotting script for the demographic-
history analyses performed with easySFS and Stairway Plot 2.

All individuals were analyzed as a single population under six
scenarios combining three mutation rates and two generation times.

#### `2.3.1_plot_stairway_plot.py`

Reads the six retained Stairway Plot summary files, exports a scenario
summary table, and generates the combined demographic-history figure.

Associated files are stored under:

`data/1.8.demography/`

Detailed settings, inputs, and outputs are documented in:

`bin/2.3.demographic_history/README.md`

## 2.4 SNP outlier detection

### `2.4.snp_outlier_detection/`

Contains the workflow used to identify and compare putative candidate
SNPs using BayeScan, pcadapt, and locus-specific FST analyses.

#### `2.4.1_run_bayescan_9sites.sh`

Runs BayeScan using the nine sampling sites as separate populations.

#### `2.4.2_detect_candidate_snps.R`

Processes BayeScan results, performs pcadapt and locus-specific FST
analyses, compares candidate sets, and exports summary tables and
figures.

Associated BayeScan files are stored under:

`data/1.6.snp_outlier_detection/`

Detailed instructions are provided in:

`bin/2.4.snp_outlier_detection/README.md`

## 2.5 Sequence searches and outlier annotation

### `2.5.sequence_annotation/`

Contains the scripts used to prepare representative candidate-locus
sequences, perform BLAST searches, and compile annotation tables.

The workflow also documents the manual retrieval and inspection steps
used during sequence annotation.

Representative sequences are stored under:

`data/1.7.sequence_annotation/`

Detailed instructions are provided in:

`bin/2.5.sequence_annotation/README.md`

## 2.6 Outlier allele plots

### `2.6.outlier_allele_frequencies/`

Contains the scripts used to prepare the selected coding candidate SNP
table and generate individual allele plots.

#### `2.6.1_prepare_outlier_variant_table.R`

Extracts candidate genotypes, joins variant and annotation
information, and prepares the plotting table.

#### `2.6.2_plot_outlier_allele_frequencies.R`

Generates allele-dosage and allele-frequency plots for the retained
coding candidate SNPs.

The curated input table is stored under:

`data/1.8.outlier_allele_frequencies/`

Detailed instructions are provided in:

`bin/2.6.outlier_allele_frequencies/README.md`

Generated tables, figures, logs, and intermediate outputs are written
under `results/`.
