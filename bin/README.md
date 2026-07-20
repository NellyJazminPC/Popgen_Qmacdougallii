# Analysis pipeline

This directory contains the scripts, Jupyter notebooks, configuration files, and workflow documentation used for the population genomic analyses of *Quercus macdougallii*.

The workflow is organized into numbered stages. Each stage corresponds to a script or subdirectory and follows the general order in which the analyses were performed.

## Running the R analyses

The R scripts use relative paths assuming that `bin/` is the working directory. For reproducibility, users are advised to create an RStudio project inside the `bin/` directory before running the analyses, or to manually set `bin/` as the working directory.

The local `.Rproj` file used during the original analyses is not included in the repository because it may contain machine-specific settings. 
Input files are read from `../data/` and `../metadata/`, and generated outputs are written to `../results/`.

## 1.0 Initial quality assessment

### `1.0.quality_analysis.sh`

Performs the initial quality assessment of the 79 raw single-end GBS read files using FastQC.

## 1.1 Read trimming

### `1.1.filter_trimmomatic.sh`

Processes the raw reads with Trimmomatic using three alternative trimming strategies:

* `trim01`
* `trim02`
* `trim03`

These alternative datasets were evaluated before assembly and variant calling.

## 1.2 Post-filter quality assessment

### `1.2.post-filter_quality_analysis.sh`

Runs FastQC independently on the three trimmed-read datasets to compare their quality profiles and assess their suitability for downstream analyses.

## 1.3 Assembly and variant calling

### `1.3.assembly_variant_calling_ipyrad/`

Contains the Jupyter notebooks used for assembly and variant calling with ipyrad.

Three trimming datasets were evaluated using three assembly strategies:

* de novo assembly;
* reference-based assembly using the *Quercus lobata* genome;
* reference-based assembly using the *Quercus robur* genome.

The nine analyses were executed independently and are retained as separate notebooks. The reference-based assembly using the *Q. robur* genome and the `trim01` dataset was selected for downstream analyses.

## 1.5 Population structure

### `1.5.structure/`

Preparation and analysis of population structure using multiple approaches:

* `admixture_qmacd.sh`: Runs ADMIXTURE for ancestry estimation.
* `faststructure_qmacd.sh`: Runs fastStructure for population structure inference.
* `convert2plinkformat.sh`: Converts genotype data to PLINK format for compatibility with various tools.
* `workflow_pop_structure_analysis.md`: Documentation of the population structure analysis workflow.
* `faststructure_convert.spid`: Configuration file for data conversion.

## 1.6 Outlier SNP detection

### `1.6.snps_outliers/`

Detection of candidate loci under selection using tools such as BayeScan, and preparation of files for selection analyses.

## 2.1 Genetic diversity and differentiation

### `2.1.genetic_diversity/`

Contains the R scripts used to estimate genetic diversity and
differentiation from the final *Quercus robur*-based SNP dataset.

#### `2.1.1_heterozygosity_fstatistics.R`

Calculates observed heterozygosity (Ho), expected heterozygosity (He),
inbreeding coefficients (FIS), observed and private alleles, global FST,
pairwise FST among sampling sites, and differentiation between the
northern and southern geographic zones.

#### `2.1.2_snp_diversity_tajimasD.R`

Calculates nucleotide diversity (π), Watterson's theta (θW), and
Tajima's D by sampling site, geographic zone, and across all individuals.

## 2.2 Population structure in R

### `2.2.pop_structure.R`

R script for analyzing population structure using various methods.

## 2.3 ADMIXTURE and fastStructure plots

### `2.3.admixture_faststructure_plots.R`

R script for generating plots of ADMIXTURE and fastStructure results.

## 2.4 Outlier SNP analyses

### `2.4.snps_outliers_PCAdapt_Bayescan_FST.R`

R script for identifying outlier SNPs using pcadapt, BayeScan, and FST-based approaches.

## 2.5 Sequence searches and outlier annotation

### `2.5.seq_search_blast.sh`

### `2.5.snps_outliers_sequences.R`

Scripts for searching SNP-associated sequences with BLAST and performing further analyses of outlier SNP sequences.

## 2.6 Outlier allele frequencies

### `2.6.snps_outliers_freq.R`

R script for analyzing allele frequencies of outlier SNPs.

## 2.7 Effective population size

### `2.7.Ne.R`

R script for estimating effective population size, Ne, from filtered genotype data.

## 3.1 Demographic history

### `3.1.demography.md`

Documentation and summary of the demographic analyses.

The directory combines Bash scripts, R scripts, Jupyter notebooks, configuration files, and workflow notes to document the analytical procedures used in the study. Additional details and input/output descriptions are provided in the README files associated with individual workflow stages.
