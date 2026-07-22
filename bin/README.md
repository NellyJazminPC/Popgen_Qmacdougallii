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

Contains the R scripts used to estimate genetic diversity and differentiation from the final *Quercus robur*-based SNP dataset.

#### `2.1.1_heterozygosity_fstatistics.R`

Calculates observed heterozygosity (Ho), expected heterozygosity (He), inbreeding coefficients (FIS), observed and private alleles, global FST, pairwise FST among sampling sites, and differentiation between the northern and southern geographic zones.

#### `2.1.2_snp_diversity_tajimasD.R`

Calculates nucleotide diversity (π), Watterson's theta (θW), and Tajima's D by sampling site, geographic zone, and across all individuals.

## 2.2 Population structure

### `2.2.population_structure/`

Contains the scripts used to prepare genotype files and evaluate population structure using multivariate, network-based, and model-based approaches.

#### `2.2.1_pca.R`

Performs principal component analysis, exports PCA scores and explained variance, and regenerates Figure 4. PC1 and PC2 explain 7.31% of the total genomic variation.

#### `2.2.2_dapc.R`

Performs discriminant analysis of principal components using the nine sampling sites as a priori groups. The optional `find.clusters` exploration is retained as a disabled block. Cross-validation is conducted first across 5–50 PCs and then across 5–30 PCs. The focused analysis uses 100 serial replicates, a training proportion of 0.90, eight discriminant axes, and a random seed of 999. Seven PCs are retained based on the lowest root mean squared error.

#### `2.2.3_msn.R`

Calculates Nei's genetic distances among individuals and generates a minimum spanning network using a Kamada-Kawai layout.

#### `2.2.4_prepare_plink.sh`

Generates additive/dominance and binary PLINK files from the PED/MAP files exported with TASSEL.

#### `2.2.5_run_admixture.sh`

Runs ADMIXTURE for K = 1–10 using cross-validation.

#### `2.2.6_run_faststructure.sh`

Runs fastStructure for K = 1–10 using the simple and logistic prior models.

#### `2.2.7_plot_admixture_faststructure.R`

Generates the ADMIXTURE cross-validation plot and the ancestry-proportion plots used in the manuscript.

Conversion of the VCF to a `genlight` object excludes 41 multiallelic loci, leaving 5,385 biallelic SNPs for PCA, DAPC, and the minimum spanning network. Each script validates sample identifiers against the public metadata before running the analysis.

The retained ADMIXTURE `.Q` files, fastStructure `.meanQ` files, and model-selection summaries are stored in:

```text
data/1.4.population_structure/
```

The geographic panel of the population-structure figure is not regenerated by the public plotting script because precise coordinates for this threatened microendemic species are not included in the public metadata.

## 2.4 SNP outlier detection

### `2.4.snp_outlier_detection/`

Contains the reproducible workflow used to identify and compare putative candidate SNPs with BayeScan, pcadapt, and locus-specific FST analyses.

#### `2.4.1_run_bayescan_9sites.sh`

Runs the retained BayeScan analysis using the nine sampling sites as separate populations. The script resolves repository paths automatically, selects the appropriate bundled executable for macOS or Linux, and supports a dry-run mode for command validation.

#### `2.4.2_detect_candidate_snps.R`

Processes the nine-site BayeScan results, performs pcadapt with Bonferroni correction, calculates locus-specific FST under three population-grouping scenarios, and consolidates the candidate SNP sets.

The BayeScan input and retained output are stored in:

```text
data/1.6.snp_outlier_detection/
```

Detailed execution instructions, required inputs, software dependencies, and main outputs are documented in:

```text
bin/2.4.snp_outlier_detection/README.md
```

## 2.5 Sequence searches and outlier annotation

### `2.5.sequence_annotation/`

Contains the workflow used to prepare representative sequences for
candidate outlier loci, document the NCBI BLAST searches, compare
alternative taxonomic search scopes, and compile the annotation tables.

The workflow includes manual steps for retrieving candidate sequences
from the ipyrad `.loci` output and inspecting selected alignments in
Geneious.

Magnoliopsida was retained for the final annotation and biological
interpretation.

Detailed inputs, scripts, manual steps, and expected outputs are
documented in:

`bin/2.5.sequence_annotation/README.md`

## 2.6 Outlier allele frequencies

### `2.6.snps_outliers_freq.R`

Analyzes allele frequencies of candidate outlier SNPs.

## 2.7 Effective population size

### `2.7.Ne.R`

Estimates effective population size (Ne) from filtered genotype data.

## 3.1 Demographic history

### `3.1.demography.md`

Documents and summarizes the demographic-history analyses.

Additional input/output details and stage-specific instructions are provided in the README files associated with individual workflow sections.
