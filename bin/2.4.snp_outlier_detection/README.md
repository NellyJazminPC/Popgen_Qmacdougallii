# SNP outlier detection

This directory contains the workflow used to identify putative candidate SNPs in *Quercus macdougallii* using BayeScan, pcadapt, and locus-specific FST analyses.

## Scripts

### `2.4.1_run_bayescan_9sites.sh`

Runs BayeScan for the nine sampling sites.

The script:

- resolves repository paths automatically;
- selects the bundled macOS or Linux BayeScan executable;
- detects the available number of processors;
- accepts the `THREADS` environment variable to override the default;
- supports `DRY_RUN=1` to validate paths and commands without executing BayeScan.

Run from the repository root with:

```bash
bin/2.4.snp_outlier_detection/2.4.1_run_bayescan_9sites.sh
```

Validate without executing BayeScan:

```bash
DRY_RUN=1 \
  bin/2.4.snp_outlier_detection/2.4.1_run_bayescan_9sites.sh
```

The prepared BayeScan input file is stored in:

```text
data/1.6.snp_outlier_detection/bayescan_input/
```

BayeScan output is written to:

```text
data/1.6.snp_outlier_detection/bayescan_output/
```

The BayeScan input was generated from the filtered VCF using PGDSpider during the original analysis. The prepared input file is retained to reproduce the BayeScan run.

Exploratory BayeScan analyses using North–South and PZ-versus-remaining-site groupings were not retained in the final candidate-locus workflow.

### `2.4.2_detect_candidate_snps.R`

Processes BayeScan results and performs the remaining candidate-locus analyses.

The script includes:

1. BayeScan evaluation for the nine sampling sites.
2. pcadapt analysis with:
   - `K = 2`;
   - minimum minor allele frequency of 0.05;
   - Bonferroni correction at alpha = 0.05.
3. Locus-specific FST analyses for:
   - the nine sampling sites;
   - the northern and southern geographic zones;
   - PZ versus the remaining sampling sites.
4. Selection of the upper 1% of the empirical locus-specific FST distribution in each grouping.
5. Comparison and consolidation of candidate SNP sets.
6. Generation of summary plots, Venn diagrams, and result tables.

Run from the repository root with:

```bash
Rscript \
  bin/2.4.snp_outlier_detection/2.4.2_detect_candidate_snps.R
```

The script resolves its location automatically and does not depend on the user's working directory.

## Required inputs

The workflow uses:

```text
data/1.6.snp_outlier_detection/
data/structure_formats/qmacd_ref_gen_rob.bed
data/structure_formats/qmacd_ref_gen_rob.bim
data/structure_formats/qmacd_ref_gen_rob.fam
data/structure_formats/qmacd_ref_gen_rob.raw
data/structure_formats/qmacd_ref_gen_rob.plk.map
metadata/Qmacdougalli_79ind_.csv
```

## Main outputs

Output files are written to `results/`.

The main manuscript Venn diagram is:

```text
results/venn_pcadapt_fst_outliers.png
```

It compares pcadapt with the three locus-specific FST scenarios. BayeScan is not included in this main figure because its two candidate SNPs did not overlap with candidates detected by pcadapt or FST.

An additional exploratory diagram including all methods is written to:

```text
results/venn_all_methods_exploratory.png
```

The final consolidated candidate table is:

```text
results/consolidated_snps_results_with_shared_info_unique.xlsx
```

The retained analyses identified:

- 18 SNPs with pcadapt;
- 54–55 SNPs within each locus-specific FST scenario;
- 122 unique SNPs detected by pcadapt and/or FST;
- 2 additional SNPs detected exclusively by BayeScan;
- 124 unique candidate SNPs in total.

## R packages

The R workflow requires:

```text
boa
coda
writexl
ggplot2
pcadapt
qvalue
snpStats
adegenet
viridis
dplyr
VennDiagram
```

Package installation is intentionally not performed inside the analysis script.
