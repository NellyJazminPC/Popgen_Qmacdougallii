# Outlier SNP allele plots

This directory contains the scripts used to extract genotype and allele
information for candidate outlier SNPs and to generate the individual
plots.

## Workflow overview

A broader set of candidate SNPs was initially inspected during the
analysis. A final set of 11 coding variants was retained after manual
inspection of coding-region placement, predicted amino-acid changes,
and sequence-alignment quality.

The retained variants comprise:

- four synonymous variants;
- six missense variants;
- one stop-gain variant.

The 11 SNPs used for plotting are listed in:

`data/1.8.outlier_allele_frequencies/manuscript_outlier_snps_for_plotting.csv`

## Scripts

### `2.6.1_prepare_outlier_variant_table.R`

Reads the filtered VCF and the Magnoliopsida-focused annotation table.

The script:

- extracts the candidate SNP genotypes;
- joins REF and ALT allele information;
- translates genotype codes into nucleotide alleles;
- generates a per-locus genotype table;
- classifies single-nucleotide substitutions as transitions or
  transversions.

Main inputs:

- `data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf`
- `results/snps_outliers_with_blast_results_magnoliopsida_highlighted.xlsx`

Main outputs:

- `results/SNPs_outliers_variants_per_locus.xlsx`
- `results/snps_outliers_with_blast_and_variants.xlsx`

### `2.6.2_plot_outlier_allele_frequencies.R`

Generates individual plots for the 11 coding SNPs retained after
manual inspection. Each SNP is exported as a separate PNG file.

Two representations are produced.

#### Dosage-scale representation

Uses `tab(..., freq = FALSE)`, which returns allele-copy counts for each
diploid individual. Mean values by sampling site are therefore shown on
a scale from 0 to 2.

This representation is written to:

`results/plots_snps_outliers_dosage_scale/`

#### Allele-frequency representation

Uses `tab(..., freq = TRUE)`, which converts the allele-copy counts to
individual allele proportions before calculating site means. Values are
therefore shown on the conventional allele-frequency scale from 0 to 1.

These plots are written to:

`results/plots_snps_outliers_frequency_scale/`

Both representations preserve the same relative differentiation
patterns among sampling sites. They differ only in the numerical scale
used for the y-axis.

## Plot information

The plot titles include either the adjusted pcadapt p-value or the
locus-specific FST value associated with each SNP.

The allele symbols correspond to the REF and ALT nucleotides recorded
in the VCF. The retained Magnoliopsida BLAST match is displayed below
each plot.

## Manual selection

The original analysis generated plots for a broader manually selected
set of candidate SNPs. The final set of 11 variants was selected after
considering coding-region placement, predicted variant effects, and
alignment quality.

## Results and repository curation

Generated tables and figures are written under `results/`, which is
excluded from version control.

Comments, filenames, and paths were standardized during repository
curation. The scripts provide a practical guide to the workflow used
for the manuscript and may require adjustment for other datasets or
computing environments.
