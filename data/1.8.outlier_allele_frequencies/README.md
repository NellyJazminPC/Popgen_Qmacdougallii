# Outlier SNPs selected for allele plots

This directory contains the coding candidate SNPs retained for
individual allele plots.

## Selection of SNPs

Individual allele plots were initially generated for a broader
manually selected set of candidate SNPs.

A final set of 11 variants located in coding regions was retained
after manual inspection of coding-region placement, predicted
amino-acid changes, and sequence-alignment quality.

The retained variants comprise:

- four synonymous variants;
- six missense variants;
- one stop-gain variant.

## File

### `manuscript_outlier_snps_for_plotting.csv`

Contains the 11 coding candidate SNPs retained for individual allele plots.

The columns are:

- `locus_name`: SNP identifier in the VCF;
- `detection_method`: analysis associated with the retained selection
  metric for the SNP;
- `shared_in`: comma-separated list of all outlier-detection analyses
  that identified the SNP;
- `selection_metric`: adjusted p-value or FST;
- `metric_value`: value of the corresponding selection metric;
- `ref_allele`: reference allele;
- `alt_allele`: alternative allele;
- `best_hit_magnoliopsida`: retained Magnoliopsida BLAST match;
- `variant_effect`: synonymous, missense, or stop-gain.

The `detection_method`, `selection_metric`, and `metric_value` columns
record the analysis and value retained for labeling and consultation.
The `shared_in` column preserves all analyses in which the SNP was
identified as a candidate.

## Related workflow

The scripts used to extract variant information and generate
individual allele plots are documented under:

`bin/2.6.outlier_allele_frequencies/`

Generated tables and figures are stored locally under `results/`, which
is excluded from version control.
