# SNP outlier detection data

This directory contains the BayeScan input and retained output used in the final SNP outlier detection workflow.

## Structure

```text
1.6.snp_outlier_detection/
├── bayescan_input/
│   └── bayescan_qmacd_ref_gen_qrob_9pop
├── bayescan_output/
│   └── bayescan_qmacd_ref_gen_qrob_9pop/
└── README.md
```

## BayeScan input

The retained input represents the nine sampling sites analyzed as separate populations.

The file was generated from the filtered SNP dataset using PGDSpider during the original analysis. Because the original conversion involved the PGDSpider graphical workflow, the prepared input file is retained here to make the BayeScan execution reproducible.

## BayeScan output

The retained output directory contains:

- the posterior locus statistics;
- the sampled Markov chain;
- acceptance-rate information;
- verification output.

These files are processed by:

```text
bin/2.4.snp_outlier_detection/2.4.2_detect_candidate_snps.R
```

At a false discovery rate threshold of 0.05, BayeScan identified two candidate SNPs in the nine-site analysis.

Exploratory BayeScan analyses based on North–South and PZ-versus-remaining-site groupings were not used in the final candidate-locus workflow and are not included in this public data directory.
