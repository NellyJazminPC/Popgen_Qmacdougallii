# Population genomics of *Quercus macdougallii*

This repository contains the scripts, Jupyter notebooks, selected
lightweight data products, and supporting documentation used for the
population genomic analyses of the endangered microendemic oak
*Quercus macdougallii*.

The study includes 79 individuals sampled across nine sites in the
Sierra Juárez of Oaxaca, Mexico.

## Data availability

The filtered Variant Call Format (VCF) dataset used in the study is
available from Zenodo:

[https://doi.org/10.5281/zenodo.20548053](https://doi.org/10.5281/zenodo.20548053)

Raw sequencing reads, reference genomes, large assembly outputs,
locally generated genotype formats, and generated analysis results are
not distributed through GitHub.

## Repository structure

| Directory | Contents |
|---|---|
| [`bin/`](bin/) | Analysis scripts, notebooks, software documentation, and workflow instructions |
| [`data/`](data/) | Selected inputs, intermediate files, retained outputs, and data documentation |
| [`metadata/`](metadata/) | Public metadata for the 79 sampled individuals |
| [`results/`](results/) | Local destination for generated tables, figures, logs, and intermediate outputs |

## Workflow overview

### Sequence processing and variant calling

1. **Initial quality assessment**  
   Evaluation of the raw genotyping-by-sequencing reads with FastQC.

2. **Read trimming**  
   Processing of the sequencing reads with alternative Trimmomatic
   parameter combinations.

3. **Post-filter quality assessment**  
   Evaluation of the trimmed-read datasets with FastQC.

4. **Assembly and variant calling**  
   Comparison of de novo and reference-based ipyrad assemblies,
   followed by variant filtering for downstream analyses.

### Population genomic analyses

1. **Genetic diversity and differentiation**  
   Estimation of heterozygosity, inbreeding, genetic differentiation,
   nucleotide diversity, Watterson's theta, and Tajima's D.

2. **Population structure**  
   Principal component analysis, discriminant analysis of principal
   components, minimum spanning networks, ADMIXTURE, and fastStructure.

3. **Demographic history**  
   Site-frequency-spectrum preparation and Stairway Plot analyses.

4. **SNP outlier detection**  
   Identification and comparison of putative candidate SNPs using
   BayeScan, pcadapt, and locus-specific FST analyses.

5. **Sequence searches and annotation**  
   Preparation of candidate-locus sequences, BLAST searches, and
   compilation of annotation information.

6. **Candidate SNP allele plots**  
   Preparation of coding candidate SNP tables and generation of
   allele-dosage and allele-frequency plots.

Detailed descriptions of the scripts, required inputs, generated
outputs, and execution steps are provided in
[`bin/README.md`](bin/README.md) and in the README for each analysis
module.

## Reproducibility notes

The repository documents the workflow used for the analyses presented
in the associated manuscript. Some stages include manual steps carried
out with external software; these steps are identified in the relevant
module documentation.

Generated files under `results/` and large local working files are
excluded from version control.

## License

Licensing information is provided in [`LICENSE`](LICENSE).
