# Candidate-locus sequence searches and annotation

This directory contains the scripts used to prepare candidate-locus
sequences, document the NCBI BLAST searches, and compile the annotation
tables used in the manuscript.

## Workflow overview

The outlier-detection analysis identified 124 candidate SNPs
representing 94 unique loci.

The workflow included the following steps:

1. identify the unique candidate loci;
2. manually retrieve representative sequences from the ipyrad `.loci`
   assembly output;
3. prepare the candidate-locus sequences for BLAST searches;
4. compare unrestricted, Viridiplantae, and Magnoliopsida BLAST searches;
5. retain Magnoliopsida for the final interpretation;
6. inspect selected sequences and alignments manually in Geneious.

The manual sequence-retrieval and Geneious inspection steps are
documented but are not performed automatically by the scripts.

## Scripts

### `2.5.1_prepare_candidate_sequences.R`

Reads the candidate-SNP table, removes the `_pos` suffix from SNP
identifiers, identifies the 94 unique loci, validates the curated FASTA
files, joins sequences to the candidate-SNP table, and records the IUPAC
ambiguity codes used during query preparation.

Main inputs:

- `results/consolidated_snps_results_with_shared_info_unique.xlsx`
- `data/1.7.sequence_annotation/candidate_locus_sequences.fasta`
- `data/1.7.sequence_annotation/candidate_locus_blast_queries.fasta`

Derived tables are written under `results/`.

### `2.5.2_run_blast.sh`

Documents the remote NCBI BLAST search retained for the final analysis.

Default settings:

- program: `blastn`
- database: `nt`
- maximum hits: 5
- taxonomic filter: `Magnoliopsida[ORGANISM]`
- output directory: `results/blast_results_Magnoliopsida/`

The settings and file paths can be adjusted through environment
variables when applying the script to another dataset.

### `2.5.3_compile_blast_results.R`

Parses the retained BLAST text files from three search scopes:

- unrestricted NCBI `nt`;
- Viridiplantae;
- Magnoliopsida.

The script extracts the best hit, E-value, and percentage identity and
generates comparative and Magnoliopsida-focused annotation tables.

## Choice of taxonomic scope

Magnoliopsida was initially used for the sequence searches. A broader
Viridiplantae search was subsequently explored following discussion
with the research team.

Viridiplantae produced more broad or nonspecific matches, whereas
Magnoliopsida produced results that were generally more interpretable
for the study system. Magnoliopsida was therefore retained for the
final annotation and biological interpretation.

## Manual inspection in Geneious

Selected query sequences and corresponding reference or best-hit
sequences were inspected manually in Geneious.

Sequence alignments were used to evaluate the position of the focal SNP
relative to coding regions and to assess potential synonymous,
nonsynonymous, and stop-gain changes.

## Repository curation

Comments, filenames, and paths were standardized during repository
curation to make the workflow easier to follow.

The BLAST searches were not rerun during this process. The retained
outputs and derived spreadsheets from the original analysis were
preserved.

These scripts are provided as a transparent record and practical guide
to the workflow used for the manuscript. They may require adjustment
when applied to other datasets or computing environments.
