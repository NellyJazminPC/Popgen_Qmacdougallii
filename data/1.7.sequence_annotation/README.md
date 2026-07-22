# Candidate-locus sequences for annotation

This directory contains the candidate-locus sequences used for BLAST
searches and downstream annotation of putative outlier SNPs.

## Background

The outlier-detection workflow identified 124 candidate SNPs. After the
position suffix beginning with `_pos` was removed from the SNP
identifiers, these candidates represented 94 unique loci.

Representative sequences for these loci were manually retrieved and
curated from the ipyrad `.loci` assembly output using the corresponding
locus identifiers.

## Files

### `candidate_locus_sequences.fasta`

Contains the manually retrieved representative sequences for the 94
unique candidate loci.

### `candidate_locus_blast_queries.fasta`

Contains the query sequences used for the NCBI BLAST searches.

The following IUPAC ambiguity-code substitutions were applied during
query preparation:

- R to A
- Y to C
- S to G
- W to T
- K to G
- M to A

The character `N` was retained.

## Related workflow

The sequence-processing and BLAST-result parsing scripts are documented
in `bin/2.5.sequence_annotation/`.

Raw BLAST outputs and derived annotation tables are stored locally under
`results/`, which is excluded from version control.
