# Assembly and variant calling with ipyrad

This directory contains the Jupyter notebooks used for the comparative assembly and variant-calling analyses of the GBS data generated for *Quercus macdougallii*.

## Assembly design

Three alternative read-trimming datasets were evaluated:

* `trim01`
* `trim02`
* `trim03`

Each dataset was analyzed using three assembly strategies:

* de novo assembly;
* reference-based assembly using the *Quercus lobata* genome;
* reference-based assembly using the *Quercus robur* genome.

This resulted in nine independent ipyrad analyses:

```text
denovo_trim01.ipynb
denovo_trim02.ipynb
denovo_trim03.ipynb

ref_gen_qlobata_trim01.ipynb
ref_gen_qlobata_trim02.ipynb
ref_gen_qlobata_trim03.ipynb

ref_gen_qrobur_trim01.ipynb
ref_gen_qrobur_trim02.ipynb
ref_gen_qrobur_trim03.ipynb
```

The analyses were kept as separate notebooks because each assembly was executed independently and required substantial computational time.

## Library type

The sequencing data were analyzed in ipyrad using:

```python
datatype = "ddrad"
```

Although the sequencing approach is described generally as genotyping-by-sequencing, the library was prepared using two restriction enzymes. The `ddrad` datatype was therefore selected for the ipyrad analyses.

## Population-assignment files

Population-assignment files were used as required during the ipyrad analyses. Separate population maps were also created during exploratory analyses with Stacks.

Several working versions were generated while alternative population groupings and filtering strategies were evaluated. Because the archived versions and their historical filenames require further verification, these files are not currently included in the public repository.

The exact file structure, delimiters, and population assignments will be documented after the original input files have been validated.

## Selected assembly

The resulting VCF files were compared among assembly strategies and trimming datasets. The reference-based assembly using the *Q. robur* genome and the `trim01` reads was selected for downstream analyses.

This assembly initially contained 7,611 SNPs. After site-level filtering and the removal of indels and minor SNP states in TASSEL v5.2.95, the final dataset contained 5,426 SNPs.

## Outputs

Large intermediate assembly files and VCF outputs are not distributed in this directory. The filtered VCF used in the study is available through the Zenodo repository associated with the manuscript.

Some relative paths and historical working filenames remain in the original notebooks. These notebooks are retained as records of the analyses that were executed and will be documented further without altering their analytical logic.
