# Population structure outputs

This directory contains the lightweight derived outputs retained from the ADMIXTURE and fastStructure analyses of *Quercus macdougallii*.

The binary PLINK input files used by both programs are located in:

```text
data/structure_formats/
```

## ADMIXTURE

The `admixture_output/` directory retains:

- `chooseK.txt`, containing the cross-validation error for K = 1–10.
- `.Q` files for K = 1–10, containing the ancestry proportions estimated for each individual.

The `.Q` files were retained because they are the direct inputs used to generate the ADMIXTURE ancestry plots.

Allele-frequency files (`.P`) and execution logs were omitted because they are auxiliary outputs and can be regenerated with the corresponding analysis script.

## fastStructure

The `faststructure_output/` directory retains:

- The `chooseK` summaries for the simple and logistic prior models.
- `.meanQ` files for K = 1–10 under both prior models.

The `.meanQ` files contain the ancestry proportions used to generate the fastStructure plots.

The `.meanP`, `.varP`, `.varQ`, and log files were omitted because they are auxiliary outputs and can be regenerated using the corresponding analysis script.

## Reproducibility

The scripts used to prepare the PLINK files, run ADMIXTURE and fastStructure, and visualize the ancestry results are documented under the population-structure section of `bin/`.
