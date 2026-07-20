# Population structure workflow

This directory contains the scripts used to prepare genotype files and evaluate population structure in *Quercus macdougallii* using PCA, DAPC, minimum spanning networks, ADMIXTURE, and fastStructure.

## Input dataset

The analyses use the final reference-based SNP dataset:

```text
data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf
```

This VCF contains 5,426 filtered SNPs obtained using the *Quercus robur* reference genome and the `trim01` read-processing strategy.

## Preparing PLINK files

The final filtered VCF was opened in TASSEL v5.2.93 and exported through the graphical interface in PLINK PED/MAP format. The conversion produced:

```text
data/structure_formats/qmacd_ref_gen_rob.plk.ped
data/structure_formats/qmacd_ref_gen_rob.plk.map
```

The script `2.2.4_prepare_plink.sh` uses these files to generate:

```text
qmacd_ref_gen_rob.raw
qmacd_ref_gen_rob.bed
qmacd_ref_gen_rob.bim
qmacd_ref_gen_rob.fam
```

The `.raw` file is an additive/dominance genotype matrix produced with PLINK `--recodeAD`. The binary PLINK files (`.bed`, `.bim`, and `.fam`) were used by ADMIXTURE, fastStructure, and pcadapt.

To run the conversion when PLINK is available in `PATH`:

```bash
bash bin/2.2.population_structure/2.2.4_prepare_plink.sh
```

A specific PLINK executable can be supplied with:

```bash
PLINK_BIN=/path/to/plink \
  bash bin/2.2.population_structure/2.2.4_prepare_plink.sh
```

## Model-based population structure

- `2.2.5_run_admixture.sh` runs ADMIXTURE for K = 1-10 using cross-validation.
- `2.2.6_run_faststructure.sh` runs fastStructure for K = 1-10 using the simple and logistic prior models.
- `2.2.7_plot_admixture_faststructure.R` visualizes the retained ancestry proportions and ADMIXTURE model-selection results.

The retained ADMIXTURE `.Q` files, fastStructure `.meanQ` files, and model-selection summaries are stored in:

```text
data/1.4.population_structure/
```

Auxiliary allele-frequency, variance, and log files are excluded from the public repository because they can be regenerated using the corresponding analysis scripts.

## Multivariate analyses

The following scripts will contain the multivariate population-structure analyses:

```text
2.2.1_pca.R
2.2.2_dapc.R
2.2.3_msn.R
```

These analyses are currently being separated from the original combined R workflow.

## Geographic panel

The geographic panel of the population-structure figure is not reproduced by the public plotting script because the coordinates of this threatened microendemic species are not included in the public metadata file.

## Exploratory conversions

A PGDSpider configuration file was tested during exploratory work but was not part of the final PLINK-based ADMIXTURE and fastStructure workflow. That file is retained only in the local project archive.
