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

The script `convert2plinkformat.sh` uses these files to generate:

```text
qmacd_ref_gen_rob.raw
qmacd_ref_gen_rob.bed
qmacd_ref_gen_rob.bim
qmacd_ref_gen_rob.fam
```

The `.raw` file is an additive/dominance genotype matrix produced with PLINK `--recodeAD`. The binary PLINK files (`.bed`, `.bim`, and `.fam`) were used by ADMIXTURE, fastStructure, and pcadapt.

To run the conversion when PLINK is available in `PATH`:

```bash
bash bin/1.5.structure/convert2plinkformat.sh
```

A specific PLINK executable can be supplied with:

```bash
PLINK_BIN=/path/to/plink \
  bash bin/1.5.structure/convert2plinkformat.sh
```

## Model-based population structure

- `admixture_qmacd.sh` runs ADMIXTURE for K = 1–10 with cross-validation.
- `faststructure_qmacd.sh` runs fastStructure for K = 1–10 using the simple and logistic prior models.
- `../2.3.admixture_faststructure_plots.R` visualizes ancestry proportions and model-selection results.

These scripts will be cleaned and renamed as part of the final repository organization.

## Exploratory conversions

A PGDSpider configuration file was tested during exploratory work but was not part of the final PLINK-based ADMIXTURE and fastStructure workflow. That file is retained only in the local project archive.
