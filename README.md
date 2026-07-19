## Workflow overview

The genomic-data preprocessing and variant-calling workflow was organized into the following stages.

### 1.0 Initial quality assessment

Raw single-end GBS reads from 79 *Quercus macdougallii* individuals were evaluated using FastQC v0.11.9.

```text
Script:
bin/1.0.quality_analysis.sh

Input:
data/raw/

Documented output:
data/1.0.quality_analysis/
```

The raw sequencing files and individual FastQC reports are not included because of their size and because the reports can be regenerated from the original reads.

### 1.1 Read trimming

The raw reads were processed with Trimmomatic v0.39 using three alternative trimming strategies:

```text
trim01
trim02
trim03
```

```text
Script:
bin/1.1.filter_trimmomatic.sh

Documented output:
data/1.1.filter/
```

The three datasets represent independent parameter combinations applied to the same raw reads.

### 1.2 Post-filter quality assessment

The quality of the three trimmed-read datasets was evaluated independently using FastQC v0.11.9.

```text
Script:
bin/1.2.post-filter_quality_analysis.sh

Input:
data/1.1.filter/

Documented output:
data/1.2.post-filter_quality_analysis/
```

The individual HTML reports are treated as reproducible intermediate outputs and are not included in the repository.

### 1.3 Assembly and variant calling

Variant discovery was performed with ipyrad using three trimming datasets and three assembly strategies:

```text
                  trim01    trim02    trim03
De novo              ✓         ✓         ✓
Q. lobata genome     ✓         ✓         ✓
Q. robur genome      ✓         ✓         ✓
```

The corresponding notebooks are available in:

```text
bin/1.3.assembly_variant_calling_ipyrad/
```

The reference-based assembly using the *Q. robur* genome and the `trim01` dataset was selected for downstream analyses. After filtering in TASSEL v5.2.93, the final dataset contained 5,426 SNPs.

Assembly summaries and documentation are provided in:

```text
data/1.3.assembly_variant_calling/
```

Large intermediate assembly files are not included in the public repository.
