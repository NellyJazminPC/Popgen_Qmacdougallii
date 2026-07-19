# Read trimming outputs

This directory corresponds to the read-trimming step of the GBS preprocessing workflow for *Quercus macdougallii*.

Following the initial quality assessment of the raw sequencing reads with FastQC v0.11.9, the single-end FASTQ files from 79 individuals were processed with Trimmomatic v0.39.

The trimming analysis was performed using:

```text
bin/1.1.filter_trimmomatic.sh
```

The script used the raw compressed FASTQ files expected in:

```text
data/raw/
```

and generated three alternative trimming datasets in:

```text
data/1.1.filter/
```

Software versions used throughout the workflow are documented in:

```text
bin/software/README.md
```

## Trimming strategies

The three trimming strategies were applied independently to the same raw sequencing files. They represent alternative parameter combinations rather than sequential filtering steps.

### `trim01`

```text
LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 CROP:140 HEADCROP:10 MINLEN:100
```

This strategy used base-quality thresholds of 20 and retained reads with a minimum length of 100 bp.

### `trim02`

```text
LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 CROP:140 HEADCROP:10 MINLEN:80
```

This strategy used base-quality thresholds of 30 and retained reads with a minimum length of 80 bp.

### `trim03`

```text
SLIDINGWINDOW:4:20 CROP:120 HEADCROP:20 MINLEN:100
```

This strategy combined sliding-window filtering with fixed cropping to retain reads of 100 bp.

## Directory contents

The local working directory used for the analysis contained 237 symbolic links, corresponding to 79 individuals and three trimming strategies.

The files followed these naming patterns:

```text
CR_01_S115.trim01.fastq.gz
CR_01_S115.trim02.fastq.gz
CR_01_S115.trim03.fastq
```

The `trim01` and `trim02` outputs were retained as compressed FASTQ files. The `trim03` outputs were initially generated as compressed files and were subsequently used in uncompressed FASTQ format during the post-filter quality assessment.

The symbolic links pointed to the corresponding trimmed FASTQ files stored outside the repository and were used only within the local working environment.

The raw and trimmed FASTQ files, as well as their symbolic links, are not included in the public repository because of file size, portability, and data-management considerations. This directory is retained through this README to document the workflow and the expected intermediate outputs.

## Related workflow steps

Initial quality assessment of raw reads:

```text
bin/1.0.quality_analysis.sh
data/1.0.quality_analysis/
```

Read trimming:

```text
bin/1.1.filter_trimmomatic.sh
data/1.1.filter/
```

Post-filter quality assessment:

```text
bin/1.2.post-filter_quality_analysis.sh
data/1.2.post-filter_quality_analysis/
```
