# Post-filter quality assessment

This directory corresponds to the quality-control assessment performed after read trimming.

Raw single-end GBS reads were processed with Trimmomatic v0.39 using three alternative trimming strategies, identified as `trim01`, `trim02`, and `trim03`. The quality of the resulting reads was subsequently evaluated with FastQC v0.11.9.

## Workflow

The post-filter quality assessment was performed using:

```text
bin/1.2.post-filter_quality_analysis.sh
```

The script runs FastQC independently for each sample and for each of the three trimming strategies.

The filtered reads used as input were generated during the previous workflow step and stored in:

```text
data/1.1.filter/
```

The `trim01` and `trim02` datasets were analyzed as compressed FASTQ files (`.fastq.gz`), whereas the `trim03` dataset was analyzed in uncompressed FASTQ format (`.fastq`).

The resulting FastQC reports were written to:

```text
data/1.2.post-filter_quality_analysis/
```

FastQC was installed locally under:

```text
bin/software/FastQC/
```

Software versions used throughout the workflow are documented in:

```text
bin/software/README.md
```

## Outputs

FastQC generated one HTML quality report per sample and trimming strategy, with filenames such as:

```text
CR_01_S115.trim01_fastqc.html
CR_01_S115.trim02_fastqc.html
CR_01_S115.trim03_fastqc.html
```

These reports were used to compare the effects of the three trimming strategies on sequence quality and to assess whether the processed reads were suitable for downstream assembly and variant-calling analyses.

The individual FastQC HTML reports are not included in the public repository because they are intermediate and reproducible outputs and would add a large number of files. They are retained locally and can be regenerated using the corresponding script, provided that FastQC v0.11.9 and the filtered FASTQ files are available.

## Related workflow step

Read trimming:

```text
bin/1.1.filter_trimmomatic.sh
data/1.1.filter/
```
