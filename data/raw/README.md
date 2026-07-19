# Raw sequencing data

This directory represents the input stage of the GBS preprocessing workflow for *Quercus macdougallii*.

It originally contained the raw single-end FASTQ files generated for 79 individuals. These files were used as input for the initial quality assessment and for read trimming.

The raw FASTQ files followed sample-specific naming patterns beginning with identifiers such as:

```text
CR_01_S115
CZ_01_S100
PZ_01_S161
```

and were expected to be stored as compressed files:

```text
<sample_identifier>*.fastq.gz
```

The raw sequencing files are not included in this repository because of their size and data-management considerations. This directory is retained through this README to document the expected input location and preserve the workflow structure.

The raw reads were used by:

```text
bin/1.0.quality_analysis.sh
bin/1.1.filter_trimmomatic.sh
```
