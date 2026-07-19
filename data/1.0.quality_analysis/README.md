# Initial quality assessment

This directory corresponds to the initial quality-control assessment of the raw GBS reads for *Quercus macdougallii*.

The raw single-end FASTQ files from 79 individuals were evaluated with FastQC v0.11.9 before read trimming with Trimmomatic.

The analysis was performed using:

```text
bin/1.0.quality_analysis.sh
```

The script used the raw compressed FASTQ files expected in:

```text
data/raw/
```

and wrote the FastQC reports to:

```text
data/1.0.quality_analysis/
```

FastQC generated one quality report for each raw sequencing file. These reports were used to inspect sequence-quality profiles before read trimming.

The individual FastQC reports are not included in the public repository because they are intermediate and reproducible outputs and would add a large number of files. They can be regenerated using the corresponding script when the raw FASTQ files and FastQC v0.11.9 are available.

The software version is also documented in:

```text
bin/software/README.md
```

This directory is retained through this README to preserve the preprocessing workflow structure.

