# Filtered FASTQ symbolic links

This folder originally contained symbolic links to trimmed FASTQ files generated during the preprocessing step of the GBS dataset for *Quercus macdougallii*.

The raw and trimmed FASTQ files are not included in this repository because of file size and data management considerations. The symbolic links were used only in the local working directory to document and run the preprocessing workflow.

Files followed the naming pattern:

- `CR_01_S115.trim01.fastq`
- `CR_01_S115.trim02.fastq`
- `CR_01_S115.trim03.fastq`
- `CR_02_S127.trim01.fastq`
- `CR_02_S127.trim02.fastq`
- `CR_02_S127.trim03.fastq`

where sample identifiers indicate population/sample codes and the `trim01`, `trim02`, and `trim03` suffixes correspond to sequential trimming/filtering outputs generated during read preprocessing.
