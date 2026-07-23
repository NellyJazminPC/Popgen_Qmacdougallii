# Demographic history

This module documents the demographic-history analyses performed with
easySFS and Stairway Plot 2 for *Quercus macdougallii*.

All 79 individuals were analyzed as a single population. The
demographic reconstruction was evaluated under six scenarios combining
three mutation rates and two generation times.

## Workflow

### 1. Site-frequency spectrum

A folded one-dimensional site-frequency spectrum was generated from the
filtered VCF using easySFS.

The retained settings were:

- population grouping: all individuals as one population;
- number of sampled individuals: 79 diploid individuals;
- projection: 156 chromosomes;
- folded SFS: yes.

The population map used for this step is:

`data/structure_formats/popmap_one_pop.txt`

The filtered VCF is available from Zenodo:

`https://doi.org/10.5281/zenodo.20548053`

Selected easySFS files are stored under:

`data/1.8.demography/easy_sfs_one_pop/`

### 2. Stairway Plot

Stairway Plot 2 was evaluated using the following retained blueprint
settings:

- `nseq`: 158;
- `L`: 5,426;
- `whether_folded`: true;
- `pct_training`: 0.67;
- `nrand`: 39, 78, 117, and 156;
- `ninput`: 200.

The six scenario-specific blueprint files are stored under:

`bin/2.3.demographic_history/blueprints/`

These blueprints retain the parameter settings used to generate the
six final summary files. All six use the same folded one-dimensional
site-frequency spectrum.

Stairbuilder was used to generate the executable shell workflow for
each scenario. Generated shell files and complete working directories
are retained locally and are not distributed through GitHub.

The six scenarios combined:

- mutation rates of 1.01 × 10^-8, 4.2 × 10^-8, and 5.2 × 10^-8
  substitutions per site per generation;
- generation times of 50 and 100 years.

The retained final summary files are stored under:

`data/1.8.demography/stairway/`

Only the final summary files used for downstream visualization are
distributed through GitHub. The complete Stairway Plot working
directories are retained locally.

### 3. Combined visualization

`2.3.1_plot_stairway_plot.py` reads the six retained Stairway Plot
summary files and generates:

- a table reporting the minimum year and minimum median effective
  population size for each scenario;
- a combined demographic-history figure showing the six scenarios and
  their confidence intervals.

Generated files are written under:

`results/demography/`

## Interpretation

Because mutation rates, generation times, and the number of callable
sites are uncertain for this non-model oak, the demographic
reconstruction is treated as exploratory. Interpretation emphasizes
temporal patterns shared across scenarios rather than absolute
effective population-size estimates.

## Software

- easySFS
- Stairway Plot 2
- Python
- pandas
- matplotlib
