# General pipeline

This pipeline is organized into numbered steps, each corresponding to a script or folder in the `bin` directory. Below is a brief description of each stage, updated according to the current contents:

1. **1.0.quality_analysis.sh**  
   Initial quality assessment of raw sequencing data.

2. **1.1.filter_trimmomatic.sh**  
   Quality filtering and trimming of raw reads using Trimmomatic to remove low-quality bases and adapters.

3. **1.2.post-filter_quality_analysis.sh**  
   Quality assessment of reads after filtering to ensure data integrity for downstream analyses.

4. **1.3.assembly_variant_calling_ipyrad/**  
   De novo assembly and variant calling using ipyrad to generate loci and genotype files for each individual.

5. **1.4.populations_stacks/**  
   Population analysis with Stacks (`populations`) to calculate genetic diversity statistics and export files in various formats (VCF, Genepop, etc.).

6. **1.5.structure/**  
   Preparation and analysis of population structure using multiple approaches:
   - `admixture_qmacd.sh`: Runs ADMIXTURE for ancestry estimation.
   - `faststructure_qmacd.sh`: Runs fastStructure for population structure inference.
   - `convert2plinkformat.sh`: Converts genotype data to PLINK format for compatibility with various tools.
   - `workflow_pop_structure_analysis.md`: Documentation of the population structure analysis workflow.
   - `faststructure_convert.spid`: Configuration file for data conversion.

7. **1.6.snps_outliers/**  
   Detection of candidate loci under selection (outliers) using tools such as BayeScan, and preparation of files for selection analysis.

8. **2.1.genetic_diversity.R**  
   R script for calculating genetic diversity statistics across populations.

9. **2.2.pop_structure.R**  
   R script for analyzing population structure using various methods.

10. **2.3.admixture_faststructure_plots.R**  
    R script for generating plots of admixture and fastStructure results.

11. **2.4.snps_outliers_PCAdapt_Bayescan_FST.R**  
    R script for identifying outlier SNPs using PCAdapt, BayeScan, and FST approaches.

12. **2.5.seq_search_blast.sh / 2.5.snps_outliers_sequences.R**  
    Scripts for searching SNP sequences with BLAST and further analysis of outlier SNP sequences.

13. **2.6.snps_outliers_freq.R**  
    R script for analyzing allele frequencies of outlier SNPs.

14. **2.7.Ne.R**  
    R script for estimating effective population size (Ne) from filtered genotype data.

15. **3.1.demography.md**  
    Documentation or summary of demographic analyses.

Each of these steps is automated with bash or R scripts, ensuring reproducibility and an organized workflow for the analyses.