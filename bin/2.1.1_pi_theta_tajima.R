# ==========================================
# Tajima's D and Watterson's theta with GenoPop
# from a filtered VCF
# ==========================================
install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("Rsamtools")

library(GenoPop)
library(vcfR)
library(dplyr)

# -----------------------------
# Input files
# -----------------------------
vcf_file <- "../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf.gz"
vcf_file <- "test_500.vcf.gz"
meta_file <- "../metadata/Qmacdougalli_79ind_.csv"

# -----------------------------
# Metadata
# -----------------------------
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# Expected columns:
# ID, SITE_NAME, ZONE

# -----------------------------
# Basic VCF check
# -----------------------------
vcf_obj <- read.vcfR(vcf_file, verbose = FALSE)
vcf_samples <- colnames(vcf_obj@gt)[-1]

stopifnot(all(vcf_samples == meta$ID))

# -----------------------------
# IMPORTANT:
# seq_length should ideally be the number of callable sites
# or total evaluated sites for the dataset/group.
# Replace this when you have a better denominator.
# -----------------------------
seq_length_total <- 184703

# -----------------------------
# Total dataset
# -----------------------------
tajima_total <- TajimasD(
  vcf_path = vcf_file,
  seq_length = seq_length_total
)

theta_total <- WattersonsTheta(
  vcf_path = vcf_file,
  seq_length = seq_length_total
)

results_total <- data.frame(
  Group = "Total",
  Tajimas_D = tajima_total,
  Theta_Watterson = theta_total
)

print(results_total)

# -----------------------------
# Optional: export
# -----------------------------
write.csv(results_total, "../results/genopop_total_tajima_theta.csv", row.names = FALSE)