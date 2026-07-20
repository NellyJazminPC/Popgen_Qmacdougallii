# ============================================================
# Genetic diversity and differentiation in Quercus macdougallii
#
# Calculates:
# - Observed and expected heterozygosity
# - Inbreeding coefficients
# - Observed and private alleles by sampling site
# - Global and pairwise FST
#
# The analysis uses the 5,385 strictly biallelic SNPs retained
# automatically when converting the original VCF to genlight.
# ============================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(dartR)
  library(hierfstat)
  library(poppr)
})

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

vcf_file <- paste0(
  "../data/1.3.assembly_variant_calling/",
  "ref_gen_qrob_trim01_1_sorted.vcf"
)

metadata_file <- "../metadata/Qmacdougalli_79ind_.csv"

# ------------------------------------------------------------
# Load and validate data
# ------------------------------------------------------------

qmacd_vcf <- read.vcfR(
  vcf_file,
  verbose = FALSE
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

required_columns <- c(
  "ID",
  "SITE_NAME",
  "ZONE"
)

if (!all(required_columns %in% colnames(metadata))) {
  stop(
    "The metadata file must contain the columns: ",
    paste(required_columns, collapse = ", ")
  )
}

vcf_samples <- colnames(qmacd_vcf@gt)[-1]

if (!identical(vcf_samples, metadata$ID)) {
  stop(
    "The VCF sample order does not match the metadata file."
  )
}

# ------------------------------------------------------------
# Convert VCF to genlight
# ------------------------------------------------------------

qmacd_genlight <- vcfR2genlight(
  qmacd_vcf
)

ploidy(qmacd_genlight) <- 2

if (nLoc(qmacd_genlight) != 5385L) {
  stop(
    paste0(
      "Expected 5,385 biallelic SNPs after conversion, ",
      "but detected ",
      nLoc(qmacd_genlight),
      "."
    )
  )
}

cat(
  "\n===== SNP DATASET =====\n",
  "Individuals: ",
  nInd(qmacd_genlight),
  "\nBiallelic SNPs: ",
  nLoc(qmacd_genlight),
  "\n",
  sep = ""
)

# ------------------------------------------------------------
# Create genind objects
# ------------------------------------------------------------

site_order <- unique(
  metadata$SITE_NAME
)

zone_order <- unique(
  metadata$ZONE
)

site_genlight <- qmacd_genlight

pop(site_genlight) <- factor(
  metadata$SITE_NAME,
  levels = site_order
)

site_genind <- gl2gi(
  site_genlight,
  v = 0
)

zone_genlight <- qmacd_genlight

pop(zone_genlight) <- factor(
  metadata$ZONE,
  levels = zone_order
)

zone_genind <- gl2gi(
  zone_genlight,
  v = 0
)

# ------------------------------------------------------------
# Ho, He, and FIS by sampling site
# ------------------------------------------------------------

site_stats <- basic.stats(
  site_genind
)

site_Ho <- colMeans(
  site_stats$Ho,
  na.rm = TRUE
)

site_He <- colMeans(
  site_stats$Hs,
  na.rm = TRUE
)

site_summary <- data.frame(
  Scale = "Sampling site",
  Group = site_order,
  n_individuals = as.integer(
    table(
      factor(
        metadata$SITE_NAME,
        levels = site_order
      )
    )
  ),
  Ho = as.numeric(
    site_Ho[site_order]
  ),
  He = as.numeric(
    site_He[site_order]
  ),
  Fis = as.numeric(
    1 -
      site_Ho[site_order] /
      site_He[site_order]
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Observed and private alleles by sampling site
# ------------------------------------------------------------

site_genclone <- poppr::as.genclone(
  site_genind
)

private_allele_data <- poppr::private_alleles(
  site_genclone,
  report = "data.frame",
  level = "population"
)

allele_table <- tab(
  site_genind,
  NA.method = "zero"
)

site_summary$observed_alleles <- sapply(
  site_order,
  function(site) {

    individuals <- pop(site_genind) == site

    sum(
      colSums(
        allele_table[
          individuals,
          ,
          drop = FALSE
        ]
      ) > 0
    )
  }
)

site_summary$distinct_private_alleles <- sapply(
  site_order,
  function(site) {

    sum(
      private_allele_data$population == site &
        private_allele_data$count > 0
    )
  }
)

site_summary$private_allele_copies <- sapply(
  site_order,
  function(site) {

    selected <- (
      private_allele_data$population == site &
        private_allele_data$count > 0
    )

    sum(
      private_allele_data$count[selected]
    )
  }
)

site_summary$Fst <- NA_real_

# ------------------------------------------------------------
# Ho, He, and FIS by geographic zone
# ------------------------------------------------------------

zone_stats <- basic.stats(
  zone_genind
)

zone_Ho <- colMeans(
  zone_stats$Ho,
  na.rm = TRUE
)

zone_He <- colMeans(
  zone_stats$Hs,
  na.rm = TRUE
)

zone_summary <- data.frame(
  Scale = "Geographic zone",
  Group = zone_order,
  n_individuals = as.integer(
    table(
      factor(
        metadata$ZONE,
        levels = zone_order
      )
    )
  ),
  Ho = as.numeric(
    zone_Ho[zone_order]
  ),
  He = as.numeric(
    zone_He[zone_order]
  ),
  Fis = as.numeric(
    1 -
      zone_Ho[zone_order] /
      zone_He[zone_order]
  ),
  observed_alleles = NA_integer_,
  distinct_private_alleles = NA_integer_,
  private_allele_copies = NA_integer_,
  Fst = NA_real_,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Species-wide statistics
# ------------------------------------------------------------

overall_stats <- site_stats$overall

species_summary <- data.frame(
  Scale = "Species",
  Group = "All individuals",
  n_individuals = nrow(metadata),
  Ho = unname(
    overall_stats["Ho"]
  ),
  He = unname(
    overall_stats["Hs"]
  ),
  Fis = unname(
    overall_stats["Fis"]
  ),
  observed_alleles = sum(
    colSums(allele_table) > 0
  ),
  distinct_private_alleles = NA_integer_,
  private_allele_copies = NA_integer_,
  Fst = unname(
    overall_stats["Fst"]
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------
# Pairwise FST
# ------------------------------------------------------------

site_pairwise_Fst <- as.matrix(
  genet.dist(
    site_genind,
    method = "WC84"
  )
)

zone_pairwise_Fst <- as.matrix(
  genet.dist(
    zone_genind,
    method = "WC84"
  )
)

# ------------------------------------------------------------
# Combine and print summaries
# ------------------------------------------------------------

diversity_summary <- rbind(
  site_summary,
  zone_summary,
  species_summary
)

rownames(diversity_summary) <- NULL

cat(
  "\n===== GENETIC DIVERSITY SUMMARY =====\n"
)

print(
  diversity_summary,
  digits = 7,
  row.names = FALSE
)

cat(
  "\n===== GLOBAL FST =====\n"
)

print(
  unname(
    overall_stats["Fst"]
  ),
  digits = 7
)

cat(
  "\n===== PAIRWISE FST AMONG SAMPLING SITES =====\n"
)

print(
  round(
    site_pairwise_Fst,
    6
  )
)

cat(
  "\n===== NORTH-SOUTH FST =====\n"
)

print(
  round(
    zone_pairwise_Fst,
    6
  )
)

# ------------------------------------------------------------
# Export one consolidated results table
# ------------------------------------------------------------

dir.create(
  "../results",
  showWarnings = FALSE
)

write.csv(
  diversity_summary,
  "../results/genetic_diversity_summary.csv",
  row.names = FALSE
)
