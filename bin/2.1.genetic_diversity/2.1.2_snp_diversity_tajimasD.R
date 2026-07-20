# ============================================================
# SNP-based diversity statistics for Quercus macdougallii
#
# Calculates nucleotide diversity (pi), Watterson's theta,
# and Tajima's D from strictly biallelic SNPs.
# ============================================================

library(vcfR)

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

vcf_file <- "../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf"
metadata_file <- "../metadata/Qmacdougalli_79ind_.csv"

# ------------------------------------------------------------
# Load data
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
# Extract and standardize genotypes
# ------------------------------------------------------------

gt <- extract.gt(
  qmacd_vcf,
  element = "GT"
)

gt <- sub(":.*$", "", gt)
gt <- gsub("\\|", "/", gt)

# ------------------------------------------------------------
# Retain strictly biallelic SNPs
# ------------------------------------------------------------

# Multiallelic ALT definitions in the VCF.
multiallelic_alt <- grepl(
  ",",
  qmacd_vcf@fix[, "ALT"]
)

# Genotypes containing alleles numbered 2 or higher.
multiallelic_gt <- apply(
  gt,
  1,
  function(x) {
    x <- x[!is.na(x)]

    any(
      grepl(
        "(^|/)([2-9]|[1-9][0-9]+)($|/)",
        x
      )
    )
  }
)

multiallelic_loci <- (
  multiallelic_alt |
    multiallelic_gt
)

biallelic_loci <- !multiallelic_loci

gt_biallelic <- gt[
  biallelic_loci,
  ,
  drop = FALSE
]

filter_summary <- data.frame(
  total_SNPs = nrow(gt),
  excluded_multiallelic_SNPs = sum(multiallelic_loci),
  retained_biallelic_SNPs = nrow(gt_biallelic)
)

print(filter_summary)

if (
  nrow(gt) != 5426L ||
    sum(multiallelic_loci) != 41L ||
    nrow(gt_biallelic) != 5385L
) {
  stop(
    paste0(
      "Unexpected SNP counts. Expected 5,426 total SNPs, ",
      "41 multiallelic SNPs, and 5,385 retained ",
      "biallelic SNPs."
    )
  )
}

# ------------------------------------------------------------
# Convert genotypes to alternative-allele counts
# ------------------------------------------------------------

alt_counts <- matrix(
  NA_real_,
  nrow = nrow(gt_biallelic),
  ncol = ncol(gt_biallelic),
  dimnames = dimnames(gt_biallelic)
)

alt_counts[gt_biallelic == "0/0"] <- 0

alt_counts[
  gt_biallelic %in% c("0/1", "1/0")
] <- 1

alt_counts[gt_biallelic == "1/1"] <- 2

if (any(is.na(alt_counts))) {
  stop(
    paste0(
      "Missing or unrecognized genotypes remain in the ",
      "biallelic SNP matrix."
    )
  )
}

# ------------------------------------------------------------
# Diversity-statistics function
# ------------------------------------------------------------

calculate_snp_diversity <- function(
  mat,
  min_individuals_for_D = 4
) {

  n_individuals <- ncol(mat)
  n_sequences <- 2 * n_individuals
  loci_evaluated <- nrow(mat)

  p_alt <- rowSums(mat) / (2 * n_individuals)

  segregating <- (
    p_alt > 0 &
      p_alt < 1
  )

  segregating_sites <- sum(segregating)

  # Sample-size-corrected nucleotide diversity.
  pi_per_locus <- (
    n_sequences /
      (n_sequences - 1)
  ) * 2 * p_alt * (1 - p_alt)

  pi <- mean(pi_per_locus)

  sequence_indices <- seq_len(
    n_sequences - 1
  )

  a1 <- sum(
    1 / sequence_indices
  )

  a2 <- sum(
    1 / (sequence_indices^2)
  )

  thetaW <- (
    segregating_sites /
      a1
  ) / loci_evaluated

  b1 <- (
    n_sequences + 1
  ) / (
    3 * (n_sequences - 1)
  )

  b2 <- (
    2 *
      (
        n_sequences^2 +
          n_sequences +
          3
      )
  ) / (
    9 *
      n_sequences *
      (n_sequences - 1)
  )

  c1 <- b1 - (1 / a1)

  c2 <- (
    b2 -
      (
        (n_sequences + 2) /
          (a1 * n_sequences)
      ) +
      (
        a2 /
          (a1^2)
      )
  )

  e1 <- c1 / a1

  e2 <- c2 / (
    (a1^2) + a2
  )

  nucleotide_differences <- (
    pi * loci_evaluated
  )

  if (
    segregating_sites > 1 &&
      n_individuals >= min_individuals_for_D
  ) {
    tajimas_D <- (
      nucleotide_differences -
        (segregating_sites / a1)
    ) / sqrt(
      (e1 * segregating_sites) +
        (
          e2 *
            segregating_sites *
            (segregating_sites - 1)
        )
    )
  } else {
    tajimas_D <- NA_real_
  }

  data.frame(
    n_individuals = n_individuals,
    n_sequences = n_sequences,
    loci_evaluated = loci_evaluated,
    segregating_sites = segregating_sites,
    pi = pi,
    thetaW = thetaW,
    Tajimas_D = tajimas_D
  )
}

# ------------------------------------------------------------
# Sampling-site statistics
# ------------------------------------------------------------

site_order <- unique(
  metadata$SITE_NAME
)

site_results <- do.call(
  rbind,
  lapply(
    site_order,
    function(site) {

      ids <- metadata$ID[
        metadata$SITE_NAME == site
      ]

      result <- calculate_snp_diversity(
        alt_counts[
          ,
          ids,
          drop = FALSE
        ]
      )

      result$Scale <- "Sampling site"
      result$Group <- site

      result
    }
  )
)

# ------------------------------------------------------------
# Geographic-zone statistics
# ------------------------------------------------------------

zone_order <- unique(
  metadata$ZONE
)

zone_results <- do.call(
  rbind,
  lapply(
    zone_order,
    function(zone) {

      ids <- metadata$ID[
        metadata$ZONE == zone
      ]

      result <- calculate_snp_diversity(
        alt_counts[
          ,
          ids,
          drop = FALSE
        ]
      )

      result$Scale <- "Geographic zone"
      result$Group <- zone

      result
    }
  )
)

# ------------------------------------------------------------
# Species-wide statistics
# ------------------------------------------------------------

species_results <- calculate_snp_diversity(
  alt_counts
)

species_results$Scale <- "Species"
species_results$Group <- "All individuals"

# ------------------------------------------------------------
# Combine, print, and export results
# ------------------------------------------------------------

diversity_results <- rbind(
  site_results,
  zone_results,
  species_results
)

diversity_results <- diversity_results[
  ,
  c(
    "Scale",
    "Group",
    "n_individuals",
    "n_sequences",
    "loci_evaluated",
    "segregating_sites",
    "pi",
    "thetaW",
    "Tajimas_D"
  )
]

rownames(diversity_results) <- NULL

print(diversity_results)

dir.create(
  "../results",
  showWarnings = FALSE
)

write.csv(
  diversity_results,
  "../results/snp_diversity_biallelic_5385.csv",
  row.names = FALSE
)
