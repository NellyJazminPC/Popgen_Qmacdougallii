# -------------------------------------------------------------------------
# Script: 2.6.2_plot_outlier_allele_frequencies.R
#
# Purpose:
#   Generate individual allele plots for the 11 coding candidate SNPs.
#
# Two representations are generated:
#
#   1. Mean allele dosage using tab(..., freq = FALSE), on a 0-2 scale.
#   2. Allele frequency using tab(..., freq = TRUE), on a 0-1 scale.
#
# The calculation follows the original workflow based on:
#
#   temp <- seploc(qmacd_genind)
#   tab(temp[[snp_name]])
#   tapply(..., pop(qmacd_genind), mean)
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(dartR)
})

vcf_file <- file.path(
  "data",
  "1.3.assembly_variant_calling",
  "ref_gen_qrob_trim01_1_sorted.vcf"
)

metadata_file <- file.path(
  "metadata",
  "Qmacdougalli_79ind_.csv"
)

selected_snps_file <- file.path(
  "data",
  "1.8.outlier_allele_frequencies",
  "manuscript_outlier_snps_for_plotting.csv"
)

dosage_output_dir <- file.path(
  "results",
  "plots_snps_outliers_dosage_scale"
)

frequency_output_dir <- file.path(
  "results",
  "plots_snps_outliers_frequency_scale"
)

required_files <- c(
  vcf_file,
  metadata_file,
  selected_snps_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required files were not found. Run this script from the repository root:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

dir.create(
  dosage_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  frequency_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

qmacd_vcf <- read.vcfR(
  vcf_file,
  verbose = FALSE
)

metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

selected_snps <- read.csv(
  selected_snps_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "locus_name",
  "selection_metric",
  "metric_value",
  "ref_allele",
  "alt_allele",
  "best_hit_magnoliopsida"
)

missing_columns <- setdiff(
  required_columns,
  names(selected_snps)
)

if (length(missing_columns) > 0) {
  stop(
    "The selected-SNP table is missing columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(selected_snps) != 11) {
  stop(
    "Expected 11 selected SNPs, found ",
    nrow(selected_snps),
    ".",
    call. = FALSE
  )
}

vcf_sample_ids <- colnames(qmacd_vcf@gt)[-1]

if (!identical(vcf_sample_ids, metadata$ID)) {
  stop(
    "The VCF sample order does not match the metadata ID column.",
    call. = FALSE
  )
}

# Convert the VCF to genlight. The 41 multiallelic loci are omitted
# because genlight supports only biallelic loci.
qmacd_genlight <- vcfR2genlight(
  qmacd_vcf
)

ploidy(qmacd_genlight) <- 2

# SITE_NAME contains the nine ordered sampling-site labels:
# 1CZ, 2MT, 3MC, 4MB, 5CY, 6LS, 7PZ, 8CR, and 9IT.
pop(qmacd_genlight) <- metadata$SITE_NAME

qmacd_genind <- gl2gi(
  qmacd_genlight,
  v = 0
)

pop(qmacd_genind) <- metadata$SITE_NAME

# Preserve the original workflow by splitting the complete genind object
# into one genind object per locus.
temp <- seploc(
  qmacd_genind
)

missing_loci <- setdiff(
  selected_snps$locus_name,
  names(temp)
)

if (length(missing_loci) > 0) {
  stop(
    "Selected SNPs were not found in temp: ",
    paste(missing_loci, collapse = ", "),
    call. = FALSE
  )
}

site_axis_labels <- c(
  "CZ", "MT", "MC", "MB", "CY",
  "LS", "PZ", "CR", "IT"
)

calculate_site_means <- function(snp_name, use_frequency) {
  # This preserves the original two-step calculation.
  snp_data <- tab(
    temp[[snp_name]],
    freq = use_frequency
  )

  site_means <- apply(
    snp_data,
    2,
    function(allele_values) {
      tapply(
        allele_values,
        pop(qmacd_genind),
        mean,
        na.rm = TRUE
      )
    }
  )

  site_means <- as.matrix(
    site_means
  )

  if (nrow(site_means) != length(site_axis_labels)) {
    stop(
      "Unexpected number of sampling sites for SNP ",
      snp_name,
      ": ",
      nrow(site_means),
      ".",
      call. = FALSE
    )
  }

  if (all(is.na(site_means))) {
    stop(
      "All calculated values are missing for SNP ",
      snp_name,
      ".",
      call. = FALSE
    )
  }

  site_means
}

save_individual_plot <- function(
  plot_values,
  output_file,
  plot_title,
  allele_symbols,
  best_hit,
  y_axis_label,
  best_hit_cex
) {
  png(
    filename = output_file,
    width = 2200,
    height = 1200,
    res = 300
  )

  par(
    mar = c(6, 5, 4, 2) + 0.1
  )

  matplot(
    plot_values,
    type = "b",
    pch = allele_symbols,
    xlab = "SITE",
    ylab = y_axis_label,
    main = plot_title,
    xaxt = "n",
    cex = 1.1
  )

  axis(
    side = 1,
    at = seq_along(site_axis_labels),
    labels = site_axis_labels
  )

  mtext(
    text = paste(
      "Best hit (Magnoliopsida):",
      best_hit
    ),
    side = 1,
    line = 5,
    cex = best_hit_cex,
    col = "blue"
  )

  dev.off()
}

for (row_index in seq_len(nrow(selected_snps))) {
  snp <- selected_snps[row_index, ]

  snp_name <- snp$locus_name
  snp_value <- round(
    as.numeric(snp$metric_value),
    6
  )

  allele_symbols <- c(
    snp$ref_allele,
    snp$alt_allele
  )

  if (snp$selection_metric == "adjusted_p_value") {
    plot_prefix <- "allele_frequency_pcadapt_"

    plot_title <- paste(
      "p-value adj",
      as.character(snp_value),
      snp_name
    )

    best_hit_cex <- 0.6
  } else if (snp$selection_metric == "fst") {
    plot_prefix <- "allele_frequency_fst_"

    plot_title <- paste(
      "FST value",
      as.character(snp_value),
      snp_name
    )

    best_hit_cex <- 0.5
  } else {
    stop(
      "Unknown selection metric for ",
      snp_name,
      ": ",
      snp$selection_metric,
      call. = FALSE
    )
  }

  # Original representation: mean number of allele copies per
  # diploid individual.
  dosage_values <- calculate_site_means(
    snp_name,
    use_frequency = FALSE
  )

  # Conventional population allele-frequency representation.
  frequency_values <- calculate_site_means(
    snp_name,
    use_frequency = TRUE
  )

  # For diploid data, frequencies should equal dosage values divided by 2.
  if (
    !isTRUE(
      all.equal(
        unname(frequency_values),
        unname(dosage_values / 2),
        tolerance = 1e-8
      )
    )
  ) {
    stop(
      "Dosage and frequency values are inconsistent for SNP ",
      snp_name,
      ".",
      call. = FALSE
    )
  }

  dosage_output_file <- file.path(
    dosage_output_dir,
    paste0(
      plot_prefix,
      snp_name,
      ".png"
    )
  )

  frequency_output_file <- file.path(
    frequency_output_dir,
    paste0(
      plot_prefix,
      snp_name,
      ".png"
    )
  )

  save_individual_plot(
    plot_values = dosage_values,
    output_file = dosage_output_file,
    plot_title = plot_title,
    allele_symbols = allele_symbols,
    best_hit = snp$best_hit_magnoliopsida,
    y_axis_label = "Mean allele dosage",
    best_hit_cex = best_hit_cex
  )

  save_individual_plot(
    plot_values = frequency_values,
    output_file = frequency_output_file,
    plot_title = plot_title,
    allele_symbols = allele_symbols,
    best_hit = snp$best_hit_magnoliopsida,
    y_axis_label = "Allele frequency",
    best_hit_cex = best_hit_cex
  )

  cat(
    "Generated dosage and frequency plots for:",
    snp_name,
    "\n"
  )
}

cat(
  "\nDosage-scale plots:\n",
  dosage_output_dir,
  "\n"
)

cat(
  "Frequency-scale plots:\n",
  frequency_output_dir,
  "\n"
)
