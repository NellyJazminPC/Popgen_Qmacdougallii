# -------------------------------------------------------------------------
# Script: 2.6.1_prepare_outlier_variant_table.R
#
# Purpose:
#   Extract genotypes and allele information for candidate outlier SNPs,
#   translate genotype codes into nucleotide alleles, and add REF, ALT,
#   and mutation-type information to the retained annotation table.
#
# This script documents the table-preparation step used before generating
# the individual allele plots.
# -------------------------------------------------------------------------

required_packages <- c(
  "vcfR",
  "dplyr",
  "readxl",
  "writexl"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(vcfR)
  library(dplyr)
  library(readxl)
  library(writexl)
})

get_script_path <- function() {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(file_argument) != 1) {
    stop("Run this script with Rscript.", call. = FALSE)
  }

  normalizePath(
    sub("^--file=", "", file_argument),
    winslash = "/",
    mustWork = TRUE
  )
}

translate_genotype <- function(genotype, ref, alt) {
  if (
    is.na(genotype) ||
    genotype %in% c(".", "./.", ".|.")
  ) {
    return(NA_character_)
  }

  alternative_alleles <- strsplit(
    alt,
    split = ",",
    fixed = TRUE
  )[[1]]

  allele_lookup <- c(
    ref,
    alternative_alleles
  )

  allele_indices <- strsplit(
    genotype,
    split = "[/|]"
  )[[1]]

  if (any(allele_indices == ".")) {
    return(NA_character_)
  }

  allele_indices <- suppressWarnings(
    as.integer(allele_indices) + 1L
  )

  if (
    any(is.na(allele_indices)) ||
    any(allele_indices > length(allele_lookup))
  ) {
    return(NA_character_)
  }

  separator <- if (grepl("\\|", genotype)) "|" else "/"

  paste(
    allele_lookup[allele_indices],
    collapse = separator
  )
}

script_path <- get_script_path()
script_dir <- dirname(script_path)

repo_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(
  repo_root,
  "results"
)

vcf_file <- file.path(
  repo_root,
  "data",
  "1.3.assembly_variant_calling",
  "ref_gen_qrob_trim01_1_sorted.vcf"
)

annotation_file <- file.path(
  results_dir,
  "snps_outliers_with_blast_results_magnoliopsida_highlighted.xlsx"
)

genotype_output_file <- file.path(
  results_dir,
  "SNPs_outliers_variants_per_locus.xlsx"
)

annotation_output_file <- file.path(
  results_dir,
  "snps_outliers_with_blast_and_variants.xlsx"
)

required_files <- c(
  vcf_file,
  annotation_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Required files were not found:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

vcf <- read.vcfR(
  vcf_file,
  verbose = FALSE
)

annotation_table <- read_xlsx(
  annotation_file
)

if (!"locus_name" %in% names(annotation_table)) {
  stop(
    "The annotation table does not contain a locus_name column.",
    call. = FALSE
  )
}

loci_of_interest <- unique(
  annotation_table$locus_name
)

variant_information <- as.data.frame(
  vcf@fix,
  stringsAsFactors = FALSE
)

genotype_matrix <- extract.gt(
  vcf,
  element = "GT"
)

genotype_table <- data.frame(
  locus_name = rownames(genotype_matrix),
  genotype_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sample_columns <- colnames(
  genotype_matrix
)

variant_table <- genotype_table %>%
  left_join(
    variant_information %>%
      select(ID, REF, ALT),
    by = c("locus_name" = "ID")
  ) %>%
  filter(
    locus_name %in% loci_of_interest
  )

cat(
  "Candidate SNPs retained from the VCF:",
  nrow(variant_table),
  "\n"
)

variant_table <- variant_table %>%
  mutate(
    across(
      all_of(sample_columns),
      ~ mapply(
        translate_genotype,
        .,
        REF,
        ALT,
        USE.NAMES = FALSE
      )
    )
  ) %>%
  relocate(
    REF,
    ALT,
    .after = last_col()
  )

write_xlsx(
  variant_table,
  genotype_output_file
)

annotation_with_variants <- annotation_table %>%
  left_join(
    variant_table %>%
      select(
        locus_name,
        REF,
        ALT
      ),
    by = "locus_name"
  ) %>%
  mutate(
    mutation_type = case_when(
      paste0(REF, ALT) %in%
        c("AG", "GA", "CT", "TC") ~ "TRANSITION",
      nchar(REF) == 1 &
        nchar(ALT) == 1 ~ "TRANSVERSION",
      TRUE ~ NA_character_
    )
  )

write_xlsx(
  annotation_with_variants,
  annotation_output_file
)

cat(
  "Variant table written to:\n",
  genotype_output_file,
  "\n"
)

cat(
  "Annotation table with REF, ALT, and mutation type written to:\n",
  annotation_output_file,
  "\n"
)
