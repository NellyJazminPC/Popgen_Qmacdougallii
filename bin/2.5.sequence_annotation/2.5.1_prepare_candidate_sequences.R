# -------------------------------------------------------------------------
# Script: 2.5.1_prepare_candidate_sequences.R
#
# Purpose:
#   Prepare and validate the manually curated sequences associated with
#   candidate outlier SNPs.
#
# Workflow context:
#   The outlier workflow identified 124 candidate SNPs representing
#   94 unique loci. Representative locus sequences were manually retrieved
#   from the ipyrad .loci assembly output and retained as a FASTA file.
#
# This script does not automate the manual sequence-retrieval step.
# It validates the curated FASTA files and joins the sequences to the
# candidate-SNP table.
# -------------------------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
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
  library(readxl)
  library(dplyr)
  library(stringr)
  library(writexl)
})

get_script_path <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) != 1) {
    stop(
      "The script path could not be determined. Run this file with Rscript.",
      call. = FALSE
    )
  }

  normalizePath(
    sub("^--file=", "", file_arg),
    winslash = "/",
    mustWork = TRUE
  )
}

read_fasta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]

  header_positions <- which(startsWith(lines, ">"))

  if (length(header_positions) == 0) {
    stop("No FASTA headers were found in: ", path, call. = FALSE)
  }

  sequence_end_positions <- c(
    header_positions[-1] - 1,
    length(lines)
  )

  sequences <- vapply(
    seq_along(header_positions),
    function(index) {
      start_position <- header_positions[index] + 1
      end_position <- sequence_end_positions[index]

      if (start_position > end_position) {
        return("")
      }

      paste0(lines[start_position:end_position], collapse = "")
    },
    character(1)
  )

  locus_names <- sub(
    "^>\\s*",
    "",
    lines[header_positions]
  )

  fasta_table <- data.frame(
    locus_name = locus_names,
    sequence = toupper(sequences),
    stringsAsFactors = FALSE
  )

  if (anyDuplicated(fasta_table$locus_name)) {
    duplicated_names <- unique(
      fasta_table$locus_name[duplicated(fasta_table$locus_name)]
    )

    stop(
      "Duplicated FASTA identifiers were found: ",
      paste(duplicated_names, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(!nzchar(fasta_table$sequence))) {
    stop(
      "One or more FASTA records contain an empty sequence.",
      call. = FALSE
    )
  }

  fasta_table
}

script_path <- get_script_path()
script_dir <- dirname(script_path)

repo_root <- normalizePath(
  file.path(script_dir, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(repo_root, "results")
data_dir <- file.path(
  repo_root,
  "data",
  "1.7.sequence_annotation"
)

candidate_snp_file <- file.path(
  results_dir,
  "consolidated_snps_results_with_shared_info_unique.xlsx"
)

candidate_sequence_fasta <- file.path(
  data_dir,
  "candidate_locus_sequences.fasta"
)

retained_blast_query_fasta <- file.path(
  data_dir,
  "candidate_locus_blast_queries.fasta"
)

sequence_output_file <- file.path(
  results_dir,
  "consolidated_snps_with_sequences.xlsx"
)

ambiguity_output_file <- file.path(
  results_dir,
  "consolidated_snps_with_sequences_and_ambiguities.xlsx"
)

candidate_locus_list_file <- file.path(
  results_dir,
  "candidate_locus_names.txt"
)

required_files <- c(
  candidate_snp_file,
  candidate_sequence_fasta,
  retained_blast_query_fasta
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required files were not found:\n",
    paste0("  - ", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

candidate_snps <- read_xlsx(candidate_snp_file)

if (!"locus_name" %in% names(candidate_snps)) {
  stop(
    "The candidate-SNP table does not contain a locus_name column.",
    call. = FALSE
  )
}

candidate_snps <- candidate_snps %>%
  mutate(
    locus_name_clean = sub("_pos.*", "", locus_name)
  )

candidate_loci <- unique(candidate_snps$locus_name_clean)
duplicate_snp_records <- sum(duplicated(candidate_snps$locus_name_clean))

cat("Candidate SNP records:", nrow(candidate_snps), "\n")
cat("Unique candidate loci:", length(candidate_loci), "\n")
cat(
  "Repeated locus records after removing _pos suffixes:",
  duplicate_snp_records,
  "\n"
)

if (nrow(candidate_snps) != 124) {
  warning(
    "The candidate-SNP table contains ",
    nrow(candidate_snps),
    " records instead of the 124 records in the retained analysis."
  )
}

if (length(candidate_loci) != 94) {
  warning(
    "The candidate-SNP table contains ",
    length(candidate_loci),
    " unique loci instead of the 94 loci in the retained analysis."
  )
}

writeLines(
  candidate_loci,
  candidate_locus_list_file
)

candidate_sequences <- read_fasta(candidate_sequence_fasta)
retained_blast_queries <- read_fasta(retained_blast_query_fasta)

cat(
  "Curated candidate sequences:",
  nrow(candidate_sequences),
  "\n"
)

cat(
  "Retained BLAST query sequences:",
  nrow(retained_blast_queries),
  "\n"
)

missing_candidate_sequences <- setdiff(
  candidate_loci,
  candidate_sequences$locus_name
)

unexpected_candidate_sequences <- setdiff(
  candidate_sequences$locus_name,
  candidate_loci
)

if (length(missing_candidate_sequences) > 0) {
  stop(
    "Candidate loci without a curated sequence: ",
    paste(missing_candidate_sequences, collapse = ", "),
    call. = FALSE
  )
}

if (length(unexpected_candidate_sequences) > 0) {
  stop(
    "Curated FASTA identifiers not found in the candidate table: ",
    paste(unexpected_candidate_sequences, collapse = ", "),
    call. = FALSE
  )
}

candidate_snps_with_sequences <- candidate_snps %>%
  left_join(
    candidate_sequences,
    by = c("locus_name_clean" = "locus_name")
  )

if (any(is.na(candidate_snps_with_sequences$sequence))) {
  stop(
    "One or more candidate SNPs could not be matched to a sequence.",
    call. = FALSE
  )
}

write_xlsx(
  candidate_snps_with_sequences,
  sequence_output_file
)

candidate_snps_with_ambiguities <- candidate_snps_with_sequences %>%
  mutate(
    amb_R = str_count(sequence, "R"),
    amb_Y = str_count(sequence, "Y"),
    amb_S = str_count(sequence, "S"),
    amb_W = str_count(sequence, "W"),
    amb_K = str_count(sequence, "K"),
    amb_M = str_count(sequence, "M"),
    total_ambiguities = amb_R + amb_Y + amb_S + amb_W + amb_K + amb_M
  )

write_xlsx(
  candidate_snps_with_ambiguities,
  ambiguity_output_file
)

generated_blast_queries <- candidate_sequences %>%
  mutate(
    sequence = sequence %>%
      str_replace_all("R", "A") %>%
      str_replace_all("Y", "C") %>%
      str_replace_all("S", "G") %>%
      str_replace_all("W", "T") %>%
      str_replace_all("K", "G") %>%
      str_replace_all("M", "A")
  )

blast_query_comparison <- generated_blast_queries %>%
  rename(generated_sequence = sequence) %>%
  full_join(
    retained_blast_queries %>%
      rename(retained_sequence = sequence),
    by = "locus_name"
  ) %>%
  mutate(
    sequences_match =
      generated_sequence == retained_sequence
  )

if (
  any(is.na(blast_query_comparison$sequences_match)) ||
  any(!blast_query_comparison$sequences_match)
) {
  mismatched_loci <- blast_query_comparison %>%
    filter(is.na(sequences_match) | !sequences_match) %>%
    pull(locus_name)

  stop(
    "The regenerated BLAST-query sequences do not match the retained ",
    "retained BLAST-query FASTA for the following loci: ",
    paste(mismatched_loci, collapse = ", "),
    call. = FALSE
  )
}

cat(
  "OK: regenerated BLAST-query sequences match the retained BLAST-query FASTA.\n"
)

cat("Candidate-locus list written to:\n", candidate_locus_list_file, "\n")
cat("Sequence table written to:\n", sequence_output_file, "\n")
cat("Ambiguity table written to:\n", ambiguity_output_file, "\n")
