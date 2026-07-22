# -------------------------------------------------------------------------
# Script: 2.5.3_compile_blast_results.R
#
# Purpose:
#   Parse the retained NCBI BLAST text outputs, compare the three search
#   scopes explored during the analysis workflow, and generate the
#   consolidated annotation tables.
#
# Search scopes:
#   1. NCBI nt without a taxonomic restriction
#   2. Viridiplantae
#   3. Magnoliopsida
#
# Magnoliopsida was retained for final interpretation because it produced
# more taxonomically relevant matches and fewer nonspecific results.
# -------------------------------------------------------------------------

required_packages <- c(
  "readxl",
  "dplyr",
  "stringr",
  "tidyr",
  "writexl",
  "openxlsx"
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
  library(tidyr)
  library(writexl)
  library(openxlsx)
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

parse_blast_file <- function(path, source_name) {
  blast_content <- readLines(path, warn = FALSE)

  locus_name <- sub(
    "_blast\\.txt$",
    "",
    basename(path)
  )

  alignments_start <- grep(
    "^ALIGNMENTS",
    blast_content
  )

  if (length(alignments_start) == 0) {
    return(
      data.frame(
        locus_name = locus_name,
        best_hit = NA_character_,
        e_value = NA_character_,
        perc_identity = NA_character_,
        source = source_name,
        stringsAsFactors = FALSE
      )
    )
  }

  best_hit_index <- alignments_start[1] + 1

  best_hit <- if (best_hit_index <= length(blast_content)) {
    blast_content[best_hit_index] %>%
      str_remove("^>") %>%
      str_trim()
  } else {
    NA_character_
  }

  e_value_lines <- grep(
    "Expect =",
    blast_content,
    value = TRUE
  )

  e_value <- if (length(e_value_lines) > 0) {
    str_extract(
      e_value_lines[1],
      "(?<=Expect = )\\S+"
    )
  } else {
    NA_character_
  }

  identity_lines <- grep(
    "Identities =",
    blast_content,
    value = TRUE
  )

  perc_identity <- if (length(identity_lines) > 0) {
    str_extract(
      identity_lines[1],
      "(?<=\\()\\d+%"
    )
  } else {
    NA_character_
  }

  data.frame(
    locus_name = locus_name,
    best_hit = best_hit,
    e_value = e_value,
    perc_identity = perc_identity,
    source = source_name,
    stringsAsFactors = FALSE
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

sequence_table_file <- file.path(
  results_dir,
  "consolidated_snps_with_sequences_and_ambiguities.xlsx"
)

blast_dirs <- c(
  all_db = file.path(
    results_dir,
    "blast_results_all_db"
  ),
  Magnoliopsida = file.path(
    results_dir,
    "blast_results_Magnoliopsida"
  ),
  Viridiplantae = file.path(
    results_dir,
    "blast_results_Viridiplantae"
  )
)

combined_output_file <- file.path(
  results_dir,
  "snps_outliers_with_blast_results_combined.xlsx"
)

magnoliopsida_output_file <- file.path(
  results_dir,
  "snps_outliers_with_blast_results_magnoliopsida_highlighted.xlsx"
)

if (!file.exists(sequence_table_file)) {
  stop(
    "The sequence table was not found:\n",
    sequence_table_file,
    call. = FALSE
  )
}

missing_blast_dirs <- blast_dirs[
  !dir.exists(blast_dirs)
]

if (length(missing_blast_dirs) > 0) {
  stop(
    "The following BLAST result directories were not found:\n",
    paste0("  - ", missing_blast_dirs, collapse = "\n"),
    call. = FALSE
  )
}

candidate_snps_with_sequences <- read_xlsx(
  sequence_table_file
)

required_columns <- c(
  "locus_name",
  "locus_name_clean"
)

missing_columns <- setdiff(
  required_columns,
  names(candidate_snps_with_sequences)
)

if (length(missing_columns) > 0) {
  stop(
    "The sequence table is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

blast_tables <- lapply(
  names(blast_dirs),
  function(source_name) {
    blast_dir <- blast_dirs[[source_name]]

    blast_files <- list.files(
      blast_dir,
      pattern = "_blast\\.txt$",
      full.names = TRUE
    )

    cat(
      source_name,
      "BLAST files:",
      length(blast_files),
      "\n"
    )

    if (length(blast_files) == 0) {
      stop(
        "No BLAST files were found in: ",
        blast_dir,
        call. = FALSE
      )
    }

    if (length(blast_files) != 94) {
      warning(
        source_name,
        " contains ",
        length(blast_files),
        " BLAST files instead of the 94 files in the retained analysis."
      )
    }

    parsed_table <- bind_rows(
      lapply(
        blast_files,
        parse_blast_file,
        source_name = source_name
      )
    )

    if (anyDuplicated(parsed_table$locus_name)) {
      duplicated_loci <- unique(
        parsed_table$locus_name[
          duplicated(parsed_table$locus_name)
        ]
      )

      stop(
        "Duplicated BLAST results were found for ",
        source_name,
        ": ",
        paste(duplicated_loci, collapse = ", "),
        call. = FALSE
      )
    }

    parsed_table
  }
)

blast_summary <- bind_rows(
  blast_tables
)

missing_hits <- blast_summary %>%
  filter(is.na(e_value)) %>%
  select(source, locus_name)

cat(
  "BLAST records without a parsed hit:",
  nrow(missing_hits),
  "\n"
)

if (nrow(missing_hits) > 0) {
  print(missing_hits)
}

blast_summary_wide <- blast_summary %>%
  pivot_wider(
    names_from = source,
    values_from = c(
      best_hit,
      e_value,
      perc_identity
    ),
    names_sep = "_"
  ) %>%
  mutate(
    is_consistent = if_else(
      !is.na(best_hit_all_db) &
        !is.na(best_hit_Magnoliopsida) &
        !is.na(best_hit_Viridiplantae) &
        best_hit_all_db == best_hit_Magnoliopsida &
        best_hit_all_db == best_hit_Viridiplantae,
      "Yes",
      "No"
    )
  )

candidate_snps_with_blast <- candidate_snps_with_sequences %>%
  left_join(
    blast_summary_wide,
    by = c(
      "locus_name_clean" = "locus_name"
    )
  )

write_xlsx(
  candidate_snps_with_blast,
  combined_output_file
)

magnoliopsida_results <- candidate_snps_with_blast %>%
  select(
    -any_of(
      c(
        "best_hit_all_db",
        "best_hit_Viridiplantae",
        "e_value_all_db",
        "e_value_Viridiplantae",
        "perc_identity_all_db",
        "perc_identity_Viridiplantae",
        "is_consistent"
      )
    )
  ) %>%
  mutate(
    e_value_Magnoliopsida = suppressWarnings(
      as.numeric(e_value_Magnoliopsida)
    )
  )

workbook <- createWorkbook()

addWorksheet(
  workbook,
  "Magnoliopsida Results"
)

writeData(
  workbook,
  "Magnoliopsida Results",
  magnoliopsida_results
)

freezePane(
  workbook,
  "Magnoliopsida Results",
  firstRow = TRUE
)

green_style <- createStyle(
  fontColour = "#006400",
  bgFill = "#C6EFCE"
)

e_value_column <- which(
  names(magnoliopsida_results) ==
    "e_value_Magnoliopsida"
)

if (
  length(e_value_column) == 1 &&
  nrow(magnoliopsida_results) > 0
) {
  conditionalFormatting(
    workbook,
    sheet = "Magnoliopsida Results",
    cols = e_value_column,
    rows = 2:(nrow(magnoliopsida_results) + 1),
    rule = "<=0.05",
    style = green_style
  )
}

saveWorkbook(
  workbook,
  magnoliopsida_output_file,
  overwrite = TRUE
)

cat(
  "Combined BLAST table written to:\n",
  combined_output_file,
  "\n"
)

cat(
  "Magnoliopsida annotation table written to:\n",
  magnoliopsida_output_file,
  "\n"
)
