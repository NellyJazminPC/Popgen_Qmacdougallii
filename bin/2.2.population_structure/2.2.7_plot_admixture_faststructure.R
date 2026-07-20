#!/usr/bin/env Rscript

# Plot ADMIXTURE and fastStructure ancestry results for Quercus macdougallii.
#
# Public outputs reproduced by this script:
#   1. ADMIXTURE cross-validation error for K = 1-10.
#   2. ADMIXTURE ancestry proportions for K = 1 and K = 2.
#   3. fastStructure ancestry proportions for K = 1 under the simple prior.
#   4. fastStructure ancestry proportions for K = 2 under the logistic prior.
#   5. A combined three-panel figure corresponding to panels A-C of Figure 6.
#
# The geographic map in panel D is not generated here because the coordinates
# of this threatened microendemic species are not included in the public
# metadata file.
#
# Run from the repository root or any subdirectory:
#
#   Rscript bin/2.2.population_structure/2.2.7_plot_admixture_faststructure.R

required_packages <- c("ggplot2", "dplyr", "tidyr", "patchwork")

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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

find_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    has_git <- dir.exists(file.path(current, ".git"))
    has_data <- dir.exists(file.path(current, "data"))
    has_metadata <- dir.exists(file.path(current, "metadata"))

    if (has_git && has_data && has_metadata) {
      return(current)
    }

    parent <- dirname(current)

    if (identical(parent, current)) {
      stop(
        "Repository root not found. Run this script from the repository ",
        "root or one of its subdirectories.",
        call. = FALSE
      )
    }

    current <- parent
  }
}

read_ancestry <- function(path, fam_ids, cluster_names) {
  if (!file.exists(path)) {
    stop("Missing ancestry file: ", path, call. = FALSE)
  }

  ancestry <- read.table(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(ancestry) != length(fam_ids)) {
    stop(
      "Unexpected number of rows in ", basename(path), ": ",
      nrow(ancestry), " instead of ", length(fam_ids), ".",
      call. = FALSE
    )
  }

  if (ncol(ancestry) != length(cluster_names)) {
    stop(
      "Unexpected number of ancestry columns in ", basename(path), ": ",
      ncol(ancestry), " instead of ", length(cluster_names), ".",
      call. = FALSE
    )
  }

  colnames(ancestry) <- cluster_names
  ancestry$ID <- fam_ids

  ancestry
}

prepare_long_ancestry <- function(ancestry, metadata, cluster_names) {
  joined <- ancestry |>
    left_join(metadata, by = "ID")

  if (anyNA(joined$SITE_NAME) || anyNA(joined$NUM_SAMPLE)) {
    missing_ids <- joined$ID[
      is.na(joined$SITE_NAME) | is.na(joined$NUM_SAMPLE)
    ]

    stop(
      "Metadata could not be matched for: ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }

  joined |>
    mutate(
      NUM_SAMPLE = as.numeric(NUM_SAMPLE),
      individual = factor(ID, levels = ID[order(NUM_SAMPLE)])
    ) |>
    pivot_longer(
      cols = all_of(cluster_names),
      names_to = "cluster",
      values_to = "ancestry"
    ) |>
    mutate(
      cluster = factor(cluster, levels = cluster_names),
      SITE_NAME = factor(
        SITE_NAME,
        levels = unique(SITE_NAME[order(NUM_SAMPLE)])
      ),
      ZONE = factor(
        ZONE,
        levels = c("NORTH", "SOUTH")
      )
    )
}

make_ancestry_plot <- function(data, title, cluster_colors) {
  ggplot(
    data,
    aes(x = individual, y = ancestry, fill = cluster)
  ) +
    geom_col(width = 1) +
    facet_grid(
      cols = vars(ZONE, SITE_NAME),
      scales = "free_x",
      space = "free_x"
    ) +
    scale_fill_manual(
      values = cluster_colors,
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      x = NULL,
      y = "Ancestry proportion",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.spacing.x = grid::unit(0.08, "lines"),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(size = 8),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
}

repo_root <- find_repo_root()

metadata_path <- file.path(
  repo_root,
  "metadata",
  "Qmacdougalli_79ind_.csv"
)

fam_path <- file.path(
  repo_root,
  "data",
  "structure_formats",
  "qmacd_ref_gen_rob.fam"
)

admixture_dir <- file.path(
  repo_root,
  "data",
  "1.4.population_structure",
  "admixture_output"
)

faststructure_dir <- file.path(
  repo_root,
  "data",
  "1.4.population_structure",
  "faststructure_output"
)

output_dir <- file.path(
  repo_root,
  "results",
  "population_structure"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(metadata_path)) {
  stop("Missing metadata file: ", metadata_path, call. = FALSE)
}

if (!file.exists(fam_path)) {
  stop("Missing PLINK FAM file: ", fam_path, call. = FALSE)
}

metadata <- read.csv(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_columns <- c(
  "ID",
  "SITE_NAME",
  "ZONE",
  "NUM_SAMPLE"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(metadata$ID)) {
  stop("The metadata ID column contains duplicates.", call. = FALSE)
}

fam <- read.table(
  fam_path,
  header = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(fam) < 2) {
  stop("Unexpected PLINK FAM format.", call. = FALSE)
}

fam_ids <- fam[[2]]

if (anyDuplicated(fam_ids)) {
  stop("The PLINK FAM file contains duplicate sample IDs.", call. = FALSE)
}

if (!setequal(fam_ids, metadata$ID)) {
  only_fam <- setdiff(fam_ids, metadata$ID)
  only_metadata <- setdiff(metadata$ID, fam_ids)

  stop(
    "Sample IDs differ between the PLINK FAM file and metadata.\n",
    "Only in FAM: ", paste(only_fam, collapse = ", "), "\n",
    "Only in metadata: ", paste(only_metadata, collapse = ", "),
    call. = FALSE
  )
}

# ADMIXTURE cross-validation error

choose_k_path <- file.path(admixture_dir, "chooseK.txt")

if (!file.exists(choose_k_path)) {
  stop("Missing ADMIXTURE chooseK file: ", choose_k_path, call. = FALSE)
}

choose_k_lines <- readLines(choose_k_path, warn = FALSE)

choose_k <- data.frame(
  K = as.integer(
    sub(
      ".*CV error \\(K=([0-9]+)\\).*",
      "\\1",
      choose_k_lines
    )
  ),
  CV_error = as.numeric(
    sub(
      ".*:\\s*([0-9.]+)\\s*$",
      "\\1",
      choose_k_lines
    )
  )
) |>
  arrange(K)

if (anyNA(choose_k$K) || anyNA(choose_k$CV_error)) {
  stop("Could not parse ADMIXTURE cross-validation values.", call. = FALSE)
}

cv_plot <- ggplot(choose_k, aes(x = K, y = CV_error)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = choose_k$K) +
  labs(
    x = "K",
    y = "Cross-validation error"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = file.path(output_dir, "admixture_cv_error.png"),
  plot = cv_plot,
  width = 7,
  height = 5,
  dpi = 300
)

# Cluster labels are arbitrary. The column assignments below preserve the
# orientation used in the manuscript figures.

cluster_colors <- c(
  K1 = "#0072B2",
  K2 = "#E69F00"
)

admixture_k1 <- read_ancestry(
  file.path(
    admixture_dir,
    "qmacd_ref_gen_rob.1.Q"
  ),
  fam_ids,
  cluster_names = "K1"
)

admixture_k2 <- read_ancestry(
  file.path(
    admixture_dir,
    "qmacd_ref_gen_rob.2.Q"
  ),
  fam_ids,
  cluster_names = c("K2", "K1")
)

faststructure_simple_k1 <- read_ancestry(
  file.path(
    faststructure_dir,
    "qmacd_ref_gen_rob.simple.1.meanQ"
  ),
  fam_ids,
  cluster_names = "K1"
)

faststructure_logistic_k2 <- read_ancestry(
  file.path(
    faststructure_dir,
    "qmacd_ref_gen_rob.logistic.2.meanQ"
  ),
  fam_ids,
  cluster_names = c("K2", "K1")
)

admixture_k1_long <- prepare_long_ancestry(
  admixture_k1,
  metadata,
  cluster_names = "K1"
)

admixture_k2_long <- prepare_long_ancestry(
  admixture_k2,
  metadata,
  cluster_names = c("K1", "K2")
)

faststructure_simple_k1_long <- prepare_long_ancestry(
  faststructure_simple_k1,
  metadata,
  cluster_names = "K1"
)

faststructure_logistic_k2_long <- prepare_long_ancestry(
  faststructure_logistic_k2,
  metadata,
  cluster_names = c("K1", "K2")
)

plot_admixture_k1 <- make_ancestry_plot(
  admixture_k1_long,
  "ADMIXTURE, K = 1",
  cluster_colors["K1"]
)

plot_faststructure_k1 <- make_ancestry_plot(
  faststructure_simple_k1_long,
  "fastStructure, K = 1, simple prior",
  cluster_colors["K1"]
)

plot_admixture_k2 <- make_ancestry_plot(
  admixture_k2_long,
  "ADMIXTURE, K = 2",
  cluster_colors
)

plot_faststructure_k2 <- make_ancestry_plot(
  faststructure_logistic_k2_long,
  "fastStructure, K = 2, logistic prior",
  cluster_colors
)

panel_a <- (
  plot_admixture_k1 +
    labs(title = "ADMIXTURE")
) / (
  plot_faststructure_k1 +
    labs(title = "fastStructure, simple prior")
) +
  plot_annotation(title = "A) K = 1")

plot_admixture_k2_panel <- plot_admixture_k2 +
  labs(title = "B) ADMIXTURE, K = 2")

plot_faststructure_k2_panel <- plot_faststructure_k2 +
  labs(title = "C) fastStructure, K = 2, logistic prior")

combined_figure <- panel_a /
  plot_admixture_k2_panel /
  plot_faststructure_k2_panel +
  plot_layout(heights = c(2, 1, 1))

ggsave(
  filename = file.path(
    output_dir,
    "figure6_population_structure_panels_A-C.png"
  ),
  plot = combined_figure,
  width = 12,
  height = 14,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "admixture_K1.png"),
  plot = plot_admixture_k1,
  width = 12,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "faststructure_simple_K1.png"),
  plot = plot_faststructure_k1,
  width = 12,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "admixture_K2.png"),
  plot = plot_admixture_k2,
  width = 12,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "faststructure_logistic_K2.png"),
  plot = plot_faststructure_k2,
  width = 12,
  height = 4,
  dpi = 300
)

message("Population-structure plots were written to:")
message(output_dir)
