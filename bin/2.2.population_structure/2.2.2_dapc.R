# Discriminant analysis of principal components (DAPC) for
# Quercus macdougallii using the nine sampling sites as a priori groups.
#
# Workflow:
#   1. Optional exploratory clustering with find.clusters().
#   2. Broad DAPC cross-validation across 5-50 retained PCs.
#   3. Focused DAPC cross-validation across 5-30 retained PCs.
#   4. Final DAPC using the number of PCs with the lowest RMSE.
#
# The broad run reproduces the exploratory supplementary analysis. The
# focused run uses 100 serial replicates to obtain a stable and reproducible
# PC selection. In the final dataset, this procedure retained 7 PCs.
#
# Input:
#   data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf
#   metadata/Qmacdougalli_79ind_.csv
#
# Outputs:
#   results/population_structure/figure5_dapc.png
#   results/population_structure/figure5_dapc.tiff
#   results/population_structure/supplementary_dapc_xval_broad.png
#   results/population_structure/supplementary_dapc_xval_focused.png
#   results/population_structure/dapc_scores.csv
#   results/population_structure/dapc_membership_probabilities.csv
#   results/population_structure/dapc_assignment_table.csv
#   results/population_structure/dapc_cross_validation_summary.csv
#   results/population_structure/dapc_cross_validation_summary.txt
#   results/population_structure/dapc_xval_broad.rds
#   results/population_structure/dapc_xval_focused.rds
#
# Run from the repository root or any subdirectory:
#
#   Rscript bin/2.2.population_structure/2.2.2_dapc.R

required_packages <- c(
  "vcfR",
  "adegenet",
  "dartR",
  "poppr",
  "ggplot2"
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
  library(adegenet)
  library(dartR)
  library(poppr)
  library(ggplot2)
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

repo_root <- find_repo_root()

vcf_path <- file.path(
  repo_root,
  "data",
  "1.3.assembly_variant_calling",
  "ref_gen_qrob_trim01_1_sorted.vcf"
)

metadata_path <- file.path(
  repo_root,
  "metadata",
  "Qmacdougalli_79ind_.csv"
)

output_dir <- file.path(
  repo_root,
  "results",
  "population_structure"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(vcf_path)) {
  stop("Missing VCF file: ", vcf_path, call. = FALSE)
}

if (!file.exists(metadata_path)) {
  stop("Missing metadata file: ", metadata_path, call. = FALSE)
}

metadata <- read.csv(
  metadata_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_metadata_columns <- c("ID", "SITE_NAME", "ZONE")

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

qmacd_vcf <- read.vcfR(vcf_path, verbose = FALSE)
vcf_ids <- colnames(qmacd_vcf@gt)[-1]

if (anyDuplicated(vcf_ids)) {
  stop("The VCF contains duplicate sample IDs.", call. = FALSE)
}

if (!setequal(vcf_ids, metadata$ID)) {
  only_vcf <- setdiff(vcf_ids, metadata$ID)
  only_metadata <- setdiff(metadata$ID, vcf_ids)

  stop(
    "Sample IDs differ between the VCF and metadata.\n",
    "Only in VCF: ", paste(only_vcf, collapse = ", "), "\n",
    "Only in metadata: ", paste(only_metadata, collapse = ", "),
    call. = FALSE
  )
}

metadata <- metadata[match(vcf_ids, metadata$ID), , drop = FALSE]
metadata$SITE_CODE <- sub("^[0-9]+", "", metadata$SITE_NAME)
metadata$ZONE <- toupper(metadata$ZONE)

if (!all(metadata$ZONE %in% c("NORTH", "SOUTH"))) {
  stop(
    "ZONE must contain only NORTH or SOUTH after conversion to uppercase.",
    call. = FALSE
  )
}

site_order <- unique(metadata$SITE_CODE)

qmacd_genlight <- vcfR2genlight(qmacd_vcf)
ploidy(qmacd_genlight) <- 2
pop(qmacd_genlight) <- factor(
  metadata$SITE_NAME,
  levels = unique(metadata$SITE_NAME)
)

message(
  sprintf(
    "Converted VCF to genlight: %d individuals and %d biallelic SNPs retained.",
    nInd(qmacd_genlight),
    nLoc(qmacd_genlight)
  )
)

qmacd_genind <- gl2gi(qmacd_genlight, v = 0)
qmacd_genclone <- as.genclone(qmacd_genind)

genotype_matrix <- tab(
  qmacd_genclone,
  NA.method = "mean"
)

site_groups <- pop(qmacd_genclone)

if (length(site_groups) != nrow(genotype_matrix)) {
  stop(
    "The number of DAPC groups does not match the genotype matrix.",
    call. = FALSE
  )
}

# -------------------------------------------------------------------------
# Optional exploratory clustering
# -------------------------------------------------------------------------
#
# This step was used to explore possible genetic clusters. It is disabled
# because it is interactive and did not define the final nine-site DAPC.

RUN_FIND_CLUSTERS <- FALSE

if (RUN_FIND_CLUSTERS) {
  set.seed(999)

  exploratory_clusters <- find.clusters(
    qmacd_genind,
    max.n.clust = 8
  )

  saveRDS(
    exploratory_clusters,
    file.path(output_dir, "dapc_find_clusters_exploratory.rds")
  )
}

# Parameters shared by both cross-validation stages.
XVAL_SEED <- 999
TRAINING_SET <- 0.9
N_DA <- nlevels(site_groups) - 1L

# -------------------------------------------------------------------------
# Broad cross-validation: 5-50 PCs
# -------------------------------------------------------------------------

set.seed(XVAL_SEED)

png(
  filename = file.path(
    output_dir,
    "supplementary_dapc_xval_broad.png"
  ),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

dapc_xval_broad <- xvalDapc(
  genotype_matrix,
  site_groups,
  n.pca = seq(5, 50, by = 5),
  n.rep = 30,
  n.da = N_DA,
  training.set = TRAINING_SET,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  parallel = "no"
)

dev.off()

saveRDS(
  dapc_xval_broad,
  file.path(output_dir, "dapc_xval_broad.rds")
)

# -------------------------------------------------------------------------
# Focused cross-validation: 5-30 PCs, 100 serial replicates
# -------------------------------------------------------------------------

set.seed(XVAL_SEED)

message(
  "Running focused DAPC cross-validation with 100 serial replicates."
)

png(
  filename = file.path(
    output_dir,
    "supplementary_dapc_xval_focused.png"
  ),
  width = 8,
  height = 6,
  units = "in",
  res = 300
)

dapc_xval_focused <- xvalDapc(
  genotype_matrix,
  site_groups,
  n.pca = 5:30,
  n.rep = 100,
  n.da = N_DA,
  training.set = TRAINING_SET,
  result = "groupMean",
  center = TRUE,
  scale = FALSE,
  parallel = "no"
)

dev.off()

saveRDS(
  dapc_xval_focused,
  file.path(output_dir, "dapc_xval_focused.rds")
)

selected_n_pca <- as.integer(
  dapc_xval_focused[["Number of PCs Achieving Lowest MSE"]]
)

if (
  length(selected_n_pca) != 1L ||
  is.na(selected_n_pca) ||
  selected_n_pca < 1L
) {
  stop(
    "The number of PCs with the lowest RMSE could not be determined.",
    call. = FALSE
  )
}

focused_rmse <- dapc_xval_focused[[
  "Root Mean Squared Error by Number of PCs of PCA"
]]

focused_success <- dapc_xval_focused[[
  "Mean Successful Assignment by Number of PCs of PCA"
]]

cross_validation_summary <- data.frame(
  n_pca = as.integer(names(focused_rmse)),
  mean_success = as.numeric(focused_success),
  rmse = as.numeric(focused_rmse)
)

write.csv(
  cross_validation_summary,
  file.path(
    output_dir,
    "dapc_cross_validation_summary.csv"
  ),
  row.names = FALSE
)

summary_lines <- c(
  "DAPC cross-validation",
  "=====================",
  "",
  paste("Individuals:", nrow(genotype_matrix)),
  paste("Biallelic SNPs:", nLoc(qmacd_genlight)),
  paste("A priori groups:", nlevels(site_groups)),
  paste("Random seed:", XVAL_SEED),
  paste("Training proportion:", TRAINING_SET),
  paste("Discriminant axes:", N_DA),
  "",
  "Broad exploratory run",
  "---------------------",
  "Candidate PCs: 5, 10, ..., 50",
  "Replicates: 30",
  paste(
    "PCs with lowest RMSE:",
    dapc_xval_broad[["Number of PCs Achieving Lowest MSE"]]
  ),
  "",
  "Focused confirmatory run",
  "------------------------",
  "Candidate PCs: 5-30",
  "Replicates: 100",
  paste("PCs with lowest RMSE:", selected_n_pca),
  paste(
    "PCs with highest mean success:",
    dapc_xval_focused[[
      "Number of PCs Achieving Highest Mean Success"
    ]]
  )
)

writeLines(
  summary_lines,
  file.path(
    output_dir,
    "dapc_cross_validation_summary.txt"
  )
)

# The final DAPC is the model returned by the focused cross-validation.
dapc_result <- dapc_xval_focused$DAPC

if (as.integer(dapc_result$n.pca) != selected_n_pca) {
  stop(
    "The DAPC object does not use the selected number of PCs.",
    call. = FALSE
  )
}

dapc_scores <- as.data.frame(dapc_result$ind.coord)
dapc_scores$ID <- rownames(dapc_result$ind.coord)

if (!setequal(dapc_scores$ID, metadata$ID)) {
  stop(
    "DAPC score identifiers do not match the metadata IDs.",
    call. = FALSE
  )
}

dapc_scores <- merge(
  dapc_scores,
  metadata,
  by = "ID",
  sort = FALSE
)

dapc_scores <- dapc_scores[
  match(rownames(dapc_result$ind.coord), dapc_scores$ID),
  ,
  drop = FALSE
]

dapc_scores$SITE_CODE <- factor(
  dapc_scores$SITE_CODE,
  levels = site_order
)

dapc_scores$ZONE <- factor(
  dapc_scores$ZONE,
  levels = c("NORTH", "SOUTH")
)

write.csv(
  dapc_scores,
  file.path(output_dir, "dapc_scores.csv"),
  row.names = FALSE
)

membership_probabilities <- as.data.frame(
  dapc_result$posterior
)

membership_probabilities$ID <- rownames(
  dapc_result$posterior
)

membership_probabilities <- merge(
  membership_probabilities,
  metadata[, c("ID", "SITE_NAME", "SITE_CODE", "ZONE")],
  by = "ID",
  sort = FALSE
)

membership_probabilities <- membership_probabilities[
  match(
    rownames(dapc_result$posterior),
    membership_probabilities$ID
  ),
  ,
  drop = FALSE
]

write.csv(
  membership_probabilities,
  file.path(
    output_dir,
    "dapc_membership_probabilities.csv"
  ),
  row.names = FALSE
)

assignment_table <- as.data.frame.matrix(
  table(
    observed_site = sub(
      "^[0-9]+",
      "",
      as.character(dapc_result$grp)
    ),
    assigned_site = sub(
      "^[0-9]+",
      "",
      as.character(dapc_result$assign)
    )
  )
)

assignment_table$observed_site <- rownames(
  assignment_table
)

rownames(assignment_table) <- NULL

assignment_table <- assignment_table[
  ,
  c(
    "observed_site",
    setdiff(
      colnames(assignment_table),
      "observed_site"
    )
  ),
  drop = FALSE
]

write.csv(
  assignment_table,
  file.path(output_dir, "dapc_assignment_table.csv"),
  row.names = FALSE
)

zone_colors <- c(
  NORTH = "#E69F00",
  SOUTH = "#0072B2"
)

site_shapes <- c(
  CZ = 21,
  MT = 22,
  MC = 23,
  MB = 24,
  CY = 25,
  LS = 21,
  PZ = 22,
  CR = 23,
  IT = 24
)

missing_site_shapes <- setdiff(
  levels(dapc_scores$SITE_CODE),
  names(site_shapes)
)

if (length(missing_site_shapes) > 0) {
  stop(
    "No DAPC shape was defined for: ",
    paste(missing_site_shapes, collapse = ", "),
    call. = FALSE
  )
}

if (!all(c("LD1", "LD2") %in% colnames(dapc_scores))) {
  stop(
    "The DAPC result does not contain both LD1 and LD2.",
    call. = FALSE
  )
}

ellipse_data <- dapc_scores[, c("LD1", "LD2", "ZONE")]

dapc_plot <- ggplot(
  dapc_scores,
  aes(
    x = LD1,
    y = LD2,
    fill = ZONE,
    shape = SITE_CODE
  )
) +
  stat_ellipse(
    data = ellipse_data,
    aes(
      x = LD1,
      y = LD2,
      group = ZONE,
      colour = ZONE
    ),
    inherit.aes = FALSE,
    type = "norm",
    level = 0.95,
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linewidth = 0.4
  ) +
  geom_point(
    size = 5,
    alpha = 0.65,
    colour = "black",
    stroke = 0.5
  ) +
  scale_fill_manual(
    values = zone_colors,
    labels = c(
      NORTH = "North",
      SOUTH = "South"
    ),
    drop = FALSE
  ) +
  scale_colour_manual(
    values = zone_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = site_shapes,
    drop = FALSE
  ) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        colour = "black",
        alpha = 1,
        size = 4
      )
    ),
    shape = guide_legend(
      order = 2,
      override.aes = list(
        fill = "white",
        colour = "black",
        alpha = 1,
        size = 4
      )
    )
  ) +
  labs(
    x = "LD1",
    y = "LD2",
    fill = "Geographic zone",
    shape = "Sampling site"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )

ggsave(
  filename = file.path(output_dir, "figure5_dapc.png"),
  plot = dapc_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "figure5_dapc.tiff"),
  plot = dapc_plot,
  width = 10,
  height = 8,
  dpi = 300,
  compression = "lzw"
)

message(
  sprintf(
    "DAPC completed using nine sampling sites and %d retained PCs.",
    selected_n_pca
  )
)

message("DAPC outputs were written to:")
message(output_dir)
