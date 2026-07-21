# Principal component analysis of the final Quercus macdougallii SNP dataset.
#
# Input:
#   data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf
#   metadata/Qmacdougalli_79ind_.csv
#
# Outputs:
#   results/population_structure/figure4_pca.png
#   results/population_structure/figure4_pca.tiff
#   results/population_structure/pca_scores.csv
#   results/population_structure/pca_explained_variance.csv
#
# Run from the repository root or any subdirectory:
#
#   Rscript bin/2.2.population_structure/2.2.1_pca.R

required_packages <- c("vcfR", "adegenet", "ggplot2")

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

# Match metadata explicitly to the sample order in the VCF.
metadata <- metadata[match(vcf_ids, metadata$ID), , drop = FALSE]

qmacd_genlight <- vcfR2genlight(qmacd_vcf)
ploidy(qmacd_genlight) <- 2
pop(qmacd_genlight) <- metadata$SITE_NAME

message(
  sprintf(
    "Converted VCF to genlight: %d individuals and %d biallelic SNPs retained.",
    nInd(qmacd_genlight),
    nLoc(qmacd_genlight)
  )
)

max_components <- min(
  nInd(qmacd_genlight) - 1L,
  nLoc(qmacd_genlight)
)

if (max_components < 2L) {
  stop("At least two principal components are required.", call. = FALSE)
}

# Retain all mathematically available axes so percentages are calculated
# relative to the complete eigenvalue spectrum.
qmacd_pca <- glPca(
  qmacd_genlight,
  nf = max_components,
  parallel = FALSE
)

explained_variance <- 100 * qmacd_pca$eig / sum(qmacd_pca$eig)

pca_variance <- data.frame(
  component = seq_along(explained_variance),
  explained_variance_percent = explained_variance,
  cumulative_variance_percent = cumsum(explained_variance)
)

write.csv(
  pca_variance,
  file.path(output_dir, "pca_explained_variance.csv"),
  row.names = FALSE
)

pca_scores <- as.data.frame(qmacd_pca$scores)
pca_scores$ID <- vcf_ids

pca_scores <- merge(
  pca_scores,
  metadata,
  by = "ID",
  sort = FALSE
)

# Restore the original VCF sample order after the explicit metadata join.
pca_scores <- pca_scores[
  match(vcf_ids, pca_scores$ID),
  ,
  drop = FALSE
]

# SITE_NAME includes a numeric prefix (for example, 1CZ and 2MT).
# Create a short site code for plotting while retaining the original field.
pca_scores$SITE_CODE <- sub("^[0-9]+", "", pca_scores$SITE_NAME)

# Preserve the site order already used in the metadata and VCF:
# CZ, MT, MC, MB, CY, LS, PZ, CR, and IT.
site_order <- unique(
  sub("^[0-9]+", "", metadata$SITE_NAME)
)

zone_order <- c("NORTH", "SOUTH")

pca_scores$SITE_CODE <- factor(
  pca_scores$SITE_CODE,
  levels = site_order
)

pca_scores$ZONE <- factor(
  toupper(pca_scores$ZONE),
  levels = zone_order
)

if (anyNA(pca_scores$SITE_CODE)) {
  stop("One or more sampling-site codes could not be assigned.", call. = FALSE)
}

if (anyNA(pca_scores$ZONE)) {
  stop(
    "ZONE must contain only NORTH or SOUTH after conversion to uppercase.",
    call. = FALSE
  )
}

write.csv(
  pca_scores,
  file.path(output_dir, "pca_scores.csv"),
  row.names = FALSE
)

site_colors <- c(
  CZ = "#D55E00",
  MT = "#000000",
  MC = "#56B4E9",
  MB = "#009E73",
  CY = "#F0E442",
  LS = "#999999",
  PZ = "#0072B2",
  CR = "#CC79A7",
  IT = "#E69F00"
)

missing_site_colors <- setdiff(
  levels(pca_scores$SITE_CODE),
  names(site_colors)
)

if (length(missing_site_colors) > 0) {
  stop(
    "No PCA color was defined for: ",
    paste(missing_site_colors, collapse = ", "),
    call. = FALSE
  )
}

zone_shapes <- c(
  NORTH = 21,
  SOUTH = 24
)

pca_plot <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    fill = SITE_CODE,
    shape = ZONE
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(
    size = 5,
    alpha = 0.65,
    colour = "black",
    stroke = 0.5
  ) +
  scale_fill_manual(
    values = site_colors,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = zone_shapes,
    labels = c(
      NORTH = "North",
      SOUTH = "South"
    ),
    drop = FALSE
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      override.aes = list(
        fill = "white",
        colour = "black",
        alpha = 1,
        size = 4
      )
    ),
    fill = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,
        colour = "black",
        alpha = 1,
        size = 4
      )
    )
  ) +
  labs(
    x = sprintf("PC1 (%.2f%%)", explained_variance[1]),
    y = sprintf("PC2 (%.2f%%)", explained_variance[2]),
    fill = "Sampling site",
    shape = "Geographic zone"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 11)
  )

ggsave(
  filename = file.path(output_dir, "figure4_pca.png"),
  plot = pca_plot,
  width = 10,
  height = 8,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "figure4_pca.tiff"),
  plot = pca_plot,
  width = 10,
  height = 8,
  dpi = 300,
  compression = "lzw"
)

message(
  sprintf(
    "PCA completed: PC1 = %.2f%%, PC2 = %.2f%%, cumulative = %.2f%%.",
    explained_variance[1],
    explained_variance[2],
    sum(explained_variance[1:2])
  )
)

message("PCA outputs were written to:")
message(output_dir)
