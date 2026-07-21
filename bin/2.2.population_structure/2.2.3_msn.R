# Minimum spanning network (MSN) for Quercus macdougallii
# based on Nei's genetic distances among individuals.
#
# Input:
#   data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf
#   metadata/Qmacdougalli_79ind_.csv
#
# Outputs:
#   results/population_structure/msn_nei_distances.png
#   results/population_structure/msn_nei_distances.tiff
#   results/population_structure/msn_nei_distance_matrix.csv
#
# Run from the repository root or any subdirectory:
#
#   Rscript bin/2.2.population_structure/2.2.3_msn.R

required_packages <- c(
  "vcfR",
  "adegenet",
  "dartR",
  "poppr",
  "igraph"
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
  library(igraph)
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

plot_msn <- function(
  genclone_object,
  msn_object,
  palette,
  seed = 69
) {
  set.seed(seed)

  plot_poppr_msn(
    genclone_object,
    msn_object,
    mlg = FALSE,
    inds = character(0),
    gadj = 25,
    nodescale = 51,
    palette = palette,
    cutoff = NULL,
    quantiles = FALSE,
    beforecut = TRUE,
    pop.leg = TRUE,
    size.leg = FALSE,
    scale.leg = TRUE,
    scale.leg.title = "Nei distances",
    layfun = igraph::layout_with_kk
  )
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
  "ZONE"
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

metadata <- metadata[match(vcf_ids, metadata$ID), , drop = FALSE]
metadata$SITE_CODE <- sub("^[0-9]+", "", metadata$SITE_NAME)
metadata$ZONE <- toupper(metadata$ZONE)

if (!all(metadata$ZONE %in% c("NORTH", "SOUTH"))) {
  stop(
    "ZONE must contain only NORTH or SOUTH after conversion to uppercase.",
    call. = FALSE
  )
}

qmacd_genlight <- vcfR2genlight(qmacd_vcf)
ploidy(qmacd_genlight) <- 2
pop(qmacd_genlight) <- factor(
  metadata$SITE_CODE,
  levels = unique(metadata$SITE_CODE)
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

if (!identical(indNames(qmacd_genclone), vcf_ids)) {
  stop(
    "Individual order changed during conversion to genclone.",
    call. = FALSE
  )
}

# Calculate Nei's genetic distance among individuals.
nei_distance <- nei.dist(
  qmacd_genclone,
  warning = TRUE
)

if (!inherits(nei_distance, "dist")) {
  stop("nei.dist did not return a dist object.", call. = FALSE)
}

if (attr(nei_distance, "Size") != nInd(qmacd_genclone)) {
  stop(
    "The Nei distance matrix does not match the number of individuals.",
    call. = FALSE
  )
}

nei_matrix <- as.matrix(nei_distance)
rownames(nei_matrix) <- vcf_ids
colnames(nei_matrix) <- vcf_ids

write.csv(
  nei_matrix,
  file.path(output_dir, "msn_nei_distance_matrix.csv"),
  row.names = TRUE
)

# Generate the MSN and retain ties between equally short connections.
qmacd_msn <- poppr.msn(
  qmacd_genclone,
  nei_distance,
  showplot = FALSE,
  include.ties = TRUE
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

population_levels <- levels(pop(qmacd_genclone))
population_codes <- sub("^[0-9]+", "", population_levels)

missing_colors <- setdiff(
  population_codes,
  names(site_colors)
)

if (length(missing_colors) > 0) {
  stop(
    "No MSN color was defined for: ",
    paste(missing_colors, collapse = ", "),
    call. = FALSE
  )
}

msn_palette <- adjustcolor(
  site_colors[population_codes],
  alpha.f = 0.5
)

names(msn_palette) <- population_levels

png(
  filename = file.path(
    output_dir,
    "msn_nei_distances.png"
  ),
  width = 10,
  height = 8,
  units = "in",
  res = 300
)

plot_msn(
  qmacd_genclone,
  qmacd_msn,
  palette = msn_palette
)

dev.off()

tiff_arguments <- list(
  filename = file.path(
    output_dir,
    "msn_nei_distances.tiff"
  ),
  width = 10,
  height = 8,
  units = "in",
  res = 300
)

if (capabilities("cairo")) {
  tiff_arguments$type <- "cairo"
  tiff_arguments$compression <- "lzw"
}

do.call(tiff, tiff_arguments)

plot_msn(
  qmacd_genclone,
  qmacd_msn,
  palette = msn_palette
)

dev.off()

message(
  sprintf(
    "MSN completed for %d individuals using Nei's genetic distance.",
    nInd(qmacd_genclone)
  )
)

message("MSN outputs were written to:")
message(output_dir)
