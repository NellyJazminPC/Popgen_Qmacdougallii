# -------------------------------
# SNP-based nucleotide diversity (pi) from VCF
# -------------------------------

library(vcfR)
library(dplyr)

# 1. Load VCF
qmacd_vcf <- read.vcfR("../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf")

# 2. Load metadata
pop.metadata <- read.csv("../metadata/Qmacdougalli_79ind_.csv", stringsAsFactors = FALSE)

# Check that IDs match the VCF sample order
stopifnot(all(colnames(qmacd_vcf@gt)[-1] == pop.metadata$ID))

# 3. Extract genotype matrix (remove FORMAT column)
gt <- extract.gt(qmacd_vcf, element = "GT")

# 4. Function to convert diploid genotypes to alt allele counts
gt_to_altcount <- function(x) {
  x[x %in% c("./.", ".|.")] <- NA
  out <- rep(NA_real_, length(x))
  
  out[x %in% c("0/0", "0|0")] <- 0
  out[x %in% c("0/1", "1/0", "0|1", "1|0")] <- 1
  out[x %in% c("1/1", "1|1")] <- 2
  
  return(out)
}

# Apply per locus
alt_counts <- apply(gt, 1, gt_to_altcount)  # loci x individuals? depends on apply
alt_counts <- t(alt_counts)                 # rows = loci, cols = individuals

colnames(alt_counts) <- pop.metadata$ID
rownames(alt_counts) <- rownames(gt)

# 5. Function to calculate SNP-based pi
# For biallelic diploid SNPs: pi = 2p(1-p), averaged across polymorphic loci in the VCF
calc_pi <- function(mat) {
  # mat: rows = loci, cols = individuals, values = 0/1/2 alt allele counts
  p_alt <- rowMeans(mat / 2, na.rm = TRUE)
  pi_locus <- 2 * p_alt * (1 - p_alt)
  mean(pi_locus, na.rm = TRUE)
}

# 6. Total pi
pi_total <- calc_pi(alt_counts)

# 7. Pi by site
site_list <- split(pop.metadata$ID, pop.metadata$SITE_NAME)

pi_by_site <- data.frame(
  Site = names(site_list),
  pi_snp = sapply(site_list, function(ids) {
    calc_pi(alt_counts[, ids, drop = FALSE])
  }),
  row.names = NULL
)

# 8. Pi by zone
zone_list <- split(pop.metadata$ID, pop.metadata$ZONE)

pi_by_zone <- data.frame(
  Zone = names(zone_list),
  pi_snp = sapply(zone_list, function(ids) {
    calc_pi(alt_counts[, ids, drop = FALSE])
  }),
  row.names = NULL
)

# 9. Add total
pi_total_df <- data.frame(
  Group = "Total",
  pi_snp = pi_total
)

# 10. Print results
print(pi_by_site)
print(pi_by_zone)
print(pi_total_df)

# 11. Optional: export
write.csv(pi_by_site, "../results/pi_by_site_snp_based.csv", row.names = FALSE)
write.csv(pi_by_zone, "../results/pi_by_zone_snp_based.csv", row.names = FALSE)
write.csv(pi_total_df, "../results/pi_total_snp_based.csv", row.names = FALSE)





# mat = rows are loci, columns are individuals
# values = 0, 1, 2 alt allele counts

calc_snp_diversity <- function(mat, n_mode = c("chromosomes", "taxa")) {
  
  n_mode <- match.arg(n_mode)
  
  # Number of individuals genotyped per locus
  n_ind_locus <- rowSums(!is.na(mat))
  
  # This function assumes no missing data or constant sample size across loci
  if (length(unique(n_ind_locus)) > 1) {
    warning("Sample size varies among loci. Tajima's D formula assumes constant n.")
  }
  
  n_ind <- min(n_ind_locus)
  
  if (n_mode == "chromosomes") {
    n <- 2 * n_ind
  } else {
    n <- n_ind
  }
  
  # Allele frequency of the alternative allele per locus
  p_alt <- rowSums(mat, na.rm = TRUE) / (2 * n_ind_locus)
  
  # Segregating sites within the group
  segregating <- p_alt > 0 & p_alt < 1
  S <- sum(segregating, na.rm = TRUE)
  
  # Number of loci evaluated
  L <- length(p_alt)
  
  # Sample-corrected nucleotide diversity per locus
  pi_locus <- (n / (n - 1)) * 2 * p_alt * (1 - p_alt)
  pi <- mean(pi_locus, na.rm = TRUE)
  
  # Watterson's theta
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  
  thetaW <- (S / a1) / L
  
  # Tajima's D
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
  
  e1 <- c1 / a1
  e2 <- c2 / ((a1^2) + a2)
  
  k <- pi * L
  
  if (S > 1) {
    tajimaD <- (k - (S / a1)) / sqrt((e1 * S) + (e2 * S * (S - 1)))
  } else {
    tajimaD <- NA
  }
  
  data.frame(
    n_individuals = n_ind,
    n_sequences = n,
    segregating_sites = S,
    loci_evaluated = L,
    pi = pi,
    thetaW = thetaW,
    Tajimas_D = tajimaD
  )
}


# By site
site_list <- split(pop.metadata$ID, pop.metadata$SITE_NAME)

div_by_site <- do.call(rbind, lapply(names(site_list), function(site) {
  ids <- site_list[[site]]
  out <- calc_snp_diversity(alt_counts[, ids, drop = FALSE], n_mode = "chromosomes")
  out$Group <- site
  out
}))

# By zone
zone_list <- split(pop.metadata$ID, pop.metadata$ZONE)

div_by_zone <- do.call(rbind, lapply(names(zone_list), function(zone) {
  ids <- zone_list[[zone]]
  out <- calc_snp_diversity(alt_counts[, ids, drop = FALSE], n_mode = "chromosomes")
  out$Group <- zone
  out
}))

# Species-wide
div_total <- calc_snp_diversity(alt_counts, n_mode = "chromosomes")
div_total$Group <- "Total"

# Combine
div_all <- rbind(div_by_site, div_by_zone, div_total)

div_all <- div_all[, c(
  "Group",
  "n_individuals",
  "n_sequences",
  "segregating_sites",
  "loci_evaluated",
  "pi",
  "thetaW",
  "Tajimas_D"
)]

print(div_all)




calc_snp_diversity(alt_counts, n_mode = "taxa")
calc_snp_diversity(alt_counts, n_mode = "chromosomes")






gt <- extract.gt(qmacd_vcf, element = "GT")

gt_simple <- sub(":.*$", "", gt)
gt_simple <- gsub("\\|", "/", gt_simple)

sort(table(as.vector(gt_simple), useNA = "ifany"), decreasing = TRUE)



sum(is.na(gt_simple))
mean(is.na(gt_simple))

row_missing <- rowSums(is.na(gt_simple))
summary(row_missing)

col_missing <- colSums(is.na(gt_simple))
summary(col_missing)


summary(rowSums(!is.na(alt_counts)))
summary(colSums(is.na(alt_counts)))
###############################################
#############################################
# Prueba para SNPs bialélicos


# Genotypes already simplified
gt <- extract.gt(qmacd_vcf, element = "GT")
gt_simple <- sub(":.*$", "", gt)
gt_simple <- gsub("\\|", "/", gt_simple)

# Loci with multiallelic ALT field
multi_alt_fix <- grepl(",", qmacd_vcf@fix[, "ALT"])

# Loci where genotypes contain allele 2
multi_gt <- apply(gt_simple, 1, function(x) any(grepl("2", x)))

# Summary
table(multi_alt_fix, multi_gt)

sum(multi_alt_fix)
sum(multi_gt)
sum(multi_alt_fix | multi_gt)

#### Solo loci bialélicos

biallelic_loci <- !(multi_alt_fix | multi_gt)

gt_bi <- gt_simple[biallelic_loci, , drop = FALSE]

# Check genotype classes after filtering
sort(table(as.vector(gt_bi), useNA = "ifany"), decreasing = TRUE)

#### Conteos alélicos

alt_counts_bi <- matrix(
  NA_real_,
  nrow = nrow(gt_bi),
  ncol = ncol(gt_bi),
  dimnames = dimnames(gt_bi)
)

alt_counts_bi[gt_bi == "0/0"] <- 0
alt_counts_bi[gt_bi %in% c("0/1", "1/0")] <- 1
alt_counts_bi[gt_bi == "1/1"] <- 2

# This should be TRUE
stopifnot(!any(is.na(alt_counts_bi)))

#### Recalcular pi, theta y Tajima

calc_snp_diversity_biallelic <- function(mat, min_ind_for_D = 4) {
  
  n_ind <- ncol(mat)
  n <- 2 * n_ind
  L <- nrow(mat)
  
  p_alt <- rowSums(mat) / (2 * n_ind)
  
  segregating <- p_alt > 0 & p_alt < 1
  S <- sum(segregating)
  
  # Sample-corrected nucleotide diversity
  pi_locus <- (n / (n - 1)) * 2 * p_alt * (1 - p_alt)
  pi <- mean(pi_locus)
  
  a1 <- sum(1 / (1:(n - 1)))
  a2 <- sum(1 / ((1:(n - 1))^2))
  
  thetaW <- (S / a1) / L
  
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  
  c1 <- b1 - (1 / a1)
  c2 <- b2 - ((n + 2) / (a1 * n)) + (a2 / (a1^2))
  
  e1 <- c1 / a1
  e2 <- c2 / ((a1^2) + a2)
  
  k <- pi * L
  
  if (S > 1 && n_ind >= min_ind_for_D) {
    tajimaD <- (k - (S / a1)) / sqrt((e1 * S) + (e2 * S * (S - 1)))
  } else {
    tajimaD <- NA
  }
  
  data.frame(
    n_individuals = n_ind,
    n_sequences = n,
    loci_evaluated = L,
    segregating_sites = S,
    pi = pi,
    thetaW = thetaW,
    Tajimas_D = tajimaD
  )
}

#### Por sitio y zona

# By site
site_list <- split(pop.metadata$ID, pop.metadata$SITE_NAME)

div_by_site <- do.call(rbind, lapply(names(site_list), function(site) {
  ids <- site_list[[site]]
  out <- calc_snp_diversity_biallelic(alt_counts_bi[, ids, drop = FALSE])
  out$Group <- site
  out
}))

# By zone
zone_list <- split(pop.metadata$ID, pop.metadata$ZONE)

div_by_zone <- do.call(rbind, lapply(names(zone_list), function(zone) {
  ids <- zone_list[[zone]]
  out <- calc_snp_diversity_biallelic(alt_counts_bi[, ids, drop = FALSE])
  out$Group <- zone
  out
}))

# Species-wide
div_total <- calc_snp_diversity_biallelic(alt_counts_bi)
div_total$Group <- "Total"

# Combine
div_all <- rbind(div_by_site, div_by_zone, div_total)

div_all <- div_all[, c(
  "Group",
  "n_individuals",
  "n_sequences",
  "loci_evaluated",
  "segregating_sites",
  "pi",
  "thetaW",
  "Tajimas_D"
)]

print(div_all)
