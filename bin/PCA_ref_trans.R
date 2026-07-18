library(vcfR)
library(dartR)
library(adegenet)
library(tidyverse)
library(janitor)
library(forcats)
library(ggrepel)

okabe_ito <- c("#D55E00","#000000","#56B4E9","#009E73","#F0E442",
               "#999999","#0072B2","#CC79A7","#E69F00")

run_pca_for_vcf <- function(vcf_path, ref_label, out_prefix, meta_path) {
  # --- Lectura VCF -> genlight ---
  genl <- read.vcfR(vcf_path) |> vcfR2genlight()
  
  # --- Metadata ---
  meta <- read.csv(meta_path) |> clean_names()
  stopifnot(all(colnames(genl) == meta$id))
  
  # --- Poblaciones y ploidía ---
  ploidy(genl) <- 2
  pop(genl) <- meta$site_name
  genl_zone <- genl; pop(genl_zone) <- meta$pop
  
  # --- PCA ---
  set.seed(12345)
  pca <- glPca(genl, nf = 80)
  var_exp <- 100 * pca$eig / sum(pca$eig)
  pc1_lab <- sprintf("PC1 (%.2f%%)", var_exp[1])
  pc2_lab <- sprintf("PC2 (%.2f%%)", var_exp[2])
  
  # --- Scores + etiquetas cortas ---
  scores <- as.data.frame(pca$scores) |>
    mutate(
      id       = meta$id,
      site_raw = meta$site_name,
      zone_raw = meta$pop,
      # limpiar el prefijo numérico del sitio:  "1_CZ" -> "CZ"
      site_short = sub("^[0-9]+_?", "", site_raw),
      # si tus IDs ya incluyen número de réplica, úsalo; si no, asigna índice
      indiv_num  = row_number(),
      label_short = paste0(site_short, "_", indiv_num),
      pop_num  = as.numeric(sub("^([0-9]+).*", "\\1", site_raw)),
      pop_lab  = fct_reorder(site_short, pop_num),
      zone     = factor(recode(zone_raw, "POP1"="North", "POP2"="South"),
                        levels = c("North","South"))
    )
  
  # Colores
  pop_levels <- levels(scores$pop_lab)
  col_map <- setNames(rep(okabe_ito, length.out = length(pop_levels)), pop_levels)
  
  # --- Plot con etiquetas cortas ---
  p <- ggplot(scores, aes(PC1, PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(aes(fill = pop_lab, shape = zone), size = 3.8,
               color = "black", alpha = 0.75, stroke = 0.5) +
    geom_text_repel(aes(label = label_short), size = 3, max.overlaps = 40) +
    scale_fill_manual(values = col_map, name = "Populations") +
    scale_shape_manual(values = c("North"=21, "South"=24), name = "Zones") +
    labs(
      x = pc1_lab, y = pc2_lab,
      title = sprintf("PCA of SNPs (reference: %s)", ref_label),
      subtitle = "Individuals labeled with short site codes; fill = population, shape = zone"
    ) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.text  = element_text(size = 12),
      axis.title   = element_text(size = 14),
      axis.text    = element_text(size = 12),
      plot.title   = element_text(face = "bold")
    )
  
  # Guardar
  ggsave(paste0(out_prefix, "_with_short_labels.png"),  p, width = 10, height = 7, dpi = 300)
  ggsave(paste0(out_prefix, "_with_short_labels.tiff"), p, width = 10, height = 7, dpi = 300, compression = "lzw")
  
  return(p)
}

# ---- Ejecutar para ambos VCF ----
p_qlob <- run_pca_for_vcf(
  vcf_path   = "../data/1.3.assembly_variant_calling/ref_trans_qlob_trim01_1_sorted_names.vcf",
  ref_label  = "Q. lobata transcriptome",
  out_prefix = "../results/qmacd_PCA_qlob_daltonic_friendly",
  meta_path  = "../metadata/Qmacdougalli_79ind_.csv"
)
p_qlob

p_qrob <- run_pca_for_vcf(
  vcf_path   = "../data/1.3.assembly_variant_calling/ref_trans_qrob_trim01_1_sorted_names.vcf",
  ref_label  = "Q. robur transcriptome",
  out_prefix = "../results/qmacd_PCA_qrob_daltonic_friendly_names",
  meta_path  = "../metadata/Qmacdougalli_79ind_.csv"
)
p_qrob













# Paquetes básicos
library(vcfR)
library(adegenet)  # glPca
library(dartR)     # vcfR2genlight
library(tidyverse)
library(ggrepel)   # para evitar traslapes en etiquetas (opcional)

# --- Función PCA básico tomando labels del VCF ---
run_basic_pca_from_vcf <- function(vcf_path, out_prefix = NULL, nf = 50, seed = 123) {
  # 1) Leer VCF y extraer nombres de muestra del VCF
  vcf <- read.vcfR(vcf_path)
  sample_ids <- colnames(vcf@gt)[-1]  # nombres tal cual en el VCF
  
  # 2) VCF -> genlight
  genl <- vcfR2genlight(vcf)
  ploidy(genl) <- 2
  
  # 3) Verificación de correspondencia (importante)
  #   Debe ser TRUE para asegurar que los labels van con los puntos correctos.
  message("IDs iguales entre VCF y genlight: ",
          all(indNames(genl) == sample_ids))
  
  # 4) PCA
  set.seed(seed)
  pca <- glPca(genl, nf = nf)
  
  # 5) Porcentaje de varianza
  var_exp <- 100 * pca$eig / sum(pca$eig)
  pc1_lab <- sprintf("PC1 (%.2f%%)", var_exp[1])
  pc2_lab <- sprintf("PC2 (%.2f%%)", var_exp[2])
  
  # 6) Data frame para graficar (usa labels del VCF)
  scr <- as.data.frame(pca$scores) %>%
    mutate(sample_id = sample_ids)
  
  # 7) Gráfico
  gp <- ggplot(scr, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.2, color = "grey70") +
    geom_point(size = 3.2, color = "black") +
    geom_text_repel(aes(label = sample_id), size = 3, max.overlaps = 50) +
    labs(x = pc1_lab, y = pc2_lab,
         title = "PCA of SNPs",
         subtitle = "Labels taken directly from VCF sample names") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold"))
  
  print(gp)
  
  # 8) Guardar (opcional)
  if (!is.null(out_prefix)) {
    ggsave(paste0(out_prefix, ".png"),  gp, width = 10, height = 8, dpi = 300)
    ggsave(paste0(out_prefix, ".tiff"), gp, width = 10, height = 8, dpi = 300, compression = "lzw")
  }
  
  invisible(list(plot = gp, pca = pca, scores = scr))
}

# --- Ejemplos de uso ---
# 1) Con Q. lobata
res_qlob <- run_basic_pca_from_vcf(
  vcf_path   = "../data/1.3.assembly_variant_calling/ref_trans_qlob_trim01_1_sorted.vcf",
  out_prefix = "../results/pca_basic_qlob"
)

# 2) Con Q. robur
res_qrob <- run_basic_pca_from_vcf(
  vcf_path   = "../data/1.3.assembly_variant_calling/ref_trans_qrob_trim01_1_sorted.vcf",
  out_prefix = "../results/pca_basic_qrob"
)
