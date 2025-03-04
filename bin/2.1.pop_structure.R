# Population genomics of Quercus macdougallii
## PCA, DAPC, and MNS analysis for population structure
## Visualization of fastStructure and ADMIXTURE results
## Nelly J. Pacheco Cruz
## August 2024

# Load libraries
library(vcfR)
library(dartR)
library(ggplot2)

# Load the VCF file into R
qmacd_vcf <- read.vcfR("../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf")

# Convert the VCF to a genlight object
qmacd_genlight <- vcfR2genlight(qmacd_vcf)

# Load the metadata file
pop.metadata <- read.csv("../metadata/Qmacdougalli_79ind_.csv")
head(pop.metadata)

# Verify that the individual names in the VCF match the metadata
all(colnames(qmacd_vcf@gt)[-1] == pop.metadata$ID)

# Set ploidy
ploidy(qmacd_genlight) <- 2
qmacd_genlight@ploidy

# Add population levels using SITE and ZONE metadata
pop(qmacd_genlight) <- pop.metadata$SITE  # Assign SITE as population
qmacd_genlight_pop <- qmacd_genlight  # Create a copy to use POP
pop(qmacd_genlight_pop) <- pop.metadata$POP  # Assign POP as population

# Convert from genlight to genind and genclone (DartR)
qmacd_genind <- gl2gi(qmacd_genlight, v=1)
qmacd_genclone <- as.genclone(qmacd_genind)

# For POP
qmacd_genind_pop <- gl2gi(qmacd_genlight_pop, v=1)
qmacd_genclone_pop <- as.genclone(qmacd_genind_pop)

# Color assignment
cols <- c("#7570B3", "#075277", "#00B1E8", "#1FC944",
          "#E6AB02", "#E7298A", "#E07E34", "#F15858", "red")

# PCA 
qmacd_pca <- glPca(qmacd_genlight, nf=80)

# Summary of eigenvalues
sum(100 * qmacd_pca$eig / sum(qmacd_pca$eig))
(qmacd_pca$eig[1] / sum(qmacd_pca$eig)) * 100
(qmacd_pca$eig[2] / sum(qmacd_pca$eig)) * 100
(qmacd_pca$eig[3] / sum(qmacd_pca$eig)) * 100
(qmacd_pca$eig[4] / sum(qmacd_pca$eig)) * 100

# Barplot of eigenvalues
barplot(100 * qmacd_pca$eig / sum(qmacd_pca$eig), main="PCA Eigenvalues")
title(ylab="Percentage of explained variance", line=2)
title(xlab="Eigenvalues", line=1)

# Adjust margins
par(mar=c(4, 4, 4, 4))

# PCA scores
qmacd_pca_scores <- as.data.frame(qmacd_pca$scores)
qmacd_pca_scores$pop <- pop(qmacd_genlight)  # Use SITE as population
qmacd_pca_scores$zone <- pop(qmacd_genlight_pop)  # Use POP as zone

# Replace POP1 and POP2 with North and South in the zone column
qmacd_pca_scores$zone <- ifelse(qmacd_pca_scores$zone == "POP1", "North", 
                                ifelse(qmacd_pca_scores$zone == "POP2", "South", qmacd_pca_scores$zone))

# Crear un vector con los nombres de las poblaciones sin el número inicial
pop_labels <- gsub("^[0-9]+", "", unique(qmacd_pca_scores$pop))

# Asegurarse de que los nombres estén en el orden numérico original
pop_labels <- pop_labels[order(unique(qmacd_pca_scores$pop))]

# PCA PLOT
# Define colors
cols <- c("#00B1E8", "#075277", "#E7298A", "#E07E34", "#F15858", "#E6AB02", "#7570B3", "#1FC944", "red")

# Create the plot
set.seed(12345)
qmacd_PCA_plot <- ggplot(qmacd_pca_scores, aes(x=PC1, y=PC2, colour=as.factor(pop), shape=as.factor(zone), fill=as.factor(pop))) + 
  geom_point(size=4, alpha=0.7) + 
  scale_color_manual(values=cols, name="Populations", labels=pop_labels) +  # Usar etiquetas editadas
  scale_fill_manual(values=cols, name="Populations", labels=pop_labels) +  # Usar etiquetas editadas
  scale_shape_manual(name="Zones", values=c(21, 24, 25), labels=c("North", "South")) +  # Cambiar etiquetas de la leyenda para zone
  geom_hline(yintercept=0) + 
  geom_vline(xintercept=0) + 
  theme_bw() +
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=16),
        axis.title.x=element_text(size=18), 
        axis.text.x=element_text(size=14),
        axis.title.y=element_text(size=18), 
        axis.text.y=element_text(size=16)) + 
  xlab("PC1 %4.62") + 
  ylab("PC2 %2.76")

# Display the plot
print(qmacd_PCA_plot)

# Save the plot in TIFF format (high resolution, widely accepted)
ggsave("../results/qmacd_PCA_plot.tiff", qmacd_PCA_plot, width=10, height=8, dpi=300, compression="lzw")

#################
# DAPC Analysis
################


