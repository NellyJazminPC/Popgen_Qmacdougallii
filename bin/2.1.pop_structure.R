# Population genomics of Quercus macdougallii
## PCA, DAPC, and MNS analysis for population structure
## Visualization of fastStructure and ADMIXTURE results
## Nelly J. Pacheco Cruz
## August 2024

# Load libraries
library(vcfR)       # For handling VCF files
library(dartR)      # For working with genomic data in genlight format
library(ggplot2)    # For high-quality graphics
library(adegenet)   # For DAPC analysis and other genetic functions
library(parallel)   # For parallel processing (used in xvalDapc)

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

#################
# PCA #
#################
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
# With the find.clusters function of the adegenet package
# These functions implement the clustering procedure used in Discriminant Analysis of Principal Components (DAPC, Jombart et al. 2010).
# This procedure consists in running successive K-means with an increasing number of clusters (k), after transforming data using a principal component analysis (PCA). For each model, a statistical measure of goodness of fit (by default, BIC; Bayesian Information Criterion) is computed, which allows to choose the optimal k.

grp <- find.clusters(qmacd_genind, max.n.clust=8)

# To obtain a graph with the proportions belonging to each of the two possible groups
table.value(table(pop(qmacd_genind), grp$grp), col.lab=paste("inf", 1:9),
            row.lab=paste("ori", 1:9))

# On an exploratory basis, we use dapc from the adegenet package
dapc1 <- dapc(qmacd_genind, grp$grp)

# Plot
scatter(dapc1)
summary(dapc1)

# Graph of group assignment per individual
assignplot(dapc1, subset=1:79)
dapc1$assign

# dapc1 plot
scatter(dapc1, col = cols, cex = 4, cell=0, cstar = 1, mstree = TRUE, 
        lwd = 2, lty = 2, legend = T, clabel = T, 
        posi.leg = "topleft", txt.leg=c("1", "2", "3", "4"),
        scree.da = T, posi.da = "topright", posi.pca = "topright", 
        scree.pca = T, pch=20, solid = 0.4)

# Using other parameters with dapc:
pnw.dapc <- dapc(qmacd_genlight, var.contrib = TRUE, scale = FALSE, n.pca = NULL, n.da = nPop(qmacd_genlight)-1)

# mstree = Minimal Spanning Tree
# It can be added with True or False
scatter(pnw.dapc, cell = 0, pch = 18:25, 
        cstar = 0, mstree = TRUE, lwd = 2, legend = T, 
        lty = 2)

# Now we reset the plotting parameters to default
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, las = 0)

### xvalDapc is Cross-validation for Discriminant Analysis of Principal Components
set.seed(999)
pramx <- xvalDapc(tab(qmacd_genclone, NA.method = "mean"), pop(qmacd_genclone))

# Then adjust the test range of number of PC, in this case from 5 to 25
# n.pca is the number of different number of PCA axes to be retained for the cross-validation
set.seed(999)
system.time(pramx <- xvalDapc(tab(qmacd_genclone, NA.method = "mean"), 
                              pop(qmacd_genclone), n.pca = 5:25, n.rep = 30, 
                              parallel = "multicore", ncpus = 6))

names(pramx) # The first element are all the samples
pramx[2:6] # To detect the best number of PC retained based on the cross-validation error

# Plot the results
# clabel= site names
scatter(pramx$DAPC, cex = 2, col = cols, cell=1.3, cstar = 0, legend = F, mstree = TRUE, lwd = 2, lty = 2,
        clabel = T, posi.leg = "topleft", scree.pca = F, scree.da = F,
        posi.pca = "topright", posi.da = "bottomleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1, pch=19)

# Open a TIFF device to save the DAPC plot with high resolution
tiff("../results/dapc_plot.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Remove the first number from the group names in pramx$DAPC$grp
pramx$DAPC$grp <- gsub("^[0-9]+", "", pramx$DAPC$grp)

# Convert pramx$DAPC$grp to a factor
pramx$DAPC$grp <- as.factor(pramx$DAPC$grp)

# Plot the DAPC results with small scree plots inside the main plot
scatter(pramx$DAPC, cex = 2, col = cols, cell = 1.3, cstar = 0, legend = F, mstree = TRUE, 
        lwd = 2, lty = 2, clabel = T, posi.leg = "topright", scree.pca = T, scree.da = T,
        posi.pca = "topright", posi.da = "topleft", cleg = 0.2, xax = 1, yax = 2, inset.solid = 1, pch = 19,
        ratio.pca = 0.2, ratio.da = 0.2)  # Adjust the size of scree plots

# Close the TIFF device
dev.off()




# Analyze variable contributions
# Extract the contribution of each SNP to the discriminant functions
var_contrib <- loadingplot(pramx$DAPC$var.contr, axis = 1, lab.jitter = 1, main = "Contributions of SNPs to DA1")

# Save the variable contributions to a CSV file
write.csv(var_contrib, "../results/dapc_variable_contributions.csv", row.names = FALSE)

# Analyze membership probabilities
# Extract the membership probabilities for each individual
membership_probs <- pramx$DAPC$posterior
colnames(membership_probs) <- paste("Cluster", 1:ncol(membership_probs), sep = "_")

# Add individual names and population information
membership_probs <- data.frame(
  Individual = rownames(membership_probs),
  Population = pop(qmacd_genind),
  membership_probs
)

# Save the membership probabilities to a CSV file
write.csv(membership_probs, "../results/dapc_membership_probabilities.csv", row.names = FALSE)

# Plot membership probabilities
# Create a barplot of membership probabilities


# Open a TIFF device to save the membership probabilities plot
tiff("../results/dapc_membership_probabilities_plot.tiff", width = 12, height = 6, units = "in", res = 300, compression = "lzw")

# Adjust plot margins to make space for the rotated labels
par(mar = c(10, 4, 4, 4))  # Bottom, Left, Top, Right margins

# Create the barplot
bp <- barplot(t(as.matrix(membership_probs[, -c(1, 2)])), 
              col = cols, 
              main = "Membership Probabilities of Individuals",
              xlab = "",  # Remove default x-axis label
              ylab = "Membership Probability",
              las = 2,  # Rotate x-axis labels 90 degrees
              cex.names = 0.7)  # Reduce the size of x-axis labels

# Add a manual legend outside the plot area
legend("topright", legend = colnames(membership_probs[, -c(1, 2)]), 
       fill = cols, bty = "n", title = "Clusters", 
       xpd = TRUE, inset = c(-0.2, 0))

# Close the TIFF device
dev.off()


# Summary of DAPC results
# Print a summary of the DAPC analysis
summary(pramx$DAPC)

# Print the number of individuals assigned to each cluster
table(pramx$DAPC$assign)

# Print the discriminant functions and their eigenvalues
print(pramx$DAPC$eig)

# Print the proportion of variance explained by each discriminant function
print(pramx$DAPC$eig / sum(pramx$DAPC$eig) * 100)



#######
# Minimum Spanning Networks #
#######

# Load the igraph package (if not already loaded)
library(igraph)

# Calculate genetic distance
qmacd_dist <- bitwise.dist(qmacd_genclone)

# Generate the Minimum Spanning Network (MSN)
qmacd_msn <- poppr.msn(qmacd_genclone, qmacd_dist, showplot = FALSE, include.ties = TRUE)

# Adjust node sizes using igraph functions
node.size <- rep(2, times = nInd(qmacd_genclone))
names(node.size) <- indNames(qmacd_genclone)
V(qmacd_msn$graph)$size <- node.size  # Use V() from igraph to set vertex attributes

# Plot the MSN
set.seed(12345)
plot_poppr_msn(qmacd_genclone, qmacd_msn, 
               palette = cols,
               gadj = 500)

# Interactive mode (optional)
# imsn()

# Subset the data (if needed)
qmacd_genclone_sub <- popsub(qmacd_genclone, exclude = character(0))

# Handle missing data by imputing with mean
qmacd_genclone_nomiss <- missingno(qmacd_genclone, type = 'mean')

# Calculate Nei's genetic distance
qmacd_genclone_dist <- nei.dist(qmacd_genclone_nomiss, warning = TRUE)

# Generate another MSN with the subsetted data
min_span_net <- poppr.msn(qmacd_genclone_sub, qmacd_genclone_dist, showplot = T, include.ties = TRUE)

# Open a TIFF device to save the MSN plot with high resolution
tiff("../results/msn_plot.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Plot the MSN with Kamada-Kawai layout
set.seed(69)
plot_poppr_msn(qmacd_genclone,
               min_span_net,
               inds = c("CR_01", "CR_02", "IT_01", "IT_02", "IT_03", 
                        "CY_02", "CY_08", "MT_05", "MB_03", "MC_05", 
                        "MT_06", "LS_01", "LS_02", "LS_03", "LS_04"),
               mlg = FALSE,
               gadj = 25,
               nodescale = 51,
               palette = cols,
               cutoff = NULL,  # No aplicar cutoff
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = FALSE,  # Ocultar leyenda de poblaciones
               size.leg = FALSE,  # Ocultar leyenda de sample/node
               scale.leg = TRUE,
               layfun = igraph::layout_with_kk)  # Usar Kamada-Kawai layout

# Close the TIFF device
dev.off()







