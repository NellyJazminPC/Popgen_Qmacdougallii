# Population genomics of Quercus macdougallii
## PCA, DAPC, and MNS analysis for population structure
## Visualization of fastStructure and ADMIXTURE results
## Nelly J. Pacheco Cruz
## January 2025

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

# Define a consistent color palette
cols <- c("#00B1E8", "#075277", "#E7298A", "#E07E34", "#F15858", "#E6AB02", "#7570B3", "#1FC944", "red")

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

# Create a vector with the population names without the initial number
pop_labels <- gsub("^[0-9]+", "", unique(qmacd_pca_scores$pop))

# Ensure that the names are in the original numerical order
pop_labels <- pop_labels[order(unique(qmacd_pca_scores$pop))]

# PCA PLOT
# Create the plot
set.seed(12345)
qmacd_PCA_plot <- ggplot(qmacd_pca_scores, aes(x=PC1, y=PC2, colour=as.factor(pop), shape=as.factor(zone), fill=as.factor(pop))) + 
  geom_point(size=4, alpha=0.7) + 
  scale_color_manual(values=cols, name="Populations", labels=pop_labels) +  # Use edited labels
  scale_fill_manual(values=cols, name="Populations", labels=pop_labels) +  # Use edited labels
  scale_shape_manual(name="Zones", values=c(21, 24, 25), labels=c("North", "South")) +  # Change legend labels for zone
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

######
# DAPC plot scatter - 9 sites
#####
# clabel= site names
scatter(pramx$DAPC, cex = 2, col = cols, cell=1.3, cstar = 0, legend = F, mstree = TRUE, lwd = 2, lty = 2,
        clabel = T, posi.leg = "topleft", scree.pca = F, scree.da = F,
        posi.pca = "topright", posi.da = "bottomleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1, pch=19)

# Open a TIFF device to save the DAPC plot with high resolution
tiff("../results/dapc_plot_9sites.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Remove the first number from the group names in pramx$DAPC$grp
pramx$DAPC$grp <- gsub("^[0-9]+", "", pramx$DAPC$grp)

# Convert pramx$DAPC$grp to a factor
pramx$DAPC$grp <- factor(pramx$DAPC$grp, levels = unique(pramx$DAPC$grp))

# Plot the DAPC results with small scree plots inside the main plot
scatter(pramx$DAPC, cex = 2, col = cols, cell = 1.3, cstar = 0, legend = F, mstree = TRUE, 
        lwd = 2, lty = 2, clabel = T, posi.leg = "topright", scree.pca = T, scree.da = T,
        posi.pca = "topright", posi.da = "topleft", cleg = 0.2, xax = 1, yax = 2, inset.solid = 1, pch = 19,
        ratio.pca = 0.2, ratio.da = 0.2)  # Adjust the size of scree plots

# Close the TIFF device
dev.off()


#######
# DAPC plot en scatter - 2 pop - North and South
#######
# Define a consistent color palette for the 9 sites
site_cols <- c("CZ" = "#FFA500",  # Orange for CZ
               "MT" = "#FFA500",  # Orange for MT
               "MC" = "#FFA500",  # Orange for MC
               "MB" = "#FFA500",  # Orange for MB
               "CY" = "#FFA500",  # Orange for CY
               "LS" = "#008000",  # Green for LS
               "PZ" = "#008000",  # Green for PZ
               "CR" = "#008000",  # Green for CR
               "IT" = "#008000")  # Green for IT

# Ensure pramx$DAPC$grp is a factor
pramx$DAPC$grp <- factor(pramx$DAPC$grp)

# Map colors to the levels of pramx$DAPC$grp
level_colors <- site_cols[levels(pramx$DAPC$grp)]

# Debugging: Check the mapping of levels to colors
print("Mapping of levels to colors:")
print(levels(pramx$DAPC$grp))  # Verify the levels
print(level_colors)  # Verify the color assignments

# Open a TIFF device to save the DAPC plot with high resolution
tiff("../results/dapc_plot_2pop.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Plot the DAPC results without mstree, centroids, and labels, and add a legend for sites
scatter(pramx$DAPC, 
        cex = 2, 
        col = level_colors,  # Use the manually assigned colors
        cell = 0,            # Remove background shading
        cstar = 0,           # Remove centroids (circles)
        legend = F,       # Add a legend for the sites
        mstree = FALSE,      # Remove the minimum spanning tree
        lwd = 2, 
        lty = 2, 
        clabel = F,      # Remove labels for the points
        posi.leg = "topright",  # Position the legend in the top-right corner
        scree.pca = T, 
        scree.da = T,
        posi.pca = "topright", 
        posi.da = "topleft", 
        cleg = 0.75,         # Adjust the size of the legend
        xax = 1, 
        yax = 2, 
        inset.solid = 1, 
        pch = 19,            # Use solid points
        ratio.pca = 0.2, 
        ratio.da = 0.2)
# Close the TIFF device
dev.off()

#####
# DAPC - circunferencias de distancias máximas - ggplot
######
# Load ggplot2 and ggforce
library(ggplot2)
library(ggforce)  # For geom_circle

# Extract individual coordinates from the DAPC
dapc_coords <- as.data.frame(pramx$DAPC$ind.coord)  # Extract individual coordinates
dapc_coords$group <- ifelse(pramx$DAPC$grp %in% c("CZ", "MT", "MC", "MB", "CY"), "South", "North")  # Assign groups
dapc_coords$site <- pramx$DAPC$grp  # Assign site information (e.g., CZ, MC, MB)

# Calculate centroids for each group
centroids <- aggregate(. ~ group, data = dapc_coords, FUN = mean)

# Calculate the radius for each group as the maximum distance from the centroid
dapc_coords <- merge(dapc_coords, centroids, by = "group", suffixes = c("", "_centroid"))
dapc_coords$distance <- sqrt((dapc_coords$LD1 - dapc_coords$LD1_centroid)^2 + 
                               (dapc_coords$LD2 - dapc_coords$LD2_centroid)^2)
radius <- aggregate(distance ~ group, data = dapc_coords, FUN = max)
centroids <- merge(centroids, radius, by = "group")  # Add radius to centroids

# Define shapes for each site (filled shapes: 21-25)
site_shapes <- c("CZ" = 21, "MT" = 22, "MC" = 23, "MB" = 24, "CY" = 25, 
                 "LS" = 21, "PZ" = 22, "CR" = 23, "IT" = 24)

# Create the DAPC plot with ggplot2
dapc_plot <- ggplot(dapc_coords, aes(x = LD1, y = LD2, color = group, fill = group, shape = site)) +
  # Add circunferences (borders only, no fill)
  geom_circle(data = centroids, aes(x0 = LD1, y0 = LD2, r = distance, color = group), 
              inherit.aes = FALSE, alpha = 1, fill = NA, size = 0.8) +  # Thinner border
  # Plot individual points with shapes based on site
  geom_point(size = 8, alpha = 0.5) +  # Larger and more translucent points
  # Define colors for groups
  scale_color_manual(values = c("South" = "#FFA500", "North" = "#008000"), 
                     labels = c("South", "North")) +  # Update legend labels
  scale_fill_manual(values = c("South" = "#FFA500", "North" = "#008000")) +  # Fill for shapes
  # Define shapes for sites
  scale_shape_manual(values = site_shapes) +
  # Highlight X and Y axes at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) +  # Highlight Y axis
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +  # Highlight X axis
  # Customize the theme
  theme_bw() +  # Apply theme_bw()
  theme(panel.grid = element_blank(),  # Remove internal grid lines
        legend.title = element_blank(), 
        legend.position = "right",  # Move the legend to the right
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        plot.title = element_blank()) +  # Remove the title
  labs(x = "LD1", y = "LD2")  # Add axis labels

# Print the plot
print(dapc_plot)

# Save the DAPC plot in TIFF format (high resolution, widely accepted)
ggsave("../results/qmacd_DAPC_plot_2pop_ggplot.tiff", dapc_plot, width = 10, height = 8, dpi = 300, compression = "lzw")

# Save the DAPC plot in PNG format (high resolution, widely supported)
ggsave("../results/qmacd_DAPC_plot_2pop_ggplot.png", dapc_plot, width = 10, height = 8, dpi = 300)

#####
# DAPC con Elipse de confianza (Confidence Ellipse) en ggplot
#####
# Load ggplot2
library(ggplot2)

# Extract individual coordinates from the DAPC
dapc_coords <- as.data.frame(pramx$DAPC$ind.coord)  # Extract individual coordinates
dapc_coords$group <- ifelse(pramx$DAPC$grp %in% c("CZ", "MT", "MC", "MB", "CY"), "South", "North")  # Assign groups
dapc_coords$site <- pramx$DAPC$grp  # Assign site information (e.g., CZ, MC, MB)

# Define shapes for each site (filled shapes: 21-25)
site_shapes <- c("CZ" = 21, "MT" = 22, "MC" = 23, "MB" = 24, "CY" = 25, 
                 "LS" = 21, "PZ" = 22, "CR" = 23, "IT" = 24)

# Create the DAPC plot with ggplot2
dapc_plot <- ggplot(dapc_coords, aes(x = LD1, y = LD2, color = group, fill = group)) +
  # Add confidence ellipses with no fill and marked borders
  stat_ellipse(aes(group = group), type = "norm", level = 0.95, 
               geom = "path", size = 1, linetype = "solid") +  # No fill, only border
  # Plot individual points with shapes based on site
  geom_point(aes(shape = site), size = 4, alpha = 0.5) +  # Larger and more translucent points
  # Define colors for groups
  scale_color_manual(values = c("South" = "#FFA500", "North" = "#008000"), 
                     labels = c("South", "North")) +  # Update legend labels
  scale_fill_manual(values = c("South" = "#FFA500", "North" = "#008000")) +  # Fill for points
  # Define shapes for sites
  scale_shape_manual(values = site_shapes) +
  # Highlight X and Y axes at 0
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 1) +  # Highlight Y axis
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1) +  # Highlight X axis
  # Customize the theme
  theme_bw() +  # Apply theme_bw()
  theme(panel.grid = element_blank(),  # Remove internal grid lines
        legend.title = element_blank(), 
        legend.position = "right",  # Move the legend to the right
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        plot.title = element_blank()) +  # Remove the title
  labs(x = "LD1", y = "LD2")  # Add axis labels

# Print the plot
print(dapc_plot)

# Save the DAPC plot in TIFF format (high resolution, widely accepted)
ggsave("../results/qmacd_DAPC_plot_with_ellipses_2pop_ggplot.tiff", dapc_plot, width = 10, height = 8, dpi = 300, compression = "lzw")

# Save the DAPC plot in PNG format (high resolution, widely supported)
ggsave("../results/qmacd_DAPC_plot_with_ellipses_2pop_ggplot.png", dapc_plot, width = 10, height = 8, dpi = 300)


#####
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
               cutoff = NULL,  # Do not apply cutoff
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = FALSE,  # Hide population legend
               size.leg = FALSE,  # Hide sample/node legend
               scale.leg = TRUE,
               layfun = igraph::layout_with_kk)  # Use Kamada-Kawai layout

# Close the TIFF device
dev.off()

### TEST different colours

# Define a consistent color palette with transparency
cols_transparent <- adjustcolor(cols, alpha.f = 0.5)  # 50% transparency

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
               palette = cols_transparent,
               gadj = 500,
               title = "Nei distances")  # Change the title to "Nei distances"

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
tiff("../results/msn_plot_nei_distances.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Plot the MSN with Kamada-Kawai layout and custom scale title
#               inds = c("CR_01", "CR_02", "IT_01", 
#"CY_02", "CY_08", "MT_05", "MB_03", "MC_05", 
#"LS_01", "LS_02", "LS_03", "LS_04"),
set.seed(69)
plot_poppr_msn(qmacd_genclone,
               min_span_net,
               mlg = FALSE,
               inds = character(0),
               gadj = 25,
               nodescale = 51,
               palette = cols_transparent,
               cutoff = NULL,  # Do not apply cutoff
               quantiles = FALSE,
               beforecut = TRUE,
               pop.leg = FALSE,  # Hide population legend
               size.leg = FALSE,  # Hide sample/node legend
               scale.leg = TRUE,  # Show the scale legend
               scale.leg.title = "Nei distances",  # Change the scale title
               layfun = igraph::layout_with_kk)  # Use Kamada-Kawai layout

# Close the TIFF device
dev.off()

# Confirmation message
cat("MSN plot with Nei distances saved to ../results/msn_plot_nei_distances.tiff.\n")



####

#####
# MSN pruebas copilot
#####
# Load required libraries
library(igraph)
library(poppr)

# Calculate genetic distance
qmacd_dist <- bitwise.dist(qmacd_genclone)

# Generate the Minimum Spanning Network (MSN)
qmacd_msn <- poppr.msn(qmacd_genclone, qmacd_dist, showplot = FALSE, include.ties = TRUE)

# Adjust node sizes using igraph functions
node.size <- rep(5, times = nInd(qmacd_genclone))  # Increase node size for better visibility
names(node.size) <- indNames(qmacd_genclone)
V(qmacd_msn$graph)$size <- node.size  # Set node sizes

# Adjust node colors to be more transparent
cols_transparent <- adjustcolor(cols, alpha.f = 0.7)  # 70% (opacidad) transparency
V(qmacd_msn$graph)$color <- cols_transparent[pop(qmacd_genclone)]  # Assign colors based on population

# Adjust edge colors and widths
E(qmacd_msn$graph)$color <- "gray70"  # Light gray for edges
E(qmacd_msn$graph)$width <- 1  # Thin edges for better clarity

# Open a TIFF device to save the MSN plot with high resolution
tiff("../results/msn_plot_9sites_custom.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Plot the MSN with customizations
set.seed(12345)
plot(qmacd_msn$graph, 
     layout = layout_with_kk,  # Kamada-Kawai layout for better spacing
     vertex.label = NA,        # Remove node labels for a cleaner plot
     main = "")  # Add a title

# Close the TIFF device
dev.off()

# Add a legend for the sites
legend("topright", 
       legend = unique(pop(qmacd_genclone)),  # Names of the sites
       col = unique(V(qmacd_msn$graph)$color),  # Colors corresponding to the sites
       pch = 19,              # Use solid circles for the legend
       pt.cex = 1.5,          # Size of the points in the legend
       cex = 0.8,             # Size of the text in the legend
       bty = "n",             # Remove the box around the legend
       title = "Sites")       # Title of the legend


# Plot the MSN with individual labels
set.seed(12345)
plot(qmacd_msn$graph, 
     layout = layout_with_kk,  # Kamada-Kawai layout for better spacing
     vertex.label = indNames(qmacd_genclone),  # Add individual labels
     vertex.label.cex = 0.7,   # Adjust the size of the labels
     vertex.label.color = "black",  # Set the color of the labels
     vertex.label.dist = 1,    # Distance of the labels from the nodes
     vertex.size = V(qmacd_msn$graph)$size,  # Use the node sizes already defined
     vertex.color = V(qmacd_msn$graph)$color,  # Use the node colors already defined
     edge.color = E(qmacd_msn$graph)$color,  # Use the edge colors already defined
     edge.width = E(qmacd_msn$graph)$width,  # Use the edge widths already defined
     main = "Minimum Spanning Network (MSN)")  # Add a title


# Export the MSN to Cytoscape-compatible formats

# 1. Export as GraphML
write_graph(qmacd_msn$graph, file = "../results/qmacd_msn.graphml", format = "graphml")
cat("MSN exported to ../results/qmacd_msn.graphml (GraphML format).\n")

# 2. Export as Edge List
edge_list <- as_data_frame(qmacd_msn$graph, what = "edges")
write.csv(edge_list, "../results/qmacd_msn_edgelist.csv", row.names = FALSE)
cat("MSN exported to ../results/qmacd_msn_edgelist.csv (Edge List format).\n")

# 3. Export node attributes (optional)
node_attributes <- as_data_frame(qmacd_msn$graph, what = "vertices")
write.csv(node_attributes, "../results/qmacd_msn_node_attributes.csv", row.names = FALSE)
cat("Node attributes exported to ../results/qmacd_msn_node_attributes.csv.\n")

# Save the customized MSN plot as a TIFF file
tiff("../results/msn_plot_customized.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")
plot(qmacd_msn$graph, 
     layout = layout_with_kk, 
     vertex.label = NA, 
     main = "")
dev.off()
cat("MSN plot saved to ../results/msn_plot_customized.tiff.\n")


# Define the individuals to label
individuals_to_label <- c("CZ_02", "CR_01", "MB_03")

# Create a vector for vertex labels
vertex_labels <- ifelse(indNames(qmacd_genclone) %in% individuals_to_label, 
                        indNames(qmacd_genclone), 
                        NA)  # Only label the specified individuals

# Open a TIFF device to save the MSN plot with high resolution
tiff("../results/msn_plot_selected_labels.tiff", width = 10, height = 8, units = "in", res = 300, compression = "lzw")

# Plot the MSN with customizations
set.seed(12345)
plot(qmacd_msn$graph, 
     layout = layout_with_kk,  # Kamada-Kawai layout for better spacing
     vertex.label = vertex_labels,  # Add labels only for the specified individuals
     vertex.label.cex = 0.8,        # Adjust the size of the labels
     vertex.label.color = "black",  # Set the color of the labels
     vertex.label.dist = 1,         # Distance of the labels from the nodes
     vertex.size = V(qmacd_msn$graph)$size,  # Use the node sizes already defined
     vertex.color = V(qmacd_msn$graph)$color,  # Use the node colors already defined
     edge.color = E(qmacd_msn$graph)$color,  # Use the edge colors already defined
     edge.width = E(qmacd_msn$graph)$width,  # Use the edge widths already defined
     main = "Minimum Spanning Network (MSN)")  # Add a title

# Close the TIFF device
dev.off()

# Confirmation message
cat("MSN plot with selected labels saved to ../results/msn_plot_selected_labels.tiff.\n")





#### FST ####  

# Cargar las librerías necesarias
library(hierfstat)
library(vegan)
library(ape)
library(dartR)

# Asegúrate de que las poblaciones estén asignadas al objeto genlight
pop(qmacd_genlight) <- as.factor(pop.metadata$SITE)  # Usar la columna SITE como población

# Convertir el objeto genlight a genind (necesario para hierfstat)
qmacd_genind <- gl2gi(qmacd_genlight, v = 1)  # Convertir genlight a genind

# -------------------------------
# 1. Pairwise FST
# -------------------------------
# Calcular FST por pares entre las poblaciones
pairwise_fst <- genet.dist(qmacd_genind, method = "WC84")  # Usa el método de Weir & Cockerham (1984)

# Imprimir los resultados de FST por pares
print(pairwise_fst)


### Gráfico
# Convertir la matriz de pairwise FST a un data frame
pairwise_fst_df <- as.data.frame(as.table(as.matrix(pairwise_fst)))

# Renombrar las columnas para claridad
colnames(pairwise_fst_df) <- c("Population1", "Population2", "FST")

# Filtrar para mantener solo la mitad triangular inferior (sin duplicados)
pairwise_fst_df <- pairwise_fst_df[as.numeric(pairwise_fst_df$Population1) > as.numeric(pairwise_fst_df$Population2), ]

# Crear el heatmap triangular con ggplot2
heatmap_plot <- ggplot(pairwise_fst_df, aes(x = Population1, y = Population2, fill = FST)) +
  geom_tile(color = "white") +  # Crear las celdas del heatmap
  geom_text(aes(label = sprintf("%.2f", FST)), size = 4, color = "black") +  # Agregar valores con 2 decimales
  scale_fill_gradient(low = "white", high = "red", name = "FST") +  # Gradiente de colores
  scale_x_discrete(labels = c("MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT")) +  # Etiquetas personalizadas para X
  scale_y_discrete(labels = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR")) +  # Etiquetas personalizadas para Y
  theme_minimal() +  # Tema minimalista
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotar etiquetas del eje X
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  labs(title = "Pairwise FST (9 Sites)")

# Mostrar el heatmap
print(heatmap_plot)

# Guardar el heatmap en un archivo
ggsave("../results/pairwise_fst_heatmap.tiff", heatmap_plot, width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave("../results/pairwise_fst_heatmap.png", heatmap_plot, width = 8, height = 6, dpi = 300)


# -------------------------------
# 2. Global FST - 9 sites
# -------------------------------

# Calcular estadísticas básicas, incluyendo el FST global
fst_results <- basic.stats(qmacd_genind)

# Extraer el FST global
fst_global <- fst_results$overall["Fst"]

# Imprimir el FST global
cat("Global FST:", fst_global, "\n")

# -------------------------------
# 2. Pairwise FST entre la zona norte y la zona sur
# -------------------------------
# Asegúrate de que las zonas estén asignadas correctamente
pop.metadata$POP <- as.factor(pop.metadata$POP)  # Convertir ZONE a factor
pop(qmacd_genind) <- pop.metadata$POP  # Asignar ZONE como población en el objeto genind

# Calcular FST por pares entre las zonas
pairwise_fst_zones <- genet.dist(qmacd_genind, method = "WC84")  # Método de Weir & Cockerham (1984)

# Imprimir los resultados de FST por pares entre las zonas
print(pairwise_fst_zones)


# -------------------------------
# Exportar datos
# -------------------------------

# Preparar los datos para exportar
# 1. Global FST
global_fst_df <- data.frame(Metric = "Global FST", Value = fst_global)

# 2. Pairwise FST entre zonas
pairwise_fst_zones_df <- as.data.frame(as.table(as.matrix(pairwise_fst_zones)))
colnames(pairwise_fst_zones_df) <- c("Zone1", "Zone2", "FST")

# 3. Pairwise FST entre sitios
pairwise_fst_df <- as.data.frame(as.table(as.matrix(pairwise_fst)))
colnames(pairwise_fst_df) <- c("Site1", "Site2", "FST")
# Crear un archivo Excel con múltiples hojas
writexl::write_xlsx(
  list(
    "Global FST" = global_fst_df,
    "Pairwise FST Zones" = pairwise_fst_zones_df,
    "Pairwise FST Sites" = pairwise_fst_df
  ),
  path = "../results/fst_results.xlsx"
)


# -------------------------------
# 3. PERMANOVA
# -------------------------------
# Asegúrate de que las POP_ASIGN estén asignadas correctamente
pop.metadata$POP_ASIG <- as.factor(pop.metadata$POP_ASIG)  # Convertir ZONE a factor
pop(qmacd_genind) <- pop.metadata$POP_ASIG  # Asignar ZONE como población en el objeto genind

# Calcular distancias genéticas entre individuos
genetic_dist <- dist(tab(qmacd_genind, NA.method = "mean"))  # Matriz de distancias genéticas

# Realizar PERMANOVA para evaluar la diferenciación entre poblaciones
permanova_results <- adonis2(genetic_dist ~ pop(qmacd_genind), data = pop.metadata, permutations = 999)

# Imprimir los resultados de PERMANOVA
print(permanova_results)

# Guardar los resultados de PERMANOVA en un archivo de texto
capture.output(permanova_results, file = "../results/permanova_results.txt")
cat("PERMANOVA results saved to ../results/permanova_results.txt\n")

#--------------------------------#
# 
#--------------------------------#
