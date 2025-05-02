#
#
#





# -------------------------------
# Bayescan; 9 sites
# -------------------------------

##Plots and results visualization
library(boa)
library(coda)
source("software/bayescan_distributed_2.01/R functions/plot_R.r")
# For the nine sites
# Load results
sel <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_9pop/bayescan_qmacd_ref_gen_qrob_.sel",
                  colClasses ="numeric")
bayes <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_9pop/bayescan_qmacd_ref_gen_qrob__fst.txt")
loci <- read.delim("../data/structure_formats/qmacd_ref_gen_rob.bim", header = F)

##Verificar que convergieron las cadenas
chain <- mcmc(sel,thin=10)

par("mar")
par(mar=c(4,4,4,4))

plot(chain)
summary(chain)
autocorr.diag(chain)
autocorr.plot(chain)
effectiveSize(chain)
geweke.diag(chain, frac1=0.1, frac2=0.5)

# Exportar summary
# Cargar la librería
library(writexl)

# Obtener el resumen de la cadena
summary_results <- summary(chain)

# Convertir los resultados a un data frame
summary_df <- as.data.frame(summary_results$statistics)

# Agregar los nombres de las filas como una columna adicional
summary_df <- cbind(Parameter = rownames(summary_df), summary_df)

# Exportar a un archivo Excel
write_xlsx(summary_df, "../results/summary_bayescan_snps_detected.xlsx")

# Mensaje de confirmación
cat("El resumen de Bayescan se ha exportado a '../results/summary_bayescan_snps_detected_9sites.xlsx'\n")

###Identificar outliers
#Unir información de SNP con resultados de Bayescan
bayes <- cbind(bayes, loci[,c(1,2,4)])
colnames(bayes) <- c(colnames(bayes)[1:5], "CHR", "SNP", "Pos")

# Definir el nivel de significancia
alpha <- 0.05

# Generar el gráfico de BayeScan
plot_bayescan_results <- plot_bayescan(
  "../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_9pop/bayescan_qmacd_ref_gen_qrob__fst.txt", 
  FDR = alpha
)

# Guardar el gráfico generado por plot_bayescan
jpeg("../results/plot_bayescan_results.jpg", width = 800, height = 600, res = 300)
plot_bayescan_results
dev.off()

# Identificar outliers según el nivel de significancia
outliers <- bayes[bayes$qval < alpha, ]

# Exportar los outliers a un archivo de texto
write.table(outliers, "../results/bayes_outliers_9sites.txt", row.names = FALSE, quote = FALSE)


# -------------------------------
# BayeScan; 2 sites (Zona norte y zona sur)
# -------------------------------

# Cargar funciones de BayeScan
source("software/bayescan_distributed_2.01/R functions/plot_R.r")

# Cargar resultados de BayeScan para 2 sitios
sel_2pop <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_2pop/bayescan_qmacd_ref_gen_qrob_.sel",
                       colClasses = "numeric")
bayes_2pop <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_2pop/bayescan_qmacd_ref_gen_qrob__fst.txt")
loci_2pop <- read.delim("../data/structure_formats/qmacd_ref_gen_rob.bim", header = FALSE)

# Verificar convergencia de las cadenas
chain_2pop <- mcmc(sel_2pop, thin = 10)

# Graficar y analizar la convergencia
par(mar = c(4, 4, 4, 4))
plot(chain_2pop)
summary(chain_2pop)
autocorr.diag(chain_2pop)
autocorr.plot(chain_2pop)
effectiveSize(chain_2pop)
geweke.diag(chain_2pop, frac1 = 0.1, frac2 = 0.5)

# Exportar el resumen de la cadena a un archivo Excel
summary_results_2pop <- summary(chain_2pop)
summary_df_2pop <- as.data.frame(summary_results_2pop$statistics)
summary_df_2pop <- cbind(Parameter = rownames(summary_df_2pop), summary_df_2pop)
write_xlsx(summary_df_2pop, "../results/summary_bayescan_snps_detected_2pop.xlsx")
cat("El resumen de BayeScan para 2 sitios se ha exportado a '../results/summary_bayescan_snps_detected_2pop.xlsx'\n")

# Unir información de SNP con resultados de BayeScan
bayes_2pop <- cbind(bayes_2pop, loci_2pop[, c(1, 2, 4)])
colnames(bayes_2pop) <- c(colnames(bayes_2pop)[1:5], "CHR", "SNP", "Pos")

# Definir el nivel de significancia
alpha <- 0.05

# Generar el gráfico de BayeScan
plot_bayescan_results_2pop <- plot_bayescan(
  "../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_2pop/bayescan_qmacd_ref_gen_qrob__fst.txt", 
  FDR = alpha
)

# Guardar el gráfico generado por BayeScan
jpeg("../results/plot_bayescan_results_2pop.jpg", width = 800, height = 600, res = 300)
plot_bayescan_results_2pop
dev.off()

# Identificar outliers según el nivel de significancia
outliers_2pop <- bayes_2pop[bayes_2pop$qval < alpha, ]

# Exportar los outliers a un archivo de texto
write.table(outliers_2pop, "../results/bayes_outliers_2pop.txt", row.names = FALSE, quote = FALSE)
cat("Los outliers de BayeScan para 2 sitios se han exportado a '../results/bayes_outliers_2pop.txt'\n")


# -------------------------------
# PCAdapt
# -------------------------------

library(pcadapt)
library(qvalue)
library(ggplot2)

# Cargar los datos genómicos en formato PLINK
qmacd_pcadapt <- read.pcadapt("../data/structure_formats/qmacd_ref_gen_rob.bed", type = "bed")

# Cargar los metadatos de la especie
pop.metadata <- read.csv("../metadata/Qmacdougalli_79ind_.csv")
head(pop.metadata)  # Verificar las primeras filas de los metadatos

# Cargar información de los SNPs
snps <- read.delim("../data/structure_formats/qmacd_ref_gen_rob.bim", sep = "\t", header = FALSE)

# Realizar PCA considerando 20 componentes principales y un filtro de frecuencia alélica mínima (MAF)
pca <- pcadapt(qmacd_pcadapt, K = 20, min.maf = 0.05)

# Visualizar cuántos componentes principales usar (gráfico de screeplot)
plot(pca, option = "screeplot")

# Graficar los scores del PCA, coloreados por el nombre del sitio
plot(pca, option = "scores", pop = pop.metadata$SITE_NAME)

# Repetir el análisis de PCAdapt considerando solo los dos primeros eigenvectores
pcadapt <- pcadapt(qmacd_pcadapt, K = 2, min.maf = 0.05)
str(pcadapt)  # Inspeccionar la estructura del objeto resultante

# Graficar los p-values obtenidos
hist(pcadapt$pvalues, xlab = "p-values", main = NULL, breaks = 100, col = "orange")

# Graficar el Manhattan plot y el QQ plot
plot(pcadapt, option = "manhattan")
plot(pcadapt, option = "qqplot")

# Aplicar corrección de Bonferroni para ajustar los p-values
padj <- p.adjust(pcadapt$pvalues, method = "bonferroni")

# Combinar la información de los SNPs con los p-values y los p-values ajustados
snps <- cbind(snps[, c(1, 2, 4)], pcadapt$pvalues, padj)
colnames(snps) <- c("CHR", "SNP", "Pos", "pval", "padj")

# Identificar los SNPs outliers con un nivel de significancia alfa = 0.05
alfa <- 0.05
outliers <- na.omit(snps[snps$padj < alfa, ])
nrow(outliers)  # Número de SNPs outliers identificados

# Graficar el Manhattan plot con los outliers resaltados en rojo
library(ggplot2)
manhattan_plot <- plot(pcadapt, option = "manhattan") + 
  geom_point(data = outliers, aes(x = as.numeric(rownames(outliers)), y = -log10(pval)), color = "red")
manhattan_plot
# Exportar el gráfico Manhattan con outliers a un archivo PNG
ggsave("../results/pcadapt_manhattan_outliers.png", plot = manhattan_plot, width = 10, height = 6, dpi = 300)

# Exportar los SNPs outliers a un archivo de texto
write.table(outliers, "../results/pcadapt_outliers.txt", row.names = TRUE, quote = FALSE)

cat("El análisis de PCAdapt se completó. Los resultados se han exportado a '../results/'.\n")

# -------------------------------
# SNPs outliers with the highest FST
# -------------------------------

# Instalar y cargar el paquete BiocManager si no está disponible
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Instalar y cargar el paquete snpStats
BiocManager::install("snpStats")
library(snpStats)

# Definir el directorio de trabajo y las rutas de los archivos PLINK
WD <- "/Users/nelly/bioinfo/Popgen_Qmacdougallii/data/structure_formats/"
bed <- paste(WD, "/qmacd_ref_gen_rob.bed", sep = "")
bim <- paste(WD, "/qmacd_ref_gen_rob.bim", sep = "")
fam <- paste(WD, "/qmacd_ref_gen_rob.fam", sep = "")
raw <- paste(WD, "/qmacd_ref_gen_rob.raw", sep = "")
map <- paste(WD, "/qmacd_ref_gen_rob.plk.map", sep = "")

# Cargar los datos de muestra
subject.support <- read.csv(paste("/Users/nelly/bioinfo/Popgen_Qmacdougallii/metadata/Qmacdougalli_79ind_.csv", sep = ""), header = TRUE)

# Integrar datos SNP para el manejo con snpStats
quercus <- read.plink(bed, bim, fam)
snps <- quercus$genotype
snp.support <- quercus$map

# Integrar datos SNP para el manejo con Adegenet
quercusgl <- read.PLINK(file = raw, mapfile = map)

# Eliminar SNPs duplicados
glmatrix <- as.matrix(quercusgl)
glmatrix <- glmatrix[, -grep("HET", colnames(glmatrix))]
quercusgl <- as.genlight(glmatrix)

# Calcular FST para cada SNP entre los nueve sitios
fpopSITE <- Fst(snps, subject.support[, 'SITE_NAME'], pairwise = FALSE)
names(fpopSITE$Fst) <- colnames(snps)
fpopSITE$Fst <- fpopSITE$Fst[!is.na(fpopSITE$Fst)]
cat("Weighted mean FST (9 sites):", weighted.mean(fpopSITE$Fst), "\n")

# Calcular FST para cada SNP entre las zonas (Norte y Sur)
fpopZONE <- Fst(snps, subject.support[, 'ZONE'], pairwise = FALSE)
names(fpopZONE$Fst) <- colnames(snps)
fpopZONE$Fst <- fpopZONE$Fst[!is.na(fpopZONE$Fst)]
cat("Weighted mean FST (zones):", weighted.mean(fpopZONE$Fst), "\n")

# Calcular FST para cada SNP considerando PZ como un clúster aparte
fpop <- Fst(snps, subject.support[, 'POP_ASIG'], pairwise = FALSE)
names(fpop$Fst) <- colnames(snps)
fpop$Fst <- fpop$Fst[!is.na(fpop$Fst)]
cat("Weighted mean FST (PZ cluster):", weighted.mean(fpop$Fst), "\n")

# Función para identificar SNPs con valores FST más altos
fst.alta <- function(snpsfst, percent) {
  alta <- quantile(snpsfst, probs = percent)  # Umbral del percentil deseado
  fst.snp <- data.frame(
    SNPs = names(snpsfst[snpsfst >= 0]),
    Fst = snpsfst[snpsfst >= 0],
    Dif = snpsfst[snpsfst >= 0] > alta
  )
  grafica <- ggplot(fst.snp, aes(x = SNPs, y = Fst, colour = Dif)) +
    geom_point(shape = 19, size = 4, alpha = 0.6) +
    scale_colour_brewer(palette = "Set1") +
    guides(colour = FALSE) +
    theme_bw() +
    ylab(expression(paste("F"[ST], sep = "")))
  snps.alta <- snpsfst[snpsfst > alta]
  return(list(Fst.value = alta, snp.info = snps.alta, fst.snp = fst.snp, grafica = grafica))
}

# Identificar SNPs con FST más altos para PZ como clúster
pop.99 <- fst.alta(fpop$Fst, .99)
num_true <- sum(pop.99$fst.snp$Dif)
cat("Número de SNPs con TRUE en la columna Dif (PZ cluster):", num_true, "\n")

# Exportar los resultados y la gráfica
write.csv(pop.99$fst.snp[pop.99$fst.snp$Dif, ], "../results/fst_outliers_2pop_PZ.csv", row.names = FALSE)
ggsave("../results/fst_outliers_2pop_PZ_plot.png", plot = pop.99$grafica, width = 10, height = 6, dpi = 300)

# Identificar SNPs con FST más altos para las zonas (Norte y Sur)
pop2.99 <- fst.alta(fpopZONE$Fst, .99)
num_true <- sum(pop2.99$fst.snp$Dif)
cat("Número de SNPs con TRUE en la columna Dif (zones):", num_true, "\n")

# Exportar los resultados y la gráfica
write.csv(pop2.99$fst.snp[pop2.99$fst.snp$Dif, ], "../results/fst_outliers_2pop_NS.csv", row.names = FALSE)
ggsave("../results/fst_outliers_2pop_NS_plot.png", plot = pop2.99$grafica, width = 10, height = 6, dpi = 300)

# Identificar SNPs con FST más altos para los nueve sitios
pop3.99 <- fst.alta(fpopSITE$Fst, .99)
num_true <- sum(pop3.99$fst.snp$Dif)
cat("Número de SNPs con TRUE en la columna Dif (9 sites):", num_true, "\n")

# Exportar los resultados y la gráfica
write.csv(pop3.99$fst.snp[pop3.99$fst.snp$Dif, ], "../results/fst_outliers_9sites.csv", row.names = FALSE)
ggsave("../results/fst_outliers_9sites_plot.png", plot = pop3.99$grafica, width = 10, height = 6, dpi = 300)

#-----------------------------
# Frecuencias alélicas de los SNPs outliers - VCF file - genind format
#-----------------------------

# Cargar las librerías
library(adegenet)  # Manejo de objetos genéticos y estadísticas de diversidad
library(vcfR)      # Lectura y conversión de archivos VCF
library(dartR)     # Conversión de genlight a genind y genclone

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

# Convert from genlight to genind and genclone (DartR)
qmacd_genind <- gl2gi(qmacd_genlight, v=1)
# Add population levels using SITE and ZONE metadata
pop(qmacd_genind) <- pop.metadata$SITE_NAME  # Assign SITE as population


#para quitar los que tengan NAs
temp <- seploc(qmacd_genind) 

# Extraer los nombres de los locus de los SNPs outliers
# Leer el archivo y extraer la columna de los locus
outliers <- read.delim(
  "../results/pcadapt_outliers.txt", 
  sep = "",  # Separa por cualquier espacio/tab
  header = FALSE,
  stringsAsFactors = FALSE
)

locus_names <- outliers$V3  # La tercera columna contiene los nombres
# Eliminar el primer elemento "Pos"
locus_names <- locus_names[-1]

# Verificar los nombres actualizados
print(locus_names)

#-----------------------------
# Graficar SNPs outliers de PCAdapt
#-----------------------------

library(adegenet)

# Lista de SNPs con sus valores de pch y p-value adj
snp_info <- list(
  list(name = "loc1548_pos36", pval = "4.72e-9", pch = c("c", "t")), #
  list(name = "loc1479_pos25", pval = "2.503e-5", pch = c("c", "t")),##
  list(name = "loc1479_pos27", pval = "2.674e-5", pch = c("g", "a")),##
  list(name = "loc1237_pos66", pval = "1.569e-6", pch = c("t", "c")),#
  list(name = "loc1302_pos23", pval = "1.881e-8", pch = c("a", "g")),#
  list(name = "loc756_pos43", pval = "1.993e-6", pch = c("g", "a")),#
  list(name = "loc771_pos122", pval = "3.997e-6", pch = c("c", "t")),#
  list(name = "loc891_pos102", pval = "9.2284e-8", pch = c("g", "a")),#
  list(name = "loc1171_pos126", pval = "1.08e-5", pch = c("g", "a")),#
  list(name = "loc656_pos101", pval = "1.865e-5", pch = c("t", "c")),#
  list(name = "loc1192_pos113", pval = "1.442e-5", pch = c("t", "c")),#
  list(name = "loc1192_pos123", pval = "1.442e-5", pch = c("g", "t")),#
  list(name = "loc809_pos115", pval = "2.592e-5", pch = c("a", "t")),#
  list(name = "loc126_pos84", pval = "2.674e-5", pch = c("t", "a")),#
  list(name = "loc127_pos199", pval = "7.292e-13", pch = c("c", "t")),#
  list(name = "loc1380_pos57", pval = "8.989e-7", pch = c("g", "t")),#
  list(name = "loc421_pos96", pval = "2.297e-6", pch = c("t", "a")),#
  list(name = "loc428_pos105", pval = "1.951e-5", pch = c("g", "a"))#
)

# Iterar sobre cada SNP en la lista
for (snp in snp_info) {
  # Extraer el nombre del SNP, el p-value ajustado y los alelos
  snp_name <- snp$name
  snp_pval <- snp$pval
  snp_pch <- snp$pch
  
  # Extraer las frecuencias alélicas del SNP
  snp_data <- tab(temp[[snp_name]])
  freq_snp <- apply(snp_data, 2, function(e) tapply(e, pop(qmacd_genind), mean, na.rm = TRUE))
  
  # Configurar el gráfico
  par(mar = c(5, 5, 4, 2) + 0.1)  # Márgenes ajustados
  matplot(freq_snp, type = "b", pch = snp_pch,
          xlab = "SITE", ylab = "Allele frequency", 
          main = paste("p-value adj", snp_pval, snp_name),
          xaxt = "n", cex = 1.5)
  axis(side = 1, at = 1:9, lab = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT"))
  
  # Exportar el gráfico como PNG
  png_filename <- paste0("../results/allele_frequency_", snp_name, ".png")
  png(png_filename, width = 2000, height = 1200, res = 300)  # Tamaño y resolución ajustados
  par(mar = c(5, 5, 4, 2) + 0.1)  # Márgenes ajustados para el archivo exportado
  matplot(freq_snp, type = "b", pch = snp_pch,
          xlab = "SITE", ylab = "Allele frequency", 
          main = paste("p-value adj", snp_pval, snp_name),
          xaxt = "n", cex = 1.5)
  axis(side = 1, at = 1:9, lab = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT"))
  dev.off()  # Cerrar el dispositivo gráfico
  
  # Mensaje de confirmación
  cat("Gráfico guardado para SNP:", snp_name, "en", png_filename, "\n")
}


#-----------------------------
# Graficar SNPs outliers de Bayescan
#-----------------------------

# Lista de SNPs con sus valores de pch y p-value adj
snp_info <- list(
  list(name = "loc880_pos42", qval = "0.033", pch = c("c", "t")), #
  list(name = "loc379_pos13", qval = "0.046", pch = c("c", "g"))
)

# Iterar sobre cada SNP en la lista
for (snp in snp_info) {
  # Extraer el nombre del SNP, el p-value ajustado y los alelos
  snp_name <- snp$name
  snp_qval <- snp$qval
  snp_pch <- snp$pch
  
  # Extraer las frecuencias alélicas del SNP
  snp_data <- tab(temp[[snp_name]])
  freq_snp <- apply(snp_data, 2, function(e) tapply(e, pop(qmacd_genind), mean, na.rm = TRUE))
  
  # Configurar el gráfico
  par(mar = c(5, 5, 4, 2) + 0.1)  # Márgenes ajustados
  matplot(freq_snp, type = "b", pch = snp_pch,
          xlab = "SITE", ylab = "Allele frequency", 
          main = paste("q-value", snp_qval, snp_name),
          xaxt = "n", cex = 1.5)
  axis(side = 1, at = 1:9, lab = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT"))
  
  # Exportar el gráfico como PNG
  png_filename <- paste0("../results/allele_frequency_bayescan_", snp_name, ".png")
  png(png_filename, width = 2000, height = 1200, res = 300)  # Tamaño y resolución ajustados
  par(mar = c(5, 5, 4, 2) + 0.1)  # Márgenes ajustados para el archivo exportado
  matplot(freq_snp, type = "b", pch = snp_pch,
          xlab = "SITE", ylab = "Allele frequency", 
          main = paste("p-value adj", snp_qval, snp_name),
          xaxt = "n", cex = 1.5)
  axis(side = 1, at = 1:9, lab = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT"))
  dev.off()  # Cerrar el dispositivo gráfico
  
  # Mensaje de confirmación
  cat("Gráfico guardado para SNP:", snp_name, "en", png_filename, "\n")
}

# -------------------------------
# SNPs identificados con PCAdapt y Bayescan
# -------------------------------



# -------------------------------
# SNPs identificados con PCAdapt y Bayescan
# -------------------------------

library(VennDiagram)

#Cargar información de los genes
genes <- read.delim("../data/annotation/Pvulgaris.v2.1.gene.bed", header = F)
genes[,2] <- as.numeric(as.character(genes[,2]))
genes[,3] <- as.numeric(as.character(genes[,3]))
genes[,1] <- as.character(genes[,1])

############################################################################
#Cargar datos
outl.pca <- read.delim("../out/pcadapt_outliers", sep = " ")
outl.pca$Pos <- as.numeric(as.character(outl.pca$Pos))
outl.pca$CHR <- as.character(outl.pca$CHR)

##Identify the genes
genes.pca <- numeric(0)

for (i in 1:nrow(outl.pca)) {
  if (nrow(genes[genes[,1] == outl.pca[i,1] & genes[,2] <= outl.pca[i,3] & genes[,3] >= outl.pca[i,3],])>0) {
    a<- cbind(genes[genes[,1] == outl.pca[i,1] & genes[,2] <= outl.pca[i,3] & genes[,3] >= outl.pca[i,3],],
              outl.pca[i,c(2:5)])
    genes.pca <- rbind(genes.pca, a)
  }
}

#IDs de genes
genes.id.pca <- droplevels(unique(genes.pca[,4]))

#######################################################################################
################ Outliers identified with Bayescan ##################
outl.bay <- read.delim("../out/bayes_outliers", sep = " ")
outl.bay$Pos <- as.numeric(as.character(outl.bay$Pos))
outl.bay$CHR <- as.character(outl.bay$CHR)

##Identificar los genes
genes.bay <- numeric(0)

for (i in 1:nrow(outl.bay)) {
  if (nrow(genes[genes[,1] == outl.bay[i,6] & genes[,2] <= outl.bay[i,8] & genes[,3] >= outl.bay[i,8],])>0) {
    a<- cbind(genes[genes[,1] == outl.bay[i,6] & genes[,2] <= outl.bay[i,8] & genes[,3] >= outl.bay[i,8],],
              outl.bay[i,c(3,4,7,8)])
    genes.bay <- rbind(genes.bay, a)
  }
}

#IDs de genes
genes.id.bay <- droplevels(unique(genes.bay[,4]))

##################################################################
### Outliers identified by the two methods

#NOTA: Si hubieramos encontrado más genes candidatos podríamos buscar cuáles encontraron ambos programas
outliers <- intersect(genes.id.bay, genes.id.pca)

venn.diagram(list(genes.id.bay, genes.id.pca), 
             category.names = c("BayeScan" , "PCAdapt"),
             filename = "../out/venn_genes.png",
             output = TRUE ,
             imagetype="png" ,
             resolution = 600,
             compression = "lzw",
             lwd = 1,
             lty = "blank",
             fill = c("#5BACED", "#F5AC2E"),
             cex = 1.25,
             fontfamily = "times",
             cat.cex = .8,
             cat.default.pos = "outer",
             cat.pos = c(-15, 15),
             #cat.dist = c(0.055, 0.055, 0.055),
             cat.fontfamily = "times"
)


outl.bay.8sites <- read.delim("../results/databases/bayes_outliers_8sites", sep = " ")
outl.bay.8sites <- outl.bay.8sites[,7]

outl.bay.2pop <- read.delim("../results/databases/bayes_outliers_2pop", sep = " ")
outl.bay.2pop <- outl.bay.2pop[,7]


outl.pcadapt <- read.delim("../results/databases/pcadapt_outliers", sep = " ")
outl.pcadapt <- outl.pcadapt[,2]


outl.snpstats.8sites <- read.csv("../results/databases/pop.99.fst.snps.79.2_8sitios.csv")
outl.snpstats.8sites <- outl.snpstats.8sites[,1]

outl.snpstats.2pop <- read.csv("../results/databases/pop2.99.fst.snps.79.2_2poblaciones.csv")
outl.snpstats.2pop <- outl.snpstats.2pop[,1]
class(outl.snpstats.2pop)

outliers <- intersect(outl.bay.2pop, outl.bay.8sites, outl.snpstats.2pop, outl.snpstats.8sites, outl.pcadapt)
class(outliers)

outliers_intersect <- Reduce(intersect, list(outl.bay.2pop, outl.bay.8sites, outl.snpstats.2pop, outl.snpstats.8sites, outl.pcadapt))
outliers_intersect

venn.diagram( x = list(outl.bay.2pop, outl.bay.8sites, outl.snpstats.2pop, outl.snpstats.8sites, outl.pcadapt), 
              category.names = c("Bayescan2pop" ,"Bayescan8sites", "SnpStat2pop", "SnpStat8sites", "PCAdapt"),
              filename = "../results/figures/venn.outl.bay.snpstat.pcadapt.2pop.8sites.png",
              output = TRUE ,
              imagetype="png" ,
              resolution = 600,
              compression = "lzw",
              lwd = 1,
              lty = "blank",
              fill = c("blue", "red", "yellow", "green","pink"),
              cex = 1,
              fontfamily = "times",
              cat.cex = .8,
              cat.default.pos = "outer",
              #cat.pos = c(-15, 15),
              #cat.dist = c(0.055, 0.055, 0.055),
              # cat.fontfamily = "times"
)


#-----------------------------------------
# SNPs outliers FST values
#-----------------------------------------


## Load required libraries
library(ade4)
library (adegenet)
library(poppr)
library(ape)
library(NMF)
library(snpStats)
library(survival)
library(Matrix)
library(ggplot2)
library(gdsfmt)
library(SNPRelate)
library(grid)
library(gridExtra)
library(maptools)
library(ggmap)

##################################################################################################
## Load SNP data
"/media/nel_pc/n311_pc/Campos_Project_003_Analysis/data_metadata_5/CGUM_79/var.79.2.1.sorted.plk.ped"
bed<-paste(WD,"/var.79.2.1.sorted.bed",sep="")
bim<-paste(WD,"/var.79.2.1.sorted.bim",sep="")
fam<-paste(WD,"/var.79.2.1.sorted.fam",sep="")
ped<-paste(WD,"/var.79.2.1.sorted.plk.ped",sep="")
raw<-paste(WD,"/var.79.2.1.sorted.raw",sep="")
map <-paste(WD, "/var.79.2.1.sorted.plk.map",sep="")

## Load sample data
subject.support<-read.csv(paste("/media/nel_pc/n311_pc/Campos_Project_003_Analysis/data_metadata_5/Qmacdougalli_79ind_stacks.csv",sep=""),header=TRUE) 
# Load information about domestication and improvement SNPs
#snpsextract<-read.csv(paste(WD,"/snp.support.csv",sep=""),header=TRUE)

## Integrate SNP data for SNPStats package handling
quercus<-read.plink(bed,bim,fam)


snps<-quercus$genotype
snp.support<-quercus$map
## Integrate SNP data for Adegenet package handling
quercusgl <- read.PLINK(file=raw,mapfile=map)
## Remove duplicated SNPs
glmatrix <- as.matrix(quercusgl)
glmatrix <- glmatrix[,-grep("HET",colnames(glmatrix))]
quercusgl<-as.genlight(glmatrix)


###################################################################

###################################################################

###################################################################
###### Integrate data for SNPRelate package handling
genofile <- snpgdsOpen(paste(WD,"maiz_listo.gds",sep=""))

snpset <- snpgdsLDpruning(genofile,ld.threshold = 1, autosome.only=F) ## Remove monomorphic SNPs
snpset.id <- unlist(snpset) #get all selected snp id

#### Subset SNP matrix keeping only Domestication and improvement SNPs
glmatrixdom <- glmatrix[, pmatch(snpsextract[,2],colnames(glmatrix))]
glmatrixdom <- glmatrixdom[, colSums(is.na(glmatrixdom)) != nrow(glmatrixdom)]

###################################################################

###################################################################

###################################################################


Fst
## Fst measurement for every SNP among POP, keeping SNP name and removing NAs
fpop2 <- Fst (snps,subject.support[,'SITIO'], pairwise = F)
names(fpop2$Fst)<-colnames(snps)
fpop2$Fst<-fpop2$Fst[!is.na(fpop2$Fst)]

weighted.mean(fpop2$Fst)


#Vamos a cambiarlos por sitios
fpop <- Fst (snps,subject.support[,'SITIO'], pairwise = F)
names(fpop$Fst)<-colnames(snps)
fpop$Fst<-fpop$Fst[!is.na(fpop$Fst)]
###
weighted.mean(fpop$Fst)

fmun <- Fst (snps,subject.support[,'MUN'])
names(fmun$Fst)<-colnames(snps)
fmun$Fst<-fmun$Fst[!is.na(fmun$Fst)]


### Function to identify SNPs with an Fst value higher than certain percent treshold. Takes as input a named vector with the Fst value for SNPs and a given treshold. 
fst.alta <- function (snpsfst, percent){ ## snpsfst=resultado de fst por snp para cierto agrupamiento, percent=.99 i.e 99%
  alta <- quantile(snpsfst,probs=percent) ## obtener valor de Fst para el percent deseado
  fst.snp <- data.frame (SNPs=c(1:length(snpsfst[snpsfst>=0])), Fst = snpsfst[snpsfst>=0], Dif = snpsfst[snpsfst>=0]>alta) # armar matriz con snps(Fst>=0), sus valores de Fst y su condicion(>/< percent)
  grafica <- ggplot(fst.snp, aes(x=SNPs, y=Fst, colour=Dif, label=)) + geom_point(shape=19, size= 4, alpha=0.6)+scale_colour_brewer(palette="Set1")+guides(colour=FALSE)+ theme_bw()+ ylab(expression(paste("F"[ST],sep=""))) # scatterplot Fst < percent
  snps.alta <- snpsfst[snpsfst>alta] #lista de snps Fst > treshold
  return(list(Fst.value=alta, snp.info=snps.alta,fst.snp=fst.snp, grafica=grafica))
}
################################################################################################################################
## Identify high Fst SNPs among races
pop.99<-fst.alta(fpop$Fst,.99)
snps_Fstalta<-data.frame(names(pop.99$snp.info))
pop.99$fst.snp
pop.99.fst.snps <- as.data.frame(pop.99$fst.snp)
write.csv(pop.99.fst.snps,'pop.99.fst.snps.79.2_8sitios.csv')

pop2.99<-fst.alta(fpop2$Fst,.99)
snps_2_Fstalta<-data.frame(names(pop2.99$snp.info))
pop2.99$fst.snp
pop2.99.fst.snps <- as.data.frame(pop2.99$fst.snp)
write.csv(pop2.99.fst.snps,'pop2.99.fst.snps.79.2_2poblaciones.csv')
###########################################################################
### Subset SNP matrix keeping only landrace high Fst  SNPs
snpsfst99 <- names(pop.99$snp.info)
glmatrixfst <- glmatrix[, pmatch(snpsfst99, colnames(glmatrix))]
glmatrixfst <- glmatrixfst[, colSums(is.na(glmatrixfst)) != nrow(glmatrixfst)]

snpsfst99_mun <- names(mun.99$snp.info)
glmatrixfst_mun <- glmatrix[, pmatch(snpsfst99_mun, colnames(glmatrix))]
glmatrixfst_mun <- glmatrixfst_mun[, colSums(is.na(glmatrixfst_mun)) != nrow(glmatrixfst_mun)]


## Figur2 2  SNP set

pop.99$grafica <-  pop.99$grafica + 
  geom_point(data=pop.99$fst.snp, aes(x=SNPs, y=Fst),size=2.5,shape=19)

postscript("Figure_2.eps",height=5,width=5)
jpeg("Figure_2.jpg",height=15,width=15, units="cm",res=300, quality=100)
pop.99$grafica  
dev.off()
####################################################################

cgum_79_1_genind <- gl2gi(cgum_79_1_vcf, v=1)
cgum_79_1_genclone <- as.genclone(cgum_79_1_genind)


temp <- seploc(cgum_79_1_genind) # seploc {adegenet} creates a list of individual loci.
#Selección de los snps con mayor fst
library(adegenet)
tab(temp)

########## 8 sitios ###############


snp35_ <- tab(temp[["SSCAFFOLD-11_43273854"]]) #A C
snp17_ <- tab(temp[["SSCAFFOLD-7_32320484"]]) # G A
snp43_ <- tab(temp[["SSCAFFOLD-7_15681531"]]) # C T
snp32_ <- tab(temp[["SSCAFFOLD-15_49540115"]]) # G A
snp22_ <- tab(temp[["SSCAFFOLD-8_7336592"]]) # T C
snp15_ <- tab(temp[["SSCAFFOLD-3_23401487"]]) # G A
snp36_ <- tab(temp[["SSCAFFOLD-14_28801306"]]) # C A
snp2227_ <- tab(temp[["SSCAFFOLD-4_7811362"]]) # G A

(freq35_ <- apply(snp35_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq17_ <- apply(snp17_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq43_ <- apply(snp43_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq32_ <- apply(snp32_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq22_ <- apply(snp22_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq15_ <- apply(snp15_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq36_ <- apply(snp36_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq2227_ <- apply(snp2227_, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))


par(mfrow = c(2, 4), mar = c(5, 4, 4, 0) + 0.1, las = 3)
matplot(freq35_, type = "b", pch = c("a", "c"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.3257",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq17_, type = "b", pch = c("g", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.3009",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq43_, type = "b", pch = c("c", "t"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2803",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq32_, type = "b", pch = c("g", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2765",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq22_, type = "b", pch = c("t", "c"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2731",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq15_, type = "b", pch = c("g", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2693",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq36_, type = "b", pch = c("c", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2645",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq2227_, type = "b", pch = c("g", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2458",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))

################################
######## 2 pop ########################################################################3


snp35_2 <- tab(temp[["SSCAFFOLD-8_7336592"]]) # t c
snp17_2 <- tab(temp[["SSCAFFOLD-1_3049986"]]) # C T
snp43_2 <- tab(temp[["SSCAFFOLD-14_28801306"]]) # c a
snp32_2 <- tab(temp[["SSCAFFOLD-15_9821005"]]) # a g 
snp22_2 <- tab(temp[["SSCAFFOLD-14_4422352"]]) # c g
snp15_2 <- tab(temp[["SSCAFFOLD-4_25297891"]]) # c t
snp36_2 <- tab(temp[["SSCAFFOLD-5_29286424"]]) # g a
snp2227_2 <- tab(temp[["SSCAFFOLD-1_40304988"]]) # A G

(freq35_2 <- apply(snp35_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq17_2 <- apply(snp17_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq43_2 <- apply(snp43_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq32_2 <- apply(snp32_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq22_2 <- apply(snp22_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq15_2 <- apply(snp15_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq36_2 <- apply(snp36_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))
(freq2227_2 <- apply(snp2227_2, 2, function(e) tapply(e, pop(cgum_79_1_genind), mean, na.rm = TRUE)))


par(mfrow = c(2, 4), mar = c(5, 4, 4, 0) + 0.1, las = 3)
matplot(freq35_2, type = "b", pch = c("t", "c"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2096",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq17_2, type = "b", pch = c("c", "t"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.2059",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq43_2, type = "b", pch = c("c", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.1991",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq32_2, type = "b", pch = c("a", "g"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.1951",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq22_2, type = "b", pch = c("c", "g"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.1853",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq15_2, type = "b", pch = c("c", "t"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.1825",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq36_2, type = "b", pch = c("g", "a"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.1796",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))
matplot(freq2227_2, type = "b", pch = c("a", "g"),
        xlab = "SITIO", ylab = "Frecuencia alélica", main = "Fst 0.1759",
        xaxt = "n", cex = 1.5)
axis(side = 1, at = 1:8, lab = c("PZ","CR","LS","MC","MB","CY","MT","CZ"))




#############
Fuente<- rownames(pop2.99$fst.snp)%in%rownames(pop.99$fst.snp) #Error por que snpsextract fue para los domesticados 

pop2.99$grafica<-  pop2.99$grafica + 
  geom_point(data=pop.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst),size=1,shape=1)+
  geom_point(data=pop2.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst,colour=Dif),size=1,shape=19)

pop2.99$grafica<-  pop.99$grafica + 
  geom_point(data=pop2.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst),size=2,shape=19,colour="black")+
  geom_point(data=pop2.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst,colour=Dif),size=1,shape=19)



pop.99$grafica + 
  geom_point(data=pop2.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst),size=3,shape=19,colour="black") +
  geom_point(data=pop2.99$fst.snp[Fuente,], aes(x=SNPs, y=Fst,colour=Dif),size=2,shape=19)



postscript("Figure_2c.eps",height=5,width=5)
jpeg("Figure_2c.jpg",height=15,width=15, units="cm",res=300, quality=100)
pop2.99$grafica  


pop2.99$fst.snp
pop2.99.fst.snps <- as.data.frame(pop2.99$fst.snp)
write.csv(pop2.99.fst.snps,'pop2-pop1.99.fst.snps.79.2.csv')


dev.off()


############################################################################################################3
## Clustering analysis

glfst <- as.genlight(glmatrixfst)
glfst_mun <- as.genlight(glmatrixfst_mun)


set.seed(12)

fstclusters <- find.clusters.genlight(glfst,n.pca = 80)
4


fstclusters_mun <- find.clusters.genlight(glfst_mun,n.pca = 48)
2

fstclusters <- find.clusters.genlight(glfst,n.pca = 48, n.clust = 5)
tab<-data.frame(Sample=rownames(glmatrix),Cluster=fstclusters$grp, Populations=subject.support[,'POP'],row.names = NULL)
tab<-tab[order(tab[,'Cluster'], tab[,'Populations']),]


fstclusters_mun <- find.clusters.genlight(glfst_mun,n.pca = 48, n.clust = 2)
tab_mun<-data.frame(Sample=rownames(glmatrix),Cluster=fstclusters_mun$grp, Populations=subject.support[,'MUN'],row.names = NULL)
tab_mun<-tab_mun[order(tab_mun[,'Cluster'], tab_mun[,'Populations']),]





### Figure_4

cairo_ps("Figure_4.eps",height=8,width=3)
jpeg("Figure_4.jpg",height=8,width=4,units="in",quality=100,res=300)
par(mfrow=c(3,1))
plot(fstclusters$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main=expression(paste("High F"[ST]," SNPs",sep="")))
plot(fstclusters_mun$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main=expression(paste("High F"[ST]," SNPs",sep="")))
plot()

dev.off()


glfst <- as.genlight(glmatrixfst,centers=3)


glfst_mun <- as.genlight(glmatrixfst_mun,centers=3)

fstclusters <- find.clusters.genlight(glfst,n.pca = 80)
3


fstclusters_mun <- find.clusters.genlight(glfst_mun,n.pca = 48)
4

cairo_ps("Figure_4.eps",height=8,width=3)
jpeg("Figure_4.jpg",height=8,width=4,units="in",quality=100,res=300)
par(mfrow=c(2,1))
plot(fstclusters$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main=expression(paste("High F"[ST]," SNPs",sep="")))
plot(fstclusters_mun$Kstat,xlab="Number of clusters", ylab="BIC",type="b",col="blue",font.main=1,main=expression(paste("High F"[ST],  "MUN", "SNPs",sep="")))

dev.off()


