# Script to analyze SNP outliers using Bayescan and PCAdapt

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

# Definir el nivel de significancia (FDR threshold)
# Nota importante:
# Este "alpha" NO corresponde al parámetro α de BayeScan (efecto del locus),
# sino al umbral de tasa de falsos descubrimientos (FDR) aplicado a los q-values.
# En el manuscrito, esto se reporta como: q-value < 0.05 (FDR = 5%).
alpha <- 0.05 # FDR

# Generar el gráfico de BayeScan
plot_bayescan_results <- plot_bayescan(
  "../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_9pop/bayescan_qmacd_ref_gen_qrob__fst.txt", 
  FDR = alpha
)
# Guardar el gráfico generado por plot_bayescan
jpeg("../results/plot_bayescan_results.jpg", width = 800, height = 600, res = 300)
plot_bayescan_results
dev.off()

#-------
# Gráfico de Bayescan con ggplot2
#-------
library(ggplot2)

# Leer los resultados de BayeScan
bayes <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_9pop/bayescan_qmacd_ref_gen_qrob__fst.txt", header = TRUE)

# Añadir columna para destacar outliers
bayes$Outlier <- bayes$qval < alpha

# Calcular el valor mínimo de log10.PO. para los SNPs outliers
limite_x <- min(bayes$log10.PO.[bayes$qval < alpha], na.rm = TRUE)

#Plot
plot_bayescan_results <- ggplot(bayes, aes(x = log10.PO., y = fst, color = Outlier)) +
  geom_point(size = 5, alpha = 0.7) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_vline(xintercept = limite_x, linetype = "dashed", color = "black", linewidth = 0.3) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),      # Tamaño de los números en los ejes
    axis.title = element_text(size = 20)      # Tamaño de los títulos de los ejes
  ) +
  labs(
    title = "",
    x = expression(log[10]~"(PO)"),
    y = expression(F[ST])
  )
print(plot_bayescan_results)
#Guardar el gráfico generado por ggplot2
ggsave(
  filename = "../results/plot_bayescan_results.png",
  plot = plot_bayescan_results,
  width = 8, height = 6, dpi = 300
)


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
# BayeScan; PZ como un clúster aparte
# -------------------------------
library(coda)
library(writexl)

# Cargar funciones de BayeScan
source("software/bayescan_distributed_2.01/R functions/plot_R.r")

# Cargar resultados de BayeScan para PZ como un clúster aparte
sel_pz <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_2pop_PZ/bayescan_qmacd_ref_gen_qrob_2po.sel",
                     colClasses = "numeric")
bayes_pz <- read.table("../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_2pop_PZ/bayescan_qmacd_ref_gen_qrob_2po_fst.txt")
loci_pz <- read.delim("../data/structure_formats/qmacd_ref_gen_rob.bim", header = FALSE)

# Verificar convergencia de las cadenas
chain_pz <- mcmc(sel_pz, thin = 10)

# Graficar y analizar la convergencia
par(mar = c(4, 4, 4, 4))
plot(chain_pz)
summary(chain_pz)
autocorr.diag(chain_pz)
autocorr.plot(chain_pz)
effectiveSize(chain_pz)
geweke.diag(chain_pz, frac1 = 0.1, frac2 = 0.5)

# Exportar el resumen de la cadena a un archivo Excel
summary_results_pz <- summary(chain_pz)
summary_df_pz <- as.data.frame(summary_results_pz$statistics)
summary_df_pz <- cbind(Parameter = rownames(summary_df_pz), summary_df_pz)
write_xlsx(summary_df_pz, "../results/summary_bayescan_snps_detected_pz.xlsx")
cat("El resumen de BayeScan para PZ como un clúster aparte se ha exportado a '../results/summary_bayescan_snps_detected_pz.xlsx'\n")

# Unir información de SNP con resultados de BayeScan
bayes_pz <- cbind(bayes_pz, loci_pz[, c(1, 2, 4)])
colnames(bayes_pz) <- c(colnames(bayes_pz)[1:5], "CHR", "SNP", "Pos")

# Contar cuántos valores son menores a un umbral (por ejemplo, alpha = 0.05)
alpha <- 0.05
outliers_count <- sum(bayes_pz$qval < alpha, na.rm = TRUE)
cat("Número de valores en bayes_pz$qval menores a", alpha, ":", outliers_count, "\n")

# Definir el nivel de significancia
alpha <- 0.05

# Generar el gráfico de BayeScan
plot_bayescan_results_pz <- plot_bayescan(
  "../data/1.6.snps_outliers/bayescan_output/bayescan_qmacd_ref_gen_qrob_2pop_PZ/bayescan_qmacd_ref_gen_qrob_2po_fst.txt", 
  FDR = alpha
)

# Guardar el gráfico generado por BayeScan
jpeg("../results/plot_bayescan_results_pz.jpg", width = 800, height = 600, res = 300)
plot_bayescan_results_pz
dev.off()

# Identificar outliers según el nivel de significancia
outliers_pz <- bayes_pz[bayes_pz$qval < alpha, ]

# Exportar los outliers a un archivo de texto
write.table(outliers_pz, "../results/bayes_outliers_pz.txt", row.names = FALSE, quote = FALSE)
cat("Los outliers de BayeScan para PZ como un clúster aparte se han exportado a '../results/bayes_outliers_pz.txt'\n")





# -------------------------------
# PCAdapt
# -------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")

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

# Personalizar el gráfico Manhattan
manhattan_plot <- plot(pcadapt, option = "manhattan") +
  geom_point(data = outliers, aes(x = as.numeric(rownames(outliers)), y = -log10(pval)), color = "red", size = 3, alpha = 0.7) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18)
  ) +
  labs(
    title = "",
    x = "SNP (with mAF > 0.05)",
    y = expression(-log[10](p-values))
  )

print(manhattan_plot)

# Exportar el gráfico Manhattan personalizado a un archivo PNG
ggsave("../results/pcadapt_manhattan_outliers.png", 
       plot = manhattan_plot, 
       width = 10, height = 8, dpi = 300)

# Exportar los SNPs outliers a un archivo de texto
write.table(outliers, "../results/pcadapt_outliers.txt", row.names = TRUE, quote = FALSE)

cat("El análisis de PCAdapt se completó. Los resultados se han exportado a '../results/'.\n")

# -------------------------------
# SNPs outliers with the highest FST
# -------------------------------

# Instalar y cargar el paquete BiocManager si no está disponible
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# Instalar y cargar el paquete snpStats
#BiocManager::install("snpStats", force = T)
# Cargar la biblioteca
library(snpStats)
library(adegenet)
library(ggplot2)
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
fst.alta <- function (snpsfst, percent){ ## snpsfst=resultado de fst por snp para cierto agrupamiento, percent=.99 i.e 99%
  alta <- quantile(snpsfst,probs=percent) ## obtener valor de Fst para el percent deseado
  fst.snp <- data.frame (SNPs=c(1:length(snpsfst[snpsfst>=0])), Fst = snpsfst[snpsfst>=0], Dif = snpsfst[snpsfst>=0]>alta) # armar matriz con snps(Fst>=0), sus valores de Fst y su condicion(>/< percent)
  grafica <- ggplot(fst.snp, aes(x=SNPs, y=Fst, colour=Dif, label=)) + geom_point(shape=19, size= 4, alpha=0.6)+scale_colour_brewer(palette="Set1")+guides(colour="none")+ theme_bw()+ ylab(expression(paste("F"[ST],sep=""))) # scatterplot Fst < percent
  snps.alta <- snpsfst[snpsfst>alta] #lista de snps Fst > treshold
  return(list(Fst.value=alta, snp.info=snps.alta,fst.snp=fst.snp, grafica=grafica))
}

# Identificar SNPs con FST más altos para PZ como clúster
pop.99 <- fst.alta(fpop$Fst,.99)
num_true <- sum(pop.99$fst.snp$Dif)
cat("Número de SNPs con TRUE en la columna Dif (PZ cluster):", num_true, "\n")

# Ajustar los nombres de los SNPs en el data frame
pop.99$fst.snp$SNPs <- rownames(pop.99$fst.snp)  # Usar directamente los nombres de las filas como nombres de SNPs

# Verificar si los nombres se asignaron correctamente
head(pop.99$fst.snp$SNPs)

# Exportar los resultados filtrados con nombres de SNPs originales
write.csv(pop.99$fst.snp[pop.99$fst.snp$Dif, ], "../results/fst_outliers_2pop_PZ.csv", row.names = FALSE)

# Exportar la gráfica generada por ggplot2
ggsave("../results/fst_outliers_2pop_PZ_plot.png", plot = pop.99$grafica, width = 10, height = 6, dpi = 300)

# Identificar SNPs con FST más altos para las zonas (Norte y Sur)
pop2.99 <- fst.alta(fpopZONE$Fst, .99)
num_true <- sum(pop2.99$fst.snp$Dif)
cat("Número de SNPs con TRUE en la columna Dif (zones):", num_true, "\n")

# Ajustar los nombres de los SNPs en el data frame
pop2.99$fst.snp$SNPs <- rownames(pop2.99$fst.snp)  # Usar directamente los nombres de las filas como nombres de SNPs

# Verificar si los nombres se asignaron correctamente
head(pop2.99$fst.snp$SNPs)

# Exportar los resultados filtrados con nombres de SNPs originales
write.csv(pop2.99$fst.snp[pop2.99$fst.snp$Dif, ], "../results/fst_outliers_2pop_NS.csv", row.names = FALSE)

# Exportar la gráfica generada por ggplot2
ggsave("../results/fst_outliers_2pop_NS_plot.png", plot = pop2.99$grafica, width = 10, height = 6, dpi = 300)

# Identificar SNPs con FST más altos para los nueve sitios
pop3.99 <- fst.alta(fpopSITE$Fst, .99)
num_true <- sum(pop3.99$fst.snp$Dif)
cat("Número de SNPs con TRUE en la columna Dif (9 sites):", num_true, "\n")

# Ajustar los nombres de los SNPs en el data frame
pop3.99$fst.snp$SNPs <- rownames(pop3.99$fst.snp)  # Usar directamente los nombres de las filas como nombres de SNPs

# Verificar si los nombres se asignaron correctamente
head(pop3.99$fst.snp$SNPs)

# Exportar los resultados filtrados con nombres de SNPs originales
write.csv(pop3.99$fst.snp[pop3.99$fst.snp$Dif, ], "../results/fst_outliers_9sites.csv", row.names = FALSE)

# Exportar la gráfica generada por ggplot2
ggsave("../results/fst_outliers_9sites_plot.png", plot = pop3.99$grafica, width = 10, height = 6, dpi = 300)

#-----------------------------
# FST outliers en un solo grafico
#-----------------------------
library(ggplot2)
# Combina los data frames de los tres análisis, agregando una columna "grupo"
fst_9sites <- pop3.99$fst.snp
fst_9sites$grupo <- "9 sitios"

fst_ns <- pop2.99$fst.snp
fst_ns$grupo <- "Norte-Sur"

fst_pz <- pop.99$fst.snp
fst_pz$grupo <- "PZ"

# Añade un índice de fila para cada data frame antes de combinar
fst_9sites$index <- seq_len(nrow(fst_9sites))
fst_ns$index     <- seq_len(nrow(fst_ns))
fst_pz$index     <- seq_len(nrow(fst_pz))

# Unir todos en un solo data frame
fst_todos <- rbind(fst_9sites, fst_ns, fst_pz)



#-------
# Prueba con colores azul, verde y naranja
#-------
library(ggplot2)

# Clasificar loci según el criterio real:
# TRUE = top 1% del FST empírico en cada escenario
# FALSE = loci no candidatos
fst_todos$color_fst <- ifelse(fst_todos$Dif, fst_todos$grupo, "Non-outliers")


#------------
# Gráfico combinado de las pruebas de FST con paleta de colores de Viridis
#------------

# Instala viridis si no lo tienes
#if (!requireNamespace("viridis", quietly = TRUE)) {
#  install.packages("viridis")
#}

# Cargar bibliotecas
library(viridis)
library(dplyr)

# Renombrar los grupos para la leyenda
fst_todos$color_fst <- as.character(fst_todos$color_fst)

fst_todos$color_fst <- dplyr::recode(
  fst_todos$color_fst,
  "9 sitios" = "9 sites",
  "Norte-Sur" = "North-South",
  "PZ" = "PZ",
  "Non-outliers" = "Non-outliers"
)

# Especificar el orden de la leyenda
fst_todos$color_fst <- factor(
  fst_todos$color_fst,
  levels = c("9 sites", "North-South", "PZ", "Non-outliers")
)

# Definir colores
colores_custom <- c(
  "9 sites" = "#440154",
  "North-South" = "#21908C",
  "PZ" = "#FDE725",
  "Non-outliers" = "grey70"
)

# Graficar
fst_combined_plot <- ggplot(fst_todos, aes(x = index, y = Fst)) +
  geom_point(
    aes(fill = color_fst, alpha = ifelse(color_fst == "Non-outliers", 0.5, 0.7)),
    shape = 21, size = 4, color = "black", stroke = 1
  ) +
  scale_fill_manual(values = colores_custom, name = " ") +
  scale_alpha_identity() +
  theme_minimal() +
  labs(
    title = " ",
    x = "SNP",
    y = expression(F[ST])
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20)
  )

print(fst_combined_plot)

# Exportar la figura
ggsave(
  filename = "../results/fst_outliers_combined_custom.png",
  plot = fst_combined_plot,
  width = 10, height = 8, dpi = 300
)
#-----------------------------
# Correcciones con gráfico de FST - tres escenarios
#------------------------------
# Crear una columna que identifique outliers reales (top 1%)
fst_todos$outlier <- fst_todos$Dif

# Crear una nueva columna de color basada en el criterio real
fst_todos$color_fst <- ifelse(fst_todos$outlier, fst_todos$grupo, "Below threshold")

#corregir colores
colores_custom <- c(
  "9 sites" = "#440154",
  "North-South" = "#21908C",
  "PZ" = "#FDE725",
  "Below threshold" = "#31688E"
)

# factor leyenda
fst_todos$color_fst <- factor(
  fst_todos$color_fst,
  levels = c("9 sites", "North-South", "PZ", "Below threshold")
)

fst_combined_plot <- ggplot(fst_todos, aes(x = index, y = Fst)) +
  geom_point(
    aes(fill = color_fst, alpha = ifelse(color_fst == "Below threshold", 0.5, 0.8)),
    shape = 21, size = 4, color = "black", stroke = 1
  ) +
  scale_fill_manual(values = colores_custom, name = " ") +
  scale_alpha_identity() +
  theme_minimal() +
  labs(
    title = " ",
    x = "SNP",
    y = expression(F[ST])
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  )

print(fst_combined_plot)
#-----------------------------
# Venn diagramas Bayescan + PCAdapt + FST
#------------------------------

# Cargar la librería necesaria
library(VennDiagram)

# Cargar los nombres de los SNPs identificados en cada análisis
# Asegúrate de que los archivos contengan los nombres de los SNPs en una columna específica
bayescan_snps <- read.table("../results/bayes_outliers_9sites.txt", header = TRUE)$SNP
pcadapt_snps <- read.table("../results/pcadapt_outliers.txt", header = TRUE)$SNP
fst_pz_snps <- read.csv("../results/fst_outliers_2pop_PZ.csv")$SNPs
fst_ns_snps <- read.csv("../results/fst_outliers_2pop_NS.csv")$SNPs
fst_9sites_snps <- read.csv("../results/fst_outliers_9sites.csv")$SNPs

# Crear una lista con los conjuntos de SNPs
snps_list <- list(
  BayeScan = bayescan_snps,
  PCAdapt = pcadapt_snps,
  FST_PZ = fst_pz_snps,
  FST_NS = fst_ns_snps,
  FST_9Sites = fst_9sites_snps
)

# Generar el diagrama de Venn
venn_plot <- venn.diagram(
  x = snps_list,
  category.names = c("BayeScan", "PCAdapt", "FST_PZ", "FST_NS", "FST_9Sites"),
  filename = "../results/venn_snps_analysis.png",
  output = TRUE,
  imagetype = "png",
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  lty = "blank",
  fill = c("#5BACED", "#F5AC2E", "#8BC34A", "#FF5722", "#9C27B0"),
  cex = 1.25,
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontfamily = "sans",
  cat.default.pos = "outer"
)

# Identificar los SNPs compartidos entre todos los análisis
shared_snps <- Reduce(intersect, snps_list)

# Exportar los SNPs compartidos a un archivo de texto
write.table(shared_snps, "../results/shared_snps_across_analyses.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Mensaje de confirmación
cat("El diagrama de Venn se ha guardado en '../results/venn_snps_analysis.png'.\n")
cat("Los SNPs compartidos entre todos los análisis se han exportado a '../results/shared_snps_across_analyses.txt'.\n")

### FST outliers
# Crear una lista con los conjuntos de SNPs
fst_snps_list <- list(
  FST_PZ = fst_pz_snps,
  FST_NS = fst_ns_snps,
  FST_9Sites = fst_9sites_snps
)

# Generar el diagrama de Venn
venn_plot_fst <- venn.diagram(
  x = fst_snps_list,
  category.names = c("FST_PZ", "FST_NS", "FST_9Sites"),
  filename = "../results/venn_fst_outliers.png",
  output = TRUE,
  imagetype = "png",
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  lty = "blank",
  fill = c("#8BC34A", "#FF5722", "#9C27B0"),
  cex = 1.25,
  fontfamily = "sans",
  cat.cex = 1,
  cat.fontfamily = "sans",
  cat.default.pos = "outer"
)

# Identificar los SNPs compartidos entre los análisis de FST
shared_fst_snps <- Reduce(intersect, fst_snps_list)

# Crear un data frame con los SNPs compartidos y sus valores de FST
shared_fst_values <- data.frame(
  SNPs = shared_fst_snps,
  FST_PZ = fpop$Fst[shared_fst_snps],
  FST_NS = fpopZONE$Fst[shared_fst_snps],
  FST_9Sites = fpopSITE$Fst[shared_fst_snps]
)

# Exportar los SNPs compartidos y sus valores de FST a un archivo CSV
write.csv(shared_fst_values, "../results/shared_fst_snps_with_values.csv", row.names = FALSE)

# Mensaje de confirmación
cat("Los SNPs compartidos y sus valores de FST se han exportado a '../results/shared_fst_snps_with_values.csv'.\n")#-----------------------------

#--------
# FST SNPs outliers and PCAdapt SNPs
#--------

# Crear una lista con los conjuntos de SNPs
snps_list <- list(
  PCAdapt = pcadapt_snps,
  FST_PZ = fst_pz_snps,
  FST_NS = fst_ns_snps,
  FST_9Sites = fst_9sites_snps
)


# Paleta Okabe-Ito (daltónicos friendly)
venn_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442") # naranja, azul, verde, amarillo

venn_plot <- venn.diagram(
  x = snps_list,
  category.names = c("PCAdapt", "FST_PZ", "FST_NS", "FST_9Sites"),
  filename = "../results/venn_pcadapt_fst_outliers.png",
  output = TRUE,
  imagetype = "png",
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  lty = "blank",
  fill = venn_colors,
  cex = 2.5,                # Tamaño del número dentro de los círculos
  fontfamily = "sans",
  cat.cex = 1.5,            # Tamaño de las etiquetas de los grupos
  cat.fontfamily = "sans",
  cat.default.pos = "outer"
)

# Identificar los SNPs compartidos entre PCAdapt y los análisis de FST
shared_snps <- Reduce(intersect, snps_list)

# Crear un data frame con los SNPs compartidos, sus valores de FST y p-values
shared_snps_values <- data.frame(
  SNPs = shared_snps,
  FST_PZ = fpop$Fst[shared_snps],
  FST_NS = fpopZONE$Fst[shared_snps],
  FST_9Sites = fpopSITE$Fst[shared_snps],
  PCAdapt_pval = pcadapt$pvalues[shared_snps]
)

# Ordenar los SNPs compartidos por el p-value de PCAdapt (de mayor a menor)
shared_snps_values <- shared_snps_values[order(-shared_snps_values$PCAdapt_pval), ]

# Exportar los SNPs compartidos y sus valores a un archivo CSV
write.csv(shared_snps_values, "../results/shared_pcadapt_fst_snps_with_values.csv", row.names = FALSE)

# Exportar los SNPs compartidos a un archivo de texto
write.table(shared_snps, "../results/shared_pcadapt_fst_snps.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Mensaje de confirmación
cat("El diagrama de Venn se ha guardado en '../results/venn_pcadapt_fst_outliers.png'.\n")
cat("Los SNPs compartidos entre PCAdapt y los análisis de FST se han exportado a '../results/shared_pcadapt_fst_snps.txt'.\n")
cat("Los SNPs compartidos con sus valores de FST y p-values se han exportado a '../results/shared_pcadapt_fst_snps_with_values.csv'.\n")

#-----------------------------
# SNPs outliers identificados por BayeScan, PCAdapt y FST. Compartidos y unicos
#-----------------------------
# Crear listas con los datos de cada análisis
# BayeScan (9 sitios)
bayescan_snps_9sites <- read.table("../results/bayes_outliers_9sites.txt", header = TRUE)
bayescan_snps_9sites$metodo_identificacion <- "bayescan_9sites"
bayescan_snps_9sites$parametro <- "qval"
bayescan_snps_9sites <- bayescan_snps_9sites[, c("SNP", "metodo_identificacion", "parametro", "qval")]
colnames(bayescan_snps_9sites) <- c("locus_name", "metodo_identificacion", "parametro", "value")

# PCAdapt
pcadapt_snps <- read.table("../results/pcadapt_outliers.txt", header = TRUE)
pcadapt_snps$metodo_identificacion <- "pcadapt"
pcadapt_snps$parametro <- "pval"
pcadapt_snps <- pcadapt_snps[, c("SNP", "metodo_identificacion", "parametro", "pval")]
colnames(pcadapt_snps) <- c("locus_name", "metodo_identificacion", "parametro", "value")

# FST Outliers (PZ)
fst_pz_snps <- read.csv("../results/fst_outliers_2pop_PZ.csv")
fst_pz_snps$metodo_identificacion <- "fst_pz"
fst_pz_snps$parametro <- "fst"
fst_pz_snps <- fst_pz_snps[, c("SNPs", "metodo_identificacion", "parametro", "Fst")]
colnames(fst_pz_snps) <- c("locus_name", "metodo_identificacion", "parametro", "value")

# FST Outliers (Norte-Sur)
fst_ns_snps <- read.csv("../results/fst_outliers_2pop_NS.csv")
fst_ns_snps$metodo_identificacion <- "fst_ns"
fst_ns_snps$parametro <- "fst"
fst_ns_snps <- fst_ns_snps[, c("SNPs", "metodo_identificacion", "parametro", "Fst")]
colnames(fst_ns_snps) <- c("locus_name", "metodo_identificacion", "parametro", "value")

# FST Outliers (9 sitios)
fst_9sites_snps <- read.csv("../results/fst_outliers_9sites.csv")
fst_9sites_snps$metodo_identificacion <- "fst_9sites"
fst_9sites_snps$parametro <- "fst"
fst_9sites_snps <- fst_9sites_snps[, c("SNPs", "metodo_identificacion", "parametro", "Fst")]
colnames(fst_9sites_snps) <- c("locus_name", "metodo_identificacion", "parametro", "value")

# Combinar todos los resultados en un único data frame
all_snps <- rbind(
  bayescan_snps_9sites,
  pcadapt_snps,
  fst_pz_snps,
  fst_ns_snps,
  fst_9sites_snps
)

# Identificar SNPs compartidos (shared) a partir de los análisis de Venn
shared_snps <- read.table("../results/shared_snps_across_analyses.txt", header = FALSE, col.names = "locus_name")
all_snps$shared <- ifelse(all_snps$locus_name %in% shared_snps$locus_name, "yes", "no")

# EXTRA
# Crear una lista con los conjuntos de SNPs
snps_list <- list(
  bayescan_9sites = bayescan_snps_9sites$locus_name,
  pcadapt = pcadapt_snps$locus_name,
  fst_pz = fst_pz_snps$locus_name,
  fst_ns = fst_ns_snps$locus_name,
  fst_9sites = fst_9sites_snps$locus_name
)

# Crear una matriz binaria que indique la presencia de cada SNP en cada análisis
presence_matrix <- sapply(snps_list, function(set) all_snps$locus_name %in% set)

# Contar en cuántos análisis está presente cada SNP
all_snps$shared_count <- rowSums(presence_matrix)

# Crear una columna que indique en qué análisis está presente cada SNP
all_snps$shared_in <- apply(presence_matrix, 1, function(row) {
  paste(names(snps_list)[row], collapse = ", ")
})

# Eliminar duplicados basados en la columna locus_name
unique_snps <- all_snps[!duplicated(all_snps$locus_name), ]

# Contar cuántos duplicados fueron eliminados
num_duplicates <- nrow(all_snps) - nrow(unique_snps)
cat("Número de duplicados eliminados:", num_duplicates, "\n")

# Exportar el data frame sin duplicados como archivo Excel
write_xlsx(unique_snps, "../results/consolidated_snps_results_with_shared_info_unique.xlsx")

# Mensaje de confirmación
cat("El data frame consolidado sin duplicados se ha exportado a '../results/consolidated_snps_results_with_shared_info_unique.xlsx'.\n")


