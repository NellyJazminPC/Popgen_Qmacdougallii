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
pop(qmacd_genlight) <- pop.metadata$SITE_NAME  # Assign SITE as population
qmacd_genlight_pop <- qmacd_genlight  # Create a copy to use POP
pop(qmacd_genlight_pop) <- pop.metadata$ZONE  # Assign POP as population

# Convert from genlight to genind and genclone (DartR)
qmacd_genind <- gl2gi(qmacd_genlight, v=1)
qmacd_genclone <- as.genclone(qmacd_genind)

# For POP
qmacd_genind_pop <- gl2gi(qmacd_genlight_pop, v=1)
qmacd_genclone_pop <- as.genclone(qmacd_genind_pop)


# Calcular estadísticas de diversidad genética

# ajustar

snps_divers <- summary(qmacd_genind)
names(snps_divers)

plot(snps_divers$n.by.pop, snps_divers$pop.n.all)

barplot(snps_divers$loc.n.all)

barplot(snps_divers$Hexp-snps_divers$Hobs, main="Heterozygosity: expected-observed",
        ylab="Hexp - Hobs")
barplot(snps_divers$pop.eff, main="Sample sizes per population",
        ylab="Number of genotypes",las=3)

bartlett.test(list(snps_divers$Hexp,snps_divers$Hobs))

t.test(snps_divers$Hexp,snps_divers$Hobs,pair=T,var.equal=TRUE,alter="greater")


#----------------
# FIS (inbreeding coefficient) and Heterozygosis
#----------------
# Cargar librerías
library(adegenet) # Para manejar datos genéticos
library(hierfstat) # Para cálculos de estadísticos F
library(ggplot2) # Para visualización (opcional)

### Entre la zona Norte y Sur

# Calcular las heterocigosidades observadas (Ho) y esperadas (He) por población
stats <- basic.stats(qmacd_genind_pop)

print(stats) # Para ver el FIS total

# Extraer Ho y He por población
Ho_poblaciones <- apply(stats$Ho, 2, mean, na.rm = TRUE)  # Heterocigosidad observada promedio
He_poblaciones <- apply(stats$Hs, 2, mean, na.rm = TRUE)  # Heterocigosidad esperada promedio

### 
# Crear tabla de resultados
resultados <- data.frame(
  Poblacion = names(Ho_poblaciones),
  Ho_promedio = Ho_poblaciones,
  He_promedio = He_poblaciones,
  FIS = 1 - (Ho_poblaciones/He_poblaciones)
)
print(resultados)
# Calcular FIS por población: FIS = 1 - (Ho / He)
FIS_poblaciones <- 1 - (Ho_poblaciones / He_poblaciones)

# Mostrar resultados
print(FIS_poblaciones)

# crear gráfico
plot_FIS <- ggplot(data.frame(
  Población = names(FIS_poblaciones),
  FIS = FIS_poblaciones,
  N = as.numeric(table(pop(qmacd_genind_pop)))),
  aes(x = Población, y = FIS, size = N, color = FIS)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  labs(title = "FIS vs. Tamaño poblacional", 
       y = "FIS (Endogamia si > 0)", 
       x = "Población", 
       size = "Nº individuos")
plot_FIS



### Entre los 9 sitios

# Calcular las heterocigosidades observadas (Ho) y esperadas (He) por población
stats <- basic.stats(qmacd_genind)

print(stats) # Para ver el FIS total

# Extraer Ho y He por población
Ho_poblaciones <- apply(stats$Ho, 2, mean, na.rm = TRUE)  # Heterocigosidad observada promedio
He_poblaciones <- apply(stats$Hs, 2, mean, na.rm = TRUE)  # Heterocigosidad esperada promedio

### 
# Crear tabla de resultados
resultados <- data.frame(
  Poblacion = names(Ho_poblaciones),
  Ho_promedio = Ho_poblaciones,
  He_promedio = He_poblaciones,
  FIS = 1 - (Ho_poblaciones/He_poblaciones)
)
print(resultados)
# Calcular FIS por población: FIS = 1 - (Ho / He)
FIS_poblaciones <- 1 - (Ho_poblaciones / He_poblaciones)

# Mostrar resultados
print(FIS_poblaciones)

# crear gráfico
plot_FIS <- ggplot(data.frame(
  Población = names(FIS_poblaciones),
  FIS = FIS_poblaciones,
  N = as.numeric(table(pop(qmacd_genind)))),
  aes(x = Población, y = FIS, size = N, color = FIS)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_gradient2(low = "blue", mid = "gray", high = "red", midpoint = 0) +
  labs(title = "FIS vs. Tamaño poblacional", 
       y = "FIS (Endogamia si > 0)", 
       x = "Población", 
       size = "Nº individuos")
plot_FIS

pop(qmacd_vcf) <- pop.metadata$SITE_NAME


#-------------- 
# Alelos privados y totales por sitio 
#--------------

library(ggplot2)
library(dplyr)
library(adegenet)

# Alelos privados por población (sitio)
priv_al <- private_alleles(qmacd_genclone_pop, report = "data.frame", level = "population")

# Guardar tabla de alelos privados
#write.csv(priv_al, "private_allele_ref.gen.qrob.csv", row.names = FALSE)

# Visualización de alelos privados
plot_priv_ale <- ggplot(priv_al) + 
  geom_tile(aes(x = population, y = allele, fill = count)) +
  labs(title = "Alelos privados por población", x = "Población", y = "Alelo")
print(plot_priv_ale)

# Resumen: número de alelos privados por sitio
private_alleles_por_sitio <- priv_al %>%
  group_by(population) %>%
  summarise(private_alleles = sum(count))
print(private_alleles_por_sitio)

# Alelos totales por sitio
tab_alleles <- tab(qmacd_genind_pop, NA.method = "zero")
pops <- pop(qmacd_genind_pop)
total_alleles_por_sitio <- data.frame(
  population = levels(pops),
  total_alleles = sapply(levels(pops), function(pop_name) {
    inds <- which(pops == pop_name)
    sum(colSums(tab_alleles[inds, , drop = FALSE]) > 0)
  })
)
print(total_alleles_por_sitio)


# Obtener la matriz de alelos por individuo
tab_alleles <- tab(qmacd_genind_pop, NA.method = "zero")

# Contar el número de alelos presentes en toda la muestra
total_alleles_global <- sum(colSums(tab_alleles) > 0)

cat("Número total de alelos en toda la muestra:", total_alleles_global, "\n")


# -------------------------------
# Estadisticos de FST
# -------------------------------

#### FST ####  

# Cargar las librerías necesarias
library(hierfstat)
library(vegan)
library(ape)
library(dartR)

# Asegúrate de que las poblaciones estén asignadas al objeto genlight
#pop(qmacd_genlight) <- as.factor(pop.metadata$SITE_NAME)  # Usar la columna SITE como población

# Convertir el objeto genlight a genind (necesario para hierfstat)
#qmacd_genind <- gl2gi(qmacd_genlight, v = 1)  # Convertir genlight a genind

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
library(viridis)

heatmap_plot <- ggplot(pairwise_fst_df, aes(x = Population1, y = Population2, fill = FST)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", FST)), size = 4, color = "white") +
  scale_fill_viridis_c(option = "viridis", name = "FST", direction = 1) +  # Paleta viridis
  scale_x_discrete(labels = c("MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT")) +
  scale_y_discrete(labels = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
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


