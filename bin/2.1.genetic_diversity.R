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

pop(qmacd_genind) <- pop.metadata$SITE_NAME

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
# FIS
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



# Diversidad equivalente a π (Hs = heterocigosidad esperada)
pi_por_poblacion <- apply(stats$Hs, 2, mean, na.rm = TRUE)
print(pi_por_poblacion)



# Función aproximada de Tajima's D para datos genotípicos
tajimaD_approx <- function(genind_obj) {
  # Requiere el paquete pegas
  if (!require("pegas")) install.packages("pegas")
  library(pegas)
  
  # Convertir a formato DNAbin (aproximación)
  dna <- genind2DNAbin(genind_obj)
  
  # Calcular Tajima's D (puede dar advertencias)
  tajima.test(dna)
}

# Ejecutar con precaución (puede no ser perfecto para SNPs)
tajima_result <- tajimaD_approx(qmacd_genind)
print(tajima_result)