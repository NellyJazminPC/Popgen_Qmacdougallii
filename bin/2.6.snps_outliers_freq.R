# Cargar las librerías necesarias
library(vcfR)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)

# Ruta al archivo VCF
vcf_file <- "../data/1.3.assembly_variant_calling/ref_gen_qrob_trim01_1_sorted.vcf"

# Ruta al archivo con los loci de interés
blast_results_file <- "../results/consolidated_snps_with_blast_results.xlsx"

# Leer el archivo VCF
vcf <- read.vcfR(vcf_file)

# Leer el archivo con los loci de interés
blast_results <- read_xlsx(blast_results_file)

# Extraer los loci de interés
loci_of_interest <- blast_results$locus_name  # Ajustar el nombre de la columna si es necesario

# Extraer la información de las variantes
variants <- as.data.frame(vcf@fix)  # Contiene CHROM, POS, ID, REF, ALT, etc.
genotypes <- extract.gt(vcf)        # Matriz de genotipos

# Verificar las primeras filas de las columnas relevantes
cat("Primeros valores de la columna ID en variants:\n")
print(head(variants$ID))
cat("Primeros valores de loci_of_interest:\n")
print(head(loci_of_interest))

# Convertir la matriz de genotipos en un data frame
genotypes_df <- as.data.frame(genotypes)
genotypes_df <- tibble::rownames_to_column(genotypes_df, var = "locus_name")  # Agregar locus_name como columna

# Unir la información de variantes (REF y ALT) con los genotipos
genotypes_with_variants <- genotypes_df %>%
  left_join(variants %>% select(ID, REF, ALT), by = c("locus_name" = "ID"))

# Filtrar para conservar solo los loci de interés
genotypes_with_variants <- genotypes_with_variants %>%
  filter(locus_name %in% loci_of_interest)

# Verificar si hay datos después del filtro
cat("Número de filas después del filtro:\n")
print(nrow(genotypes_with_variants))

# Función para traducir genotipos a nucleótidos
translate_genotype <- function(genotype, ref, alt) {
  if (is.na(genotype)) {
    return(NA)
  }
  alleles <- strsplit(genotype, "/")[[1]]  # Separar los alelos (por ejemplo, "0/1" -> c("0", "1"))
  nucleotides <- sapply(alleles, function(allele) {
    if (allele == "0") {
      return(ref)  # Alelo de referencia
    } else if (allele == "1") {
      return(alt)  # Alelo alternativo
    } else {
      return(NA)   # Manejar casos inesperados
    }
  })
  return(paste(nucleotides, collapse = "/"))  # Combinar los nucleótidos en un formato como "A/T"
}

# Traducir los genotipos a nucleótidos para cada individuo
genotypes_with_variants <- genotypes_with_variants %>%
  mutate(across(
    starts_with(c("CZ_", "MT_", "MC_", "MB_", "CY_", "LS_", "PZ_", "CR_", "IT_")),  # Ajustar prefijos
    ~ mapply(translate_genotype, ., REF, ALT)
  ))

# Convertir el data frame a formato largo (individuos como filas)
variants_long <- genotypes_with_variants %>%
  pivot_longer(
    cols = -c(locus_name, REF, ALT),  # Seleccionar todas las columnas excepto locus_name, REF y ALT
    names_to = "individual",         # Nombres de las columnas se convierten en "individual"
    values_to = "nucleotide"         # Valores de las celdas se convierten en "nucleotide"
  )


# Convertir el data frame a formato ancho (un locus por fila, individuos como columnas)
variants_wide <- variants_long %>%
  pivot_wider(
    names_from = individual,         # Los nombres de los individuos se convierten en columnas
    values_from = nucleotide         # Los valores de las celdas son los nucleótidos
  ) %>%
  relocate(REF, ALT, .after = last_col())  # Mover REF y ALT al final del data frame

# Exportar el data frame a un archivo Excel
write_xlsx(variants_wide, "../results/SNPs_outliers_variants_per_locus.xlsx")

# Mensaje de confirmación
cat("El archivo con las variantes por locus (en formato ancho) se ha exportado a '../results/SNPs_out_variants_per_locus.xlsx'.\n")

# ---------------------------------------------------------------
# Agregar REF y ALT al archivo general de locus
# Descripción: Este script agrega las columnas REF y ALT al archivo de resultados SNPs outliers + BLAST
#              para saber las variantes antes de hacer los gráficos de frecuencias
#              y calcula el tipo de mutación (TRANSITION o TRANSVERSION).
# ---------------------------------------------------------------

# Cargar las librerías necesarias
library(dplyr)
library(readxl)
library(writexl)

# Ruta al archivo de resultados BLAST
blast_results_file <- "../results/consolidated_snps_with_blast_results.xlsx"

# Leer el archivo de resultados BLAST
blast_results <- read_xlsx(blast_results_file)

# Agregar las columnas REF y ALT al archivo de resultados BLAST
blast_results_updated <- blast_results %>%
  left_join(variants_wide %>% select(locus_name, REF, ALT), by = c("locus_name" = "locus_name")) %>%
  mutate(
    mutation_type = case_when(
      (REF == "A" & ALT == "G") | (REF == "G" & ALT == "A") ~ "TRANSITION",
      (REF == "C" & ALT == "T") | (REF == "T" & ALT == "C") ~ "TRANSITION",
      TRUE ~ "TRANSVERSION"
    )
  )

# Exportar el archivo actualizado
write_xlsx(blast_results_updated, "../results/snps_outliers_with_blast_and_variants.xlsx")

# Mensaje de confirmación
cat("El archivo actualizado con REF, ALT y el tipo de mutación se ha exportado a '../results/snps_outliers_with_blast_and_variants.xlsx'.\n")

# ---------------------------------------------------------------
# Plots de frecuencias de nucleótidos
# Descripción: Este script genera gráficos de frecuencias de nucleótidos para cada locus
#              y los guarda en archivos PNG.
# ---------------------------------------------------------------
#--------------------------------
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

#------------------------------
# Graficar SNPs outliers de PCAdapt
#-----------------------------
# Cargar las librerías necesarias
library(dplyr)
library(readxl)
library(dplyr)
library(purrr)  # Cargar purrr para usar map2
library(readxl)

# Ruta al archivo Supplementary_2.xlsx
supplementary_file <- "../doc/Supplementary_2.xlsx"

# Leer la hoja "TS3" del archivo Supplementary_2.xlsx
supplementary_data <- read_xlsx(supplementary_file, sheet = "TS3")

# Filtrar las filas no repetidas por la columna "locus_name_clean"
supplementary_data <- supplementary_data %>%
  distinct(locus_name_clean, .keep_all = TRUE)

# Filtrar las filas que tienen "Yes" en la columna "plot"
filtered_data <- supplementary_data %>%
  filter(plot == "Yes")

# Filtrar las filas donde "parametro" es "pval" y redondear "value" a 6 decimales
filtered_pval <- filtered_data %>%
  filter(parametro == "pval") %>%
  mutate(value = round(value, 6))

# Crear la lista snp_info dinámicamente desde filtered_pval
snp_info <- filtered_pval %>%
  mutate(
    pch = map2(REF, ALT, ~ c(.x, .y))  # Crear el formato para pch como una lista con REF y ALT
  ) %>%
  select(name = locus_name, pval = value, pch) %>%  # Seleccionar las columnas necesarias
  mutate(
    pval = as.character(pval)  # Asegurarse de que pval sea un carácter
  ) %>%
  pmap(~ list(name = ..1, pval = ..2, pch = ..3))  # Crear una lista de listas

# Verificar la estructura de snp_info
str(snp_info)

#-----------------------------
# Graficar SNPs outliers de PCAdapt
#-----------------------------
library(adegenet)

# Lista de SNPs con sus valores de pch y p-value adj
snp_info 

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
  png(png_filename, width = 2400, height = 1200, res = 300)  # Tamaño y resolución ajustados
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

#------------------------------
# Graficar SNPs outliers de FST
#-----------------------------
# Cargar las librerías necesarias
library(dplyr)
library(readxl)
library(purrr)  # Cargar purrr para usar map2
library(adegenet)

# Ruta al archivo Supplementary_2.xlsx
supplementary_file <- "../doc/Supplementary_2.xlsx"

# Leer la hoja "TS3" del archivo Supplementary_2.xlsx
supplementary_data <- read_xlsx(supplementary_file, sheet = "TS3")

# Filtrar las filas no repetidas por la columna "locus_name_clean"
supplementary_data <- supplementary_data %>%
  distinct(locus_name_clean, .keep_all = TRUE)

# Filtrar las filas que tienen "Yes" en la columna "plot"
filtered_data <- supplementary_data %>%
  filter(plot == "Yes")

# Filtrar las filas donde "parametro" es "fst" y redondear "value" a 6 decimales
filtered_fst <- filtered_data %>%
  filter(parametro == "fst") %>%
  mutate(value = round(value, 6))

# Crear la lista snp_info dinámicamente desde filtered_fst
snp_info <- filtered_fst %>%
  mutate(
    pch = map2(REF, ALT, ~ c(.x, .y))  # Crear el formato para pch como una lista con REF y ALT
  ) %>%
  select(name = locus_name, fst = value, pch) %>%  # Seleccionar las columnas necesarias
  mutate(
    fst = as.character(fst)  # Asegurarse de que fst sea un carácter
  ) %>%
  pmap(~ list(name = ..1, fst = ..2, pch = ..3))  # Crear una lista de listas

# Verificar la estructura de snp_info
str(snp_info)

# Iterar sobre cada SNP en la lista
for (snp in snp_info) {
  # Extraer el nombre del SNP, el FST ajustado y los alelos
  snp_name <- snp$name
  snp_fst <- snp$fst
  snp_pch <- snp$pch
  
  # Extraer las frecuencias alélicas del SNP
  snp_data <- tab(temp[[snp_name]])
  freq_snp <- apply(snp_data, 2, function(e) tapply(e, pop(qmacd_genind), mean, na.rm = TRUE))
  
  # Configurar el gráfico
  par(mar = c(5, 5, 4, 2) + 0.1)  # Márgenes ajustados
  matplot(freq_snp, type = "b", pch = snp_pch,
          xlab = "SITE", ylab = "Allele frequency", 
          main = paste("FST value", snp_fst, snp_name),
          xaxt = "n", cex = 1.5)
  axis(side = 1, at = 1:9, lab = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT"))
  
  # Exportar el gráfico como PNG
  png_filename <- paste0("../results/allele_frequency_fst_", snp_name, ".png")
  png(png_filename, width = 2400, height = 1200, res = 300)  # Tamaño y resolución ajustados
  par(mar = c(5, 5, 4, 2) + 0.1)  # Márgenes ajustados para el archivo exportado
  matplot(freq_snp, type = "b", pch = snp_pch,
          xlab = "SITE", ylab = "Allele frequency", 
          main = paste("FST value", snp_fst, snp_name),
          xaxt = "n", cex = 1.5)
  axis(side = 1, at = 1:9, lab = c("CZ", "MT", "MC", "MB", "CY", "LS", "PZ", "CR", "IT"))
  dev.off()  # Cerrar el dispositivo gráfico
  
  # Mensaje de confirmación
  cat("Gráfico guardado para SNP:", snp_name, "en", png_filename, "\n")
}