# ---------------------------------------------------------------
# Script: 2.5.snps_outliers_sequences.R
# Descripción: Este script procesa el archivo generado en el script 2.4,
#              extrae los nombres de los loci únicos, limpia los caracteres
#              adicionales en los nombres, los formatea en estilo FASTA y
#              los exporta como un archivo .txt para su posterior uso.
# ---------------------------------------------------------------

# Cargar la librería necesaria
library(readxl)

# Cargar el archivo generado en el script 2.4
unique_snps <- read_xlsx("../results/consolidated_snps_results_with_shared_info_unique.xlsx")

# Extraer solo la columna locus_name
locus_names <- unique_snps$locus_name

# Quitar los caracteres a partir de "_pos" en cada línea
locus_names_clean <- sub("_pos.*", "", locus_names)

# Contar cuántos duplicados hay
num_duplicates <- sum(duplicated(locus_names_clean))
cat("Número de duplicados encontrados:", num_duplicates, "\n")

# Eliminar duplicados
locus_names_unique <- unique(locus_names_clean)

# Continuar con el procesamiento
# Agregar un ">" al inicio de cada línea
locus_names_fasta <- paste0(">", locus_names_unique)

# Exportar los resultados como un archivo .txt
writeLines(locus_names_fasta, "../results/locus_names_for_sequences.txt")

# Mensaje de confirmación
cat("El archivo con los nombres de los loci únicos procesados se ha exportado a '../results/locus_names_for_sequences.txt'.\n")

# ---------------------------------------------------------------
# Descripción: Agrega una columna con los nombres de loci sin "_pos*",
#              y empata las secuencias desde un archivo FASTA con los loci.
# ---------------------------------------------------------------

# Cargar las librerías necesarias
library(readxl)
library(dplyr)
library(stringr)
library(writexl)

# Cargar el archivo generado en el script 2.4
#unique_snps <- read_xlsx("../results/consolidated_snps_results_with_shared_info_unique.xlsx")

# Agregar una columna con los nombres de loci sin "_pos*"
unique_snps <- unique_snps %>%
  mutate(locus_name_clean = sub("_pos.*", "", locus_name))

# Cargar el archivo de secuencias en formato FASTA
fasta_lines <- readLines("../results/locus_names_for_sequences.txt")

# Procesar el archivo FASTA para obtener un data frame con locus_name y secuencia
fasta_df <- data.frame(
  locus_name = sub(">", "", fasta_lines[grep("^>", fasta_lines)]),  # Extraer nombres de loci
  sequence = fasta_lines[!grepl("^>", fasta_lines)],               # Extraer secuencias
  stringsAsFactors = FALSE
)

# Unir las secuencias al data frame unique_snps
unique_snps_with_sequences <- unique_snps %>%
  left_join(fasta_df, by = c("locus_name_clean" = "locus_name"))

# Exportar el data frame actualizado con las secuencias como archivo Excel
write_xlsx(unique_snps_with_sequences, "../results/consolidated_snps_with_sequences.xlsx")

# Mensaje de confirmación
cat("El archivo con las secuencias incorporadas se ha exportado a '../results/consolidated_snps_with_sequences.xlsx'.\n")
