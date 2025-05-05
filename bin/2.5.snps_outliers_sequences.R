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

# ---------------------------------------------------------------
# Descripción: Contar caracteres ambiguos en cada secuencia y exportar el resultado.
# ---------------------------------------------------------------

# Contar caracteres ambiguos en cada secuencia
unique_snps_with_sequences <- unique_snps_with_sequences %>%
  mutate(
    amb_R = str_count(sequence, "R"),  # Contar R
    amb_Y = str_count(sequence, "Y"),  # Contar Y
    amb_S = str_count(sequence, "S"),  # Contar S
    amb_W = str_count(sequence, "W"),  # Contar W
    amb_K = str_count(sequence, "K"),  # Contar K
    amb_M = str_count(sequence, "M"),  # Contar M
    total_ambiguities = amb_R + amb_Y + amb_S + amb_W + amb_K + amb_M  # Sumar todas las ambigüedades
  )

# Exportar el data frame actualizado con las secuencias y los conteos como archivo Excel
write_xlsx(unique_snps_with_sequences, "../results/consolidated_snps_with_sequences_and_ambiguities.xlsx")

# Mensaje de confirmación
cat("El archivo con las secuencias, los conteos de caracteres ambiguos y el total de ambigüedades se ha exportado a '../results/consolidated_snps_with_sequences_and_ambiguities.xlsx'.\n")

# ---------------------------------------------------------------
# Script: Sustituir ambigüedades en las secuencias del archivo FASTA
# Descripción: Este script procesa el archivo FASTA, sustituye los caracteres
#              ambiguos por los nucleótidos especificados y guarda un nuevo archivo.
# ---------------------------------------------------------------

# Cargar el archivo de secuencias en formato FASTA
fasta_lines <- readLines("../results/locus_names_for_sequences.txt")

# Identificar líneas con secuencias (no encabezados)
sequence_lines <- !grepl("^>", fasta_lines)

# Sustituir las ambigüedades en las secuencias
fasta_lines[sequence_lines] <- fasta_lines[sequence_lines] %>%
  str_replace_all("R", "A") %>%
  str_replace_all("Y", "C") %>%
  str_replace_all("S", "G") %>%
  str_replace_all("W", "T") %>%
  str_replace_all("K", "G") %>%
  str_replace_all("M", "A")

# Guardar el archivo con las secuencias corregidas
writeLines(fasta_lines, "../results/locus_names_for_sequences_wo_amb.txt")

# Mensaje de confirmación
cat("El archivo con las secuencias corregidas se ha exportado a '../results/locus_names_for_sequences_wo_amb.txt'.\n")


# ---------------------------------------------------------------
# Script: 2.6.blast_results.R
# Descripción: Este script procesa los resultados de BLAST, extrae la información
#              relevante y la une al data frame de SNPs únicos con secuencias.
# ---------------------------------------------------------------

# Cargar las librerías necesarias
library(dplyr)
library(stringr)
library(readr)

# Ruta a la carpeta con los resultados BLAST
blast_results_dir <- "../results/blast_results/"

# Leer los archivos de resultados BLAST
blast_files <- list.files(blast_results_dir, pattern = "*.txt", full.names = TRUE)

# Crear un data frame vacío para almacenar los resultados procesados
blast_summary <- data.frame(
  locus_name = character(),
  best_hit = character(),
  e_value = character(),
  perc_identity = character(),
  stringsAsFactors = FALSE
)

# Procesar cada archivo de resultados BLAST
for (file in blast_files) {
  # Extraer el nombre del locus desde el nombre del archivo
  locus_name <- gsub("_blast.txt", "", basename(file))
  
  # Leer el contenido del archivo
  blast_content <- readLines(file)
  
  # Buscar la sección ALIGNMENTS
  alignments_start <- grep("^ALIGNMENTS", blast_content)
  if (length(alignments_start) > 0) {
    # Extraer la línea del mejor hit (primera línea después de ALIGNMENTS)
    best_hit_line <- blast_content[alignments_start + 1]
    best_hit <- str_remove(best_hit_line, "^>")  # Quitar el símbolo ">"
    best_hit <- str_trim(best_hit)  # Eliminar espacios en blanco
    
    # Buscar el E-value en la sección ALIGNMENTS
    e_value_line <- grep("Expect =", blast_content, value = TRUE)
    e_value <- str_extract(e_value_line[1], "(?<=Expect = )\\S+")  # Extraer el E-value
    
    # Buscar el porcentaje de identidad en la sección ALIGNMENTS
    perc_identity_line <- grep("Identities =", blast_content, value = TRUE)
    perc_identity <- str_extract(perc_identity_line[1], "(?<=\\()\\d+%")  # Extraer el porcentaje de identidad
    
    # Agregar los datos al resumen
    blast_summary <- rbind(blast_summary, data.frame(
      locus_name = locus_name,
      best_hit = best_hit,
      e_value = e_value,
      perc_identity = perc_identity,
      stringsAsFactors = FALSE
    ))
  } else {
    # Si no hay hits, agregar un registro vacío
    blast_summary <- rbind(blast_summary, data.frame(
      locus_name = locus_name,
      best_hit = NA,
      e_value = NA,
      perc_identity = NA,
      stringsAsFactors = FALSE
    ))
  }
}

# Unir los resultados BLAST al data frame unique_snps_with_sequences
unique_snps_with_blast <- unique_snps_with_sequences %>%
  left_join(blast_summary, by = c("locus_name_clean" = "locus_name"))


# Revisar si hay NAs
na_files <- blast_summary %>% filter(is.na(e_value)) %>% pull(locus_name)
print(na_files)

# Exportar el data frame actualizado con los resultados BLAST
write_xlsx(unique_snps_with_blast, "../results/consolidated_snps_with_blast_results.xlsx")

# Mensaje de confirmación
cat("El archivo con los resultados BLAST incorporados se ha exportado a '../results/consolidated_snps_with_blast_results.xlsx'.\n")

