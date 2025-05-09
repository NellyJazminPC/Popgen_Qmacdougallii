#!/bin/bash

# Archivo con las secuencias en formato FASTA
FASTA_FILE="../results/locus_names_for_sequences_wo_amb.txt"

# Carpeta para guardar los resultados
OUTPUT_DIR="../results/blast_results"
mkdir -p "$OUTPUT_DIR"

# Base de datos de NCBI para la búsqueda (por ejemplo, nt)
DB="nt"

# Número máximo de hits por secuencia
MAX_HITS=5

# Filtro para restringir la búsqueda a Magnoliopsida
ENTREZ_QUERY="Viridiplantae[ORGANISM]"

# Arreglo para registrar las secuencias que fallaron
FAILED_SEQUENCES=()

# Contador de secuencias procesadas
TOTAL_SEQUENCES=0

# Iterar sobre cada secuencia en el archivo FASTA
while read -r line; do
  if [[ $line == ">"* ]]; then
    # Extraer el nombre de la secuencia
    SEQ_NAME=$(echo $line | sed 's/>//')
    ((TOTAL_SEQUENCES++))  # Incrementar el contador de secuencias
  else
    # Realizar la búsqueda BLAST en línea usando la API de NCBI
    echo "Realizando BLAST para la secuencia: $SEQ_NAME"

    # Paso 1: Enviar la consulta y obtener el RID
    RID=$(curl -s -X POST "https://blast.ncbi.nlm.nih.gov/Blast.cgi" \
      -d "CMD=Put" \
      -d "PROGRAM=blastn" \
      -d "DATABASE=$DB" \
      -d "QUERY=$line" \
      -d "HITLIST_SIZE=$MAX_HITS" \
      -d "ENTREZ_QUERY=$ENTREZ_QUERY" \
      | grep "RID =" | sed 's/.*RID = \(.*\)/\1/')

    if [[ -z "$RID" ]]; then
      echo "Error: No se pudo obtener un RID para la secuencia $SEQ_NAME"
      FAILED_SEQUENCES+=("$SEQ_NAME")
      continue
    fi

    echo "RID obtenido: $RID. Esperando resultados..."

    # Paso 2: Esperar y recuperar los resultados
    STATUS="WAITING"
    while [[ "$STATUS" == "WAITING" ]]; do
      sleep 60  # Esperar 60 segundos antes de verificar el estado
      STATUS=$(curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_OBJECT=SearchInfo" \
        | grep "Status=" | sed 's/.*Status=\(.*\)/\1/')
    done

    if [[ "$STATUS" == "READY" ]]; then
      # Descargar los resultados en formato tabular
      curl -s "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_TYPE=Text" \
        -o "$OUTPUT_DIR/${SEQ_NAME}_blast.txt"
      echo "Resultados guardados en $OUTPUT_DIR/${SEQ_NAME}_blast.txt"
    else
      echo "Error: La búsqueda para $SEQ_NAME falló con estado $STATUS"
      FAILED_SEQUENCES+=("$SEQ_NAME")
    fi
  fi
done < "$FASTA_FILE"

# Archivo para guardar el resumen
SUMMARY_FILE="../results/blast_summary.txt"

# Limpiar el archivo de resumen si ya existe
> "$SUMMARY_FILE"

# Imprimir resumen de secuencias que fallaron
if [[ ${#FAILED_SEQUENCES[@]} -gt 0 ]]; then
  echo "Resumen de errores:" | tee -a "$SUMMARY_FILE"
  echo "Las siguientes secuencias fallaron:" | tee -a "$SUMMARY_FILE"
  for SEQ in "${FAILED_SEQUENCES[@]}"; do
    echo "- $SEQ" | tee -a "$SUMMARY_FILE"
  done
  echo "Número total de secuencias que fallaron: ${#FAILED_SEQUENCES[@]} del total de secuencias: $TOTAL_SEQUENCES" | tee -a "$SUMMARY_FILE"
else
  echo "Todas las búsquedas BLAST se completaron exitosamente." | tee -a "$SUMMARY_FILE"
  echo "Número total de secuencias procesadas: $TOTAL_SEQUENCES" | tee -a "$SUMMARY_FILE"
fi

echo "Las búsquedas BLAST se han completado. Los resultados están en la carpeta '$OUTPUT_DIR'." | tee -a "$SUMMARY_FILE"