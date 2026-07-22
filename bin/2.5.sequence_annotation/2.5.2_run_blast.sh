#!/usr/bin/env bash

# Run remote NCBI BLAST searches for the candidate-locus sequences.
#
# The retained analysis restricts searches to Magnoliopsida because this
# taxonomic scope produced more biologically relevant matches and fewer
# nonspecific results than the broader Viridiplantae search.
#
# The FASTA file used here contains one sequence per line.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

FASTA_FILE="${FASTA_FILE:-$REPO_ROOT/data/1.7.sequence_annotation/candidate_locus_blast_queries.fasta}"
OUTPUT_DIR="${OUTPUT_DIR:-$REPO_ROOT/results/blast_results_Magnoliopsida}"
SUMMARY_FILE="${SUMMARY_FILE:-$REPO_ROOT/results/blast_summary_Magnoliopsida.txt}"

DB="${DB:-nt}"
MAX_HITS="${MAX_HITS:-5}"
ENTREZ_QUERY="${ENTREZ_QUERY:-Magnoliopsida[ORGANISM]}"

FAILED_SEQUENCES=()
TOTAL_SEQUENCES=0
SEQ_NAME=""

if ! command -v curl >/dev/null 2>&1; then
  echo "ERROR: curl is required but was not found." >&2
  exit 1
fi

if [[ ! -f "$FASTA_FILE" ]]; then
  echo "ERROR: FASTA input file was not found:" >&2
  echo "$FASTA_FILE" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

echo "Input FASTA:       $FASTA_FILE"
echo "Output directory:  $OUTPUT_DIR"
echo "Database:          $DB"
echo "Taxonomic filter:  $ENTREZ_QUERY"
echo "Maximum hits:      $MAX_HITS"
echo

while IFS= read -r line || [[ -n "$line" ]]; do
  [[ -z "$line" ]] && continue

  if [[ "$line" == ">"* ]]; then
    SEQ_NAME="${line#>}"
    ((TOTAL_SEQUENCES++))
    continue
  fi

  if [[ -z "$SEQ_NAME" ]]; then
    echo "ERROR: a sequence was found before its FASTA header." >&2
    exit 1
  fi

  echo "Submitting BLAST search for: $SEQ_NAME"

  RID=$(
    curl -sS -X POST "https://blast.ncbi.nlm.nih.gov/Blast.cgi" \
      -d "CMD=Put" \
      -d "PROGRAM=blastn" \
      -d "DATABASE=$DB" \
      -d "QUERY=$line" \
      -d "HITLIST_SIZE=$MAX_HITS" \
      -d "ENTREZ_QUERY=$ENTREZ_QUERY" \
      | grep "RID =" \
      | sed 's/.*RID = \(.*\)/\1/'
  )

  if [[ -z "$RID" ]]; then
    echo "ERROR: no RID was obtained for $SEQ_NAME."
    FAILED_SEQUENCES+=("$SEQ_NAME")
    continue
  fi

  echo "RID obtained: $RID. Waiting for results..."

  STATUS="WAITING"

  while [[ "$STATUS" == "WAITING" ]]; do
    sleep 60

    STATUS=$(
      curl -sS \
        "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_OBJECT=SearchInfo" \
        | grep "Status=" \
        | sed 's/.*Status=\(.*\)/\1/'
    )
  done

  if [[ "$STATUS" == "READY" ]]; then
    output_file="$OUTPUT_DIR/${SEQ_NAME}_blast.txt"

    curl -sS \
      "https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=$RID&FORMAT_TYPE=Text" \
      -o "$output_file"

    echo "Results saved to: $output_file"
  else
    echo "ERROR: the search for $SEQ_NAME ended with status $STATUS."
    FAILED_SEQUENCES+=("$SEQ_NAME")
  fi

  echo
done < "$FASTA_FILE"

: > "$SUMMARY_FILE"

if [[ ${#FAILED_SEQUENCES[@]} -gt 0 ]]; then
  echo "BLAST search error summary:" | tee -a "$SUMMARY_FILE"
  echo "The following sequences failed:" | tee -a "$SUMMARY_FILE"

  for sequence_name in "${FAILED_SEQUENCES[@]}"; do
    echo "- $sequence_name" | tee -a "$SUMMARY_FILE"
  done

  echo \
    "Failed sequences: ${#FAILED_SEQUENCES[@]} of $TOTAL_SEQUENCES." \
    | tee -a "$SUMMARY_FILE"
else
  echo "All BLAST searches completed successfully." \
    | tee -a "$SUMMARY_FILE"

  echo "Total sequences processed: $TOTAL_SEQUENCES." \
    | tee -a "$SUMMARY_FILE"
fi

echo "BLAST results directory: $OUTPUT_DIR" \
  | tee -a "$SUMMARY_FILE"
