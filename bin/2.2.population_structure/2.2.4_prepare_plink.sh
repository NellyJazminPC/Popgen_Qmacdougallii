#!/usr/bin/env bash

set -euo pipefail

# Prepare PLINK files for population-structure and outlier analyses.
#
# Prerequisite
# ------------
# The final filtered VCF was exported with the TASSEL v5.2.95 graphical
# interface in PLINK PED/MAP format, producing:
#
#   qmacd_ref_gen_rob.plk.ped
#   qmacd_ref_gen_rob.plk.map
#
# This script generates:
#
#   qmacd_ref_gen_rob.raw
#       Additive/dominance genotype matrix generated with --recodeAD.
#
#   qmacd_ref_gen_rob.bed
#   qmacd_ref_gen_rob.bim
#   qmacd_ref_gen_rob.fam
#       Binary PLINK files used by ADMIXTURE, fastStructure, and pcadapt.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
DATA_DIR="${REPO_ROOT}/data/structure_formats"

INPUT_PREFIX="${DATA_DIR}/qmacd_ref_gen_rob.plk"
OUTPUT_PREFIX="${DATA_DIR}/qmacd_ref_gen_rob"

# By default, use a PLINK executable available in PATH.
# A different executable can be supplied through PLINK_BIN.
PLINK_BIN="${PLINK_BIN:-plink}"

for extension in ped map; do
    input_file="${INPUT_PREFIX}.${extension}"

    if [[ ! -f "$input_file" ]]; then
        echo "Missing TASSEL-generated input file: $input_file" >&2
        exit 1
    fi
done

if [[ "$PLINK_BIN" == */* ]]; then
    if [[ ! -x "$PLINK_BIN" ]]; then
        echo "PLINK executable not found or not executable: $PLINK_BIN" >&2
        exit 1
    fi
elif ! command -v "$PLINK_BIN" >/dev/null 2>&1; then
    echo "PLINK was not found in PATH." >&2
    echo "Set PLINK_BIN to the path of the PLINK executable." >&2
    exit 1
fi

echo "Creating additive/dominance genotype matrix..."

"$PLINK_BIN" \
    --file "$INPUT_PREFIX" \
    --noweb \
    --recodeAD \
    --out "$OUTPUT_PREFIX"

echo "Creating binary PLINK files..."

"$PLINK_BIN" \
    --file "$INPUT_PREFIX" \
    --noweb \
    --make-bed \
    --out "$OUTPUT_PREFIX"

echo "PLINK preparation completed."
echo "Output prefix: $OUTPUT_PREFIX"
