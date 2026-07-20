#!/usr/bin/env bash

set -euo pipefail

# Run ADMIXTURE with cross-validation for K = 1–10.
#
# Input:
#   data/structure_formats/qmacd_ref_gen_rob.bed
#   data/structure_formats/qmacd_ref_gen_rob.bim
#   data/structure_formats/qmacd_ref_gen_rob.fam
#
# Retained outputs:
#   data/1.4.population_structure/admixture_output/
#
# Default parameters:
#   K range: 1–10
#   Threads: 8
#   Random seed: 12345
#
# Example with custom parameters:
#
#   THREADS=4 K_MIN=1 K_MAX=5 \
#     bash bin/2.2.population_structure/2.2.5_run_admixture.sh
#
# A specific ADMIXTURE executable can be provided through ADMIXTURE_BIN.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PLINK_PREFIX="${REPO_ROOT}/data/structure_formats/qmacd_ref_gen_rob"
OUTPUT_DIR="${REPO_ROOT}/data/1.4.population_structure/admixture_output"

ADMIXTURE_BIN="${ADMIXTURE_BIN:-admixture}"
THREADS="${THREADS:-8}"
SEED="${SEED:-12345}"
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-10}"

for variable_name in THREADS SEED K_MIN K_MAX; do
    value="${!variable_name}"

    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        echo "${variable_name} must be a non-negative integer." >&2
        exit 1
    fi
done

if (( THREADS < 1 )); then
    echo "THREADS must be at least 1." >&2
    exit 1
fi

if (( K_MIN < 1 || K_MAX < K_MIN )); then
    echo "Invalid K range: K_MIN=${K_MIN}, K_MAX=${K_MAX}." >&2
    exit 1
fi

for extension in bed bim fam; do
    input_file="${PLINK_PREFIX}.${extension}"

    if [[ ! -f "$input_file" ]]; then
        echo "Missing PLINK input file: $input_file" >&2
        exit 1
    fi
done

if [[ "$ADMIXTURE_BIN" == */* ]]; then
    if [[ "$ADMIXTURE_BIN" != /* ]]; then
        ADMIXTURE_BIN="${REPO_ROOT}/${ADMIXTURE_BIN#./}"
    fi

    if [[ ! -x "$ADMIXTURE_BIN" ]]; then
        echo "ADMIXTURE executable not found: $ADMIXTURE_BIN" >&2
        exit 1
    fi
elif ! command -v "$ADMIXTURE_BIN" >/dev/null 2>&1; then
    echo "ADMIXTURE was not found in PATH." >&2
    echo "Set ADMIXTURE_BIN to the executable path." >&2
    exit 1
else
    ADMIXTURE_BIN="$(command -v "$ADMIXTURE_BIN")"
fi

mkdir -p "$OUTPUT_DIR"

CHOOSE_K_FILE="${OUTPUT_DIR}/chooseK.txt"
: > "$CHOOSE_K_FILE"

cd "$OUTPUT_DIR"

echo "Running ADMIXTURE for K = ${K_MIN}–${K_MAX}"
echo "Threads: $THREADS"
echo "Random seed: $SEED"

for (( K = K_MIN; K <= K_MAX; K++ )); do
    LOG_FILE="log${K}.out"

    echo
    echo "Running K=${K}..."

    "$ADMIXTURE_BIN" \
        --cv \
        "-j${THREADS}" \
        "-s${SEED}" \
        "${PLINK_PREFIX}.bed" \
        "$K" \
        | tee "$LOG_FILE"

    grep "CV error" "$LOG_FILE" >> "$CHOOSE_K_FILE"
done

echo
echo "ADMIXTURE analysis completed."
cat "$CHOOSE_K_FILE"
