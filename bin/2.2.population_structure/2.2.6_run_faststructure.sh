#!/usr/bin/env bash

set -euo pipefail

# Run fastStructure for K = 1–10 using the simple and logistic priors.
#
# Input:
#   data/structure_formats/qmacd_ref_gen_rob.bed
#   data/structure_formats/qmacd_ref_gen_rob.bim
#   data/structure_formats/qmacd_ref_gen_rob.fam
#
# Output:
#   data/1.4.population_structure/faststructure_output/
#
# Default parameters:
#   K range: 1–10
#   Random seed: 20
#   Priors: simple and logistic
#
# The fastStructure installation directory must contain:
#   structure.py
#   chooseK.py
#
# Example:
#
#   FASTSTRUCTURE_DIR=/path/to/fastStructure \
#   PYTHON_BIN=python \
#     bash bin/2.2.population_structure/2.2.6_run_faststructure.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

PLINK_PREFIX="${REPO_ROOT}/data/structure_formats/qmacd_ref_gen_rob"
OUTPUT_DIR="${REPO_ROOT}/data/1.4.population_structure/faststructure_output"

FASTSTRUCTURE_DIR="${FASTSTRUCTURE_DIR:-}"
PYTHON_BIN="${PYTHON_BIN:-python}"
SEED="${SEED:-20}"
K_MIN="${K_MIN:-1}"
K_MAX="${K_MAX:-10}"

for variable_name in SEED K_MIN K_MAX; do
    value="${!variable_name}"

    if [[ ! "$value" =~ ^[0-9]+$ ]]; then
        echo "${variable_name} must be a non-negative integer." >&2
        exit 1
    fi
done

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

if [[ -z "$FASTSTRUCTURE_DIR" ]]; then
    echo "FASTSTRUCTURE_DIR is not set." >&2
    echo "Set it to the directory containing structure.py and chooseK.py." >&2
    exit 1
fi

if [[ "$FASTSTRUCTURE_DIR" != /* ]]; then
    FASTSTRUCTURE_DIR="${REPO_ROOT}/${FASTSTRUCTURE_DIR#./}"
fi

STRUCTURE_SCRIPT="${FASTSTRUCTURE_DIR}/structure.py"
CHOOSE_K_SCRIPT="${FASTSTRUCTURE_DIR}/chooseK.py"

for script in "$STRUCTURE_SCRIPT" "$CHOOSE_K_SCRIPT"; do
    if [[ ! -f "$script" ]]; then
        echo "Required fastStructure script not found: $script" >&2
        exit 1
    fi
done

if [[ "$PYTHON_BIN" == */* ]]; then
    if [[ ! -x "$PYTHON_BIN" ]]; then
        echo "Python executable not found: $PYTHON_BIN" >&2
        exit 1
    fi
elif ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
    echo "Python was not found in PATH." >&2
    echo "Set PYTHON_BIN to the required Python executable." >&2
    exit 1
else
    PYTHON_BIN="$(command -v "$PYTHON_BIN")"
fi

mkdir -p "$OUTPUT_DIR"

echo "Running fastStructure"
echo "Input prefix: $PLINK_PREFIX"
echo "Output directory: $OUTPUT_DIR"
echo "K range: ${K_MIN}–${K_MAX}"
echo "Random seed: $SEED"

for prior in simple logistic; do
    OUTPUT_PREFIX="${OUTPUT_DIR}/qmacd_ref_gen_rob.${prior}"

    echo
    echo "Prior: $prior"

    for (( K = K_MIN; K <= K_MAX; K++ )); do
        echo "Running K=${K}..."

        "$PYTHON_BIN" "$STRUCTURE_SCRIPT" \
            -K "$K" \
            --input="$PLINK_PREFIX" \
            --output="$OUTPUT_PREFIX" \
            --full \
            --seed="$SEED" \
            --prior="$prior" \
            --format=bed
    done

    CHOOSE_K_FILE="${OUTPUT_DIR}/chooseK_qmacd_ref_gen_rob.${prior}.txt"

    "$PYTHON_BIN" "$CHOOSE_K_SCRIPT" \
        --input="$OUTPUT_PREFIX" \
        > "$CHOOSE_K_FILE"

    echo
    echo "Model-selection summary for prior '${prior}':"
    cat "$CHOOSE_K_FILE"
done

echo
echo "fastStructure analysis completed."
