#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

INPUT_FILE="${REPO_ROOT}/data/1.6.snp_outlier_detection/bayescan_input/bayescan_qmacd_ref_gen_qrob_9pop"
OUTPUT_DIR="${REPO_ROOT}/data/1.6.snp_outlier_detection/bayescan_output/bayescan_qmacd_ref_gen_qrob_9pop"
SOFTWARE_DIR="${REPO_ROOT}/bin/software/bayescan_distributed_2.01/binaries"

detect_threads() {
  case "$(uname -s)" in
    Darwin)
      sysctl -n hw.ncpu
      ;;
    Linux)
      getconf _NPROCESSORS_ONLN
      ;;
    *)
      echo 1
      ;;
  esac
}

case "$(uname -s)" in
  Darwin)
    DEFAULT_BAYESCAN_BIN="${SOFTWARE_DIR}/BayeScan2.0_macos64bits"
    ;;
  Linux)
    DEFAULT_BAYESCAN_BIN="${SOFTWARE_DIR}/BayeScan2.0_linux64bits"
    ;;
  *)
    echo "Unsupported operating system: $(uname -s)" >&2
    echo "Set BAYESCAN_BIN to a compatible executable." >&2
    exit 1
    ;;
esac

BAYESCAN_BIN="${BAYESCAN_BIN:-${DEFAULT_BAYESCAN_BIN}}"
THREADS="${THREADS:-$(detect_threads)}"

if [[ ! -s "${INPUT_FILE}" ]]; then
  echo "BayeScan input file not found or empty:" >&2
  echo "  ${INPUT_FILE}" >&2
  exit 1
fi

if [[ ! -f "${BAYESCAN_BIN}" ]]; then
  echo "BayeScan executable not found:" >&2
  echo "  ${BAYESCAN_BIN}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

command=(
  "${BAYESCAN_BIN}"
  "${INPUT_FILE}"
  -od "${OUTPUT_DIR}"
  -threads "${THREADS}"
)

echo "BayeScan nine-site analysis"
echo "  Input:   ${INPUT_FILE}"
echo "  Output:  ${OUTPUT_DIR}"
echo "  Threads: ${THREADS}"

printf '  Command:'
printf ' %q' "${command[@]}"
printf '\n'

if [[ "${DRY_RUN:-0}" == "1" ]]; then
  echo "Dry run completed. BayeScan was not executed."
  exit 0
fi

"${command[@]}"

echo "BayeScan nine-site analysis completed."
