#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# run_bench.sh — build the bamstrom benchmark image and run bench.py
#
# Usage:
#   ./bench/run_bench.sh <bam> <bai> [-- <extra bench.py args>]
#
# Examples:
#   ./bench/run_bench.sh /data/full.bam /data/full.bam.bai
#   ./bench/run_bench.sh /data/full.bam /data/full.bam.bai -- --csv /data/results.csv
#   ./bench/run_bench.sh /data/full.bam /data/full.bam.bai -- --no-drop-cache --repeats 1
#
# The BAM/BAI directory is mounted read-only at /data inside the container.
# If --csv points inside /data the CSV file will be written to your host.
# ---------------------------------------------------------------------------

BAM="${1:?Usage: $0 <bam> <bai> [-- extra args]}"
BAI="${2:?Usage: $0 <bam> <bai> [-- extra args]}"
shift 2

# Collect any extra args after --
EXTRA_ARGS=()
if [[ "${1:-}" == "--" ]]; then
    shift
    EXTRA_ARGS=("$@")
fi

IMAGE="bamstrom-bench"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR="$(cd "$(dirname "$BAM")" && pwd)"
BAM_FILE="/data/$(basename "$BAM")"
BAI_FILE="/data/$(basename "$BAI")"

echo "=== Building image: $IMAGE ==="
docker build -t "$IMAGE" "$REPO_ROOT"

echo ""
echo "=== Running benchmark ==="
echo "    BAM : $BAM_FILE (host: $BAM)"
echo "    BAI : $BAI_FILE (host: $BAI)"
echo "    Data: $DATA_DIR → /data (ro)"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo "    Args: ${EXTRA_ARGS[*]}"
echo ""

docker run --rm \
    --privileged \
    -v "$DATA_DIR:/data:ro" \
    "$IMAGE" \
    python3 /app/bench.py "$BAM_FILE" "$BAI_FILE" "${EXTRA_ARGS[@]}"
