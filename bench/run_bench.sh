#!/usr/bin/env bash
set -euo pipefail

# ---------------------------------------------------------------------------
# run_bench.sh — build the bamstorm benchmark image and run bench.py
#
# Usage:
#   ./bench/run_bench.sh <bam> <bai> [--csv <file>] [-- <extra bench.py args>]
#
# Options:
#   --csv <file>   Write results CSV to <file> on the host (any path).
#                  The directory containing <file> is mounted rw into the container.
#
# Examples:
#   ./bench/run_bench.sh /data/full.bam /data/full.bam.bai
#   ./bench/run_bench.sh /data/full.bam /data/full.bam.bai --csv ./results.csv
#   ./bench/run_bench.sh /data/full.bam /data/full.bam.bai --csv /results/out.csv -- --repeats 1
# ---------------------------------------------------------------------------

BAM="${1:?Usage: $0 <bam> <bai> [--csv <file>] [-- extra args]}"
BAI="${2:?Usage: $0 <bam> <bai> [--csv <file>] [-- extra args]}"
shift 2

# Parse --csv and -- from remaining args
CSV_HOST=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "${1}" in
        --csv)
            CSV_HOST="${2:?--csv requires a filename}"
            shift 2
            ;;
        --)
            shift
            EXTRA_ARGS=("$@")
            break
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

IMAGE="bamstorm-bench"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DATA_DIR="$(cd "$(dirname "$BAM")" && pwd)"
BAM_FILE="/data/$(basename "$BAM")"
BAI_FILE="/data/$(basename "$BAI")"

# Build CSV mount and container path if requested
CSV_MOUNT_ARGS=()
CSV_BENCH_ARG=()
if [[ -n "$CSV_HOST" ]]; then
    CSV_HOST="$(realpath -m "$CSV_HOST")"   # resolve to absolute path
    CSV_OUT_DIR="$(dirname "$CSV_HOST")"
    mkdir -p "$CSV_OUT_DIR"
    CSV_BENCH_ARG=("--csv" "/out/$(basename "$CSV_HOST")")
    CSV_MOUNT_ARGS=(-v "$CSV_OUT_DIR:/out")
fi

echo "=== Building image: $IMAGE ==="
docker build -t "$IMAGE" "$REPO_ROOT"

echo ""
echo "=== Running benchmark ==="
echo "    BAM : $BAM_FILE (host: $BAM)"
echo "    BAI : $BAI_FILE (host: $BAI)"
echo "    Data: $DATA_DIR → /data (ro)"
[[ -n "$CSV_HOST" ]] && echo "    CSV : $CSV_HOST"
[[ ${#EXTRA_ARGS[@]} -gt 0 ]] && echo "    Args: ${EXTRA_ARGS[*]}"
echo ""

docker run --rm \
    --privileged \
    -v "$DATA_DIR:/data:ro" \
    "${CSV_MOUNT_ARGS[@]}" \
    "$IMAGE" \
    python3 /app/bench.py "$BAM_FILE" "$BAI_FILE" \
        "${CSV_BENCH_ARG[@]}" \
        "${EXTRA_ARGS[@]}"
