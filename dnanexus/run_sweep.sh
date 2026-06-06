#!/bin/bash
# Submit one benchmark job per instance type and collect job IDs.
#
# Usage:
#   ./run_sweep.sh \
#     --bam  file-xxxx \
#     --bai  file-yyyy \
#     --image docker.io/yourname/bamstorm-bench:latest \
#     [--threads 2,4,8,16,32,64] \
#     [--repeats 3]
set -euo pipefail

BAM=""
BAI=""
IMAGE=""
THREADS="1,2,4,8,16,32,64,96"
REPEATS=3

# Instance types to sweep: name → vCPUs (for reference)
# Adjust this list to fit your budget and the thread counts you want to test.
INSTANCES=(
    "mem2_ssd1_v2_x96"
)

usage() {
    cat <<EOF
Usage: $0 --bam <file-id> --bai <file-id> --image <docker-image> [OPTIONS]

Required:
  --bam    FILE-ID   DNAnexus file ID of the BAM file
  --bai    FILE-ID   DNAnexus file ID of the BAI file
  --image  IMAGE     Docker image reference

Optional:
  --threads  LIST    Comma-separated thread counts  (default: 1,2,4,8,16,32,64,96)
  --repeats  N       Repetitions per config         (default: $REPEATS)
  -h                 Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam)    BAM="$2";     shift 2 ;;
        --bai)    BAI="$2";     shift 2 ;;
        --image)  IMAGE="$2";   shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --repeats) REPEATS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$BAM" || -z "$BAI" || -z "$IMAGE" ]]; then
    echo "Error: --bam, --bai, and --image are required."
    usage
    exit 1
fi

echo "Submitting sweep: ${#INSTANCES[@]} instance types"
echo "BAM     : $BAM"
echo "Threads : $THREADS"
echo "Repeats : $REPEATS"
echo ""

JOB_IDS=()

for INST in "${INSTANCES[@]}"; do
    JOB_ID=$(dx run bamstorm_bench \
        -i bam_file="$BAM" \
        -i bai_file="$BAI" \
        -i docker_image="$IMAGE" \
        -i threads="$THREADS" \
        -i repeats="$REPEATS" \
        --instance-type "$INST" \
        --name "bamstorm_bench_${INST}" \
        --brief \
        --yes)
    JOB_IDS+=("$JOB_ID")
    echo "Submitted  $INST  →  $JOB_ID"
done

echo ""
echo "Monitor all jobs:"
echo "  dx watch ${JOB_IDS[*]}"
echo ""
echo "Download results when complete:"
for JOB_ID in "${JOB_IDS[@]}"; do
    echo "  dx download ${JOB_ID}:results_csv"
done
