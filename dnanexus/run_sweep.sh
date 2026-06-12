#!/bin/bash
# Submit one benchmark job per thread count.
#
# Each job downloads 4 BAM copies (same data, different file-IDs) so that
# each tool reads a fresh file with no cross-tool page-cache contamination:
#   bam1 → bamstorm
#   bam2 → samtools
#   bam3 → rabbitbam
#   bam4 → pysam
#
# Usage:
#   ./run_sweep.sh \
#     --bam1 file-xxxx --bai1 file-xxxx \
#     --bam2 file-xxxx --bai2 file-xxxx \
#     --bam3 file-xxxx --bai3 file-xxxx \
#     --bam4 file-xxxx --bai4 file-xxxx \
#     --image docker.io/yourname/bamstorm-bench:latest \
#     [--threads 1,2,4,8,16,32,64,96] \
#     [--repeats 3] \
#     [--instance mem2_ssd1_v2_x96]
set -euo pipefail

BAM1="file-J8Vpyv80PZKGf3zGP2BqkkQ8" BAI1="file-J8Vpyp00PZKJpx9Q3k8pvjq3"
BAM2="file-J8ZjXB00PZK5z997GZ7K4yyb" BAI2="file-J8ZkvG80PZKGXkYkpjB6Q6B7"
BAM3="file-J8Zp2yQ0PZK1Z3qF0fv26YvG" BAI3="file-J8ZpfZQ0PZK6g7y49Kfp6JVg"
BAM4="file-J8b2Vz80PZKPkKpV5QXFJ9qP" BAI4="file-J8Zz1Z80PZKBY4GxfVJ3qQx1"
IMAGE=""
THREADS="1,2,4,8,16,32,64,96"
REPEATS=3
INSTANCE="mem2_ssd1_v2_x96"

usage() {
    cat <<EOF
Usage: $0 --bam1 <id> --bai1 <id> --bam2 <id> --bai2 <id> \
          --bam3 <id> --bai3 <id> --bam4 <id> --bai4 <id> \
          --image <docker-image> [OPTIONS]

Required:
  --image     IMAGE     Docker image reference

Optional (pre-filled with test_15gb_1..4.bam defaults):
  --bam1..4   FILE-ID   DNAnexus file IDs for the 4 BAM copies
  --bai1..4   FILE-ID   DNAnexus file IDs for their indexes

Optional:
  --threads   LIST      Comma-separated thread counts (default: $THREADS)
  --repeats   N         Repetitions per tool              (default: $REPEATS)
  --instance  TYPE      DNAnexus instance type            (default: $INSTANCE)
  -h                    Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bam1)    BAM1="$2";     shift 2 ;;
        --bai1)    BAI1="$2";     shift 2 ;;
        --bam2)    BAM2="$2";     shift 2 ;;
        --bai2)    BAI2="$2";     shift 2 ;;
        --bam3)    BAM3="$2";     shift 2 ;;
        --bai3)    BAI3="$2";     shift 2 ;;
        --bam4)    BAM4="$2";     shift 2 ;;
        --bai4)    BAI4="$2";     shift 2 ;;
        --image)   IMAGE="$2";    shift 2 ;;
        --threads) THREADS="$2";  shift 2 ;;
        --repeats) REPEATS="$2";  shift 2 ;;
        --instance) INSTANCE="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

if [[ -z "$IMAGE" ]]; then
    echo "Error: --image is required."
    usage
    exit 1
fi

IFS=',' read -ra THREAD_LIST <<< "$THREADS"

echo "Submitting ${#THREAD_LIST[@]} jobs (one per thread count)"
echo "Instance : $INSTANCE"
echo "Repeats  : $REPEATS"
echo "BAM 1    : $BAM1 (bamstorm)"
echo "BAM 2    : $BAM2 (samtools)"
echo "BAM 3    : $BAM3 (rabbitbam)"
echo "BAM 4    : $BAM4 (pysam)"
echo ""

JOB_IDS=()

for T in "${THREAD_LIST[@]}"; do
    JOB_ID=$(dx run bamstorm_bench \
        -i bam_file_1="$BAM1" -i bai_file_1="$BAI1" \
        -i bam_file_2="$BAM2" -i bai_file_2="$BAI2" \
        -i bam_file_3="$BAM3" -i bai_file_3="$BAI3" \
        -i bam_file_4="$BAM4" -i bai_file_4="$BAI4" \
        -i docker_image="$IMAGE" \
        -i threads="$T" \
        -i repeats="$REPEATS" \
        --instance-type "$INSTANCE" \
        --name "bamstorm_bench_t${T}" \
        --brief \
        --yes)
    JOB_IDS+=("$JOB_ID")
    echo "Submitted  threads=$T  →  $JOB_ID"
done

echo ""
echo "Monitor all jobs:"
echo "  dx watch ${JOB_IDS[*]}"
echo ""
echo "Download results when complete:"
for JOB_ID in "${JOB_IDS[@]}"; do
    echo "  dx download ${JOB_ID}:results_csv"
done
