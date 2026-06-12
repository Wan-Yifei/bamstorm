#!/bin/bash
# DNAnexus entry point for the bamstorm benchmark applet.
#
# Inputs (set by the platform after dx-download-all-inputs):
#   $bam_file_N / $bai_file_N  - DNAnexus file IDs for BAM N and its index
#   $docker_image              - Docker image reference (string)
#   $threads                   - single thread count for this job (int)
#   $repeats                   - number of timed repetitions per tool (int)
#
# Tool-to-BAM mapping (cold-cache isolation):
#   bam 1 → bamstorm
#   bam 2 → samtools
#   bam 3 → rabbitbam
#   bam 4 → pysam
set -euo pipefail

main() {
    echo "=== Bamstorm Benchmark on DNAnexus ==="
    echo "CPUs    : $(nproc)"
    echo "RAM     : $(free -h | awk '/^Mem/{print $2}')"
    echo "Image   : $docker_image"
    echo "Threads : $threads"
    echo "Repeats : $repeats"
    echo ""

    # ── scratch directory (DNAnexus mounts NVMe RAID at /) ───────────────────
    echo "=== Storage layout ==="
    lsblk -d -o NAME,ROTA,TYPE,SIZE,MOUNTPOINT
    df -h /
    echo ""
    mkdir -p /mnt/work

    # ── ensure fio is available on the worker ─────────────────────────────────
    if ! command -v fio &>/dev/null; then
        DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends fio > /dev/null
    fi

    echo "=== fio pre-flight (sequential read, 10s) ==="
    fio --name=preflight --rw=read --bs=1m --size=512m \
        --filename=/mnt/work/fio_preflight.dat \
        --runtime=10 --time_based --direct=1 \
        --output-format=json > /tmp/fio_preflight.json
    BANDWIDTH_MB=$(python3 -c "
import json
d = json.load(open('/tmp/fio_preflight.json'))
bw = d['jobs'][0]['read']['bw'] / 1024
print(f'{bw:.0f}')
")
    rm -f /mnt/work/fio_preflight.dat
    echo "Disk sequential read: ${BANDWIDTH_MB} MB/s"

    if [ "${BANDWIDTH_MB}" -lt 400 ]; then
        echo "ERROR: disk bandwidth ${BANDWIDTH_MB} MB/s < 400 MB/s threshold"
        echo "       Likely a slow network mount — aborting."
        exit 1
    fi
    echo "Storage OK (${BANDWIDTH_MB} MB/s — local NVMe RAID with dm-crypt overhead expected)"
    echo ""

    # ── download all 4 BAM/BAI pairs directly to local NVMe ─────────────────
    mkdir -p /mnt/work/data /mnt/work/results

    for N in 1 2 3 4; do
        BAM_VAR="bam_file_${N}"
        BAI_VAR="bai_file_${N}"
        echo "Downloading BAM ${N}..."
        dx download "${!BAM_VAR}" -o "/mnt/work/data/input${N}.bam" --no-progress
        dx download "${!BAI_VAR}" -o "/mnt/work/data/input${N}.bam.bai" --no-progress
        echo "  input${N}.bam : $(du -h /mnt/work/data/input${N}.bam | cut -f1)"
    done
    echo "All downloads complete."
    echo ""

    # ── generate bench.toml ───────────────────────────────────────────────────
    cat > /mnt/work/bench.toml << TOML
[benchmark]
threads = [$threads]
repeats = $repeats
drop_cache = false

[fio]
enabled = true
size = "1g"
runtime = 20
numjobs_parallel = 0
tmpdir = "/tmp"
TOML

    # ── pull Docker image ─────────────────────────────────────────────────────
    echo "Pulling image..."
    docker pull "$docker_image"

    # ── run benchmark ─────────────────────────────────────────────────────────
    echo "Starting benchmark (threads=$threads)..."
    docker run --rm \
        -v /mnt/work/data:/data \
        -v /mnt/work/results:/results \
        -v /mnt/work/bench.toml:/app/bench.toml:ro \
        "$docker_image" \
        python3 /app/bench.py \
            /data/input1.bam \
            /data/input1.bam.bai \
            --bam2 /data/input2.bam --bai2 /data/input2.bam.bai \
            --bam3 /data/input3.bam --bai3 /data/input3.bam.bai \
            --bam4 /data/input4.bam --bai4 /data/input4.bam.bai \
            --config /app/bench.toml \
            --csv /results/benchmark.csv \
        2>&1 | tee /mnt/work/results/bench.log

    # ── upload outputs ────────────────────────────────────────────────────────
    mkdir -p ~/out/results_csv ~/out/bench_log
    cp /mnt/work/results/benchmark.csv ~/out/results_csv/benchmark.csv
    cp /mnt/work/results/bench.log     ~/out/bench_log/bench.log

    dx-upload-all-outputs --parallel
    echo "Done."
}
