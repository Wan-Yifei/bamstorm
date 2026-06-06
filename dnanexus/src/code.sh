#!/bin/bash
# DNAnexus entry point for the bamstorm benchmark applet.
# Inputs (set by dx-toolkit after dx-download-all-inputs):
#   $bam_file_path   - local path to downloaded BAM file
#   $bai_file_path   - local path to downloaded BAI file
#   $docker_image    - Docker image reference (string)
#   $threads         - comma-separated thread counts, e.g. "2,4,8,16"
#   $repeats         - number of timed repetitions (int)
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

    # ── download inputs directly to local NVMe ───────────────────────────────
    # Skip dx-download-all-inputs (which lands on EBS root) and write straight
    # to /mnt/work so the BAM never touches the slower network-backed root disk.
    mkdir -p /mnt/work/data /mnt/work/results

    echo "Downloading BAM..."
    dx download "$bam_file" -o /mnt/work/data/input.bam --no-progress
    dx download "$bai_file" -o /mnt/work/data/input.bam.bai --no-progress
    echo "Download complete: $(du -h /mnt/work/data/input.bam | cut -f1)"

    # ── generate bench.toml from job parameters ───────────────────────────────
    # Convert comma-separated string "2,4,8" → TOML array [2, 4, 8]
    threads_toml=$(echo "$threads" | sed 's/,/, /g')

    cat > /mnt/work/bench.toml << TOML
[benchmark]
threads = [$threads_toml]
repeats = $repeats
drop_cache = true

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

    # ── run benchmark (repeat loop here; drop_cache on host before each run) ──
    echo "Starting benchmark..."
    for i in $(seq 1 "$repeats"); do
        echo ""
        echo "=== Repeat ${i}/${repeats}: dropping page cache ==="
        sync && echo 3 > /proc/sys/vm/drop_caches

        if [ "$i" -eq 1 ]; then
            append_arg=""
        else
            append_arg="--append"
        fi

        docker run --rm \
            -v /mnt/work/data:/data \
            -v /mnt/work/results:/results \
            -v /mnt/work/bench.toml:/app/bench.toml:ro \
            "$docker_image" \
            python3 /app/bench.py \
                /data/input.bam \
                /data/input.bam.bai \
                --config /app/bench.toml \
                --no-drop-cache \
                --repeats 1 \
                --repeat-index "$i" \
                --csv /results/benchmark.csv \
                $append_arg \
            2>&1 | tee -a /mnt/work/results/bench.log
    done

    # ── upload outputs ────────────────────────────────────────────────────────
    mkdir -p ~/out/results_csv ~/out/bench_log
    cp /mnt/work/results/benchmark.csv ~/out/results_csv/benchmark.csv
    cp /mnt/work/results/bench.log     ~/out/bench_log/bench.log

    dx-upload-all-outputs --parallel
    echo "Done."
}
