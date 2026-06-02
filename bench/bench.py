#!/usr/bin/env python3
"""
BAM reader benchmark: bamstrom vs samtools vs pysam.

Metrics per run:
  - Wall-clock elapsed time (seconds)
  - Disk IO throughput: BAM file size / elapsed  (MB/s)
  - Record count (sanity check)

Thread scaling: bamstrom is run at 1, 2, 4, 8, and max-CPU threads.
"""

import argparse
import os
import subprocess
import sys
import time

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

BENCH_COUNT_BIN = "/app/bench_count"


# ── helpers ──────────────────────────────────────────────────────────────────

def file_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


def fmt_row(label: str, threads: int | str, elapsed: float, mb: float,
            count: int, bam_mb: float) -> str:
    throughput = bam_mb / elapsed if elapsed > 0 else float("inf")
    return (
        f"  {label:<28}  threads={str(threads):<4}  "
        f"{elapsed:7.3f}s  {throughput:8.1f} MB/s  records={count}"
    )


# ── bamstrom ──────────────────────────────────────────────────────────────────

def run_bamstrom(bam: str, bai: str, threads: int) -> tuple[float, int]:
    """Return (elapsed_s, record_count)."""
    cmd = [BENCH_COUNT_BIN, "--threads", str(threads), bam, bai]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    elapsed = time.perf_counter() - t0
    count = int(result.stdout.strip())
    return elapsed, count


# ── samtools ──────────────────────────────────────────────────────────────────

def run_samtools(bam: str, threads: int) -> tuple[float, int]:
    """Return (elapsed_s, record_count). Uses samtools view -c."""
    cmd = ["samtools", "view", "-c", "--threads", str(threads), bam]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    elapsed = time.perf_counter() - t0
    count = int(result.stdout.strip())
    return elapsed, count


# ── pysam ─────────────────────────────────────────────────────────────────────

def run_pysam(bam: str) -> tuple[float, int]:
    """Return (elapsed_s, record_count)."""
    t0 = time.perf_counter()
    n = 0
    with pysam.AlignmentFile(bam, "rb", check_sq=False) as f:
        for _ in f.fetch(until_eof=True):
            n += 1
    elapsed = time.perf_counter() - t0
    return elapsed, n


# ── warm-up ───────────────────────────────────────────────────────────────────

def drop_caches() -> None:
    """Drop OS page cache if running as root (Linux only)."""
    try:
        with open("/proc/sys/vm/drop_caches", "w") as fh:
            fh.write("1\n")
    except (PermissionError, FileNotFoundError):
        pass  # not root or not Linux — skip silently


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="BAM benchmark")
    parser.add_argument("bam", help="Path to BAM file")
    parser.add_argument("bai", help="Path to BAI index file")
    parser.add_argument(
        "--repeats", type=int, default=3,
        help="Number of timed repetitions per configuration (default: 3)",
    )
    parser.add_argument(
        "--no-drop-cache", action="store_true",
        help="Skip dropping OS page cache between runs",
    )
    args = parser.parse_args()

    bam_mb = file_mb(args.bam)
    max_cpus = os.cpu_count() or 1
    thread_counts = sorted(set([1, 2, 4, 8, max_cpus]))

    print(f"\nBAM file : {args.bam}")
    print(f"BAM size : {bam_mb:.1f} MB")
    print(f"CPU cores: {max_cpus}")
    print(f"Repeats  : {args.repeats}")
    print()
    print(f"  {'Tool':<28}  {'':10}  {'elapsed':>9}  {'throughput':>10}  records")
    print("  " + "-" * 75)

    results: list[str] = []

    def timed_best(fn, *fn_args) -> tuple[float, int]:
        """Run fn(*fn_args) args.repeats times; return (min_elapsed, last_count)."""
        best, count = float("inf"), 0
        for _ in range(args.repeats):
            if not args.no_drop_cache:
                drop_caches()
            elapsed, count = fn(*fn_args)
            if elapsed < best:
                best = elapsed
        return best, count

    # bamstrom — thread scaling
    print("  [bamstrom]")
    for t in thread_counts:
        elapsed, count = timed_best(run_bamstrom, args.bam, args.bai, t)
        row = fmt_row("bamstrom", t, elapsed, bam_mb, count, bam_mb)
        print(row)
        results.append(row)

    # samtools — 1 thread and max threads
    print()
    print("  [samtools]")
    for t in [1, max_cpus]:
        try:
            elapsed, count = timed_best(run_samtools, args.bam, t)
            row = fmt_row("samtools view -c", t, elapsed, bam_mb, count, bam_mb)
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            row = f"  {'samtools view -c':<28}  threads={t:<4}  ERROR: {e}"
        print(row)
        results.append(row)

    # pysam — single-threaded only
    print()
    print("  [pysam]")
    if HAS_PYSAM:
        elapsed, count = timed_best(run_pysam, args.bam)
        row = fmt_row("pysam fetch(until_eof)", 1, elapsed, bam_mb, count, bam_mb)
    else:
        row = "  pysam not installed — skipped"
    print(row)
    results.append(row)

    # summary table
    print()
    print("=" * 79)
    print("Summary (best of {} runs each):".format(args.repeats))
    print("=" * 79)
    print(f"  {'Tool':<28}  {'':10}  {'elapsed':>9}  {'throughput':>10}  records")
    print("  " + "-" * 75)
    for row in results:
        print(row)
    print()


if __name__ == "__main__":
    main()
