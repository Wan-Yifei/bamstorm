#!/usr/bin/env python3
"""
Standalone coverage benchmark: bamstorm AlignmentFile.coverage() vs pysam pileup().

Compares per-base coverage via Rust diff-array accumulation (bamstorm) against
Python PileupColumn/PileupRead object creation (pysam).  The key bottleneck for
pysam is GIL-bound Python object allocation: O(positions × depth) PileupRead
objects per contig.  bamstorm uses an O(reads) diff-array so it is unaffected
by sequencing depth.

Usage (Docker):
    docker run --rm --privileged -v /data:/data bamstorm-bench \\
        python3 /app/bench_coverage.py /data/sample.bam /data/sample.bam.bai \\
        --contig chr22 --csv /data/coverage.csv

Usage (local, no sudo):
    python3 bench/bench_coverage.py sample.bam sample.bam.bai \\
        --contig chr22 --no-drop-cache
"""

import argparse
import csv
import os
import signal
import sys
import time

DEFAULT_CONTIG   = "chr22"
DEFAULT_REPEATS  = 3
DEFAULT_TIMEOUT  = 600     # seconds; pysam chr22 can take 5-20 min at 30x

# ── page-cache eviction ───────────────────────────────────────────────────────

_drop_warned = False


def drop_caches(bam: str, bai: str) -> None:
    global _drop_warned
    try:
        with open("/proc/sys/vm/drop_caches", "w") as f:
            f.write("3\n")
    except OSError:
        if not _drop_warned:
            print(
                "[warn] /proc/sys/vm/drop_caches not writable (not root) — "
                "results reflect warm cache",
                flush=True,
            )
            _drop_warned = True


# ── runners ───────────────────────────────────────────────────────────────────

def run_bamstorm_coverage(
    bam: str, bai: str, contig: str,
    start: int | None, stop: int | None,
) -> tuple[float, int]:
    import bamstorm
    t0 = time.perf_counter()
    af = bamstorm.AlignmentFile(bam, "rb", bai_path=bai)
    cov = af.coverage(contig, start, stop)
    total = int(sum(cov))
    return time.perf_counter() - t0, total


class _PileupTimeout(Exception):
    pass


def run_pysam_pileup(
    bam: str, contig: str,
    start: int | None, stop: int | None,
    timeout: int,
) -> tuple[float | None, int | None]:
    """Returns (elapsed, total_coverage_bases), or (None, None) on timeout."""
    import pysam

    def _handler(signum, frame):
        raise _PileupTimeout()

    signal.signal(signal.SIGALRM, _handler)
    signal.alarm(timeout)
    try:
        t0 = time.perf_counter()
        with pysam.AlignmentFile(bam, "rb") as f:
            if start is not None and stop is not None:
                gen = f.pileup(contig, start, stop)
            else:
                gen = f.pileup(contig)
            total = sum(col.nsegments for col in gen)
        elapsed = time.perf_counter() - t0
        return elapsed, total
    except _PileupTimeout:
        return None, None
    finally:
        signal.alarm(0)


# ── formatting ────────────────────────────────────────────────────────────────

def fmt_time(s: float | None) -> str:
    if s is None:
        return "TIMEOUT"
    if s < 60:
        return f"{s:.2f}s"
    m, sec = divmod(s, 60)
    return f"{int(m)}m {sec:.1f}s"


def fmt_speedup(a: float | None, b: float | None) -> str:
    """Return b/a as a speedup string (pysam time / bamstorm time)."""
    if a is None or b is None or a == 0:
        return "N/A"
    return f"{b / a:.1f}×"


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Coverage benchmark: bamstorm diff-array vs pysam pileup",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", help="BAM file path")
    parser.add_argument("bai", help="BAI index path")
    parser.add_argument("--contig",  default=DEFAULT_CONTIG,
                        help="Contig/chromosome name")
    parser.add_argument("--start",   type=int, default=None, metavar="N",
                        help="0-based region start (default: full contig)")
    parser.add_argument("--stop",    type=int, default=None, metavar="N",
                        help="0-based region end   (default: full contig)")
    parser.add_argument("--repeats", type=int, default=DEFAULT_REPEATS,
                        help="Cold-run repetitions (best of N reported)")
    parser.add_argument("--warm-repeats", type=int, default=1, metavar="N",
                        help="Warm-cache repetitions after cold runs (0 = skip)")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT,
                        help="pysam pileup hard timeout in seconds")
    parser.add_argument("--no-drop-cache", action="store_true",
                        help="Skip page-cache eviction between cold runs")
    parser.add_argument("--csv", metavar="FILE",
                        help="Append results to CSV file (use '-' for stdout)")
    args = parser.parse_args()

    region = (
        f"{args.contig}:{args.start}-{args.stop}"
        if args.start is not None and args.stop is not None
        else args.contig
    )
    bam_mb = os.path.getsize(args.bam) / (1024 * 1024)

    print(f"\nCoverage benchmark — {region}")
    print(f"  BAM    : {args.bam}  ({bam_mb:,.0f} MB)")
    print(f"  BAI    : {args.bai}")
    print(f"  cold runs : {args.repeats}  (best time reported)")
    print(f"  warm runs : {args.warm_repeats}")
    print(f"  pysam timeout : {args.timeout}s")
    print()

    try:
        import bamstorm          # noqa: F401
        has_bamstorm = True
    except ImportError:
        has_bamstorm = False
        print("[warn] bamstorm not installed — bamstorm runs skipped")

    try:
        import pysam             # noqa: F401
        has_pysam = True
    except ImportError:
        has_pysam = False
        print("[warn] pysam not installed — pysam runs skipped")

    # Collects all individual (elapsed, total_cov) tuples per tool.
    all_rows: list[dict] = []

    def _record(tool: str, cache: str, repeat: int,
                elapsed: float | None, total: int | None) -> dict:
        return {
            "tool": tool, "region": region,
            "cache": cache, "repeat": repeat,
            "elapsed_s": f"{elapsed:.4f}" if elapsed is not None else "",
            "total_cov_bases": total if total is not None else "",
            "timed_out": "1" if elapsed is None else "0",
        }

    # ── cold runs ─────────────────────────────────────────────────────────────

    bs_cold: list[tuple[float, int]] = []
    if has_bamstorm:
        print(f"  [bamstorm coverage]  cold × {args.repeats}")
        for i in range(1, args.repeats + 1):
            if not args.no_drop_cache:
                drop_caches(args.bam, args.bai)
            elapsed, total = run_bamstorm_coverage(
                args.bam, args.bai, args.contig, args.start, args.stop)
            bs_cold.append((elapsed, total))
            all_rows.append(_record("bamstorm coverage", "cold", i, elapsed, total))
            print(f"    run {i}: {fmt_time(elapsed)}  total_cov={total:,}", flush=True)
        print()

    ps_cold: list[tuple[float | None, int | None]] = []
    ps_timed_out = False
    if has_pysam:
        print(f"  [pysam pileup]  cold × {args.repeats}  (timeout={args.timeout}s each)")
        for i in range(1, args.repeats + 1):
            if not args.no_drop_cache:
                drop_caches(args.bam, args.bai)
            elapsed, total = run_pysam_pileup(
                args.bam, args.contig, args.start, args.stop, args.timeout)
            ps_cold.append((elapsed, total))
            all_rows.append(_record("pysam pileup", "cold", i, elapsed, total))
            if elapsed is None:
                print(f"    run {i}: TIMEOUT (>{args.timeout}s)", flush=True)
                ps_timed_out = True
                break
            print(f"    run {i}: {fmt_time(elapsed)}  total_cov={total:,}", flush=True)
        print()

    # ── warm runs ─────────────────────────────────────────────────────────────

    bs_warm: list[tuple[float, int]] = []
    ps_warm: list[tuple[float | None, int | None]] = []

    if args.warm_repeats > 0:
        print(f"  --- warm cache (no eviction, repeats={args.warm_repeats}) ---")

        if has_bamstorm:
            print(f"\n  [bamstorm coverage]  warm × {args.warm_repeats}")
            for i in range(1, args.warm_repeats + 1):
                elapsed, total = run_bamstorm_coverage(
                    args.bam, args.bai, args.contig, args.start, args.stop)
                bs_warm.append((elapsed, total))
                all_rows.append(_record("bamstorm coverage", "warm", i, elapsed, total))
                print(f"    run {i}: {fmt_time(elapsed)}  total_cov={total:,}", flush=True)

        if has_pysam and not ps_timed_out:
            print(f"\n  [pysam pileup]  warm × {args.warm_repeats}  (timeout={args.timeout}s)")
            for i in range(1, args.warm_repeats + 1):
                elapsed, total = run_pysam_pileup(
                    args.bam, args.contig, args.start, args.stop, args.timeout)
                ps_warm.append((elapsed, total))
                all_rows.append(_record("pysam pileup", "warm", i, elapsed, total))
                if elapsed is None:
                    print(f"    run {i}: TIMEOUT", flush=True)
                    break
                print(f"    run {i}: {fmt_time(elapsed)}  total_cov={total:,}", flush=True)
        print()

    # ── summary ───────────────────────────────────────────────────────────────

    def best_time(runs: list) -> float | None:
        valid = [e for e, _ in runs if e is not None]
        return min(valid) if valid else None

    bs_best = best_time(bs_cold)
    ps_best = best_time(ps_cold)
    bs_warm_best = best_time(bs_warm)
    ps_warm_best = best_time(ps_warm)

    print("  ── Summary " + "─" * 55)
    print(f"  {'Tool':<28}  {'cold best':>12}  {'warm best':>12}  {'cold total_cov':>18}")
    print("  " + "─" * 76)

    def best_total(runs: list) -> int | None:
        valid = [(e, t) for e, t in runs if e is not None]
        if not valid:
            return None
        return min(valid, key=lambda x: x[0])[1]

    if has_bamstorm:
        tot = best_total(bs_cold)
        tot_str = f"{tot:,}" if tot is not None else "N/A"
        warm_str = fmt_time(bs_warm_best) if bs_warm else "—"
        print(f"  {'bamstorm coverage':<28}  {fmt_time(bs_best):>12}  {warm_str:>12}  {tot_str:>18}")
    if has_pysam:
        tot = best_total(ps_cold)
        tot_str = f"{tot:,}" if tot is not None else "N/A"
        warm_str = fmt_time(ps_warm_best) if ps_warm else "—"
        print(f"  {'pysam pileup':<28}  {fmt_time(ps_best):>12}  {warm_str:>12}  {tot_str:>18}")

    if has_bamstorm and has_pysam:
        print()
        print(f"  Cold speedup  (pysam / bamstorm): {fmt_speedup(bs_best, ps_best)}")
        if bs_warm and ps_warm:
            print(f"  Warm speedup  (pysam / bamstorm): {fmt_speedup(bs_warm_best, ps_warm_best)}")
    print()

    # ── CSV ───────────────────────────────────────────────────────────────────

    if args.csv and all_rows:
        fieldnames = list(all_rows[0].keys())
        if args.csv == "-":
            writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_rows)
        else:
            with open(args.csv, "w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(all_rows)
            print(f"CSV written to {args.csv}")


if __name__ == "__main__":
    main()
