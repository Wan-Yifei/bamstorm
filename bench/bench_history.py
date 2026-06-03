#!/usr/bin/env python3
"""
Cross-version bamstrom benchmark.

Builds bench_count for each git tag via temporary worktrees, runs the same
thread-scaling benchmark for every version, and collects results alongside
reference tools (samtools, rabbitbam, pysam) into a single CSV.

Usage:
    python bench_history.py <bam> <bai> --tags v0.1.0 v0.2.0 HEAD [options]

The special tag "HEAD" skips the build step and reuses the binary already at
BENCH_COUNT_BIN (i.e. whatever is currently deployed in the container).
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib          # type: ignore[no-redef]
    except ImportError:
        tomllib = None                   # type: ignore[assignment]

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

BENCH_COUNT_BIN = "/app/bench_count"
RABBITBAM_BIN   = "/opt/RabbitBAM/rabbitbam"
DEFAULT_CONFIG  = Path(__file__).parent / "bench.toml"
REPO_ROOT       = Path(__file__).parent.parent


# ── config ────────────────────────────────────────────────────────────────────

def load_config(path: Path) -> dict:
    if tomllib is None:
        return {}
    try:
        with open(path, "rb") as fh:
            return tomllib.load(fh)
    except FileNotFoundError:
        return {}


def resolve_threads(thread_list: list[int], max_cpus: int) -> list[int]:
    seen, result = set(), []
    for t in thread_list:
        t = max_cpus if t == 0 else t
        if t not in seen:
            seen.add(t)
            result.append(t)
    return result


# ── build ─────────────────────────────────────────────────────────────────────

def build_tag(tag: str, build_dir: Path) -> Path:
    """
    Check out <tag> into a temporary git worktree, build bench_count --release,
    copy the binary to <build_dir>/bench_count_<tag>, remove the worktree.
    Returns the path to the copied binary.
    """
    safe = tag.replace("/", "_").replace(".", "_")
    worktree = build_dir / f"wt_{safe}"
    binary   = build_dir / f"bench_count_{safe}"

    print(f"  [build] {tag} → {binary}")
    subprocess.run(
        ["git", "worktree", "add", "--detach", str(worktree), tag],
        cwd=REPO_ROOT, check=True, capture_output=True,
    )
    try:
        subprocess.run(
            ["cargo", "build", "--release", "--bin", "bench_count"],
            cwd=worktree, check=True,
        )
        shutil.copy(worktree / "target" / "release" / "bench_count", binary)
        binary.chmod(0o755)
    finally:
        subprocess.run(
            ["git", "worktree", "remove", "--force", str(worktree)],
            cwd=REPO_ROOT, check=True, capture_output=True,
        )
    return binary


# ── runners ───────────────────────────────────────────────────────────────────

def file_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


def drop_caches() -> None:
    try:
        with open("/proc/sys/vm/drop_caches", "w") as fh:
            fh.write("1\n")
    except OSError:
        pass


def timed_best(fn, repeats: int, drop_cache: bool, *fn_args) -> tuple[float, int]:
    best, count = float("inf"), 0
    for _ in range(repeats):
        if drop_cache:
            drop_caches()
        elapsed, count = fn(*fn_args)
        if elapsed < best:
            best = elapsed
    return best, count


def run_bench_count(binary: str, bam: str, bai: str, threads: int) -> tuple[float, int]:
    cmd = [binary, "--threads", str(threads), bam, bai]
    t0 = time.perf_counter()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"{binary} failed (exit {exc.returncode}):\n{exc.stderr.strip()}"
        ) from exc
    return time.perf_counter() - t0, int(result.stdout.strip())


def run_samtools(bam: str, threads: int) -> tuple[float, int]:
    cmd = ["samtools", "view", "-c", "--threads", str(threads), bam]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return time.perf_counter() - t0, int(result.stdout.strip())


def run_rabbitbam(bam: str, threads: int) -> tuple[float, int]:
    cmd = [RABBITBAM_BIN, "benchmark_count", "-i", bam, "-w", str(threads)]
    t0 = time.perf_counter()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"rabbitbam failed (exit {exc.returncode}):\n{exc.stderr.strip()}"
        ) from exc
    elapsed = time.perf_counter() - t0
    for line in result.stdout.splitlines():
        if line.startswith("Bam number is"):
            return elapsed, int(line.split()[-1])
    raise RuntimeError(
        f"rabbitbam output missing 'Bam number is' line:\n{result.stdout.strip()}"
    )


def run_pysam(bam: str, threads: int) -> tuple[float, int]:
    t0 = time.perf_counter()
    with pysam.AlignmentFile(bam, "rb", check_sq=False, threads=threads) as f:
        n = sum(1 for _ in f.fetch(until_eof=True))
    return time.perf_counter() - t0, n


# ── result helpers ────────────────────────────────────────────────────────────

def ok(tag: str, tool: str, threads: int, elapsed: float,
       count: int, bam_mb: float) -> dict:
    return {
        "tag":             tag,
        "tool":            tool,
        "threads":         threads,
        "elapsed_s":       round(elapsed, 3),
        "throughput_mb_s": round(bam_mb / elapsed, 1) if elapsed > 0 else float("inf"),
        "records":         count,
        "error":           "",
    }


def err(tag: str, tool: str, threads: int, error: str) -> dict:
    return {"tag": tag, "tool": tool, "threads": threads,
            "elapsed_s": "", "throughput_mb_s": "", "records": "", "error": error}


def fmt(r: dict) -> str:
    if r["error"]:
        return f"  [{r['tag']}] {r['tool']:<28}  threads={r['threads']:<4}  ERROR: {r['error']}"
    return (
        f"  [{r['tag']}] {r['tool']:<28}  threads={r['threads']:<4}  "
        f"{r['elapsed_s']:7.3f}s  {r['throughput_mb_s']:8.1f} MB/s  records={r['records']}"
    )


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Cross-version BAM benchmark")
    parser.add_argument("bam", help="Path to BAM file")
    parser.add_argument("bai", help="Path to BAI index file")
    parser.add_argument(
        "--tags", nargs="+", required=True, metavar="TAG",
        help="Git tags (or 'HEAD') to benchmark. HEAD reuses the deployed binary.",
    )
    parser.add_argument(
        "--config", default=str(DEFAULT_CONFIG),
        help=f"Path to TOML config (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--no-drop-cache", action="store_true", default=None)
    parser.add_argument(
        "--no-reference", action="store_true",
        help="Skip samtools / rabbitbam / pysam reference runs",
    )
    parser.add_argument(
        "--csv", metavar="FILE",
        help="Write results to CSV (use '-' for stdout)",
    )
    args = parser.parse_args()

    cfg       = load_config(Path(args.config))
    max_cpus  = os.cpu_count() or 1
    bam_mb    = file_mb(args.bam)

    default_threads = cfg.get("benchmark", {}).get("threads", [1, 2, 4, 8, 0])
    threads   = resolve_threads(default_threads, max_cpus)
    repeats   = args.repeats if args.repeats is not None \
                else cfg.get("benchmark", {}).get("repeats", 3)
    drop_cache = (not args.no_drop_cache) if args.no_drop_cache is not None \
                 else cfg.get("benchmark", {}).get("drop_cache", True)

    print(f"\nBAM file : {args.bam}  ({bam_mb:.1f} MB)")
    print(f"CPU cores: {max_cpus}")
    print(f"Tags     : {args.tags}")
    print(f"Threads  : {threads}")
    print(f"Repeats  : {repeats}")
    print()

    results: list[dict] = []

    # ── per-tag bamstrom benchmark ────────────────────────────────────────────
    with tempfile.TemporaryDirectory(prefix="bench_history_") as tmp:
        build_dir = Path(tmp)

        for tag in args.tags:
            print(f"\n{'─' * 60}")
            print(f"  Tag: {tag}")
            print(f"{'─' * 60}")

            if tag.upper() == "HEAD":
                binary = BENCH_COUNT_BIN
                print(f"  [build] HEAD → reusing {binary}")
            else:
                try:
                    binary = str(build_tag(tag, build_dir))
                except subprocess.CalledProcessError as e:
                    print(f"  [build] FAILED: {e}")
                    for t in threads:
                        results.append(err(tag, "bamstrom", t, f"build failed: {e}"))
                    continue

            for t in threads:
                try:
                    elapsed, count = timed_best(
                        run_bench_count, repeats, drop_cache, binary, args.bam, args.bai, t
                    )
                    r = ok(tag, "bamstrom", t, elapsed, count, bam_mb)
                except Exception as e:
                    r = err(tag, "bamstrom", t, str(e))
                print(fmt(r))
                results.append(r)

    # ── reference tools (once, tagged as "reference") ────────────────────────
    if not args.no_reference:
        ref_tag = "reference"

        print(f"\n{'─' * 60}")
        print(f"  Reference tools")
        print(f"{'─' * 60}")

        # samtools
        for t in threads:
            try:
                elapsed, count = timed_best(run_samtools, repeats, drop_cache, args.bam, t)
                r = ok(ref_tag, "samtools view -c", t, elapsed, count, bam_mb)
            except Exception as e:
                r = err(ref_tag, "samtools view -c", t, str(e))
            print(fmt(r))
            results.append(r)

        # rabbitbam
        for t in threads:
            try:
                elapsed, count = timed_best(run_rabbitbam, repeats, drop_cache, args.bam, t)
                r = ok(ref_tag, "rabbitbam benchmark_count", t, elapsed, count, bam_mb)
            except Exception as e:
                r = err(ref_tag, "rabbitbam benchmark_count", t, str(e))
            print(fmt(r))
            results.append(r)

        # pysam
        if HAS_PYSAM:
            for t in threads:
                try:
                    elapsed, count = timed_best(run_pysam, repeats, drop_cache, args.bam, t)
                    r = ok(ref_tag, "pysam fetch(until_eof)", t, elapsed, count, bam_mb)
                except Exception as e:
                    r = err(ref_tag, "pysam fetch(until_eof)", t, str(e))
                print(fmt(r))
                results.append(r)
        else:
            print("  pysam not installed — skipped")

    # ── CSV output ────────────────────────────────────────────────────────────
    if args.csv:
        fields = ["tag", "tool", "threads", "elapsed_s", "throughput_mb_s", "records", "error"]
        fh = sys.stdout if args.csv == "-" else open(args.csv, "w", newline="")
        try:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(results)
        finally:
            if fh is not sys.stdout:
                fh.close()
        if args.csv != "-":
            print(f"\nCSV written to {args.csv}")


if __name__ == "__main__":
    main()
