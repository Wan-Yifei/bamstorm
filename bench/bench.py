#!/usr/bin/env python3
"""
BAM reader benchmark: bamstorm vs samtools vs rabbitbam vs pysam.

Metrics per run:
  - Wall-clock elapsed time (seconds)
  - Disk IO throughput: BAM file size / elapsed  (MB/s)
  - Record count (sanity check)

Thread scaling is configured in bench.toml (same directory as this script).
CLI flags override config values for one-off runs.
"""

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
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

BENCH_COUNT_BIN  = "/app/bench_count"
RABBITBAM_BIN    = "/opt/RabbitBAM/rabbitbam"
DEFAULT_CONFIG   = Path(__file__).parent / "bench.toml"

# ── config loading ────────────────────────────────────────────────────────────

def load_config(path: Path) -> dict:
    if tomllib is None:
        print("[warn] tomllib not available — using built-in defaults")
        return {}
    try:
        with open(path, "rb") as fh:
            return tomllib.load(fh)
    except FileNotFoundError:
        print(f"[warn] Config not found: {path} — using built-in defaults")
        return {}


def resolve_threads(thread_list: list[int], max_cpus: int) -> list[int]:
    """Replace 0 with max_cpus and deduplicate, preserving order."""
    seen, result = set(), []
    for t in thread_list:
        t = max_cpus if t == 0 else t
        if t not in seen:
            seen.add(t)
            result.append(t)
    return result


def tool_threads(cfg: dict, tool: str, default: list[int], max_cpus: int) -> list[int]:
    """Return per-tool thread list, falling back to default."""
    override = cfg.get(tool, {}).get("threads")
    return resolve_threads(override if override is not None else default, max_cpus)


# ── helpers ───────────────────────────────────────────────────────────────────

def file_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


def fmt_row(r: dict, bam_mb: float) -> str:
    tag = "  [warm]" if r.get("cache") == "warm" else ""
    if "error" in r:
        return (
            f"  {r['tool']:<28}  threads={str(r['threads']):<4}  ERROR: {r['error']}{tag}"
        )
    throughput = bam_mb / r["elapsed"] if r["elapsed"] > 0 else float("inf")
    return (
        f"  {r['tool']:<28}  threads={str(r['threads']):<4}  "
        f"{r['elapsed']:7.3f}s  {throughput:8.1f} MB/s  records={r['records']}{tag}"
    )


# ── fio helpers ───────────────────────────────────────────────────────────────

def detect_fs(path: str) -> str:
    try:
        out = subprocess.run(
            ["df", "--output=fstype,target", path],
            capture_output=True, text=True, check=True,
        )
        parts = out.stdout.strip().splitlines()
        if len(parts) >= 2:
            fstype, mount = parts[1].split(None, 1)
            return f"{fstype} on {mount}"
    except Exception:
        pass
    return "unknown"


def run_fio(tmpfile: str, numjobs: int, size: str, runtime: int):
    """Return aggregate sequential read bandwidth in MB/s, or None on failure."""
    cmd = [
        "fio", "--name=bamstorm-io",
        "--rw=read", "--bs=1M",
        f"--size={size}", f"--numjobs={numjobs}",
        f"--runtime={runtime}", "--time_based",
        "--direct=1", "--group_reporting",
        f"--filename={tmpfile}",
        "--output-format=json",
    ]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)
        bw_kb = json.loads(out.stdout)["jobs"][0]["read"]["bw"]
        return bw_kb / 1024
    except Exception:
        return None


# ── runners ───────────────────────────────────────────────────────────────────

def run_bamstorm(bam: str, bai: str, threads: int) -> tuple[float, int]:
    cmd = [BENCH_COUNT_BIN, "--threads", str(threads), bam, bai]
    t0 = time.perf_counter()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f"bench_count failed (exit {exc.returncode}):\n{exc.stderr.strip()}"
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


def drop_caches(bam: str, bai: str) -> None:
    import ctypes, ctypes.util
    POSIX_FADV_DONTNEED = 4
    try:
        _libc = ctypes.CDLL(ctypes.util.find_library("c") or "libc.so.6", use_errno=True)
    except OSError:
        _libc = None

    for path in (bam, bai):
        tmp = path + ".tmp_nocache"
        try:
            shutil.copy2(path, tmp)
        except OSError as e:
            print(f"[error] drop_caches: copy failed for {path}: {e} — aborting benchmark", flush=True)
            sys.exit(0)
        # fsync makes dirty write-pages clean so posix_fadvise DONTNEED can evict them.
        if _libc is not None:
            try:
                wfd = os.open(tmp, os.O_RDWR)
                try:
                    os.fsync(wfd)
                    size = os.fstat(wfd).st_size
                    _libc.posix_fadvise(wfd, 0, size, POSIX_FADV_DONTNEED)
                finally:
                    os.close(wfd)
            except OSError as e:
                print(f"[warn] drop_caches: fadvise failed for {tmp}: {e}", flush=True)
        try:
            os.replace(tmp, path)
        except OSError as e:
            print(f"[error] drop_caches: replace failed for {path}: {e} — aborting benchmark", flush=True)
            sys.exit(0)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="BAM benchmark")
    parser.add_argument("bam", help="Path to BAM file (used by bamstorm; fallback for all tools)")
    parser.add_argument("bai", help="Path to BAI index file (paired with bam)")
    parser.add_argument("--bam2", default=None, metavar="FILE", help="BAM for samtools (default: bam)")
    parser.add_argument("--bai2", default=None, metavar="FILE", help="BAI for samtools")
    parser.add_argument("--bam3", default=None, metavar="FILE", help="BAM for rabbitbam (default: bam)")
    parser.add_argument("--bai3", default=None, metavar="FILE", help="BAI for rabbitbam")
    parser.add_argument("--bam4", default=None, metavar="FILE", help="BAM for pysam (default: bam)")
    parser.add_argument("--bai4", default=None, metavar="FILE", help="BAI for pysam")
    parser.add_argument(
        "--config", default=str(DEFAULT_CONFIG),
        help=f"Path to TOML config file (default: {DEFAULT_CONFIG})",
    )
    parser.add_argument(
        "--repeats", type=int, default=None,
        help="Override config: number of timed repetitions per configuration",
    )
    parser.add_argument(
        "--no-drop-cache", action="store_true", default=None,
        help="Override config: skip dropping OS page cache between runs",
    )
    parser.add_argument(
        "--csv", metavar="FILE",
        help="Write results to CSV file (use '-' for stdout)",
    )
    parser.add_argument(
        "--repeat-index", type=int, default=1, metavar="N",
        help="Repeat index label written to the CSV repeat column (default: 1)",
    )
    parser.add_argument(
        "--append", action="store_true", default=False,
        help="Append to existing CSV without writing header",
    )
    parser.add_argument(
        "--warm-repeats", type=int, default=None, metavar="N",
        help="Override config: warm-cache repetitions after cold runs (0 = skip)",
    )
    args = parser.parse_args()

    cfg = load_config(Path(args.config))
    max_cpus = os.cpu_count() or 1

    default_threads = cfg.get("benchmark", {}).get("threads", [1, 2, 4, 8, 0])
    bamstorm_threads  = tool_threads(cfg, "bamstorm",  default_threads, max_cpus)
    samtools_threads  = tool_threads(cfg, "samtools",  default_threads, max_cpus)
    rabbitbam_threads = tool_threads(cfg, "rabbitbam", default_threads, max_cpus)
    pysam_threads     = tool_threads(cfg, "pysam",     default_threads, max_cpus)

    repeats      = args.repeats if args.repeats is not None \
                   else cfg.get("benchmark", {}).get("repeats", 3)
    drop_cache   = (not args.no_drop_cache) if args.no_drop_cache is not None \
                   else cfg.get("benchmark", {}).get("drop_cache", True)
    warm_repeats = args.warm_repeats if args.warm_repeats is not None \
                   else cfg.get("benchmark", {}).get("warm_repeats", 0)

    # Per-tool BAM/BAI: fall back to the primary bam/bai if not specified.
    # Tool assignment: bamstorm=bam1, samtools=bam2, rabbitbam=bam3, pysam=bam4
    bam1, bai1 = args.bam, args.bai
    bam2, bai2 = args.bam2 or args.bam, args.bai2 or args.bai
    bam3, bai3 = args.bam3 or args.bam, args.bai3 or args.bai
    bam4, bai4 = args.bam4 or args.bam, args.bai4 or args.bai

    bam_mb = file_mb(bam1)

    print(f"\nConfig   : {args.config}")
    print(f"BAM 1 (bamstorm)  : {bam1}  ({bam_mb:.1f} MB)")
    print(f"BAM 2 (samtools)  : {bam2}")
    print(f"BAM 3 (rabbitbam) : {bam3}")
    print(f"BAM 4 (pysam)     : {bam4}")
    print(f"CPU cores: {max_cpus}")
    print(f"Repeats  : cold={repeats} (best of N), warm={warm_repeats}")
    print(f"threads  : {default_threads} (bamstorm={bamstorm_threads} samtools={samtools_threads} "
          f"rabbitbam={rabbitbam_threads} pysam={pysam_threads})")
    print()
    print(f"  {'Tool':<28}  {'':10}  {'elapsed':>9}  {'throughput':>10}  records")
    print("  " + "-" * 75)

    results: list[dict] = []

    # ── fio disk bandwidth ────────────────────────────────────────────────────
    fio_cfg     = cfg.get("fio", {})
    fio_enabled = fio_cfg.get("enabled", True)

    if fio_enabled:
        fio_size    = fio_cfg.get("size",             "1g")
        fio_runtime = fio_cfg.get("runtime",          20)
        fio_par_n   = fio_cfg.get("numjobs_parallel", 0)
        fio_par_n   = max_cpus if fio_par_n == 0 else fio_par_n
        bam_dir     = os.path.dirname(os.path.abspath(args.bam))
        fio_dir     = fio_cfg.get("tmpdir", bam_dir)

        print("  [fio disk bandwidth]")
        if shutil.which("fio"):
            print(f"  filesystem : {detect_fs(args.bam)}")
            tmpfile = os.path.join(fio_dir, ".bamstorm_fio.tmp")
            seq_bw = par_bw = None
            try:
                seq_bw = run_fio(tmpfile, 1,         fio_size, fio_runtime)
                par_bw = run_fio(tmpfile, fio_par_n, fio_size, fio_runtime)
            finally:
                try:
                    os.unlink(tmpfile)
                except OSError:
                    pass
            if seq_bw is not None:
                print(f"  sequential (1 job)          : {seq_bw:8.1f} MB/s")
                results.append({"tool": "fio-seq", "threads": 1,
                                 "elapsed": fio_runtime, "throughput": seq_bw, "records": ""})
            if par_bw is not None:
                print(f"  parallel   ({fio_par_n} jobs) : {par_bw:8.1f} MB/s")
                results.append({"tool": "fio-par", "threads": fio_par_n,
                                 "elapsed": fio_runtime, "throughput": par_bw, "records": ""})
            if seq_bw is None and par_bw is None:
                print("  fio failed — check permissions or available disk space")
        else:
            print("  fio not installed — skipped")
        print()

    def run_repeats(bam, bai, fn, *fn_args) -> list[tuple[float, int]]:
        runs = []
        for _ in range(repeats):
            if drop_cache:
                drop_caches(bam, bai)
            runs.append(fn(*fn_args))
        return runs

    def record_ok(tool: str, threads: int, runs: list[tuple[float, int]],
                  cache: str = "cold") -> dict:
        elapsed, count = min(runs, key=lambda x: x[0])
        return {"tool": tool, "threads": threads, "cache": cache, "elapsed": elapsed,
                "throughput": bam_mb / elapsed if elapsed > 0 else float("inf"),
                "records": count,
                "all_elapsed": [e for e, _ in runs]}

    def record_err(tool: str, threads: int, error: str, cache: str = "cold") -> dict:
        return {"tool": tool, "threads": threads, "cache": cache, "error": error}

    # bamstorm → bam1
    print("  [bamstorm]")
    for t in bamstorm_threads:
        try:
            runs = run_repeats(bam1, bai1, run_bamstorm, bam1, bai1, t)
            r = record_ok("bamstorm", t, runs)
        except Exception as e:
            r = record_err("bamstorm", t, str(e))
        print(fmt_row(r, bam_mb))
        results.append(r)

    # samtools → bam2
    print()
    print("  [samtools]")
    for t in samtools_threads:
        try:
            runs = run_repeats(bam2, bai2, run_samtools, bam2, t)
            r = record_ok("samtools view -c", t, runs)
        except Exception as e:
            r = record_err("samtools view -c", t, str(e))
        print(fmt_row(r, bam_mb))
        results.append(r)

    # rabbitbam → bam3
    print()
    print("  [rabbitbam]")
    for t in rabbitbam_threads:
        try:
            runs = run_repeats(bam3, bai3, run_rabbitbam, bam3, t)
            r = record_ok("rabbitbam benchmark_count", t, runs)
        except Exception as e:
            r = record_err("rabbitbam benchmark_count", t, str(e))
        print(fmt_row(r, bam_mb))
        results.append(r)

    # pysam → bam4
    print()
    print("  [pysam]")
    if HAS_PYSAM:
        for t in pysam_threads:
            try:
                runs = run_repeats(bam4, bai4, run_pysam, bam4, t)
                r = record_ok("pysam fetch(until_eof)", t, runs)
            except Exception as e:
                r = record_err("pysam fetch(until_eof)", t, str(e))
            print(fmt_row(r, bam_mb))
            results.append(r)
    else:
        print("  pysam not installed — skipped")

    # ── warm-cache runs ────────────────────────────────────────────────────────
    if warm_repeats > 0:
        print()
        print(f"  --- warm cache (no eviction, repeats={warm_repeats}) ---")

        print()
        print("  [bamstorm]  [warm]")
        for t in bamstorm_threads:
            try:
                runs = [run_bamstorm(bam1, bai1, t) for _ in range(warm_repeats)]
                r = record_ok("bamstorm", t, runs, cache="warm")
            except Exception as e:
                r = record_err("bamstorm", t, str(e), cache="warm")
            print(fmt_row(r, bam_mb))
            results.append(r)

        print()
        print("  [samtools]  [warm]")
        for t in samtools_threads:
            try:
                runs = [run_samtools(bam2, t) for _ in range(warm_repeats)]
                r = record_ok("samtools view -c", t, runs, cache="warm")
            except Exception as e:
                r = record_err("samtools view -c", t, str(e), cache="warm")
            print(fmt_row(r, bam_mb))
            results.append(r)

        print()
        print("  [rabbitbam]  [warm]")
        for t in rabbitbam_threads:
            try:
                runs = [run_rabbitbam(bam3, t) for _ in range(warm_repeats)]
                r = record_ok("rabbitbam benchmark_count", t, runs, cache="warm")
            except Exception as e:
                r = record_err("rabbitbam benchmark_count", t, str(e), cache="warm")
            print(fmt_row(r, bam_mb))
            results.append(r)

        if HAS_PYSAM:
            print()
            print("  [pysam]  [warm]")
            for t in pysam_threads:
                try:
                    runs = [run_pysam(bam4, t) for _ in range(warm_repeats)]
                    r = record_ok("pysam fetch(until_eof)", t, runs, cache="warm")
                except Exception as e:
                    r = record_err("pysam fetch(until_eof)", t, str(e), cache="warm")
                print(fmt_row(r, bam_mb))
                results.append(r)

    print()

    # CSV output
    if args.csv:
        mode = "a" if args.append else "w"
        fh = sys.stdout if args.csv == "-" else open(args.csv, mode, newline="")
        try:
            writer = csv.DictWriter(
                fh,
                fieldnames=["tool", "threads", "cache", "repeat", "elapsed_s", "throughput_mb_s", "records", "error"],
                extrasaction="ignore",
            )
            if not args.append:
                writer.writeheader()
            for r in results:
                if "error" in r:
                    writer.writerow({
                        "tool": r["tool"], "threads": r["threads"],
                        "cache": r.get("cache", "cold"), "repeat": "",
                        "elapsed_s": "", "throughput_mb_s": "", "records": "",
                        "error": r["error"],
                    })
                elif "all_elapsed" in r:
                    for i, elapsed in enumerate(r["all_elapsed"], args.repeat_index):
                        writer.writerow({
                            "tool": r["tool"], "threads": r["threads"],
                            "cache": r.get("cache", "cold"), "repeat": i,
                            "elapsed_s": elapsed,
                            "throughput_mb_s": f"{bam_mb / elapsed:.1f}" if elapsed > 0 else "",
                            "records": r["records"],
                            "error": "",
                        })
                else:
                    writer.writerow({
                        "tool": r["tool"], "threads": r["threads"],
                        "cache": r.get("cache", "cold"), "repeat": "",
                        "elapsed_s": r.get("elapsed", ""),
                        "throughput_mb_s": f"{r['throughput']:.1f}" if "throughput" in r else "",
                        "records": r.get("records", ""),
                        "error": "",
                    })
        finally:
            if fh is not sys.stdout:
                fh.close()
        if args.csv != "-":
            print(f"CSV written to {args.csv}")


if __name__ == "__main__":
    main()
