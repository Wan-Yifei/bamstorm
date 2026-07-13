#!/usr/bin/env python3
"""
Generate a pileup benchmark visualization from coverage.csv.

Usage:
    python bench/plot_coverage.py bench/result/coverage/coverage.csv
    python bench/plot_coverage.py bench/result/coverage/coverage.csv --out docs/benchmark_pileup_chr4.png
"""

import argparse
import csv
import statistics
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.patches import Patch
except Exception as e:
    sys.exit(f"matplotlib import failed: {e}")

# ── palette (colour-blind friendly, matches plot_report.py) ──────────────────
PALETTE = {
    "bamstorm base_pileup": "#E05C2A",   # main bamstorm orange
    "bamstorm coverage":    "#F5A97B",   # lighter orange
    "pysam pileup":         "#8172B2",   # pysam purple
}
HATCH = {
    "bamstorm base_pileup": "",
    "bamstorm coverage":    "",
    "pysam pileup":         "",
}
TOOL_LABEL = {
    "bamstorm base_pileup": "bamstorm\nbase_pileup()\n(parallel, 16 threads)",
    "bamstorm coverage":    "bamstorm\ncoverage()\n(single thread)",
    "pysam pileup":         "pysam\ncol.pileups\n(single thread)",
}
TOOL_ORDER = ["pysam pileup", "bamstorm coverage", "bamstorm base_pileup"]


# ── data loading ─────────────────────────────────────────────────────────────

def load(path: Path):
    cold, warm = {}, {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("timed_out", "0") == "1":
                continue
            tool    = row["tool"]
            elapsed = float(row["elapsed_s"])
            cache   = row.get("cache", "cold").strip()
            target  = warm if cache == "warm" else cold
            target.setdefault(tool, []).append(elapsed)
    return cold, warm


def best(vals):
    return min(vals)


def fmt_time(s):
    if s < 60:
        return f"{s:.1f} s"
    return f"{s/60:.1f} min"


# ── figure ───────────────────────────────────────────────────────────────────

def make_figure(cold, warm, out_path, region="chr4"):
    tools       = [t for t in TOOL_ORDER if t in cold]
    cold_best   = [best(cold[t]) for t in tools]
    warm_best   = [best(warm[t]) for t in tools if t in warm]
    has_warm    = len(warm_best) == len(tools)

    pysam_best  = best(cold["pysam pileup"])
    speedups    = [pysam_best / best(cold[t]) for t in tools]

    labels      = [TOOL_LABEL[t] for t in tools]
    colors      = [PALETTE[t] for t in tools]

    fig, (ax_time, ax_speedup) = plt.subplots(
        1, 2, figsize=(13, 4.2),
        facecolor="white",
        gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.38},
    )
    fig.suptitle(
        f"Pileup Benchmark — {region} full chromosome  "
        f"(AWS i4i.4xlarge, 16 vCPU, local NVMe)",
        fontsize=11, fontweight="bold", y=1.01,
    )

    # ── left: horizontal bar — wall-clock time (log scale) ───────────────────
    y_pos  = list(range(len(tools)))
    bar_h  = 0.42

    if has_warm:
        y_cold = [y + bar_h / 2 for y in y_pos]
        y_warm = [y - bar_h / 2 for y in y_pos]
        bars_cold = ax_time.barh(y_cold, cold_best, height=bar_h,
                                 color=colors, zorder=3, edgecolor="white")
        bars_warm = ax_time.barh(y_warm, warm_best, height=bar_h,
                                 color=colors, alpha=0.45, hatch="///",
                                 zorder=3, edgecolor="white")
        ax_time.legend(
            handles=[Patch(facecolor="#888", alpha=1.0,  label="cold cache (best of 3)"),
                     Patch(facecolor="#888", alpha=0.45, hatch="///", label="warm cache")],
            fontsize=8, loc="lower right", framealpha=0.85,
        )
        label_y = y_cold
    else:
        bars_cold = ax_time.barh(y_pos, cold_best, height=bar_h,
                                 color=colors, zorder=3, edgecolor="white")
        label_y = y_pos

    # label each cold bar with the formatted time
    x_max = max(cold_best) * (1.5 if has_warm else 1.45)
    for y, val in zip(label_y, cold_best):
        ax_time.text(val * 1.06, y, fmt_time(val),
                     va="center", ha="left", fontsize=8.5, fontweight="bold",
                     color="#333333")

    ax_time.set_xscale("log")
    ax_time.set_xlim(cold_best[-1] * 0.5, x_max)
    ax_time.set_yticks(y_pos)
    ax_time.set_yticklabels(labels, fontsize=9)
    ax_time.set_xlabel("Wall-clock time  (log scale)", fontsize=9)
    ax_time.set_title("Cold-cache elapsed time  — lower is better",
                      fontsize=10, fontweight="bold", pad=6)
    ax_time.spines["top"].set_visible(False)
    ax_time.spines["right"].set_visible(False)
    ax_time.grid(axis="x", color="#e0e0e0", linewidth=0.7, zorder=0)
    ax_time.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: fmt_time(v))
    )

    # ── right: speedup bars ───────────────────────────────────────────────────
    bars_spd = ax_speedup.bar(labels, speedups, color=colors, zorder=3,
                               width=0.55, edgecolor="white")
    for bar, spd in zip(bars_spd, speedups):
        ax_speedup.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(speedups) * 0.02,
            f"{spd:.1f}×",
            ha="center", va="bottom", fontsize=9, fontweight="bold",
            color="#333333",
        )
    ax_speedup.axhline(1.0, color="#aaaaaa", linewidth=1.0, linestyle="--", zorder=2)
    ax_speedup.set_ylim(0, max(speedups) * 1.22)
    ax_speedup.set_ylabel("Speedup  (×  vs pysam col.pileups)", fontsize=9)
    ax_speedup.set_title("Speedup vs pysam  — higher is better",
                         fontsize=10, fontweight="bold", pad=6)
    ax_speedup.tick_params(axis="x", labelsize=8)
    ax_speedup.spines["top"].set_visible(False)
    ax_speedup.spines["right"].set_visible(False)
    ax_speedup.grid(axis="y", color="#e0e0e0", linewidth=0.7, zorder=0)

    # ── footer annotation ─────────────────────────────────────────────────────
    fig.text(
        0.5, -0.04,
        f"Dataset: {region} full chromosome (~190 Mbp, 57.6M reads, ~46× coverage)  ·  "
        "bamstorm base_pileup: 16-thread parallel Rust (rayon)  ·  "
        "pysam col.pileups: single-thread Python (GIL-bound)",
        ha="center", fontsize=7.5, color="#555555",
    )

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("csv", help="coverage.csv from bench_coverage.py")
    parser.add_argument("--out", default=None)
    parser.add_argument("--region", default="chr4")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    out_path = Path(args.out) if args.out else Path("docs/benchmark_pileup_chr4.png")

    cold, warm = load(csv_path)
    if not cold:
        sys.exit(f"No valid rows in {csv_path}")

    make_figure(cold, warm, out_path, region=args.region)


if __name__ == "__main__":
    main()
