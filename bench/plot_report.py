#!/usr/bin/env python3
"""
Generate a visualization report from a benchmark CSV produced by bench.py.
Handles multi-repeat CSV (repeat column).

Usage:
    python bench/plot_report.py benchmark.csv
    python bench/plot_report.py benchmark.csv --out report.png
"""

import argparse
import csv
import random
import statistics
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    from matplotlib.gridspec import GridSpec
except Exception as _e:
    sys.exit(f"matplotlib import failed: {_e}")

# ── colour palette (colour-blind friendly) ────────────────────────────────────
PALETTE = {
    "bamstorm":                   "#E05C2A",
    "samtools view -c":           "#4C72B0",
    "rabbitbam benchmark_count":  "#55A868",
    "pysam fetch(until_eof)":     "#8172B2",
}
MARKER = {
    "bamstorm":                   "o",
    "samtools view -c":           "s",
    "rabbitbam benchmark_count":  "^",
    "pysam fetch(until_eof)":     "D",
}
DISPLAY = {
    "bamstorm":                   "bamstorm",
    "samtools view -c":           "samtools view -c",
    "rabbitbam benchmark_count":  "rabbitbam",
    "pysam fetch(until_eof)":     "pysam",
}
FIO_LINES = {
    "fio-seq": {"color": "#999999", "ls": "--"},
    "fio-par": {"color": "#555555", "ls": ":"},
}


# ── data loading ──────────────────────────────────────────────────────────────

def load(path: Path):
    """
    Returns:
      data: {tool: {threads: [throughput_mb_s, ...]}}  — list holds all repeats
      fio:  {"fio-seq": bw, "fio-par": bw}
    """
    data: dict[str, dict[int, list[float]]] = {}
    fio:  dict[str, float] = {}

    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("error"):
                continue
            tool = row["tool"]
            if tool.startswith("fio-"):
                fio[tool] = float(row["throughput_mb_s"])
                continue
            if not row.get("threads") or not row.get("throughput_mb_s"):
                continue
            t  = int(row["threads"])
            bw = float(row["throughput_mb_s"])
            data.setdefault(tool, {}).setdefault(t, []).append(bw)
    return data, fio


# ── helpers ───────────────────────────────────────────────────────────────────

def _style(ax, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.7, zorder=0)
    ax.grid(axis="x", color="#e0e0e0", linewidth=0.4, zorder=0)
    ax.legend(fontsize=8, framealpha=0.85)


def _fio_hlines(ax, fio: dict) -> None:
    for key in ("fio-seq", "fio-par"):
        if key not in fio:
            continue
        bw    = fio[key]
        style = FIO_LINES[key]
        ax.axhline(bw, color=style["color"], linestyle=style["ls"],
                   linewidth=1.4, label=f"{key}  {bw:.0f} MB/s", zorder=2)


def _thread_ticks(data: dict) -> list[int]:
    return sorted({t for td in data.values() for t in td})


def _peak_thread(td: dict[int, list[float]]) -> int:
    return max(td, key=lambda t: statistics.mean(td[t]))


# ── charts ────────────────────────────────────────────────────────────────────

def plot_throughput(ax, data: dict, fio: dict) -> None:
    """Mean line + min/max shaded band + individual repeat dots + fio reference lines."""
    for tool, td in data.items():
        xs    = sorted(td)
        means = [statistics.mean(td[t]) for t in xs]
        lo    = [min(td[t])             for t in xs]
        hi    = [max(td[t])             for t in xs]
        color  = PALETTE.get(tool, "#888888")
        marker = MARKER.get(tool, "o")
        label  = DISPLAY.get(tool, tool)
        ax.fill_between(xs, lo, hi, color=color, alpha=0.13, zorder=1)
        # individual repeat dots
        for x, vals in td.items():
            ax.scatter([x] * len(vals), vals, color=color, s=18,
                       alpha=0.45, zorder=2, edgecolors="none")
        ax.plot(xs, means, color=color, marker=marker, linewidth=2,
                markersize=6, label=label, zorder=3)
        # annotate only the peak point
        peak_x = xs[means.index(max(means))]
        peak_y = max(means)
        ax.annotate(f"{peak_y:.0f}", (peak_x, peak_y), textcoords="offset points",
                    xytext=(0, 8), ha="center", fontsize=8,
                    fontweight="bold", color=color)
    _fio_hlines(ax, fio)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(_thread_ticks(data))
    _style(ax, "Throughput vs Thread Count  (mean ± range, dots = repeats)",
           "Threads (log₂)", "Throughput  (MB/s)")


def plot_box(ax, data: dict) -> None:
    """
    Grouped box plot: x-axis = thread count, one box per tool per thread.
    Each box covers the N repeats at that (tool, threads) combination.
    Individual data points shown as jittered scatter on top.
    """
    tools      = list(data.keys())
    all_threads = _thread_ticks(data)
    n_tools    = len(tools)
    group_w    = 0.75
    box_w      = group_w / n_tools * 0.88
    offsets    = [(i - (n_tools - 1) / 2) * group_w / n_tools
                  for i in range(n_tools)]

    from matplotlib.patches import Patch
    random.seed(42)
    for j, tool in enumerate(tools):
        td     = data[tool]
        color  = PALETTE.get(tool, "#888888")
        positions = []
        box_data  = []
        for i, t in enumerate(all_threads):
            if t not in td:
                continue
            positions.append(i + offsets[j])
            box_data.append(td[t])

        bp = ax.boxplot(
            box_data,
            positions=positions,
            widths=box_w,
            patch_artist=True,
            notch=False,
            medianprops=dict(color="white", linewidth=1.8),
            whiskerprops=dict(linewidth=1.0, color=color),
            capprops=dict(linewidth=1.0, color=color),
            flierprops=dict(marker="", markersize=0),
            manage_ticks=False,
            zorder=3,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor(color)
            patch.set_alpha(0.60)

        for pos, vals in zip(positions, box_data):
            jitter = [pos + random.uniform(-box_w * 0.28, box_w * 0.28)
                      for _ in vals]
            ax.scatter(jitter, vals, color=color, s=20, alpha=0.85,
                       zorder=5, edgecolors="none")

    ax.set_xticks(range(len(all_threads)))
    ax.set_xticklabels([str(t) for t in all_threads], fontsize=9)
    ax.set_xlabel("Threads", fontsize=9)
    ax.set_ylabel("Throughput  (MB/s)", fontsize=9)
    ax.set_title("Throughput Distribution per Thread Count  (N repeats per box)",
                 fontsize=11, fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.7, zorder=0)

    handles = [Patch(facecolor=PALETTE.get(t, "#888"), alpha=0.7,
                     label=DISPLAY.get(t, t)) for t in tools]
    ax.legend(handles=handles, fontsize=8, framealpha=0.85)


def plot_speedup(ax, data: dict) -> None:
    """Speedup relative to single-thread mean, log/log axes."""
    all_threads = _thread_ticks(data)
    base = all_threads[0]
    for tool, td in data.items():
        if base not in td:
            continue
        base_mean = statistics.mean(td[base])
        xs = sorted(td)
        ys = [statistics.mean(td[t]) / base_mean for t in xs]
        color  = PALETTE.get(tool, "#888888")
        marker = MARKER.get(tool, "o")
        label  = DISPLAY.get(tool, tool)
        ax.plot(xs, ys, color=color, marker=marker, linewidth=2,
                markersize=6, label=label, zorder=3)
    ideal_x = [all_threads[0], all_threads[-1]]
    ideal_y = [1.0, all_threads[-1] / all_threads[0]]
    ax.plot(ideal_x, ideal_y, "--", color="#aaaaaa", linewidth=1.2,
            label="ideal", zorder=2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(all_threads)
    _style(ax, "Speedup vs Single Thread  (log/log)",
           "Threads (log₂)", "Speedup  (×)")


def plot_peak_bar(ax, data: dict) -> None:
    """Peak throughput bar (mean at best thread count) with min-max error bars."""
    tools  = list(data.keys())
    labels = [DISPLAY.get(t, t) for t in tools]
    colors = [PALETTE.get(t, "#888888") for t in tools]

    peak_means, err_lo, err_hi = [], [], []
    for td in data.values():
        best_t = _peak_thread(td)
        vals   = td[best_t]
        m      = statistics.mean(vals)
        peak_means.append(m)
        err_lo.append(m - min(vals))
        err_hi.append(max(vals) - m)

    bars = ax.bar(labels, peak_means, color=colors, zorder=3, width=0.55,
                  edgecolor="white", linewidth=0.8)
    ax.errorbar(range(len(labels)), peak_means,
                yerr=[err_lo, err_hi], fmt="none",
                ecolor="#333333", elinewidth=1.6, capsize=5, zorder=4)
    for bar, val in zip(bars, peak_means):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(peak_means) * 0.02,
                f"{val:.0f}", ha="center", va="bottom",
                fontsize=9, fontweight="bold")

    ax.set_ylabel("Peak Throughput  (MB/s)", fontsize=9)
    ax.set_title("Peak Throughput  (mean ± range at best thread count)",
                 fontsize=11, fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.7, zorder=0)
    ax.tick_params(axis="x", labelsize=9)
    ax.set_ylim(0, max(peak_means) * 1.2)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot bamstorm benchmark report")
    parser.add_argument("csv", help="Path to benchmark CSV")
    parser.add_argument("--out", default=None,
                        help="Output image path (default: <csv stem>_report.png)")
    args = parser.parse_args()

    csv_path = Path(args.csv)
    out_path = Path(args.out) if args.out else csv_path.with_name(
        csv_path.stem + "_report.png"
    )

    data, fio = load(csv_path)
    if not data:
        sys.exit(f"No valid rows found in {csv_path}")

    fig = plt.figure(figsize=(14, 10), facecolor="white")
    fig.suptitle(
        f"BAM Reader Benchmark — {csv_path.stem}",
        fontsize=14, fontweight="bold", y=0.98,
    )
    gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.34,
                  left=0.08, right=0.97, top=0.92, bottom=0.09,
                  height_ratios=[1.2, 1])

    plot_throughput(fig.add_subplot(gs[0, :]), data, fio)   # full-width top
    plot_speedup(fig.add_subplot(gs[1, 0]), data)
    plot_peak_bar(fig.add_subplot(gs[1, 1]), data)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Report saved: {out_path}")


if __name__ == "__main__":
    main()
