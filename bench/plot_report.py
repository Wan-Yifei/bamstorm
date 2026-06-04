#!/usr/bin/env python3
"""
Generate a visualization report from a benchmark CSV produced by bench.py or bench_history.py.

Usage:
    python bench/plot_report.py benchmark_v0.3.0.csv
    python bench/plot_report.py benchmark_v0.3.0.csv --out report_v0.3.0.png
"""

import argparse
import csv
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
    "bamstrom":                   "#E05C2A",   # orange-red  (hero)
    "samtools view -c":           "#4C72B0",   # steel blue
    "rabbitbam benchmark_count":  "#55A868",   # green
    "pysam fetch(until_eof)":     "#8172B2",   # purple
}
MARKER = {
    "bamstrom":                   "o",
    "samtools view -c":           "s",
    "rabbitbam benchmark_count":  "^",
    "pysam fetch(until_eof)":     "D",
}
DISPLAY = {
    "bamstrom":                   "bamstrom",
    "samtools view -c":           "samtools view -c",
    "rabbitbam benchmark_count":  "rabbitbam",
    "pysam fetch(until_eof)":     "pysam",
}


# ── data loading ──────────────────────────────────────────────────────────────

def load(path: Path) -> dict[str, dict[int, dict]]:
    """Return {tool: {threads: {elapsed_s, throughput_mb_s}}}."""
    data: dict[str, dict[int, dict]] = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("error"):
                continue
            tool = row["tool"]
            t    = int(row["threads"])
            data.setdefault(tool, {})[t] = {
                "elapsed_s":       float(row["elapsed_s"]),
                "throughput_mb_s": float(row["throughput_mb_s"]),
            }
    return data


# ── plot helpers ──────────────────────────────────────────────────────────────

def _style(ax, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.grid(axis="x", color="#dddddd", linewidth=0.4, zorder=0)
    ax.legend(fontsize=9, framealpha=0.8)


def plot_throughput(ax, data: dict) -> None:
    for tool, rows in data.items():
        xs = sorted(rows)
        ys = [rows[t]["throughput_mb_s"] for t in xs]
        color  = PALETTE.get(tool, "#888888")
        marker = MARKER.get(tool, "o")
        label  = DISPLAY.get(tool, tool)
        ax.plot(xs, ys, color=color, marker=marker, linewidth=2,
                markersize=7, label=label, zorder=3)
        for x, y in zip(xs, ys):
            ax.annotate(f"{y:.0f}", (x, y), textcoords="offset points",
                        xytext=(0, 7), ha="center", fontsize=7.5,
                        color=color)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(sorted({t for rows in data.values() for t in rows}))
    _style(ax, "Throughput vs Thread Count",
           "Threads (log₂ scale)", "Throughput  (MB/s)")


def plot_elapsed(ax, data: dict) -> None:
    for tool, rows in data.items():
        xs = sorted(rows)
        ys = [rows[t]["elapsed_s"] for t in xs]
        color  = PALETTE.get(tool, "#888888")
        marker = MARKER.get(tool, "o")
        label  = DISPLAY.get(tool, tool)
        ax.plot(xs, ys, color=color, marker=marker, linewidth=2,
                markersize=7, label=label, zorder=3)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(sorted({t for rows in data.values() for t in rows}))
    _style(ax, "Elapsed Time vs Thread Count",
           "Threads (log₂ scale)", "Elapsed  (s)")


def plot_speedup(ax, data: dict) -> None:
    for tool, rows in data.items():
        if 1 not in rows:
            continue
        base = rows[1]["elapsed_s"]
        xs = sorted(rows)
        ys = [base / rows[t]["elapsed_s"] for t in xs]
        color  = PALETTE.get(tool, "#888888")
        marker = MARKER.get(tool, "o")
        label  = DISPLAY.get(tool, tool)
        ax.plot(xs, ys, color=color, marker=marker, linewidth=2,
                markersize=7, label=label, zorder=3)
    # ideal line
    all_threads = sorted({t for rows in data.values() for t in rows})
    ideal_x = [all_threads[0], all_threads[-1]]
    ideal_y = [1.0, all_threads[-1] / all_threads[0]]
    ax.plot(ideal_x, ideal_y, "--", color="#aaaaaa", linewidth=1,
            label="ideal", zorder=2)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(all_threads)
    _style(ax, "Speedup vs Single Thread  (log/log)",
           "Threads (log₂)", "Speedup  (×)")


def plot_peak_bar(ax, data: dict) -> None:
    tools  = list(data.keys())
    peaks  = [max(v["throughput_mb_s"] for v in rows.values())
              for rows in data.values()]
    colors = [PALETTE.get(t, "#888888") for t in tools]
    labels = [DISPLAY.get(t, t) for t in tools]

    bars = ax.bar(labels, peaks, color=colors, zorder=3, width=0.55,
                  edgecolor="white", linewidth=0.8)
    for bar, val in zip(bars, peaks):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(peaks) * 0.015,
                f"{val:.0f}", ha="center", va="bottom",
                fontsize=9, fontweight="bold")
    ax.set_ylabel("Peak Throughput  (MB/s)", fontsize=10)
    ax.set_title("Peak Throughput per Tool", fontsize=12,
                 fontweight="bold", pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.7, zorder=0)
    ax.tick_params(axis="x", labelsize=9)
    ax.set_ylim(0, max(peaks) * 1.18)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Plot bamstrom benchmark report")
    parser.add_argument("csv", help="Path to benchmark CSV")
    parser.add_argument(
        "--out", default=None,
        help="Output image path (default: <csv stem>_report.png)",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    out_path = Path(args.out) if args.out else csv_path.with_name(
        csv_path.stem + "_report.png"
    )

    data = load(csv_path)
    if not data:
        sys.exit(f"No valid rows found in {csv_path}")

    fig = plt.figure(figsize=(14, 10), facecolor="white")
    fig.suptitle(
        f"BAM Reader Benchmark — {csv_path.stem}",
        fontsize=15, fontweight="bold", y=0.98,
    )

    gs = GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                  left=0.08, right=0.97, top=0.92, bottom=0.08)

    plot_throughput(fig.add_subplot(gs[0, 0]), data)
    plot_elapsed(fig.add_subplot(gs[0, 1]), data)
    plot_speedup(fig.add_subplot(gs[1, 0]), data)
    plot_peak_bar(fig.add_subplot(gs[1, 1]), data)

    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Report saved: {out_path}")


if __name__ == "__main__":
    main()
