#!/usr/bin/env python3
"""
fastq2tracks v3.0.4 — Chromosome-wide coverage karyogram plot
Replicates the ChIPQC "ChIP Peaks over Chromosomes" panel.

Usage:
    python3 plot_chrom_coverage.py \
        --bedgraph sample1.bg sample2.bg \
        --labels   SAMPLE1 SAMPLE2 \
        --genome   hg38 \
        --chrom-sizes /path/to/hs38n.chrom.sizes \
        --outdir   qc_post_alignment/plots/chromosome_coverage

Input:
    Per-sample bedGraph files produced by bamCoverage --binSize 100000
    OR any bedGraph/bedGraph.gz with 4 columns: chrom start end score

Output:
    <outdir>/<LABEL>_karyogram.png   — per-sample ChIPQC-style panel
    <outdir>/karyogram_all_samples.png — multi-panel (one row per sample,
                                         same chromosome layout)
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── Standard chromosome order ──────────────────────────────────────────────────
HG38_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
MM39_CHROMS = [f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY", "chrM"]

def get_chrom_order(genome):
    return HG38_CHROMS if genome.lower() in ("hg38", "hg19", "hs") else MM39_CHROMS

def read_chrom_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes

def read_bedgraph(path):
    """Return dict: chrom -> list of (start, end, score)."""
    data = defaultdict(list)
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, score = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
            data[chrom].append((start, end, score))
    return data

def plot_karyogram(data, chrom_sizes, chrom_order, label, outpath,
                   max_signal=None, color="#2166ac", bg_color="#e8e8e8",
                   figwidth=12):
    """
    Draw one karyogram panel:
      - one row per chromosome
      - x-axis: chromosomal position
      - bars: normalised signal height within each row
    """
    chroms = [c for c in chrom_order if c in chrom_sizes]
    n = len(chroms)
    if n == 0:
        print(f"  WARN: no chromosomes to plot for {label}", file=sys.stderr)
        return

    global_max = max_signal
    if global_max is None:
        global_max = 0
        for c in chroms:
            for _, _, sc in data.get(c, []):
                if sc > global_max:
                    global_max = sc
    if global_max == 0:
        global_max = 1

    genome_max = max(chrom_sizes.get(c, 1) for c in chroms)

    row_height = 0.6      # fraction of available row for signal
    row_spacing = 1.0     # total height per chromosome row
    fig_height = n * row_spacing * 0.35 + 1.5

    fig, ax = plt.subplots(figsize=(figwidth, fig_height))
    ax.set_xlim(0, genome_max)
    ax.set_ylim(0, n * row_spacing)
    ax.invert_yaxis()

    for i, chrom in enumerate(chroms):
        y_base = i * row_spacing + row_spacing * 0.5
        chrom_len = chrom_sizes.get(chrom, genome_max)

        # Background chromosome bar
        ax.barh(y_base, chrom_len, height=row_height * 0.25,
                left=0, color=bg_color, zorder=1, linewidth=0)

        # Signal bars
        bins = data.get(chrom, [])
        if bins:
            starts = np.array([b[0] for b in bins], dtype=float)
            widths = np.array([b[1] - b[0] for b in bins], dtype=float)
            scores = np.array([b[2] for b in bins], dtype=float)
            # Normalise scores to row_height
            norm_h = (scores / global_max) * row_height
            # Draw as vertical bars centred on y_base
            ax.bar(starts + widths / 2, norm_h, width=widths,
                   bottom=y_base - norm_h / 2,
                   color=color, linewidth=0, zorder=2)

    # Y-axis: chromosome labels
    ax.set_yticks([i * row_spacing + row_spacing * 0.5 for i in range(n)])
    ax.set_yticklabels(chroms, fontsize=7)
    ax.tick_params(axis="y", length=0)

    # X-axis
    ax.set_xlabel("Chromosome Size (bp)", fontsize=9)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:.1e}")
    )
    ax.tick_params(axis="x", labelsize=8)

    ax.set_title(f"ChIP Peaks over Chromosomes — {label}", fontsize=10, pad=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Chromosome name boxes (right side, mimicking ChIPQC)
    for i, chrom in enumerate(chroms):
        y = i * row_spacing + row_spacing * 0.5
        ax.text(genome_max * 1.005, y, chrom,
                va="center", ha="left", fontsize=6,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec="grey", lw=0.5))

    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {outpath}")


def plot_multi_karyogram(all_data, chrom_sizes, chrom_order, labels, outpath,
                         color="#2166ac", bg_color="#e8e8e8", figwidth=14):
    """
    Multi-sample panel: one sub-panel per sample, stacked vertically,
    sharing the same chromosome axis — allows direct comparison.
    """
    chroms = [c for c in chrom_order if c in chrom_sizes]
    n_chroms = len(chroms)
    n_samples = len(labels)
    genome_max = max(chrom_sizes.get(c, 1) for c in chroms)

    # Compute global max signal across all samples for consistent scaling
    global_max = 0
    for data in all_data:
        for c in chroms:
            for _, _, sc in data.get(c, []):
                if sc > global_max:
                    global_max = sc
    if global_max == 0:
        global_max = 1

    row_height = 0.6
    row_spacing = 1.0
    panel_height = n_chroms * row_spacing * 0.35 + 0.8
    fig_height = panel_height * n_samples + 1.0

    fig, axes = plt.subplots(n_samples, 1,
                              figsize=(figwidth, fig_height),
                              sharex=True)
    if n_samples == 1:
        axes = [axes]

    colors_list = plt.cm.tab10.colors
    for s_idx, (ax, data, label) in enumerate(zip(axes, all_data, labels)):
        col = colors_list[s_idx % len(colors_list)]
        ax.set_xlim(0, genome_max)
        ax.set_ylim(0, n_chroms * row_spacing)
        ax.invert_yaxis()

        for i, chrom in enumerate(chroms):
            y_base = i * row_spacing + row_spacing * 0.5
            chrom_len = chrom_sizes.get(chrom, genome_max)
            ax.barh(y_base, chrom_len, height=row_height * 0.25,
                    left=0, color=bg_color, zorder=1, linewidth=0)
            bins = data.get(chrom, [])
            if bins:
                starts = np.array([b[0] for b in bins], dtype=float)
                widths = np.array([b[1] - b[0] for b in bins], dtype=float)
                scores = np.array([b[2] for b in bins], dtype=float)
                norm_h = (scores / global_max) * row_height
                ax.bar(starts + widths / 2, norm_h, width=widths,
                       bottom=y_base - norm_h / 2,
                       color=col, linewidth=0, zorder=2)

        ax.set_yticks([i * row_spacing + row_spacing * 0.5 for i in range(n_chroms)])
        ax.set_yticklabels(chroms, fontsize=6)
        ax.tick_params(axis="y", length=0)
        ax.set_ylabel(label, fontsize=8, rotation=0, labelpad=60, va="center")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[-1].set_xlabel("Chromosome Size (bp)", fontsize=9)
    axes[-1].xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:.1e}")
    )
    axes[-1].tick_params(axis="x", labelsize=8)
    fig.suptitle("ChIP Peaks over Chromosomes — all samples", fontsize=11, y=1.005)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved multi-sample: {outpath}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bedgraph", nargs="+", required=True,
                    help="Per-sample 100 kb bedGraph files (space-separated)")
    ap.add_argument("--labels", nargs="+", required=True,
                    help="Sample labels (same order as --bedgraph)")
    ap.add_argument("--genome", default="hg38",
                    help="Genome build: hg38 or mm39 (default: hg38)")
    ap.add_argument("--chrom-sizes", required=True,
                    help="Path to chrom.sizes file")
    ap.add_argument("--outdir", default=".",
                    help="Output directory (default: current dir)")
    ap.add_argument("--color", default="#2166ac",
                    help="Bar colour for per-sample plots (hex, default: #2166ac)")
    args = ap.parse_args()

    if len(args.bedgraph) != len(args.labels):
        print("ERROR: --bedgraph and --labels must have the same number of items",
              file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    chrom_order = get_chrom_order(args.genome)
    chrom_sizes = read_chrom_sizes(args.chrom_sizes)

    all_data = []
    valid_labels = []
    valid_paths = []

    for bg_path, label in zip(args.bedgraph, args.labels):
        if not os.path.exists(bg_path):
            print(f"  WARN: bedGraph not found, skipping: {bg_path}", file=sys.stderr)
            continue
        print(f"  Reading: {label} ({bg_path})")
        data = read_bedgraph(bg_path)
        all_data.append(data)
        valid_labels.append(label)
        valid_paths.append(bg_path)

        out_png = os.path.join(args.outdir, f"{label}_karyogram.png")
        plot_karyogram(data, chrom_sizes, chrom_order, label, out_png,
                       color=args.color)

    if len(all_data) >= 2:
        multi_out = os.path.join(args.outdir, "karyogram_all_samples.png")
        plot_multi_karyogram(all_data, chrom_sizes, chrom_order,
                             valid_labels, multi_out)
    elif len(all_data) == 1:
        print("  Only one sample — skipping multi-sample panel")

    print("Done.")

if __name__ == "__main__":
    main()
