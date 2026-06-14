#!/usr/bin/env python3
"""
Assembly Benchmark Summary Script
==================================
Parses QUAST and BUSCO outputs from a genome assembly benchmark pipeline,
generates per-species and master summary tables, and produces publication-
quality plots.

Directory structure assumed:
  <base_dir>/
    <species>/
      STATS/
        Assembly/
          <name>_QUAST/report.tsv         (QUAST output)
          BUSCO/<name>_BUSCO/
            short_summary.specific.*.txt  (BUSCO output)

Usage:
  python assembly_benchmark_summary.py \\
      --base_dir ~/path/to/benchmark \\
      --species brunetti coccidium sarcocystis neospora theileria \\
      --out_dir ./benchmark_results \\
      --contig_thresholds 10000 50000 500000 1000000

Author: Generated for assembly benchmark analysis
"""

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Colour palettes
# ---------------------------------------------------------------------------
BUSCO_COLOURS = {
    "Complete (S)":           "#2196F3",   # blue
    "Duplicated (D)":         "#FF9800",   # orange
    "Fragmented (F)":         "#FFC107",   # amber
    "Missing (M)":            "#F44336",   # red
    "Internal stop codons":   "#9C27B0",   # purple
}

CONTIG_COLOURS = [
    "#1565C0", "#1E88E5", "#42A5F5", "#90CAF9",
    "#B3E5FC", "#E1F5FE", "#E8F5E9", "#C8E6C9",
]


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_quast(report_path: Path) -> dict:
    """
    Parse a QUAST report.tsv.
    Returns a dict of metric -> value.
    """
    data = {}
    with open(report_path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                data[parts[0].strip()] = parts[1].strip()
    return data


def parse_busco(summary_path: Path) -> dict:
    """
    Parse a BUSCO short_summary.specific.*.txt file.
    Returns a dict with keys:
        complete, single, duplicated, fragmented, missing,
        total, internal_stop_codons, complete_pct, fragmented_pct, missing_pct
    """
    result = {}
    text = summary_path.read_text()

    # Main BUSCO notation line:  C:99.1%[S:99.1%,D:0.0%],F:0.4%,M:0.4%,n:461,E:12.7%
    notation = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
        text,
    )
    if notation:
        result["complete_pct"]    = float(notation.group(1))
        result["single_pct"]      = float(notation.group(2))
        result["duplicated_pct"]  = float(notation.group(3))
        result["fragmented_pct"]  = float(notation.group(4))
        result["missing_pct"]     = float(notation.group(5))
        result["total"]           = int(notation.group(6))

    # Absolute counts
    def _grab(pattern, text, cast=int):
        m = re.search(pattern, text)
        return cast(m.group(1)) if m else None

    result["complete"]   = _grab(r"(\d+)\s+Complete BUSCOs", text)
    result["single"]     = _grab(r"(\d+)\s+Complete and single-copy BUSCOs", text)
    result["duplicated"] = _grab(r"(\d+)\s+Complete and duplicated BUSCOs", text)
    result["fragmented"] = _grab(r"(\d+)\s+Fragmented BUSCOs", text)
    result["missing"]    = _grab(r"(\d+)\s+Missing BUSCOs", text)

    # Internal stop codons – note in parentheses after Complete BUSCOs line
    isc = re.search(r"of which (\d+) contain internal stop codons", text)
    result["internal_stop_codons"] = int(isc.group(1)) if isc else 0

    return result


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------

def find_assemblies(species_base: Path):
    """
    Walk <species_base>/STATS/Assembly and return a list of dicts:
        {name, quast_report, busco_summary}
    where `name` is taken from the QUAST directory stem (minus '_QUAST').
    """
    stats_dir = species_base / "STATS" / "Assembly"
    if not stats_dir.is_dir():
        print(f"  [WARN] Stats directory not found: {stats_dir}", file=sys.stderr)
        return []

    assemblies = []

    # Collect QUAST directories
    quast_dirs = sorted(stats_dir.glob("*_QUAST"))
    if not quast_dirs:
        print(f"  [WARN] No *_QUAST directories under {stats_dir}", file=sys.stderr)

    for qdir in quast_dirs:
        name = qdir.name[: -len("_QUAST")]   # strip trailing '_QUAST'
        report = qdir / "report.tsv"
        if not report.is_file():
            print(f"  [WARN] report.tsv not found in {qdir}", file=sys.stderr)
            continue

        # Match BUSCO directory by the same assembly name
        busco_dir = stats_dir / "BUSCO" / f"{name}_BUSCO"
        busco_txt = None
        if busco_dir.is_dir():
            candidates = list(busco_dir.glob("short_summary.specific.*.txt"))
            if candidates:
                busco_txt = candidates[0]   # take the first match
            else:
                print(f"  [WARN] No short_summary file in {busco_dir}", file=sys.stderr)
        else:
            print(f"  [WARN] BUSCO directory not found: {busco_dir}", file=sys.stderr)

        assemblies.append({"name": name, "quast_report": report, "busco_summary": busco_txt})

    return assemblies


# ---------------------------------------------------------------------------
# Data compilation
# ---------------------------------------------------------------------------

CONTIG_THRESHOLDS_DEFAULT = [1_000, 5_000, 10_000, 25_000, 50_000]
QUAST_METRICS = [
    "# contigs (>= 0 bp)",
    "# contigs (>= 1000 bp)",
    "# contigs (>= 5000 bp)",
    "# contigs (>= 10000 bp)",
    "# contigs (>= 25000 bp)",
    "# contigs (>= 50000 bp)",
    "Total length (>= 0 bp)",
    "Largest contig",
    "N50",
    "N90",
    "L50",
    "Avg. coverage depth",
    "GC (%)",
]


def compile_species(species_name: str, species_base: Path, contig_thresholds: list) -> pd.DataFrame:
    """Return a tidy DataFrame with one row per assembly for a given species."""
    assemblies = find_assemblies(species_base)
    if not assemblies:
        return pd.DataFrame()

    rows = []
    for asm in assemblies:
        row = {"species": species_name, "assembly": asm["name"]}

        # QUAST
        q = parse_quast(asm["quast_report"])
        # Contig counts above each threshold
        for t in contig_thresholds:
            key = f"# contigs (>= {t} bp)"
            row[f"contigs_>={t}bp"] = int(q.get(key, 0) or 0)
        # Other assembly stats
        for metric in QUAST_METRICS:
            safe = metric.replace(" ", "_").replace("(", "").replace(")", "").replace("%", "pct").replace(">=", "ge")
            try:
                row[safe] = float(q.get(metric, "nan"))
            except ValueError:
                row[safe] = q.get(metric, "")

        # BUSCO
        if asm["busco_summary"]:
            b = parse_busco(asm["busco_summary"])
            for k, v in b.items():
                row[f"busco_{k}"] = v
        else:
            for k in ["complete", "single", "duplicated", "fragmented", "missing",
                      "internal_stop_codons", "total", "complete_pct", "fragmented_pct",
                      "missing_pct"]:
                row[f"busco_{k}"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Summary tables
# ---------------------------------------------------------------------------

SUMMARY_COLS_DISPLAY = [
    ("assembly",            "Assembly"),
    ("busco_complete_pct",  "BUSCO Complete (%)"),
    ("busco_fragmented_pct","BUSCO Fragmented (%)"),
    ("busco_missing_pct",   "BUSCO Missing (%)"),
    ("busco_internal_stop_codons", "Internal stop codons"),
    ("N50",                 "N50"),
    ("Largest_contig",      "Largest contig (bp)"),
    ("#_contigs___>=_0_bp_","# Contigs"),
    ("GC___pct_",           "GC (%)"),
    ("Avg._coverage_depth", "Avg. coverage depth"),
]


def _safe_col(col: str) -> str:
    """Normalise a column name for display lookup."""
    return col.replace(" ", "_").replace("(", "").replace(")", "").replace("%", "pct").replace(">=", "ge")


def make_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """Produce a human-readable summary table from the compiled DataFrame."""
    keep = {}
    for raw, label in SUMMARY_COLS_DISPLAY:
        # try direct match first, then normalised
        if raw in df.columns:
            keep[raw] = label
        else:
            norm = _safe_col(raw)
            for col in df.columns:
                if _safe_col(col) == norm:
                    keep[col] = label
                    break

    sub = df[list(keep.keys())].rename(columns=keep)
    # Round floats
    for col in sub.select_dtypes("float").columns:
        sub[col] = sub[col].round(2)
    return sub


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _ax_style(ax, title="", xlabel="", ylabel=""):
    ax.set_title(title, fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_contig_bars(df: pd.DataFrame, thresholds: list, out_path: Path, title: str):
    """
    Stacked horizontal bar chart: number of contigs above each size threshold.
    Each bar segment represents the *increment* between consecutive thresholds.
    """
    assemblies = df["assembly"].tolist()
    n = len(assemblies)
    if n == 0:
        return

    # Build threshold columns present in df
    th_cols = [f"contigs_>={t}bp" for t in thresholds if f"contigs_>={t}bp" in df.columns]
    th_vals = [int(t) for t in thresholds if f"contigs_>={t}bp" in df.columns]

    if not th_cols:
        return

    # Compute incremental segments (largest → smallest threshold)
    # e.g. contigs >= 50 kb, then the extra ones >=25 kb but <50 kb, etc.
    sorted_pairs = sorted(zip(th_vals, th_cols), reverse=True)
    sorted_vals  = [p[0] for p in sorted_pairs]
    sorted_cols  = [p[1] for p in sorted_pairs]

    segments     = []   # (label, array)
    prev         = np.zeros(n)
    for i, (tval, tcol) in enumerate(sorted_pairs):
        cur  = df[tcol].fillna(0).values.astype(float)
        segs = np.maximum(cur - prev, 0)
        prev = cur
        if i == 0:
            label = f"≥ {tval:,} bp"
        else:
            lo = sorted_vals[i - 1]
            label = f"{lo:,} – {tval:,} bp"
        segments.append((label, segs))

    # Reverse so smallest threshold is at the "bottom" (right-most in stacked bar)
    segments = list(reversed(segments))

    colours = CONTIG_COLOURS[: len(segments)]
    fig, ax = plt.subplots(figsize=(max(8, n * 1.2), 5))

    x      = np.arange(n)
    bottom = np.zeros(n)
    for (label, segs), colour in zip(segments, colours):
        ax.bar(x, segs, bottom=bottom, color=colour, label=label, width=0.65, edgecolor="white", linewidth=0.4)
        bottom += segs

    ax.set_xticks(x)
    ax.set_xticklabels(assemblies, rotation=30, ha="right", fontsize=8)
    ax.legend(title="Contig size", fontsize=7, title_fontsize=8,
               loc="upper right", framealpha=0.8)
    _ax_style(ax, title=title, ylabel="Number of contigs")
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_busco_bars(df: pd.DataFrame, out_path: Path, title: str):
    """
    Stacked bar chart of BUSCO components per assembly.
    Order (bottom to top): Complete(S), Duplicated, Fragmented, Missing.
    Internal stop codons shown as a hatched overlay on the Complete bar.
    """
    assemblies = df["assembly"].tolist()
    n = len(assemblies)
    if n == 0:
        return

    def _col(c):
        return df[c].fillna(0).values.astype(float) if c in df.columns else np.zeros(n)

    complete   = _col("busco_single")      # single-copy complete
    duplicated = _col("busco_duplicated")
    fragmented = _col("busco_fragmented")
    missing    = _col("busco_missing")
    isc        = _col("busco_internal_stop_codons")  # subset of complete

    x = np.arange(n)
    w = 0.65

    fig, ax = plt.subplots(figsize=(max(8, n * 1.2), 5))

    b0 = np.zeros(n)
    ax.bar(x, complete,   width=w, bottom=b0, color=BUSCO_COLOURS["Complete (S)"],
           label="Complete (S)", edgecolor="white", linewidth=0.4)
    b1 = b0 + complete
    ax.bar(x, duplicated, width=w, bottom=b1, color=BUSCO_COLOURS["Duplicated (D)"],
           label="Duplicated (D)", edgecolor="white", linewidth=0.4)
    b2 = b1 + duplicated
    ax.bar(x, fragmented, width=w, bottom=b2, color=BUSCO_COLOURS["Fragmented (F)"],
           label="Fragmented (F)", edgecolor="white", linewidth=0.4)
    b3 = b2 + fragmented
    ax.bar(x, missing,    width=w, bottom=b3, color=BUSCO_COLOURS["Missing (M)"],
           label="Missing (M)", edgecolor="white", linewidth=0.4)

    # Internal stop codons: hatched overlay starting from 0
    ax.bar(x, isc, width=w, bottom=b0,
           color="none", edgecolor=BUSCO_COLOURS["Internal stop codons"],
           hatch="///", linewidth=0, label="Internal stop codons (subset of C)")

    ax.set_xticks(x)
    ax.set_xticklabels(assemblies, rotation=30, ha="right", fontsize=8)
    ax.legend(fontsize=7, loc="lower right", framealpha=0.8)
    _ax_style(ax, title=title, ylabel="Number of BUSCOs")

    # Annotate total BUSCOs above each bar
    totals = complete + duplicated + fragmented + missing
    for xi, tot in zip(x, totals):
        if tot > 0:
            ax.text(xi, tot + totals.max() * 0.01, f"{int(tot)}", ha="center",
                    va="bottom", fontsize=7, color="#333")

    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Master (multi-species) plots  – grouped by parameter type
# ---------------------------------------------------------------------------

def _parameter_from_name(assembly_name: str) -> str:
    """
    Infer a parameter label from an assembly name.
    Expects a structure like:  <species>_<depth>X_<params>
    The species prefix (everything up to the first digit group followed by 'X')
    is stripped; what remains is treated as the parameter identifier.

    Fallback: return the full name unchanged.
    """
    # Try to strip a leading "<word>_<digits>X_" prefix
    m = re.match(r"^[A-Za-z]+_\d+X_(.+)$", assembly_name)
    if m:
        return m.group(1)
    # Try stripping a leading "<word>_" prefix only
    m2 = re.match(r"^[A-Za-z]+_(.+)$", assembly_name)
    if m2:
        return m2.group(1)
    return assembly_name


def plot_master_busco(master_df: pd.DataFrame, out_path: Path):
    """
    Multi-species BUSCO bar chart, assemblies grouped by parameter type.
    Each parameter group is separated by a small gap; species are colour-coded.
    """
    df = master_df.copy()
    df["param"] = df["assembly"].apply(_parameter_from_name)

    species_list = sorted(df["species"].unique())
    param_list   = sorted(df["param"].unique())

    sp_colours = plt.cm.tab10(np.linspace(0, 0.9, len(species_list)))
    sp_colour_map = dict(zip(species_list, sp_colours))

    # Build x-positions with gaps between parameter groups
    gap       = 0.5
    n_species = len(species_list)
    positions = []
    labels    = []
    x         = 0
    group_mid = {}
    for param in param_list:
        grp_positions = []
        for sp in species_list:
            sub = df[(df["param"] == param) & (df["species"] == sp)]
            if sub.empty:
                x += 0.8   # placeholder
                grp_positions.append(None)
                labels.append("")
            else:
                grp_positions.append(x)
                labels.append(sp)
                x += 0.8
        group_mid[param] = np.mean([p for p in grp_positions if p is not None])
        x += gap

    fig, ax = plt.subplots(figsize=(max(12, len(param_list) * n_species * 0.9), 5))

    plotted_species = set()
    x = 0
    for param in param_list:
        for sp in species_list:
            sub = df[(df["param"] == param) & (df["species"] == sp)]
            if sub.empty:
                x += 0.8
                continue
            row = sub.iloc[0]
            complete   = float(row.get("busco_single",   0) or 0)
            duplicated = float(row.get("busco_duplicated", 0) or 0)
            fragmented = float(row.get("busco_fragmented", 0) or 0)
            missing    = float(row.get("busco_missing",  0) or 0)
            isc        = float(row.get("busco_internal_stop_codons", 0) or 0)
            col = sp_colour_map[sp]
            label = sp if sp not in plotted_species else "_nolegend_"
            plotted_species.add(sp)
            ax.bar(x, complete,   width=0.7, bottom=0,
                   color=col, alpha=0.9, label=label, edgecolor="white", linewidth=0.3)
            ax.bar(x, duplicated, width=0.7, bottom=complete,
                   color=col, alpha=0.6, edgecolor="white", linewidth=0.3)
            ax.bar(x, fragmented, width=0.7, bottom=complete + duplicated,
                   color="gold", alpha=0.8, edgecolor="white", linewidth=0.3)
            ax.bar(x, missing,    width=0.7, bottom=complete + duplicated + fragmented,
                   color="tomato", alpha=0.9, edgecolor="white", linewidth=0.3)
            # ISC overlay
            ax.bar(x, isc, width=0.7, bottom=0,
                   color="none", edgecolor="purple", hatch="///",
                   linewidth=0, zorder=3)
            x += 0.8
        x += gap

    # Parameter group labels below x-axis
    for param, mid in group_mid.items():
        ax.text(mid, -ax.get_ylim()[1] * 0.07, param, ha="center", va="top",
                fontsize=8, fontweight="bold", rotation=20)

    ax.set_xticks([])
    handles = [mpatches.Patch(color=sp_colour_map[sp], label=sp) for sp in species_list]
    handles += [
        mpatches.Patch(color="gold",  label="Fragmented"),
        mpatches.Patch(color="tomato", label="Missing"),
        mpatches.Patch(facecolor="none", edgecolor="purple", hatch="///",
                       label="Internal stop codons"),
    ]
    ax.legend(handles=handles, fontsize=7, loc="lower right", framealpha=0.8, ncol=2)
    _ax_style(ax, title="Master BUSCO summary – grouped by parameter", ylabel="Number of BUSCOs")
    ax.margins(x=0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_master_contig(master_df: pd.DataFrame, thresholds: list, out_path: Path):
    """
    Multi-species contig bar chart, assemblies grouped by parameter type.
    """
    df = master_df.copy()
    df["param"] = df["assembly"].apply(_parameter_from_name)

    species_list = sorted(df["species"].unique())
    param_list   = sorted(df["param"].unique())

    sp_colours = plt.cm.tab10(np.linspace(0, 0.9, len(species_list)))
    sp_colour_map = dict(zip(species_list, sp_colours))

    th_cols = [f"contigs_>={t}bp" for t in thresholds if f"contigs_>={t}bp" in df.columns]
    th_vals = [int(t) for t in thresholds if f"contigs_>={t}bp" in df.columns]

    gap = 0.5
    fig, ax = plt.subplots(figsize=(max(12, len(param_list) * len(species_list) * 0.9), 5))

    plotted_species = set()
    x = 0
    group_mid = {}
    for param in param_list:
        xs = []
        for sp in species_list:
            sub = df[(df["param"] == param) & (df["species"] == sp)]
            if sub.empty:
                x += 0.8
                xs.append(None)
                continue
            row = sub.iloc[0]
            col = sp_colour_map[sp]
            label = sp if sp not in plotted_species else "_nolegend_"
            plotted_species.add(sp)
            # Stacked by threshold increment (largest first)
            sorted_pairs = sorted(zip(th_vals, th_cols), reverse=True)
            bottom = 0
            prev   = 0
            for i, (tval, tcol) in enumerate(sorted_pairs):
                cur = float(row.get(tcol, 0) or 0)
                seg = max(cur - prev, 0)
                prev = cur
                alpha = 0.3 + 0.7 * (i / max(len(sorted_pairs) - 1, 1))
                ax.bar(x, seg, width=0.7, bottom=bottom,
                       color=col, alpha=alpha, label=label if i == 0 else "_nolegend_",
                       edgecolor="white", linewidth=0.3)
                bottom += seg
            xs.append(x)
            x += 0.8
        group_mid[param] = np.mean([p for p in xs if p is not None]) if any(xs) else x

        x += gap

    for param, mid in group_mid.items():
        ax.text(mid, -ax.get_ylim()[1] * 0.07, param, ha="center", va="top",
                fontsize=8, fontweight="bold", rotation=20)

    ax.set_xticks([])
    handles = [mpatches.Patch(color=sp_colour_map[sp], label=sp) for sp in species_list]
    ax.legend(handles=handles, fontsize=7, loc="upper right", framealpha=0.8)
    _ax_style(ax, title="Master contig size summary – grouped by parameter",
              ylabel="Number of contigs")
    ax.margins(x=0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compile QUAST + BUSCO results and generate summary tables/plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--base_dir", required=True,
        help="Root directory that contains one sub-folder per species.",
    )
    parser.add_argument(
        "--species", nargs="+", required=True,
        help="Species folder names (space-separated).",
    )
    parser.add_argument(
        "--out_dir", default="./benchmark_results",
        help="Output directory for tables and figures (default: ./benchmark_results).",
    )
    parser.add_argument(
        "--contig_thresholds", nargs="+", type=int,
        default=CONTIG_THRESHOLDS_DEFAULT,
        help=(
            "Contig size thresholds in bp used for the stacked bar chart "
            f"(default: {' '.join(map(str, CONTIG_THRESHOLDS_DEFAULT))})."
        ),
    )
    args = parser.parse_args()

    base_dir   = Path(args.base_dir).expanduser()
    out_dir    = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    thresholds = sorted(args.contig_thresholds)

    all_dfs = []

    for sp in args.species:
        sp_dir = base_dir / sp
        print(f"\n── Processing species: {sp}  ({sp_dir})")
        if not sp_dir.is_dir():
            print(f"  [WARN] Directory not found, skipping: {sp_dir}", file=sys.stderr)
            continue

        df = compile_species(sp, sp_dir, thresholds)
        if df.empty:
            print(f"  [WARN] No data found for {sp}")
            continue
        all_dfs.append(df)

        sp_out = out_dir / sp
        sp_out.mkdir(exist_ok=True)

        # --- per-species summary table ---
        summary = make_summary_table(df)
        csv_path = sp_out / f"{sp}_summary.csv"
        summary.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")

        # --- per-species contig chart ---
        plot_contig_bars(
            df, thresholds,
            out_path=sp_out / f"{sp}_contigs.png",
            title=f"{sp} – Contig size distribution",
        )

        # --- per-species BUSCO chart ---
        plot_busco_bars(
            df,
            out_path=sp_out / f"{sp}_BUSCO.png",
            title=f"{sp} – BUSCO summary",
        )

    if not all_dfs:
        print("\n[ERROR] No data compiled. Check your --base_dir and --species arguments.")
        sys.exit(1)

    master_df = pd.concat(all_dfs, ignore_index=True)

    # --- master summary table ---
    master_summary = make_summary_table(master_df)
    master_csv = out_dir / "master_summary.csv"
    master_summary.to_csv(master_csv, index=False)
    print(f"\nSaved master summary: {master_csv}")

    # --- master plots (grouped by parameter) ---
    plot_master_busco(
        master_df,
        out_path=out_dir / "master_BUSCO.png",
    )
    plot_master_contig(
        master_df, thresholds,
        out_path=out_dir / "master_contigs.png",
    )

    print("\n✓ Done. All outputs written to:", out_dir)


if __name__ == "__main__":
    main()