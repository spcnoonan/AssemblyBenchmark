#!/usr/bin/env python3
"""
Assembly Benchmark Summary Script
==================================
Parses QUAST and BUSCO outputs from a genome assembly benchmark pipeline,
and generates per-species and master summary tables.

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
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_quast(report_path: Path) -> dict:
    """Parse a QUAST report.tsv → dict of metric -> value."""
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
        total, internal_stop_codons, complete_pct, fragmented_pct,
        missing_pct, duplicated_pct
    """
    result = {}
    text = summary_path.read_text()

    notation = re.search(
        r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)",
        text,
    )
    if notation:
        result["complete_pct"]   = float(notation.group(1))
        result["single_pct"]     = float(notation.group(2))
        result["duplicated_pct"] = float(notation.group(3))
        result["fragmented_pct"] = float(notation.group(4))
        result["missing_pct"]    = float(notation.group(5))
        result["total"]          = int(notation.group(6))

    def _grab(pattern, text, cast=int):
        m = re.search(pattern, text)
        return cast(m.group(1)) if m else None

    result["complete"]   = _grab(r"(\d+)\s+Complete BUSCOs", text)
    result["single"]     = _grab(r"(\d+)\s+Complete and single-copy BUSCOs", text)
    result["duplicated"] = _grab(r"(\d+)\s+Complete and duplicated BUSCOs", text)
    result["fragmented"] = _grab(r"(\d+)\s+Fragmented BUSCOs", text)
    result["missing"]    = _grab(r"(\d+)\s+Missing BUSCOs", text)

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
    quast_dirs = sorted(stats_dir.glob("*_QUAST"))
    if not quast_dirs:
        print(f"  [WARN] No *_QUAST directories under {stats_dir}", file=sys.stderr)

    for qdir in quast_dirs:
        name   = qdir.name[: -len("_QUAST")]
        report = qdir / "report.tsv"
        if not report.is_file():
            print(f"  [WARN] report.tsv not found in {qdir}", file=sys.stderr)
            continue

        busco_dir = stats_dir / "BUSCO" / f"{name}_BUSCO"
        busco_txt = None
        if busco_dir.is_dir():
            candidates = list(busco_dir.glob("short_summary.specific.*.txt"))
            if candidates:
                busco_txt = candidates[0]
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

        q = parse_quast(asm["quast_report"])
        for t in contig_thresholds:
            key = f"# contigs (>= {t} bp)"
            row[f"contigs_>={t}bp"] = int(q.get(key, 0) or 0)
        for metric in QUAST_METRICS:
            safe = (metric.replace(" ", "_").replace("(", "").replace(")", "")
                         .replace("%", "pct").replace(">=", "ge"))
            try:
                row[safe] = float(q.get(metric, "nan"))
            except ValueError:
                row[safe] = q.get(metric, "")

        if asm["busco_summary"]:
            b = parse_busco(asm["busco_summary"])
            for k, v in b.items():
                row[f"busco_{k}"] = v
        else:
            for k in ["complete", "single", "duplicated", "fragmented", "missing",
                      "internal_stop_codons", "total", "complete_pct", "fragmented_pct",
                      "missing_pct", "duplicated_pct"]:
                row[f"busco_{k}"] = np.nan

        rows.append(row)

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Summary tables
# ---------------------------------------------------------------------------

SUMMARY_COLS_DISPLAY = [
    ("assembly",                    "Assembly"),
    # QUAST metrics  (keys match what compile_species() stores after sanitising)
    ("#_contigs_ge_0_bp",           "# Contigs (total)"),
    ("Total_length_ge_0_bp",        "Total length (bp)"),
    ("Largest_contig",              "Largest contig (bp)"),
    ("N50",                         "N50"),
    ("N90",                         "N90"),
    ("GC_pct",                      "GC (%)"),
    ("Avg._coverage_depth",         "Avg. coverage depth"),
    # BUSCO metrics
    ("busco_complete_pct",          "BUSCO Complete (%)"),
    ("busco_duplicated_pct",        "BUSCO Duplicated (%)"),
    ("busco_fragmented_pct",        "BUSCO Fragmented (%)"),
    ("busco_missing_pct",           "BUSCO Missing (%)"),
    ("busco_internal_stop_codons",  "Internal stop codons"),
]


def _safe_col(col: str) -> str:
    return (col.replace(" ", "_").replace("(", "").replace(")", "")
               .replace("%", "pct").replace(">=", "ge").lower())


def make_summary_table(df: pd.DataFrame) -> pd.DataFrame:
    """Produce a human-readable summary table from the compiled DataFrame."""
    keep = {}
    for raw, label in SUMMARY_COLS_DISPLAY:
        if raw in df.columns:
            keep[raw] = label
        else:
            norm = _safe_col(raw)
            for col in df.columns:
                if _safe_col(col) == norm:
                    keep[col] = label
                    break

    sub = df[list(keep.keys())].rename(columns=keep).copy()
    for col in sub.select_dtypes("float").columns:
        sub[col] = sub[col].round(2)
    return sub


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Compile QUAST + BUSCO results and generate summary tables.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--base_dir", required=True,
                        help="Root directory containing one sub-folder per species.")
    parser.add_argument("--species", nargs="+", required=True,
                        help="Species folder names (space-separated).")
    parser.add_argument("--out_dir", default="./benchmark_results",
                        help="Output directory (default: ./benchmark_results).")
    parser.add_argument("--contig_thresholds", nargs="+", type=int,
                        default=CONTIG_THRESHOLDS_DEFAULT,
                        help=(f"Contig size thresholds in bp "
                              f"(default: {' '.join(map(str, CONTIG_THRESHOLDS_DEFAULT))})."))
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

        summary  = make_summary_table(df)
        csv_path = sp_out / f"{sp}_summary.csv"
        summary.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")

    if not all_dfs:
        print("\n[ERROR] No data compiled. Check --base_dir and --species.")
        sys.exit(1)

    master_df = pd.concat(all_dfs, ignore_index=True)

    master_summary = make_summary_table(master_df)
    master_csv     = out_dir / "master_summary.csv"
    master_summary.to_csv(master_csv, index=False)
    print(f"\nSaved master summary: {master_csv}")

    print("\n✓ Done. All outputs written to:", out_dir)


if __name__ == "__main__":
    main()