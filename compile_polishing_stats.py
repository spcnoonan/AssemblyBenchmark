#!/usr/bin/env python3
"""
compile_polishing_stats.py

Compile genome-assembly polishing benchmark metrics (Dorado vs Pilon) into a
single CSV.

Expected directory layout (PREFIX = assembly name):

  <results>/Stats/Dorado/<PREFIX>_DoradoBUSCO/short_summary.*.<PREFIX>_DoradoBUSCO.txt
  <results>/Stats/Dorado/<PREFIX>_DORADO_QUAST/report.tsv
  <results>/Stats/Dorado/<PREFIX>_dorado-zero_cov.txt
  <results>/Stats/Dorado/<PREFIX>_dorado-clipping.txt

  <results>/Stats/Pilon/<PREFIX>_PilonBUSCO/short_summary.*.<PREFIX>_PilonBUSCO.txt
  <results>/Stats/Pilon/<PREFIX>_PILON_QUAST/report.tsv
  <results>/Stats/Pilon/<PREFIX>_pilon-zero_cov.txt
  <results>/Stats/Pilon/<PREFIX>_pilon-clipping.txt

Usage:
  python3 compile_polishing_stats.py <results_dir> [-o output.csv]
"""

import argparse
import csv
import glob
import os
import re
import sys

# ---------------------------------------------------------------------------
# Per-tool file/dir naming conventions (case matters on Linux filesystems).
# ---------------------------------------------------------------------------
TOOLS = {
    "Dorado": {
        "busco_suffix":   "_DoradoBUSCO",
        "quast_suffix":   "_DORADO_QUAST",
        "clip_suffix":    "_dorado-clipping.txt",
        "zerocov_suffix": "_dorado-zero_cov.txt",
    },
    "Pilon": {
        "busco_suffix":   "_PilonBUSCO",
        "quast_suffix":   "_PILON_QUAST",
        "clip_suffix":    "_pilon-clipping.txt",
        "zerocov_suffix": "_pilon-zero_cov.txt",
    },
}

# Thresholds
INTERNAL_LOW = 0.15   # inclusive lower bound for an "internal" clipping event
INTERNAL_HIGH = 0.85  # inclusive upper bound for an "internal" clipping event
LARGE_ZEROCOV_MIN = 1  # a zero-cov locus is "large" when range_size > this value

NA = "NA"

# BUSCO one-line summary, e.g.:
#   C:99.1%[S:99.1%,D:0.0%],F:0.4%,M:0.4%,n:461,E:12.3%
# The E (internal stop codons) field only exists in newer BUSCO; make it optional.
BUSCO_RE = re.compile(
    r"C:(?P<C>[\d.]+)%\[S:(?P<S>[\d.]+)%,D:(?P<D>[\d.]+)%\],"
    r"F:(?P<F>[\d.]+)%,M:(?P<M>[\d.]+)%,n:(?P<n>\d+)"
    r"(?:,E:(?P<E>[\d.]+)%)?"
)


def parse_busco(path):
    """Return (percent_complete_C, percent_internal_stop_E)."""
    try:
        with open(path) as fh:
            text = fh.read()
    except OSError:
        return NA, NA
    m = BUSCO_RE.search(text)
    if not m:
        return NA, NA
    c = m.group("C")
    e = m.group("E") if m.group("E") is not None else NA
    return c, e


def parse_quast(path):
    """
    Parse a QUAST report.tsv (metric<TAB>value) and return the plain metric
    rows (not the '(>= N bp)' threshold variants).
    Returns dict with keys: total_length, n_contigs, n50, n90, largest_contig.
    """
    wanted = {
        "Total length":   "total_length",
        "# contigs":      "n_contigs",
        "N50":            "n50",
        "N90":            "n90",
        "Largest contig": "largest_contig",
    }
    out = {v: NA for v in wanted.values()}
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                metric, value = parts[0], parts[1]
                if metric in wanted:          # exact match => plain row only
                    out[wanted[metric]] = value
    except OSError:
        pass
    return out


def _field4_values(path):
    """Yield the float in field 4 (0-based index 3) for every data row."""
    with open(path) as fh:
        next(fh, None)  # skip header (row 1)
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            try:
                yield float(fields[3])
            except ValueError:
                continue


def parse_clipping(path):
    """
    Return (total_clipping_loci, internal_clipping_events).
    total = data rows (excluding header)
    internal = rows whose field-4 relative position is within [0.15, 0.85]
    """
    total = 0
    internal = 0
    try:
        for val in _field4_values(path):
            total += 1
            if INTERNAL_LOW <= val <= INTERNAL_HIGH:
                internal += 1
    except OSError:
        return NA, NA
    return total, internal


def parse_zerocov(path):
    """
    Return (total_zerocov_loci, large_zerocov_loci).
    total = data rows (excluding header)
    large = rows whose field-4 range_size is > 1
    """
    total = 0
    large = 0
    try:
        for val in _field4_values(path):
            total += 1
            if val > LARGE_ZEROCOV_MIN:
                large += 1
    except OSError:
        return NA, NA
    return total, large


def find_busco_summary(tool_dir, prefix, busco_suffix):
    """Locate the short_summary.*.txt inside the BUSCO output dir."""
    busco_dir = os.path.join(tool_dir, prefix + busco_suffix)
    hits = sorted(glob.glob(os.path.join(busco_dir, "short_summary.*.txt")))
    return hits[0] if hits else None


def discover_prefixes(tool_dir, clip_suffix):
    """Find all PREFIXes in a tool dir via the *-clipping.txt files."""
    prefixes = []
    for p in sorted(glob.glob(os.path.join(tool_dir, "*" + clip_suffix))):
        base = os.path.basename(p)
        prefixes.append(base[: -len(clip_suffix)])
    return prefixes


HEADERS = [
    "PREFIX",
    "tool",
    "total length",
    "contig number",
    "N50",
    "N90",
    "largest contig",
    "percent complete BUSCOs",
    "percent of complete BUSCOs with internal stop codons",
    "total number of clipping events",
    "number of internal clipping events",
    "total number of zero-coverage loci",
    "number of large zero-coverage loci",
]


def process(results_dir):
    stats_dir = os.path.join(results_dir, "Stats")
    if not os.path.isdir(stats_dir):
        # Allow pointing directly at the Stats dir too.
        if os.path.isdir(os.path.join(results_dir, "Dorado")) or \
           os.path.isdir(os.path.join(results_dir, "Pilon")):
            stats_dir = results_dir
        else:
            sys.exit(f"error: no 'Stats' dir found under {results_dir!r}")

    rows = []
    for tool, cfg in TOOLS.items():
        tool_dir = os.path.join(stats_dir, tool)
        if not os.path.isdir(tool_dir):
            continue
        for prefix in discover_prefixes(tool_dir, cfg["clip_suffix"]):
            busco_path = find_busco_summary(tool_dir, prefix, cfg["busco_suffix"])
            quast_path = os.path.join(tool_dir, prefix + cfg["quast_suffix"], "report.tsv")
            clip_path = os.path.join(tool_dir, prefix + cfg["clip_suffix"])
            zcov_path = os.path.join(tool_dir, prefix + cfg["zerocov_suffix"])

            c, e = parse_busco(busco_path) if busco_path else (NA, NA)
            q = parse_quast(quast_path)
            clip_total, clip_internal = parse_clipping(clip_path)
            zcov_total, zcov_large = parse_zerocov(zcov_path)

            rows.append([
                prefix,
                tool,
                q["total_length"],
                q["n_contigs"],
                q["n50"],
                q["n90"],
                q["largest_contig"],
                c,
                e,
                clip_total,
                clip_internal,
                zcov_total,
                zcov_large,
            ])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results_dir", help="Directory containing Stats/Dorado and Stats/Pilon")
    ap.add_argument("-o", "--output", default="polishing_benchmark_summary.csv",
                    help="Output CSV path (default: polishing_benchmark_summary.csv)")
    args = ap.parse_args()

    rows = process(args.results_dir)
    rows.sort(key=lambda r: (r[0], r[1]))  # by PREFIX then tool

    with open(args.output, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADERS)
        w.writerows(rows)

    print(f"Wrote {len(rows)} row(s) to {args.output}")


if __name__ == "__main__":
    main()
