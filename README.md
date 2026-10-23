# Benchmarking *Eimeria* Genome Assembly and Polishing

Scripts for benchmarking how read pre-processing, read filtering and polishing affect the quality of *Eimeria* genome assemblies built from Oxford Nanopore ligation sequencing reads, with optional PacBio HiFi reads and Illumina paired-end reads.

The benchmark runs in two stages:

1. **Assembly benchmark** (`run_benchmark.sh`): for each species, raw Nanopore reads are cleaned, filtered into six subsets that differ in target coverage and Filtlong weighting, and each subset is assembled with Flye. Every assembly is then assessed with QUAST, BUSCO and read-mapping statistics.
2. **Polishing benchmark** (`polishing_benchmark.sh`): selected assemblies from stage 1 are polished with Dorado (Nanopore reads) and then with two rounds of Pilon (Illumina reads). Each polished assembly is assessed with QUAST, BUSCO, CRAQ and Anvi'o.

Two Python scripts then collect the statistics from each stage into CSV tables for plotting and analysis.

Assemblies were also generated with MaSuRCA for comparison, but that was done in a separate script and is not part of this repository.

---

## Repository contents

| File | Role |
|---|---|
| [`run_benchmark.sh`](run_benchmark.sh) | Main SLURM pipeline for the assembly benchmark; calls the three accessory scripts below |
| [`preprocessing.sh`](preprocessing.sh) | Accessory script: adapter trimming, quality filtering and DNA CS removal |
| [`filtering.sh`](filtering.sh) | Accessory script: Filtlong read subsetting (6 subsets per read type) |
| [`flye.sh`](flye.sh) | Accessory script: Flye assembly and assembly QC for every read subset |
| [`polishing_benchmark.sh`](polishing_benchmark.sh) | SLURM pipeline for the polishing benchmark (Dorado, then 2× Pilon) |
| [`compile_benchmark_stats.py`](compile_benchmark_stats.py) | Compiles QUAST and BUSCO results from the assembly benchmark |
| [`compile_polishing_stats.py`](compile_polishing_stats.py) | Compiles QUAST, BUSCO and Anvi'o results from the polishing benchmark |

---

## Assembly benchmark: `run_benchmark.sh`

### Workflow

`run_benchmark.sh` reads a list of species data folders and runs the three accessory scripts in turn for each species. Each step's output is logged to `<species>/LOGS/`.

**1. Pre-processing** (`preprocessing.sh`)

1. The raw Nanopore FASTQ is split into chunks of 5,000,000 lines to limit memory use.
2. Each chunk is adapter-trimmed with Porechop, duplicate read IDs are renamed with `seqkit rename`, and reads are filtered with Filtlong (`--min_length 1000 --min_mean_q 10`).
3. The chunks are merged, and reads that map to the DNA control strand (DNA CS) sequence with minimap2 are removed, keeping only unmapped reads.
4. Read statistics are generated with NanoPlot and the final reads are compressed with `bgzip`.

**2. Filtering** (`filtering.sh`)

The pre-processed reads are subset with Filtlong to a target number of bases equal to **50×** and **100×** the species' genome size. Each coverage level is filtered with three weightings, giving six subsets:

| Subset name | Filtlong settings | Effect |
|---|---|---|
| `Default` | `--target_bases` only | Filtlong's default balance of read length and quality |
| `LengthWeight` | `--length_weight 10` | Strongly favours longer reads |
| `QualityWeight` | `--length_weight 0 --window_q_weight 0` | Selects on mean read quality only |

Each subset is summarised with NanoPlot. If PacBio HiFi reads are present, they are combined with the pre-processed Nanopore reads and the same six subsets are generated from the combined file.

**3. Assembly** (`flye.sh`)

Every subset is assembled with Flye (`--nano-hq` for Nanopore subsets, `--pacbio-hifi` for the combined HiFi subsets) using the species' estimated genome size. Each assembly is then assessed with:

- **QUAST**, with the subset reads supplied so that coverage depth is reported
- **BUSCO** in genome mode against the supplied lineage (`coccidia_odb12`)
- **minimap2 + samtools**: reads are mapped back to the assembly to produce `samtools flagstat` and `samtools coverage` reports

### Requirements

- SLURM (the script requests 64 CPUs and about 80 GB RAM)
- Conda/micromamba environments:
  - `assembly`: Porechop, seqkit, minimap2, samtools (with `bgzip`), NanoPlot
  - `bench_preprocess`: Filtlong, NanoPlot
  - `bench_assm`: QUAST, BUSCO, minimap2, samtools
- Flye (tested with v2.9.6) and Filtlong executables
- The DNA CS reference FASTA (Oxford Nanopore lambda control strand)
- A downloaded BUSCO lineage (`coccidia_odb12`), as BUSCO runs with `--offline`

The accessory scripts activate their own micromamba environments from hard-coded paths near the top of each script, so change these to match your system.

### Input

`DIR_PATHS` points to a text file with one species data folder per line:

```
/path/to/data/brunetti
/path/to/data/mitis
```

The **species name is taken from the folder name**, so each folder must be named after its species. Each folder should contain:

- A Nanopore FASTQ ending in `_SUP.fastq.gz` (up to one subfolder deep)
- *Optional:* a `pacbio/` subfolder containing a PacBio HiFi `*.fastq.gz` file

Genome sizes used for Filtlong and Flye are set per species inside `run_benchmark.sh`, based on initial assemblies:

| Species | Genome size (bp) |
|---|---|
| `necatrix` | 58,000,000 |
| `brunetti` | 77,000,000 |
| `mitis` | 72,000,000 |
| `lata` | 50,000,000 |
| `zaria` | 60,000,000 |

The script exits if a folder name does not match one of these species. Add a new `elif` branch to benchmark another species.

### Configuration

Edit the `SET PATHS` block at the top of `run_benchmark.sh`:

| Variable | Description |
|---|---|
| `WORKINGDIR` | Benchmark output directory |
| `SCRIPTS` | Folder containing `preprocessing.sh`, `filtering.sh` and `flye.sh` |
| `DIR_PATHS` | Text file listing species data folders (see above) |
| `DNACS` | DNA CS reference FASTA |
| `THREADS` | Should match `--cpus-per-task` |
| `LINEAGE` | BUSCO lineage directory |
| `FLYE` | Flye executable |
| `FILTLONG` | Filtlong executable |

### Usage

```bash
sbatch run_benchmark.sh
```

Species whose output folder already exists in `WORKINGDIR` are skipped, so the script can be resubmitted after adding new species to `DIR_PATHS` without re-running finished ones.

### Outputs

Assemblies and their statistics are named `<species>_<coverage>_<weighting>`, for example `brunetti_50X_LengthWeight`. HiFi assemblies have a `_HIFI` suffix.

```
WORKINGDIR/
└── <species>/
    ├── LOGS/                         # PreProcessing.log, Filtering.log, Assembly.log
    ├── ReadProcessing/
    │   └── <species>_PreProcessedReads.fastq.gz
    ├── ReadFiltering/
    │   ├── 50XREADS/                 # <species>_50X_{Default,LengthWeight,QualityWeight}.fastq.gz
    │   └── 100XREADS/                # <species>_100X_{Default,LengthWeight,QualityWeight}.fastq.gz
    ├── PacbioFiltering/              # Only if PacBio reads were found
    │   ├── <species>_NP_PB_combined_HIFI.fastq.gz
    │   ├── 50XREADS/
    │   └── 100XREADS/
    ├── Assembly/
    │   └── FLYE/
    │       ├── NANOPORE/<species>_<cov>_<weighting>/        # Flye output, incl. assembly.fasta
    │       └── HIFI/<species>_<cov>_<weighting>_HIFI/
    └── STATS/
        ├── ReadProcessing/NanoPlot/
        ├── ReadFiltering/            # <COV>_<WEIGHTING>_NANOPLOT/ (HIFI_ prefix for HiFi subsets)
        └── Assembly/
            ├── <name>_QUAST/         # QUAST report
            ├── BUSCO/<name>_BUSCO/   # BUSCO results
            ├── <name>.flagstats      # samtools flagstat
            └── <name>.cov            # samtools coverage
```

---

## Polishing benchmark: `polishing_benchmark.sh`

This script is run after the assembly benchmark on the assemblies chosen for polishing.

### Workflow

For each assembly listed in `DIR_PATHS`:

**Dorado polishing**

1. Nanopore reads are aligned to the draft assembly with `dorado aligner`, then sorted and indexed.
2. A read group recording the basecalling model (`dna_r10.4.1_e8.2_400bps_sup@v5.2.0`) is added to the BAM header, since `dorado polish` uses it to choose a polishing model.
3. The assembly is polished with `dorado polish` on the GPU (`--device cuda:0`).

**Pilon polishing** (run on the Dorado-polished assembly)

4. Illumina reads are aligned with `bwa mem`, then sorted and indexed.
5. The assembly is polished with Pilon (`--fix all`).
6. Steps 4–5 are repeated on the first Pilon output for a second iteration.

**Quality assessment** (after Dorado and again after the second Pilon iteration)

- **BUSCO** in genome mode against the supplied lineage
- **QUAST**
- **CRAQ**, using Nanopore reads mapped with minimap2 and the Illumina reads, to find regional and structural assembly errors
- **Anvi'o** `anvi-script-find-misassemblies`, using the Nanopore read alignment, to report clipping events and zero-coverage regions

### Requirements

- SLURM (the script requests 64 CPUs and about 80 GB RAM)
- A CUDA-capable GPU for `dorado polish`
- Conda environments:
  - `bench_pol`: samtools, BUSCO, QUAST, minimap2, bwa, Java, Perl
  - `anvio-dev`: Anvi'o
- Dorado (tested with v2.0.0) and a local directory of Dorado polishing models
- Pilon JAR (tested with v1.24)
- CRAQ
- A downloaded BUSCO lineage (`coccidia_odb12`)

Pilon runs with a Java heap of 80 GB (`-Xmx80G`). Change this to suit the memory you have available.

### Input

`DIR_PATHS` points to a whitespace-separated text file with one assembly per line:

```
/path/to/assembly.fasta    /path/to/nanopore_reads.fastq    /path/to/illumina_R1.fastq    /path/to/illumina_R2.fastq    Species_name
```

The last column is used to name all outputs for that assembly, so use a unique name for each line.

If your reads were basecalled with a different model, change the `DS:basecall_model=` value in the read group line of the script.

### Configuration

Edit the `SET PATHS` block at the top of `polishing_benchmark.sh`:

| Variable | Description |
|---|---|
| `WORKINGDIR` | Polishing benchmark output directory |
| `DIR_PATHS` | Text file listing assemblies and reads (see above) |
| `THREADS` | Should match `--cpus-per-task` |
| `LINEAGE` | BUSCO lineage directory |
| `DORADO` | Dorado executable |
| `MODEL` | Directory of Dorado polishing models |
| `PILON` | Pilon JAR |
| `CRAQ` | CRAQ executable |

### Usage

```bash
sbatch polishing_benchmark.sh
```

### Outputs

```
WORKINGDIR/
├── <species>/
│   ├── Dorado/
│   │   ├── <species>_dorado/consensus.fasta    # Dorado-polished assembly
│   │   ├── <species>_vs_assembly.reheadered.bam
│   │   └── nanopore_vs_<species>_dorado.bam
│   └── Pilon/
│       ├── temp/iter1/                          # First Pilon iteration
│       ├── temp/iter2/<species>.fasta           # Final Dorado + 2× Pilon assembly
│       └── nanopore_vs_<species>_pilon.bam
└── Stats/
    ├── Dorado/
    │   ├── <species>_DoradoBUSCO/
    │   ├── <species>_DORADO_QUAST/
    │   ├── <species>_CRAQ/
    │   ├── <species>_dorado-clipping.txt
    │   └── <species>_dorado-zero_cov.txt
    └── Pilon/
        ├── <species>_PilonBUSCO/
        ├── <species>_PILON_QUAST/
        ├── <species>_CRAQ/
        ├── <species>_pilon-clipping.txt
        └── <species>_pilon-zero_cov.txt
```

---

## Compiling statistics

### `compile_benchmark_stats.py`

Collects the QUAST and BUSCO results for every assembly in the assembly benchmark into summary CSV tables.

```bash
python3 compile_benchmark_stats.py \
    --base_dir /path/to/WORKINGDIR \
    --species brunetti mitis necatrix \
    --out_dir ./benchmark_results
```

| Argument | Description |
|---|---|
| `--base_dir` | The assembly benchmark `WORKINGDIR` (one subfolder per species) |
| `--species` | Species folder names to include, space-separated |
| `--out_dir` | Output directory (default: `./benchmark_results`) |

For each species, the script finds every `STATS/Assembly/<name>_QUAST/report.tsv` and its matching `STATS/Assembly/BUSCO/<name>_BUSCO/short_summary.specific.*.txt`. Assemblies with missing BUSCO results are kept, with empty BUSCO columns.

**Outputs**

- `<out_dir>/<species>/<species>_summary.csv`: one row per assembly for that species
- `<out_dir>/master_summary.csv`: all assemblies from all species

Columns: assembly name, number of contigs, total length, largest contig, N50, N90, GC (%), average coverage depth, and BUSCO complete, duplicated, fragmented and missing (%) and the number of internal stop codons.

**Requirements:** Python 3 with `pandas` and `numpy`.

### `compile_polishing_stats.py`

Collects the QUAST, BUSCO and Anvi'o results for the Dorado and Pilon assemblies into a single CSV.

```bash
python3 compile_polishing_stats.py /path/to/WORKINGDIR -o polishing_benchmark_summary.csv
```

The first argument is the polishing benchmark `WORKINGDIR` (or its `Stats/` folder directly). Assemblies are found from the `*-clipping.txt` files in `Stats/Dorado/` and `Stats/Pilon/`. Any missing value is written as `NA`.

**Output:** one row per assembly per polishing tool, sorted by assembly name and then tool, with these columns:

| Column | Source |
|---|---|
| `PREFIX`, `tool` | Assembly name and polisher (`Dorado` or `Pilon`) |
| `total length`, `contig number`, `N50`, `N90`, `largest contig` | QUAST |
| `percent complete BUSCOs` | BUSCO (C) |
| `percent of complete BUSCOs with internal stop codons` | BUSCO (E), only reported by newer BUSCO versions |
| `total number of clipping events` | Anvi'o clipping report: all rows |
| `number of internal clipping events` | Clipping events with a relative position between 0.15 and 0.85 along the contig, i.e. not near contig ends |
| `total number of zero-coverage loci` | Anvi'o zero-coverage report: all rows |
| `number of large zero-coverage loci` | Zero-coverage regions longer than 1 bp |

CRAQ results are not included in this table and should be read from the `<species>_CRAQ/` folders directly.

**Requirements:** Python 3 (standard library only).
