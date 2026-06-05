#!/bin/bash
#Script to assemble pre-processed Nanopore reads
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
filt_dir=$1 # path to the folder containing filtered reads
threads=$2 # number of threads
id=$3 # id to use for files
lineage=$4 # path to BUSCO lineage
stats=$5 # path to stats folder
pacbio=$6

set -e

# Creating sub-directories
mkdir -p FLYE
mkdir -p HIFIASM
mkdir -p $stats/BUSCO

# Setting genome size based on ID
genome_size=
echo "Genome size: ${genome_size}"
echo

# For loop to run through the results of find
# Directories etc will be extracted from file paths because they're descriptive

for reads in $(find "$filt_dir" -type f -name "*.fastq.gz"); do

    # Extracting the filtering criteria used from file path
    FILTER_PARAMS=

    # Making sub-directories
    mkdir -p FLYE/${FILTER_PARAMS}
    mkdir -p HIFIASM/${FILTER_PARAMS}

    ## Running FLYE
    echo
    echo "Running Flye for ${FILTER_PARAMS} reads..."
    echo
    flye --nano-hq $reads --genome-size $genome_size --threads $threads --out-dir FLYE/${id}_${FILTER_PARAMS}

    # Generating stats
    echo
    echo "Generating stats ..."
    echo
    assembly-stats FLYE/${id}_${FILTER_PARAMS}/assembly.fasta > $stats/Assembly/${id}_${FILTER_PARAMS}_stats.txt
    busco -m genome -l $lineage -i FLYE/${id}_${FILTER_PARAMS}/assembly.fasta -o ${id}_${FILTER_PARAMS}_BUSCO --offline
    mv ${id}_${FILTER_PARAMS}_BUSCO $stats/BUSCO

    # Running FLYE for hifi daya if $pacbio contains a file, otherwise returns a warning
    # Run using just nanopore, just pacbio, and hybrid

    echo
    echo "Running Flye for ${FILTER_PARAMS} using Hifi reads..."
    echo
    flye --pacbio-hifi $pacbio --genome-size $genome_size --threads $threads --out-dir FLYE/${id}_${FILTER_PARAMS}_HIFI

    # Generating stats
    echo
    echo "Generating stats ..."
    echo
    assembly-stats FLYE/${id}_${FILTER_PARAMS}_HIFI/assembly.fasta > $stats/Assembly/${id}_${FILTER_PARAMS}_HIFI_stats.txt
    busco -m genome -l $lineage -i FLYE/${id}_${FILTER_PARAMS}_HIFI/assembly.fasta -o ${id}_${FILTER_PARAMS}_HIFI_BUSCO --offline
    mv ${id}_${FILTER_PARAMS}_HIFI_BUSCO $stats/BUSCO

done