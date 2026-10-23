#!/bin/bash
#Script to take adaptor-trimmed nanopore reads and filtering them by different parameters
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
input=$1 # path to the nanopore fastq file to filter
id=$2 # name for output files
threads=$3 # number of threads
genome_size=$4 # genome size
filtlong=$5 # path to filtlong
stats=$6 # path to stats folder
pacbio=$7 # path to pacbio file (if one is set)

set -e

echo "activating bench_preprocess micromamba environment..."
unset PYTHONPATH
source ~/path/to/micromamba/etc/profile.d/mamba.sh
micromamba activate bench_preprocess
echo

## Creating output folders

# Subset 1: 50X reads
mkdir -p 50XREADS

# Subset 2: 100X reads
mkdir -p 100XREADS

## Setting 50X and 100X lengths based on $id value so species-specific genome sizes can be used
base_size_100x=$(( $genome_size * 100 ))
base_size_50x=$(( $genome_size * 50 ))

echo "100X genome size: ${base_size_100x}"
echo
echo "50X genome size: ${base_size_50x}"
echo

## Running filtlong

## Subset 1

# Default settings
filtlong --target_bases ${base_size_50x} ${input} | gzip > 50XREADS/${id}_50X_Default.fastq.gz
echo
echo "Filtering complete for 50X_DEFAULT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/50X_DEFAULT_NANOPLOT --plots dot --N50 --fastq 50XREADS/${id}_50X_Default.fastq.gz
echo
echo "Stats generated for 50X_DEFAULT_NANOPLOT"
echo

# Length weight
filtlong --target_bases ${base_size_50x} --length_weight 10 ${input} | gzip > 50XREADS/${id}_50X_LengthWeight.fastq.gz
echo
echo "Filtering complete for 50X_LENGTHWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/50X_LENGTHWEIGHT_NANOPLOT --plots dot --N50 --fastq 50XREADS/${id}_50X_LengthWeight.fastq.gz
echo
echo "Stats generated for 50X_LENGTHWEIGHT_NANOPLOT"
echo

# Quality weight 
filtlong --target_bases ${base_size_50x} --length_weight 0 --window_q_weight 0 ${input} | gzip > 50XREADS/${id}_50X_QualityWeight.fastq.gz
echo
echo "Filtering complete for 50X_QUALWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/50X_QUALWEIGHT_NANOPLOT --plots dot --N50 --fastq 50XREADS/${id}_50X_QualityWeight.fastq.gz
echo
echo "Stats generated for 50X_QUALWEIGHT_NANOPLOT"
echo

## Subset 2

# Default settings
filtlong --target_bases ${base_size_100x} ${input} | gzip > 100XREADS/${id}_100X_Default.fastq.gz
echo
echo "Filtering complete for 100X_DEFAULT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/100X_DEFAULT_NANOPLOT --plots dot --N50 --fastq 100XREADS/${id}_100X_Default.fastq.gz
echo
echo "Stats generated for 100X_DEFAULT_NANOPLOT"
echo

# Length weight
filtlong --target_bases ${base_size_100x} --length_weight 10 ${input} | gzip > 100XREADS/${id}_100X_LengthWeight.fastq.gz
echo
echo "Filtering complete for 100X_LENGTHWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/100X_LENGTHWEIGHT_NANOPLOT --plots dot --N50 --fastq 100XREADS/${id}_100X_LengthWeight.fastq.gz
echo
echo "Stats generated for 100X_LENGTHWEIGHT_NANOPLOT"
echo

# Quality weight 
filtlong --target_bases ${base_size_100x} --length_weight 0 --window_q_weight 0 ${input} | gzip > 100XREADS/${id}_100X_QualityWeight.fastq.gz
echo
echo "Filtering complete for 100X_QUALWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/100X_QUALWEIGHT_NANOPLOT --plots dot --N50 --fastq 100XREADS/${id}_100X_QualityWeight.fastq.gz
echo
echo "Stats generated for 100X_QUALWEIGHT_NANOPLOT"
echo

###############################################################
#           If there are PacBio reads also ...                #
###############################################################

if [[ -n "$pacbio" ]]; then

    # Creating separate output
    mkdir -p ../PacbioFiltering
    cat "$input" "$pacbio" > ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz

    # Subset 1: 50X reads
    mkdir -p ../PacbioFiltering/50XREADS

    # Subset 2: 100X reads
    mkdir -p ../PacbioFiltering/100XREADS

    # Default settings
    filtlong --target_bases ${base_size_50x} ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz | gzip > ../PacbioFiltering/50XREADS/${id}_50X_Default_HIFI.fastq.gz
    echo
    echo "Filtering complete for 50X_DEFAULT_NANOPLOT. Obtaining stats ..."
    echo
    NanoPlot -t $threads --verbose -p $id -o $stats/HIFI_50X_DEFAULT_NANOPLOT --plots dot --N50 --fastq ../PacbioFiltering/50XREADS/${id}_50X_Default_HIFI.fastq.gz
    echo
    echo "Stats generated for 50X_DEFAULT_NANOPLOT"
    echo

    # Length weight
    filtlong --target_bases ${base_size_50x} --length_weight 10 ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz | gzip > ../PacbioFiltering/50XREADS/${id}_50X_LengthWeight_HIFI.fastq.gz
    echo
    echo "Filtering complete for 50X_LENGTHWEIGHT_NANOPLOT. Obtaining stats ..."
    echo
    NanoPlot -t $threads --verbose -p $id -o $stats/HIFI_50X_LENGTHWEIGHT_NANOPLOT --plots dot --N50 --fastq ../PacbioFiltering/50XREADS/${id}_50X_LengthWeight_HIFI.fastq.gz
    echo
    echo "Stats generated for 50X_LENGTHWEIGHT_NANOPLOT"
    echo

    # Quality weight 
    filtlong --target_bases ${base_size_50x} --length_weight 0 --window_q_weight 0 ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz | gzip > ../PacbioFiltering/50XREADS/${id}_50X_QualityWeight_HIFI.fastq.gz
    echo
    echo "Filtering complete for 50X_QUALWEIGHT_NANOPLOT. Obtaining stats ..."
    echo
    NanoPlot -t $threads --verbose -p $id -o $stats/HIFI_50X_QUALWEIGHT_NANOPLOT --plots dot --N50 --fastq ../PacbioFiltering/50XREADS/${id}_50X_QualityWeight_HIFI.fastq.gz
    echo
    echo "Stats generated for 50X_QUALWEIGHT_NANOPLOT"
    echo

    ## Subset 2

    # Default settings
    filtlong --target_bases ${base_size_100x} ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz | gzip > ../PacbioFiltering/100XREADS/${id}_100X_Default_HIFI.fastq.gz
    echo
    echo "Filtering complete for 100X_DEFAULT_NANOPLOT. Obtaining stats ..."
    echo
    NanoPlot -t $threads --verbose -p $id -o $stats/HIFI_100X_DEFAULT_NANOPLOT --plots dot --N50 --fastq ../PacbioFiltering/100XREADS/${id}_100X_Default_HIFI.fastq.gz
    echo
    echo "Stats generated for 100X_DEFAULT_NANOPLOT"
    echo

    # Length weight
    filtlong --target_bases ${base_size_100x} --length_weight 10 ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz | gzip > ../PacbioFiltering/100XREADS/${id}_100X_LengthWeight_HIFI.fastq.gz
    echo
    echo "Filtering complete for 100X_LENGTHWEIGHT_NANOPLOT. Obtaining stats ..."
    echo
    NanoPlot -t $threads --verbose -p $id -o $stats/HIFI_100X_LENGTHWEIGHT_NANOPLOT --plots dot --N50 --fastq ../PacbioFiltering/100XREADS/${id}_100X_LengthWeight_HIFI.fastq.gz
    echo
    echo "Stats generated for 100X_LENGTHWEIGHT_NANOPLOT"
    echo

    # Quality weight 
    filtlong --target_bases ${base_size_100x} --length_weight 0 --window_q_weight 0 ../PacbioFiltering/${id}_NP_PB_combined_HIFI.fastq.gz | gzip > ../PacbioFiltering/100XREADS/${id}_100X_QualityWeight_HIFI.fastq.gz
    echo
    echo "Filtering complete for 100X_QUALWEIGHT_NANOPLOT. Obtaining stats ..."
    echo
    NanoPlot -t $threads --verbose -p $id -o $stats/HIFI_100X_QUALWEIGHT_NANOPLOT --plots dot --N50 --fastq ../PacbioFiltering/100XREADS/${id}_100X_QualityWeight_HIFI.fastq.gz
    echo
    echo "Stats generated for 100X_QUALWEIGHT_NANOPLOT"
    echo
else
    echo "No PacBio file found; skipping PacBio analysis"
fi
