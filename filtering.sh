#!/bin/bash
#Script to take adaptor-trimmed nanopore reads and filtering them by different parameters
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
input=$1 # path to the fastq file to filter
id=$2 # name for output files
threads=$3 # number of threads
stats=$4 # path to stats folder

set -e

## Creating output folders

# Subset 1: 50X reads
mkdir -p 50XREADS
mkdir -p 50XREADS/Default
mkdir -p 50XREADS/LengthWeight
mkdir -p 50XREADS/QualityWeight

# Subset 2: 100X reads
mkdir -p 100XREADS
mkdir -p 100XREADS/Default
mkdir -p 100XREADS/LengthWeight
mkdir -p 100XREADS/QualityWeight

## Setting 50X and 100X lengths based on $id value so species-specific genome sizes can be used
base_size_100x=
base_size_50x=

echo "100X genome size: ${base_size_100x}"
echo
echo "500X genome size: ${base_size_50x}"
echo

## Running filtlong

## Subset 1

# Default settings
filtlong --target_bases ${base_size_50x} ${input} | gzip > 50XREADS/Default/${id}_50X_Default.fastq.gz
echo
echo "Filtering complete for 50X_DEFAULT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/ReadFiltering/50X_DEFAULT_NANOPLOT --plots dot --N50 --fastq 50XREADS/Default/${id}_50X_Default.fastq.gz
echo
echo "Stats generated for 50X_DEFAULT_NANOPLOT"
echo

# Length weight
filtlong --target_bases ${base_size_50x} --length_weight 10 ${input} | gzip > 50XREADS/LengthWeight/${id}_50X_LengthWeight.fastq.gz
echo
echo "Filtering complete for 50X_LENGTHWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/ReadFiltering/50X_LENGTHWEIGHT_NANOPLOT --plots dot --N50 --fastq 50XREADS/LengthWeight/${id}_50X_LengthWeight.fastq.gz
echo
echo "Stats generated for 50X_LENGTHWEIGHT_NANOPLOT"
echo

# Quality weight 
filtlong --target_bases ${base_size_50x} --length_weight 0 --window_q_weight 0 ${input} | gzip > 50XREADS/QualityWeight/${id}_50X_QualityWeight.fastq.gz
echo
echo "Filtering complete for 50X_QUALWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/ReadFiltering/50X_QUALWEIGHT_NANOPLOT --plots dot --N50 --fastq 50XREADS/QualityWeight/${id}_50X_QualityWeight.fastq.gz
echo
echo "Stats generated for 50X_QUALWEIGHT_NANOPLOT"
echo

## Subset 2

# Default settings
filtlong --target_bases ${base_size_100x} ${input} | gzip > 100XREADS/Default/${id}_100X_Default.fastq.gz
echo
echo "Filtering complete for 100X_DEFAULT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/ReadFiltering/100X_DEFAULT_NANOPLOT --plots dot --N50 --fastq 100XREADS/Default/${id}_100X_Default.fastq.gz
echo
echo "Stats generated for 100X_DEFAULT_NANOPLOT"
echo

# Length weight
filtlong --target_bases ${base_size_100x} --length_weight 10 ${input} | gzip > 100XREADS/LengthWeight/${id}_100X_LengthWeight.fastq.gz
echo
echo "Filtering complete for 100X_LENGTHWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/ReadFiltering/100X_LENGTHWEIGHT_NANOPLOT --plots dot --N50 --fastq 100XREADS/LengthWeight/${id}_100X_LengthWeight.fastq.gz
echo
echo "Stats generated for 100X_LENGTHWEIGHT_NANOPLOT"
echo

# Quality weight 
filtlong --target_bases ${base_size_100x} --length_weight 0 --window_q_weight 0 ${input} | gzip > 100XREADS/QualityWeight/${id}_100X_QualityWeight.fastq.gz
echo
echo "Filtering complete for 100X_QUALWEIGHT_NANOPLOT. Obtaining stats ..."
echo
NanoPlot -t $threads --verbose -p $id -o $stats/ReadFiltering/100X_QUALWEIGHT_NANOPLOT --plots dot --N50 --fastq 100XREADS/QualityWeight/${id}_100X_QualityWeight.fastq.gz
echo
echo "Stats generated for 100X_QUALWEIGHT_NANOPLOT"
echo