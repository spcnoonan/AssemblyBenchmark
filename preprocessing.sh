#!/bin/bash
#Script to take input Nanopore fastq files and pre-processing them to get ready for assembly
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
input=$1 # path to the nanopore fastq files to process
id=$2 # name for output files
dnacs=$3 # path to DNA CS fasta file
threads=$4 # number of threads
filtlong=$5 # path to filtlong
stats=$6 # path to stats folder

set -e

echo "activating assembly micromamba environment..."
unset PYTHONPATH
source /storage/users/cnoonan22/micromamba/etc/profile.d/mamba.sh
micromamba activate /storage/users/cnoonan22/micromamba/envs/assembly
echo


# Read splitting for memory purposes
zcat $input | split -l 5000000 - split

# Creating a folder to store temporary files
mkdir -p temp

## Raw read processing
for file in split*; do

    # Read filtering
    echo "Working with file $file ..."
    echo
    echo "Beginning Porechop ..."
    porechop -i $file -o temp/${file}_FilteredReads1.fastq.gz --threads $threads
    echo
    echo "Renaming duplicate IDs with seqkit ..."
    echo
    seqkit rename temp/${file}_FilteredReads1.fastq.gz --threads $threads > temp/${file}_FilteredReads2.fastq
    echo
    echo "Beginning Chopper ..."
    $filtlong --min_length 1000 --min_mean_q 10 temp/${file}_FilteredReads2.fastq > temp/${file}_FilteredReads3.fastq
    echo "Quality trimming complete."
    echo "Removing file $file ..."
    rm $file temp/${file}_FilteredReads1.fastq temp/${file}_FilteredReads2.fastq
    echo

done

cat temp/*_FilteredReads3.fastq > temp/Merged_FilteredReads.fastq

# DNA CS removal
echo "Removing DNA CS sequences from reads ..."
minimap2 -ax lr:hq $dnacs temp/Merged_FilteredReads.fastq -t $threads | samtools view  -f 4 -@ $threads -b - | samtools sort -@ $threads | samtools fastq -@ $threads > ${id}_PreProcessedReads.fastq
NanoPlot -t $threads --verbose -p $id -o $stats/ReadProcessing/NanoPlot --plots dot --N50 --fastq ${id}_PreProcessedReads.fastq
bgzip ${id}_PreProcessedReads.fastq

# Removing temp files
rm -r temp/
echo