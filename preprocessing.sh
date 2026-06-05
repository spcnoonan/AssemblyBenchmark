#!/bin/bash
#Script to take input Nanopore fastq files and pre-processing them to get ready for assembly
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
input=$1 # path to the  fastq files to process
dorado=$2 # path to Dorado
id=$3 # name for output files
dnacs=$4 # path to DNA CS fasta file
threads=$5 # number of threads
stats=$6 # path to stats folder

set -e


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
    echo "Beginning Dorado ..."
    $dorado trim temp/${file}_FilteredReads1.fastq.gz | samtools fastq -@ $threads > temp/${file}_FilteredReads2.fastq
    echo "Beginning Chopper ..."
    chopper -l 1000 -q 10 -i temp/${file}_FilteredReads2.fastq | gzip > temp/${file}_FilteredReads3.fastq.gz
    echo "Quality trimming complete."
    echo "Removing file $file ..."
    rm $file temp/${file}_FilteredReads1.fastq.gz temp/${file}_FilteredReads2.fastq
    echo

done

cat temp/*_FilteredReads3.fastq.gz > temp/Merged_FilteredReads.fastq.gz

# DNA CS removal
echo "Removing DNA CS sequences from reads ..."
minimap2 -ax lr:hq $dnacs temp/Merged_FilteredReads.fastq.gz -t $threads | samtools view  -f 4 -@ $threads -b - | samtools sort -@ $threads | samtools fastq -@ $threads > ${id}_QualityFilteredReads.fastq
NanoPlot -t $threads --verbose -p $id -o $stats/ReadProcessing/NanoPlot --plots dot --N50 --fastq ${id}_QualityFilteredReads.fastq
bgzip ${id}_QualityFilteredReads.fastq

# Removing temp files
rm -r temp/
echo