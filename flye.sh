#!/bin/bash
#Script to assemble pre-processed Nanopore reads
#.................................................................................................................
#											SET PATHS BASED ON FLAGS
#.................................................................................................................
filt_dir=$1 # path to the folder containing filtered reads
id=$2 # id to use for files
threads=$3 # number of threads
lineage=$4 # path to BUSCO lineage
genome_size=$5 # genome size
stats=$6 # path to stats folder
flye=$7 # path to pacbio reads
pacbio=$8 # path to flye exe

set -e

echo "activating bench_assm micromamba environment..."
unset PYTHONPATH
source /storage/users/cnoonan22/micromamba/etc/profile.d/mamba.sh
micromamba activate /storage/users/cnoonan22/micromamba/envs/bench_assm
echo

# Creating sub-directories
mkdir -p FLYE
mkdir -p FLYE/NANOPORE
mkdir -p $stats/BUSCO

# For loop to run through the results of find
# Directories etc will be extracted from file paths because they're descriptive

for reads in $(find "$filt_dir" -type f -name "*.fastq.gz"); do

    # Extracting the filtering criteria used from file path
    FILTER_PARAMS=$(basename "$reads" | cut -d"_" -f2,3 | sed 's/.fastq.gz//g')

    ## Running FLYE
    echo
    echo "Running Flye for ${FILTER_PARAMS} reads..."
    echo
    $flye --nano-hq $reads --genome-size $genome_size --threads $threads --out-dir FLYE/NANOPORE/${id}_${FILTER_PARAMS}

    # Generating stats
    echo
    echo "Generating stats ..."
    echo
    quast.py FLYE/NANOPORE/${id}_${FILTER_PARAMS}/assembly.fasta -o $stats/${id}_${FILTER_PARAMS}_QUAST --threads $threads --nanopore $reads
    busco -m genome -l $lineage -i FLYE/NANOPORE/${id}_${FILTER_PARAMS}/assembly.fasta -o ${id}_${FILTER_PARAMS}_BUSCO --offline -c $threads
    mv ${id}_${FILTER_PARAMS}_BUSCO $stats/BUSCO
    minimap2 -ax lr:hq $reads FLYE/NANOPORE/${id}_${FILTER_PARAMS}/assembly.fasta -t $threads | samtools view -b -@ $threads - | samtools sort -@ $threads -o cov.bam
    samtools flagstat -@ $threads cov.bam > $stats/${id}_${FILTER_PARAMS}.flagstats
    samtools coverage cov.bam > $stats/${id}_${FILTER_PARAMS}.cov
    rm cov.bam
    echo
    echo "Stats obtained. Moving on to next assembly ..."
done

# Running FLYE for hifi data if $pacbio contains a file, otherwise returns a warning

if [[ -n "$pacbio" ]]; then
    echo "Found PacBio file: $pacbio"
    mkdir -p FLYE/HIFI

    for reads in $(find ../PacbioFiltering -type f -name "*.fastq.gz"); do

        # Extracting the filtering criteria used from file path
        FILTER_PARAMS=$(basename "$reads" | cut -d"_" -f2,3 | sed 's/_HIFI.fastq.gz//g')

        # Making sub-directories
        mkdir -p FLYE/HIFI/${FILTER_PARAMS}

        ## Running FLYE
        echo
        echo echo "Running Flye for ${FILTER_PARAMS} using Hifi reads..."
        echo
        $flye --pacbio-hifi $reads --genome-size $genome_size --threads $threads --out-dir FLYE/HIFI/${id}_${FILTER_PARAMS}_HIFI


        # Generating stats
        echo
        echo "Generating stats ..."
        echo
        quast.py FLYE/HIFI/${id}_${FILTER_PARAMS}_HIFI/assembly.fasta -o $stats/${id}_${FILTER_PARAMS}_QUAST --threads $threads --nanopore $reads
        busco -m genome -l $lineage -i FLYE/HIFI/${id}_${FILTER_PARAMS}_HIFI/assembly.fasta -o ${id}_${FILTER_PARAMS}_HIFI_BUSCO --offline -c $threads
        mv ${id}_${FILTER_PARAMS}_HIFI_BUSCO $stats/BUSCO
        minimap2 -ax lr:hq $reads FLYE/HIFI/${id}_${FILTER_PARAMS}_HIFI/assembly.fasta -t $threads | samtools view -b -@ $threads - | samtools sort -@ $threads -o cov.bam
        samtools flagstat -@ $threads cov.bam > $stats/${id}_${FILTER_PARAMS}_HIFI.flagstats
        samtools coverage cov.bam > $stats/${id}_${FILTER_PARAMS}_HIFI.cov
        rm cov.bam
        echo
        echo "Stats obtained. Moving on to next assembly ..."
    done
    
else
    echo "No PacBio file found; skipping PacBio analysis"
fi