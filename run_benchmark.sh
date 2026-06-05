#!/bin/bash

#SBATCH --job-name=mitis_assembly
#SBATCH --output=/storage/users/cnoonan22/experiments/EXP0011/testrun/mitis_pipeline/assembly.log
#SBATCH --mem=100000
#SBATCH --cpus-per-task=64

set -e

##################################################################################################################################################################
##																																							    
##								***	Script for the assembly benchmark of Eimeria genomes using Nanopore ligation sequencing reads and Illumina Paired-End       
##                                  Reads ***		    				                                                                                        
##															     Version 1.0																			
##																																																					##
##																																								
##																																							    
##																																							    
##								This benchmark includes the following steps:																						
##										1. Nanopore adaptor trimming and removal of DNA CS sequence																						
##										2. Filtration of reads with FiltLong using 3 different parameters																				
##										3. Assembly using Flye or Hifiasm																				
##										4. Generation of comparative assembly stats with BUSCO, assembly-stats, anvio-dev, 																
##																																																						
##																																								
##																																								
##								For the setup of the pipeline, modify the paths below accordingly. This script has been prepared for submission to a			
##								SLURM workload manager, so change the parameters above as necessary. The full pipeline requires two environments: 'assembly' 	
##								and mitohivi_env. These can be installed via conda using the .yml environment files in the GitHub repo.							
##								For best results, the data should be basecalled	using the latest super high accuracy model from Dorado.             			
##																																							    
##								This first version of the script and is subject to modifications if they have been shown to improve the final assembly.     	
##																																							    
##																																							    
##################################################################################################################################################################

#.......................................................
#					SET PATHS
#  These are paths that you will be required to change to suit your particular paths
#.......................................................
# Path to the benchmark folder
WORKINGDIR=
# Path to accessory scripts
SCRIPTS=
# Path to folder paths
DIR_PATHS=
# Path to DNA CS sequence
DNACS=/storage/users/cnoonan22/experiments/EXP0011/data/dna_cs.fasta
# Number of threads
THREADS="48"
# Model used for basecalling (must be converted to Medaka wording)
MODEL=/storage/users/cnoonan22/scratch/basecalling/r1041_e82_400bps_sup_v5.0.0_model_pt.tar.gz 
# Path to BUSCO lineage file
LINEAGE=/storage/users/cnoonan22/tool_testing/busco_test/busco_downloads/lineages/coccidia_odb10
# Path to Dorado bin
DORADO=/storage/users/cnoonan22/bin/dorado-0.7.1-linux-x64/bin/dorado
# Date
DATE=$(date +%y%m%d)

#........................................................
#			CHECKS
#........................................................

echo "activating assembly micromamba environment..."
module load apps/anaconda-4.7.12.tcl
unset PYTHONPATH


cd $WORKINGDIR
echo "Your working directory is currently ${WORKINGDIR}"
echo

STATS=${WORKINGDIR}/Stats
mkdir -p $STATS

echo "Your total stats will be stored in ${STATS}"

#................................................................
#			   BEGINNING BENCHMARK
#................................................................

while read -r data_paths; do

    ### Extract species name from end of path ###
    SPECIES=$(echo "$data_paths" | rev | cut -d"/" -f1 | rev)

    #### Check for which file types are present and set HIFI flag as positive if folder not empty####
    ### Check if output exists. Will skip species if it is (should allow for multiple script submissions)###
    ### Setting path to read files and flagging if absent###
    NANOPORE=
    PACBIO=

    ### Make output folder for everything for each species ###
    SPECIES_OUT=${WORKINGDIR}/${SPECIES}
    mkdir -p ${SPECIES_OUT}
    mkdir -p ${SPECIES_OUT}/LOGS
    mkdir -p ${SPECIES_OUT}/STATS

    cd $SPECIES_OUT


    ### Begin initial read processing
    READPROCESS=${SPECIES_OUT}/ReadProcessing
    mkdir -p $READFILTER
    mkdir -p ${SPECIES_OUT}/STATS/ReadProcessing
    cd $READPROCESS

    preprocessing="${SCRIPTS}/preprocessing.sh $NANOPORE $DORADO $SPECIES $DNACS $THREADS ${SPECIES_OUT}/STATS/ReadProcessing 2>&1 | tee ${SPECIES_OUT}/LOGS/PreProcessing.log"
    echo $preprocessing
    echo
    eval $preprocessing
    echo "Read pre-processing complete."
    echo
    echo "Reads for filtering can be found in $READPROCESS"
    echo


    ### Begin read filtering
    READFILTER=${SPECIES_OUT}/ReadFiltering
    mkdir -p $READFILTER
    mkdir -p ${SPECIES_OUT}/STATS/ReadFiltering
    cd $READFILTER

    readfilter="${SCRIPTS}/filtering.sh ${SPECIES_OUT}/ReadProcessing/${SPECIES}_QualityFilteredReads.fastq $SPECIES $THREADS ${SPECIES_OUT}/STATS/ReadFiltering2>&1 | tee ${SPECIES_OUT}/LOGS/Filtering.log"
    echo $readfilter
    echo
    eval $readfilter
    echo "Read filtering complete."
    echo
    echo "Reads for assembly can be found in $READFILTER"
    echo


    ### Generate draft assemblies for benchmark
    ASSEMBLY=${SPECIES_OUT}/Assembly
    mkdir -p $ASSEMBLY
    mkdir -p ${SPECIES_OUT}/STATS/Assembly
    cd $ASSEMBLY

    denovoassembly="${SCRIPTS}/assembly.sh ${SPECIES_OUT}/ReadFiltering/ $SPECIES $THREADS $LINEAGE ${SPECIES_OUT}/STATS/Assembly 2>&1 | tee ${SPECIES_OUT}/LOGS/Assembly.log"
    echo $denovoassembly
    echo
    eval $denovoassembly
    echo "Genome draft assembly generation complete."
    echo
    echo "Assemblies can be found in $ASSEMBLY"
    echo

done < "$DIR_PATHS"