#!/bin/bash

#SBATCH --job-name=
#SBATCH --output=
#SBATCH --mem=80000
#SBATCH --cpus-per-task=64

set -e

##################################################################################################################################################################
##																																							    
##								***	Script for the polishing benchmark of Eimeria genomes using Nanopore ligation sequencing reads and Illumina Paired-End       
##                                  Reads ***		    				                                                                                        
##															     Version 1.0																			
##																																																					##
##																																								
##																																							    
##								This benchmark includes the following steps:																						
##										1. Align Nanopore reads with Dorado align tool and sort/reheader/index BAM file																					
##										2. Polish using Nanopore reads with Dorado polish tool																			
##										3. Generate quality stats for Dorado polished assembly with BUSCO, QUAST, CRAQ and Anvi'o																		
##										4. Generation of comparative assembly stats with BUSCO, assembly-stats, anvio-dev
##										5. Align Illumina reads with BWA aligner and sort/index BAM file																					
##										6. Polish using Illumina reads with Pilon
##										7. Repeat steps 6 and 7 for a second iteration of Pilon polishing																		
##										8. Generate quality stats for Dorado polished assembly with BUSCO, QUAST, CRAQ and Anvi'o															
##																																																						
##																																								
##								This pipeline takes a text file containing the absolute paths to the assembly, Nanopore reads, Illumina reads and species name as input.
##								The text file should be formatted as follows:
##								/path/to/assembly.fasta	/path/to/nanopore_reads.fastq    /path/to/illumina_reads_1.fastq	/path/to/illumina_reads_2.fastq	Species_name
##																																
##								For the setup of the pipeline, modify the paths below accordingly. This script has been prepared for submission to a			
##								SLURM workload manager, so change the parameters above as necessary. The full pipeline requires the bench_pol environment.
#                               This can be installed via conda using the .yml environment files in the GitHub repo.							
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
WORKINGDIR=/media/lab101/4TBSSD/conor_scratch/assembly/benchmark/polishing
# Path to data paths file
DIR_PATHS=/media/lab101/4TBSSD/conor_scratch/assembly/benchmark/polishing/data_paths.txt
# Number of threads
THREADS="28"
# Path to BUSCO lineage file
LINEAGE=/media/lab101/4TBSSD/conor_scratch/assembly/data/busco_downloads/lineages/coccidia_odb12
# Path to Dorado
DORADO=/media/lab101/4TBSSD/conor_scratch/basecalling/bin/dorado-2.0.0-linux-x64/bin/dorado
# Path to DORADO model
MODEL=/media/lab101/4TBSSD/conor_scratch/assembly/data/models
# Path to Pilon
PILON=/media/lab101/4TBSSD/conor_scratch/assembly/bin/pilon-1.24.jar
# Path to CRAQ
CRAQ=/media/lab101/4TBSSD/conor_scratch/assembly/bin/CRAQ/bin/craq
# Date
DATE=$(date +%y%m%d)

#........................................................
#			CHECKS
#........................................................

cd $WORKINGDIR
echo "Your working directory is currently ${WORKINGDIR}"
echo

STATS=${WORKINGDIR}/Stats
mkdir -p $STATS
mkdir -p ${STATS}/Dorado
mkdir -p ${STATS}/Pilon

echo "Your total stats will be stored in ${STATS}"

echo "activating bench_pol conda environment..."
unset PYTHONPATH
source /home/lab101/anaconda3/etc/profile.d/conda.sh
conda activate bench_pol
echo

#........................................................
#			RUNNING BENCHMARK
#........................................................

while read -r assembly nano_data ill1 ill2 species; do

    echo "Working on $species assembly with $nano_data and $ill1, $ill2 Illumina reads"
    conda activate bench_pol

    ### Make output folder for everything for each species ###
    SPECIES_OUT=${WORKINGDIR}/${species}
    DORADO_POLISH=${WORKINGDIR}/${species}/Dorado
    PILON_POLISH=${WORKINGDIR}/${species}/Pilon
    
    mkdir -p $SPECIES_OUT $DORADO_POLISH $PILON_POLISH

    cd $SPECIES_OUT

    ### Beginning Dorado polish
    echo
    echo "Beginning Dorado polish of $assembly with $nano_data"
    echo

    cd $DORADO_POLISH
    $DORADO aligner $assembly $nano_data -t $THREADS | samtools sort -@ $THREADS -o ${species}_vs_assembly.sorted.bam
    samtools index ${species}_vs_assembly.sorted.bam

    # Adding basecall model to BAM header
    samtools view -H ${species}_vs_assembly.sorted.bam > header.sam
    echo -e '@RG\tID:1\tPL:ONT\tDS:basecall_model=dna_r10.4.1_e8.2_400bps_sup@v5.2.0' >> header.sam
    samtools reheader header.sam ${species}_vs_assembly.sorted.bam > ${species}_vs_assembly.reheadered.bam
    samtools index ${species}_vs_assembly.reheadered.bam
    rm header.sam ${species}_vs_assembly.sorted.bam
    

    $DORADO polish ${species}_vs_assembly.reheadered.bam $assembly -t $THREADS -o ${species}_dorado --device cuda:0 --models-directory $MODEL

    # Obtaining genome stats
    busco -m genome -l $LINEAGE -i ${species}_dorado/consensus.fasta -o ${species}_DoradoBUSCO --offline -c $THREADS --out_path ${STATS}/Dorado
    quast.py ${species}_dorado/consensus.fasta -o ${STATS}/Dorado/${species}_DORADO_QUAST --threads $THREADS
    minimap2 -ax lr:hq ${species}_dorado/consensus.fasta $nano_data -t $THREADS | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o nanopore_vs_${species}_dorado.bam
    samtools index nanopore_vs_${species}_dorado.bam
    perl $CRAQ -g ${species}_dorado/consensus.fasta -sms nanopore_vs_${species}_dorado.bam -ngs $ill1,$ill2 -pl T  -D ${STATS}/Dorado/${species}_CRAQ -t $THREADS

    conda activate anvio-dev
    anvi-script-find-misassemblies -b nanopore_vs_${species}_dorado.bam -o ${species}_dorado -T $THREADS --just-do-it
    mv ${species}_dorado-clipping.txt ${species}_dorado-zero_cov.txt ${STATS}/Dorado/

    ### Begin pilon correction
    cd $PILON_POLISH
    echo
    echo "Beginning Pilon correction of $assembly with $ill1 and $ill2"
    echo

    mkdir -p temp
    mkdir -p temp/iter1
    mkdir -p temp/iter2

    # Iteration one
    conda activate bench_pol

    # Carry out alignment of reads against racon/medaka output
    bwa index $DORADO_POLISH/${species}_dorado/consensus.fasta
    bwa mem -t $THREADS $DORADO_POLISH/${species}_dorado/consensus.fasta $ill1 $ill2 | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o temp/iter1/ReadsVsPolishedAssem_1.bam
    samtools index temp/iter1/ReadsVsPolishedAssem_1.bam 
    echo "bwa alignment done"
    echo

    #### pilon
    java -Xmx80G -jar $PILON --genome $DORADO_POLISH/${species}_dorado/consensus.fasta --bam temp/iter1/ReadsVsPolishedAssem_1.bam --threads $THREADS --outdir temp/iter1/ --output ${species} --tracks --fix all
    echo "iteration one done"
    echo "starting interation two"

    # Iteration two

    # Repeat alignment
    bwa index temp/iter1/${species}.fasta
    bwa mem -t $THREADS temp/iter1/${species}.fasta $ill1 $ill2 | samtools view -bS -@ $THREADS - | samtools sort -@ $THREADS -o temp/iter2/ReadsVsPolishedAssem_2.bam
    samtools index temp/iter2/ReadsVsPolishedAssem_2.bam
    echo "Second alignment done"
    echo

    #### pilon
    java -Xmx80G -jar $PILON --genome temp/iter1/${species}.fasta --bam temp/iter2/ReadsVsPolishedAssem_2.bam --threads $THREADS --outdir temp/iter2/ --output ${species} --tracks --fix all
    echo "iteration two done"

    # Obtaining genome stats
    busco -m genome -l $LINEAGE -i temp/iter2/${species}.fasta -o ${species}_PilonBUSCO --offline -c $THREADS --out_path ${STATS}/Pilon
    quast.py temp/iter2/${species}.fasta -o ${STATS}/Pilon/${species}_PILON_QUAST --threads $THREADS
    minimap2 -ax lr:hq temp/iter2/${species}.fasta $nano_data -t $THREADS | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o nanopore_vs_${species}_pilon.bam
    samtools index nanopore_vs_${species}_pilon.bam
    perl $CRAQ -g temp/iter2/${species}.fasta -sms nanopore_vs_${species}_pilon.bam -ngs $ill1,$ill2 -pl T  -D ${STATS}/Pilon/${species}_CRAQ -t $THREADS

    conda activate anvio-dev
    anvi-script-find-misassemblies -b nanopore_vs_${species}_pilon.bam -o ${species}_pilon -T $THREADS --just-do-it
    mv ${species}_pilon-clipping.txt ${species}_pilon-zero_cov.txt ${STATS}/Pilon/


done < "$DIR_PATHS"