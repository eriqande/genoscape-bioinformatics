#!/bin/bash
#$ -cwd
#$ -V
#$ -N bowtie-build-db
#$ -o bowtie-build-db.log
#$ -e bowtie-build-db.error
#$ -pe shared 1
#$ -l highp,h_data=8G,time=1:20:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from inside the "Reference" directory which should
# be inside the species directory.  It should be run under qsub.


module load bowtie


source reference-defs.sh # must have name of gzipped fasta file and short-name for reference output directory
source $PROGDEFS

mkfifo infasta

$ZCAT  $GENOME_FASTA_GZ > infasta &


bowtie-build infasta $GENOME_OUTNAME

