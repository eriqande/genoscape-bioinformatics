#!/bin/bash
#$ -cwd
#$ -V
#$ -N fastqc
#$ -o fastqc.log
#$ -e fastqc.error
#$ -pe shared 1
#$ -l highp,h_data=2G,time=1:20:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from the top level of the Plate_X directory
# using qsub.  i.e. qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh

#. /u/local/Modules/default/init/modules.sh
#module load samtools


source ~/genoscape-bioinformatics/program-defs.sh
source data-defs.sh # call this second so you can override program-defs.sh if you want to

zcat $READ1 | $FASTQC stdin


