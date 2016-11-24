#!/bin/bash
#$ -cwd
#$ -V
#$ -N sub-wifl
#$ -o sub-wifl.log
#$ -e sub-wifl.error
#$ -pe shared 1
#$ -l h_data=5G,time=10:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# run this in the directory you want the subsampled bams to go.

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java
module load samtools


for i in  /u/nobackup/klohmuel/rbay/WIFL/Plate*/bam/*; do
  j=$(basename $i);
  ~/nobackup-klohmuel/bin/bbmap/reformat.sh in=$i  out=$j  samplereadstarget=120000  mappedonly=t primaryonly=t > $j.stdout
done 