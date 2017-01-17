#!/bin/bash
#$ -cwd
#$ -V
#$ -N merge-bams
#$ -o merge-bams.log
#$ -e merge-bams.error
#$ -pe shared 1
#$ -l h_data=5G,time=24:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea


# run this in the omyV6 directory above "deduped"

source /u/local/Modules/default/init/modules.sh
PICARD_JAR=/u/nobackup/klohmuel/kruegg/bin/picard.jar

module load java
module load samtools


BAM_COMMS="$(for VAR in deduped_bams/*.bam; do printf 'I=%s ' $VAR; done)"


if [ ! -d MergedBams ]; then
 mkdir MergedBams
fi

java -Xmx2G -jar $PICARD_JAR MergeSamFiles $BAM_COMMS  OUTPUT=MergedBams/omyV6-merged.bam  SORT_ORDER=coordinate


echo "Now indexing merged bam"

samtools index MergedBams/omyV6-merged.bam

