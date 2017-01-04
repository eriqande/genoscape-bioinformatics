#!/bin/bash
#$ -cwd
#$ -V
#$ -N snp-array
#$ -o snp-array.log
#$ -e snp-array.error
#$ -pe shared 4
#$ -l h_data=3G,time=24:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-5:1
#$ -m a



source /u/local/Modules/default/init/modules.sh
module load java

GATK_JAR=/u/nobackup/klohmuel/kruegg/bin/GATK-nightly-11-29-16/GenomeAnalysisTK.jar
FASTA=../../Genome/omyV6Chr.fasta
BAM=../alignments/omyV6/MergedBams/omyV6-merged.bam
SEG=$(awk -v N=$SGE_TASK_ID '$1==N {print "-L", $3}' chromo_redos.txt)
OUT=$(awk -v N=$SGE_TASK_ID '$1==N {printf("%03d", $2)}' chromo_redos.txt)

java -Xmx4G -jar $GATK_JAR -T HaplotypeCaller \
 -R $FASTA \
 -I $BAM \
 -stand_call_conf 20.0  \
 -o $OUT.vcf --genotyping_mode DISCOVERY \
 -nct 4 \
  $SEG  > $OUT.stdout 2>$OUT.stderr





