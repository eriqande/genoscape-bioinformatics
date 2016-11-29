#!/bin/bash
#$ -cwd
#$ -V
#$ -N snp-array
#$ -o snp-array.log
#$ -e snp-array.error
#$ -l h_data=10G,time=24:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-4:1
#$ -m a


# this is to be run to get the 4 cases that errored out previously

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java

# override the one in MODULE_SOURCE
GATK_JAR=/u/nobackup/klohmuel/kruegg/bin/GATK-nightly-11-29-16/GenomeAnalysisTK.jar

FASTA=../Genome/GCA_001281735.1_ASM128173v1_genomic.fna
BAM=../MergedBams/twenty-birds-ZOLAv0-merged.bam
COMMS=comm_lines.txt

SEG=$(awk -F"\t" -v n=$SGE_TASK_ID 'NR == n {print $3}' $COMMS)

OUTN=$(awk -F"\t" -v n=$SGE_TASK_ID 'NR == n {printf("%04d", $1)}' $COMMS)
OUT=$OUTN

java -Xmx4G -jar $GATK_JAR -T HaplotypeCaller \
-R $FASTA \
-I $BAM \
-stand_call_conf 20.0 -stand_emit_conf 20.0 \
-o $OUT.vcf --genotyping_mode DISCOVERY \
 $SEG  > $OUT.stdout 2>$OUT.stderr
  
  
  