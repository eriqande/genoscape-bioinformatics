#!/bin/bash
#$ -cwd
#$ -V
#$ -N snp-array
#$ -o snp-array.log
#$ -e snp-array.error
#$ -l h_data=5G,time=24:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-187:1
#$ -m a

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java


FASTA=../Genome/GCA_001281735.1_ASM128173v1_genomic.fna
BAM=../MergedBams/twenty-birds-ZOLAv0-merged.bam
COMMS=comm_lines.txt

OUTN=$(printf '%04d' $SGE_TASK_ID)
OUT=$OUTN

SEG=$(awk -F"\t" -v n=$SGE_TASK_ID '$1 == n {print $3}' $COMMS)

java -Xmx2G -jar $GATK_JAR -T HaplotypeCaller \
-R $FASTA \
-I $BAM \
-stand_call_conf 20.0 -stand_emit_conf 20.0 \
-o $OUT.vcf --genotyping_mode DISCOVERY \
 $SEG  > $OUT.stdout 2>$OUT.stderr
  
  
  