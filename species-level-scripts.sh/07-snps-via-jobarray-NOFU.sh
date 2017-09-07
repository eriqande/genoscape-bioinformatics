#!/bin/bash
#$ -cwd
#$ -V
#$ -N snp-array
#$ -o snp-array.log
#$ -e snp-array.error
#$ -l h_data=5G,time=24:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-3:1
#$ -m a


# full run will be 286

# TO PREPARE:
# mkdir SNP-array
# cd SNP-array
# samtools view  -H ../MergedBams/NOFUv0-merged.bam | sed 's/SN://g; s/LN://g' | awk '/^@SQ/' | awk -f ~/genoscape-bioinformatics/script/assemble-scaffold-lists.awk > comm_lines.txt
#
#
#


source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java


FASTA=../Genome/Fulmarus_glacialis.fna
BAM=../MergedBams/NOFUv0-merged.bam
COMMS=comm_lines.txt

OUTN=$(printf '%04d' $SGE_TASK_ID)
OUT=$OUTN

SEG=$(awk -F"\t" -v n=$SGE_TASK_ID 'NR == n {print $2}' $COMMS)

java -Xmx2G -jar $GATK_JAR -T HaplotypeCaller \
-R $FASTA \
-I $BAM \
-stand_call_conf 20.0 -stand_emit_conf 20.0 \
-o $OUT.vcf --genotyping_mode DISCOVERY \
 $SEG  > $OUT.stdout 2>$OUT.stderr
  
  