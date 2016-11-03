#!/bin/bash
#$ -cwd
#$ -V
#$ -N radMap
#$ -o radMap.log
#$ -e radMap.error
#$ -pe shared 8
#$ -l highp,h_data=4G,time=100:00:00
#$ -M rachaelbay@gmail.com
#$ -m bea


source data-defs.sh 
source $PROGDEFS


. /u/local/Modules/default/init/modules.sh
module load samtools
module load bowtie2

PLATE=$SHORTNAME

# get absolute path to the genome database
#if [ $(uname -a | awk '{print $1}') == "Darwin" ]; then
#  GENOME_DB_ABS=$(fp $GENOME_DB)  # this is just here to test on my laptop.  Most people will not have fp. But Mac does not have readlink -f
#else
#  GENOME_DB_ABS=$(readlink -f $GENOME_DB)
#fi

GENOME_DB_ABS=$GENOME_DB   # assume it is not a relative path for now....

# get the name for the outut directory
BAMOUTDIR=alignments/$(basename $GENOME_DB)/bam
REPORTOUTDIR=alignments/$(basename $GENOME_DB)/reports
mkdir -p $BAMOUTDIR
mkdir -p $REPORTOUTDIR


###Make directory for raw bam files
cd dupfiltered


###Align each sample to genome. Note that genome reference must already be built through bowtie2-build
for sample in `ls *.1.1.fq.gz | cut -f1 -d'.'`
do
	ID="$PLATE.$sample"
	bowtie2 -x $GENOME_DB_ABS \
	--threads 8 -1 $sample.1.1.fq.gz -2 $sample.2.2.fq.gz \
	--rg-id $ID --rg SM:$sample --rg LB:$PLATE \
	--rg PU:$PLATE --rg PL:illumina  2> ../$REPORTOUTDIR/$ID  | \
	samtools view -bhS - | \
	samtools sort -  ../$BAMOUTDIR/$ID     # 	with newer version of samtools you need this:  samtools sort -  > ../$BAMOUTDIR/$ID.bam
done

