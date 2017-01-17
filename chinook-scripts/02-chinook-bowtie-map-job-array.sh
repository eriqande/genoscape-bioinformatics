#!/bin/bash
#$ -cwd
#$ -V
#$ -N radMap
#$ -o radMap.log
#$ -e radMap.error
#$ -pe shared 2
#$ -l h_data=4G,time=12:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-283:1
#$ -m a


source /u/local/Modules/default/init/modules.sh
module load samtools
module load bowtie2


GENOME_DB=/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Genome/bowtie2-omyV6/omyV6
GENOME_DB_ABS=$GENOME_DB

RADDIR=/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Prince_etal_raw/RAD_sequence

# get the name for the output directory
BAMOUTDIR=alignments/$(basename $GENOME_DB)/bam
REPORTOUTDIR=alignments/$(basename $GENOME_DB)/reports
mkdir -p $BAMOUTDIR
mkdir -p $REPORTOUTDIR


# get the relevant line from the id file and store as a bash array,
# then get the necessary parts of it
IDFILE=/u/home/k/kruegg/genoscape-bioinformatics/chinook-scripts/chinook_ids.txt
array=($(awk -v N=$SGE_TASK_ID '$1==N' $IDFILE))

theID=${array[1]}
theSM=${array[3]}
theLB=${array[5]}
thePU=${array[2]}

FQ1=$RADDIR/${theID}_R1.fastq.gz
FQ2=$RADDIR/${theID}_R2.fastq.gz

bowtie2 -x $GENOME_DB_ABS --threads 2 -1 $FQ1 -2 $FQ2 \
  --rg-id $theID --rg SM:$theSM --rg LB:$theLB --rg PU:$thePU \
  --rg PL:illumina  2> $REPORTOUTDIR/$theID  | \
  samtools view -bhS - | \
  samtools sort - $BAMOUTDIR/$theID     # 	with newer version of samtools you need this:  samtools sort -  > ../$BAMOUTDIR/$ID.bam

