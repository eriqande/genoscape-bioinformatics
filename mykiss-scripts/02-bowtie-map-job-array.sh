#!/bin/bash
#$ -cwd
#$ -V
#$ -N radMap
#$ -o radMap.log
#$ -e radMap.error
#$ -pe shared 4
#$ -l h_data=4G,time=12:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-242:1
#$ -m a


source data-defs.sh 
source $PROGDEFS


source /u/local/Modules/default/init/modules.sh
module load samtools
module load bowtie2

# get absolute path to the genome database
#if [ $(uname -a | awk '{print $1}') == "Darwin" ]; then
#  GENOME_DB_ABS=$(fp $GENOME_DB)  # this is just here to test on my laptop.  Most people will not have fp. But Mac does not have readlink -f
#else
#  GENOME_DB_ABS=$(readlink -f $GENOME_DB)
#fi

GENOME_DB=/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Genome/bowtie2-omyV6/omyV6
GENOME_DB_ABS=$GENOME_DB

RADDIR=/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Prince_etal_raw/RAD_sequence

# get the name for the output directory
BAMOUTDIR=alignments/$(basename $GENOME_DB)/bam
REPORTOUTDIR=alignments/$(basename $GENOME_DB)/reports
mkdir -p $BAMOUTDIR
mkdir -p $REPORTOUTDIR


SGE_TASK_ID=1

# get the relevant line from the id file and store as a bash array,
# then get the necessary parts of it
IDFILE=/u/home/k/kruegg/genoscape-bioinformatics/mykiss-scripts/mykiss_ids.txt
array=($(awk -v N=$SGE_TASK_ID '$1==N' $IDFILE))

theID=${array[1]}
theSM=${array[3]}
theLB=${array[5]}
thePU=${array[2]}

FQ1=$RADDIR/${theID}_R1.fastq.gz
FQ2=$RADDIR/${theID}_R2.fastq.gz

bowtie2 -x $GENOME_DB_ABS \
	--threads 4 -1 $FQ1 -2 $FQ2 \
	--rg-id $theID --rg SM:$theSM --rg LB:$theLB \
	--rg PU:$thePU --rg PL:illumina  2> ../$REPORTOUTDIR/$theID  | \
	samtools view -bhS - | \
	samtools sort -  ../$BAMOUTDIR/$theID     # 	with newer version of samtools you need this:  samtools sort -  > ../$BAMOUTDIR/$ID.bam

