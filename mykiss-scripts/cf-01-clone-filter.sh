#!/bin/bash
#$ -cwd
#$ -V
#$ -N clone_filter
#$ -o clone_filter.log
#$ -e clone_filter.error
#$ -pe shared 1
#$ -l h_data=4G,time=02:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-25:1
#$ -m a

# run this in the directory with all the fastq files, calling it with 2 arguments: 
# 1. the file that has the file name prefixes in 25 rows and 
# 2. the directory where you want the output to do
if [ $# -ne 2 ]; then
  echo "run this in the directory with all the fastq files, calling it with 2 arguments:"
  echo "   1. the file that has the file name prefixes in 25 rows and"
  echo "   2. the directory where you want the output to go"
  exit 1;
fi  

FILE_25=$1
OUTDIR=$2

CLONE_FILTER=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/clone_filter

FILES="$(awk -v N=$SGE_TASK_ID '$1==N {$1=""; print $0}' $FILE_25)"


###Filter duplicate reads from demultiplexed fastq files, creating new files within the OUTDIR directory
for sample in $FILES;
do
  $CLONE_FILTER -1 ${sample}_R1.fastq.gz -2 ${sample}_R2.fastq.gz -i gzfastq -o $OUTDIR >> $OUTDIR/clone_filter_report_${sample}.stdout 2>&1
  
  # note that this line can be removed after we go to the newer version of stacks
  gzip $OUTDIR/${sample}_R1.fastq.fil.fq_1     $OUTDIR/${sample}_R2.fastq.fil.fq_2  # gzip the results right after they are made
done

