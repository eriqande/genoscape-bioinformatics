#!/bin/bash
#$ -cwd
#$ -V
#$ -N rmdup
#$ -o rmdup.log
#$ -e rmdup.error
#$ -pe shared 1
#$ -l h_data=4G,time=1:00:00
#$ -M eric.anderson@noaa.gov
#$ -t 1-15:1
#$ -m a

# run this in the bam directory after making the bamlist.txt with about 15 bams per numbered line.

source /u/local/Modules/default/init/modules.sh
module load samtools

mkdir -p ../deduped_bams

# get the relevant line from the id file and store as a bash array,
# then get the necessary parts of it
IDFILE=bamlist.txt

LINE="$(awk -v N=$SGE_TASK_ID '$1==N {$1=""; print}' $IDFILE)"

for infile in $LINE; do 
  outfile=../deduped_bams/$infile
  report=../reports/${infile/.bam/.dedup_report}
  echo $outfile $report
  samtools rmdup $infile $outfile 2>$report
done