#!/bin/bash
#$ -cwd
#$ -V
#$ -N clone_filter
#$ -o clone_filter.log
#$ -e clone_filter.error
#$ -pe shared 1
#$ -l highp,h_data=16G,time=20:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# note that the -o and -e directives above seem important to keep the script
# from trying to write to a default output location that it is unable to write to

# This should be run from the top level of the Plate_X directory
# using qsub.  i.e. qsub ~/genoscape-bioinformatics/plate-processing-scripts/02-rad-process.sh

# this is embarrasingly parallel over the input files, but I suspect it is mostly i/o bound 
# anyway, so I will just run it all in serial...

source ~/genoscape-bioinformatics/program-defs.sh
source data-defs.sh # call this second so you can override program-defs.sh if you want to


###Move orphaned reads to their own directory
mkdir orphans
mv demultiplexed/*rem* orphans 


###Make directory for duplicate filtered reads
mkdir dupfiltered
cd demultiplexed

###Filter duplicate reads from demultiplexed fastq files, creating new files within the dupfiltered directory
for sample in `ls *.1.fq.gz | cut -f1 -d'.'`
do
  $CLONE_FILTER -1 $sample.1.fq.gz -2 $sample.2.fq.gz -i gzfastq -o ../dupfiltered >> ../dupfiltered/clone_filter_report.stdout 2>&1
done

