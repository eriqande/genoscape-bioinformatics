#!/bin/bash
#$ -cwd
#$ -V
#$ -N merge-bam-subset
#$ -o merge-bam-subset.log
#$ -e merge-bam-subset.error
#$ -pe shared 1
#$ -l h_data=5G,time=10:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java



function usage {
      echo Syntax:
      echo "  $(basename $0)  Prefix  BamPaths"
      echo "
      Prefix:  The prefix of the file that will be created inside MergedBams.
          It is recommended that this give the name of the genome used and also
          some other indication.  For example ZOLAv0-20-birds.
      BamPaths A file containg the absolute paths of the bams you want to merge.
"
}

if [ $# -eq 0 ]; then
    usage;
    exit 1;
fi;

while getopts ":h" opt; do
    case $opt in
	h    ) 
	    usage
	    exit 1
	    ;;
	#m    )  VAR=$OPTARG;
	#    ;;
	\?   )
	    usage
	    exit  1
	    ;;
    esac
done

shift $((OPTIND-1));


# uncomment to test for right number of required args
if [ $# -ne 2 ]; then
    usage;
    exit 1;
fi


PREFIX=$1
BAMPATHS=$2


BAM_COMMS=$(awk '{printf("I=%s ", $0);}' $BAMPATHS)



echo Going to process $BAM_COMMS

if [ -e MergedBams/${PREFIX}-merged.bam ]; then
  echo "File MergedBams/${PREFIX}-merged.bam already exists.  Move it somewhere else.  Exiting... "
  exit 1
fi

if [ ! -d MergedBams ]; then
 mkdir MergedBams
fi

java -Xmx2G -jar $PICARD_JAR MergeSamFiles $BAM_COMMS  OUTPUT=MergedBams/${PREFIX}-merged.bam  SORT_ORDER=coordinate


