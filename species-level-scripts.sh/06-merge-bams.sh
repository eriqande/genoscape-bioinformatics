#!/bin/bash
#$ -cwd
#$ -V
#$ -N merge-bams
#$ -o merge-bams.log
#$ -e merge-bams.error
#$ -pe shared 1
#$ -l h_data=4G,time=10:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java



function usage {
      echo Syntax:
      echo "  $(basename $0)  AlignDir  PlateOne PlateTwo PlateThree ..."
      echo "
      Alignment:  The name of the alignment directories inside each of the Plate 
          directories, inside of which is are the "bam" directories that hold all 
          the bam files. For example, ZOLAv0.
      PlateOne ... The names of the directories holding all the aligned bams from 
          each of the plates.  For example, if you have done three plates it might
          be \"Plate_1 Plate_2 PLate_3\"
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
#if [ $# -ne 4 ]; then
#    usage;
#    exit 1;
#fi



GENOME_DIR=$1
shift;

BAM_COMMS=$(while (($#)); do
    VAR=$1;
    shift;
    printf 'I=%s ' $VAR/alignments/$GENOME_DIR/bam/*.bam
done)



echo Going to process $BAM_COMMS

if [ -e MergedBams/${GENOME_DIR}-merged.bam ]; then
  echo "File MergedBams/${GENOME_DIR}-merged.bam already exists.  Move it somewhere else.  Exiting... "
  exit 1
fi

if [ ! -d MergedBams ]; then
 mkdir MergedBams
fi

java -Xmx2G -jar $PICARD_JAR MergeSamFiles $BAM_COMMS  OUTPUT=MergedBams/${GENOME_DIR}-merged.bam  SORT_ORDER=coordinate
