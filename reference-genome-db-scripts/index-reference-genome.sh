#!/bin/bash
#$ -cwd
#$ -V
#$ -N index-fasta
#$ -o index-fasta.log
#$ -e index-fasta.error
#$ -pe shared 1
#$ -l highp,h_data=5G,time=4:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from inside the "Reference" directory which should
# be inside the species directory.  It should be run under qsub.


# put default values here
# VAR=default



function usage {
      echo Syntax:
      echo
      echo "$(basename $0)  FASTA "
      echo
      echo "This script just does samtools faidx on the FASTA file and then
does a picard CreateSequenceDictionary on it.  This assumes that FASTA is named 
with the suffix .fna.  Should change that eventually."
}

if [ $# -eq 0 ]; then
    usage;
    exit 1;
fi;

# uncomment to test for right number of required args
if [ $# -ne 1 ]; then
    usage;
    exit 1;
fi


WorkDir=$(dirname $1)
FASTA=$(basename $1)

source $MODULE_SOURCE
module load samtools
module load java


cd $WorkDir

samtools faidx $FASTA


java -Xmx2G -jar $PICARD_JAR CreateSequenceDictionary \ 
      R=$FASTA \ 
      O=${FASTA/.fna/.dict}
