#!/bin/bash
#$ -cwd
#$ -V
#$ -N bowtie2-build-db
#$ -o bowtie2-build-db.log
#$ -e bowtie2-build-db.error
#$ -pe shared 1
#$ -l highp,h_data=8G,time=0:50:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from inside the "Reference" directory which should
# be inside the species directory.  It should be run under qsub.


# put default values here
# VAR=default



function usage {
      echo Syntax:
      echo
      echo "$(basename $0)  FASTA  OutputPrefix"
      echo
      echo "This script just does a simple bowtie-build from a genome is fasta format (.fa or
.fna).  It apparently cannot be gzipped, so make sure that it is not compressed. 
The script  puts the output in a in new direclty within the current working directory
bowtie2-OutputPrefix where OutputPrefix is a name passed in by the user.
The inner contents all have the name OutputPrefix with whatever extensions bowtie2 puts on it.
Inside that directory, you will also find a README.txt file that provides a log of when this
was done and what command produced it.  "
}

if [ $# -eq 0 ]; then
    usage;
    exit 1;
fi;

# uncomment to test for right number of required args
if [ $# -ne 2 ]; then
    usage;
    exit 1;
fi

FASTA=$1
OutputPrefix=$2


. /u/local/Modules/default/init/modules.sh
module load bowtie2


FASTA_ABS=$(readlink -f $FASTA)

mkdir bowtie2-$OutputPrefix
cd bowtie2-$OutputPrefix


echo "Created $(date),  by running  bowtie2-build-genome-database.sh script with arguments:
        FASTA: $FASTA
        AbsolutePath: $FASTA_ABS
        OutputPrefix: $OutputPrefix" > README.txt 


bowtie2-build $FASTA_ABS  $OutputPrefix > bowtie2-build.stdout  2> bowtie2-build.stderr



cd ../

