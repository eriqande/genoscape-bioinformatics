#!/bin/bash
#$ -cwd
#$ -V
#$ -N bowtie-build-db
#$ -o bowtie-build-db.log
#$ -e bowtie-build-db.error
#$ -pe shared 1
#$ -l highp,h_data=8G,time=1:20:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from inside the "Reference" directory which should
# be inside the species directory.  It should be run under qsub.


# put default values here
# VAR=default



function usage {
      echo Syntax:
      echo
      echo "$(basename $0)  FASTA_OR_FASTA_GZ  OutputPrefix"
      echo
      echo "This script just does a simple bowtie-build from a genome is fasta format (or 
gzipped, like fna.gz).  It puts the output in the current working directory in a new
directory called bowtie-OutputPrefix where OutputPrefix is a name passed in by the user.
The inner contents all have the name OutputPrefix with whatever extensions bowtie puts on it.
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

FASTA_OR_FASTA_GZ=$1
OutputPrefix=$2



module load bowtie

mkfifo infasta.fa

if [ ${FASTA_OR_FASTA_GZ: -3} == ".gz" ]; then
  zcat  FASTA_OR_FASTA_GZ > infasta.fa &
else 
  cat  FASTA_OR_FASTA_GZ > infasta.fa &
fi

mkdir bowtie-$OutputPrefix
cd bowtie-$OutputPrefix

echo "Created by running  bowtie-build-genome-database.sh script with arguments:
        FASTA_OR_FASTA_GZ: $FASTA_OR_FASTA_GZ
        OutputPrefix: $OutputPrefix" > README.txt 


bowtie-build ../infasta.fa $GENOME_OUTNAME



cd ../
rm infasta.fa

