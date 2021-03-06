#!/bin/bash
#$ -cwd
#$ -V
#$ -N demultiplex
#$ -o demultiplex.log
#$ -e demultiplex.error
#$ -pe shared 1
#$ -l highp,h_data=8G,time=80:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# Note.  It should take well under 30 hours for this to run on Hoffman, but 
# the trogolodytes that run the cluster kill your process if it runs over at all,
# so I just give it way more time than it might need.  Takes longer to queue, but that
# is better than completely losing a days worth of computation.

# note that the -o and -e directives above seem important to keep the script
# from trying to write to a default output location that it is unable to write to

# This should be run from the top level of the Plate_X directory
# using qsub.  i.e. qsub ~/genoscape-bioinformatics/plate-processing-scripts/02-rad-process.sh


source data-defs.sh 
source $PROGDEFS



# prepare flipped directory and output file names
mkdir flipped
OUTPATH1=${SHORTNAME}_r1
OUTPATH2=${SHORTNAME}_r2


# make a single barcodes file 
awk '{print $1}' $BC2COL > xxx_one-column-barcodes.txt

# do the flip script with FIFO's so we don't need to store any decompressed stuff around
mkfifo in1 in2
mkfifo out1 out2

$ZCAT  $READ1 > in1 &
$ZCAT  $READ2 > in2 &

cat out1 | gzip - > flipped/$OUTPATH1.fastq.gz &
cat out2 | gzip - > flipped/$OUTPATH2.fastq.gz &

###For bestRAD only: flip read pairs that are reversed and trim extra "GG"
perl ~/genoscape-bioinformatics/script/flip_trim_werrors.pl xxx_one-column-barcodes.txt in1 in2   out1 out2    true 2 > flipped/flip_script_report.txt

rm in1 in2 out1 out2 xxx_one-column-barcodes.txt


###This calls the newer g++ compiler, necessary for running the newer stacks version on Hoffman2
#source ~/scripts/mypath.sh

###Make directory for demultiplexed fastq files
mkdir demultiplexed

###Run stacks to demultiplex and trim adapter sequences
$PROC_RADTAGS -i gzfastq -1 flipped/$OUTPATH1.fastq.gz -2 flipped/$OUTPATH2.fastq.gz \
	-o ./demultiplexed -b $BC2COL -c -q -r -e sbfI --filter_illumina \
	--adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTC --adapter_2 CACTCTTTCCCTACACGACGCTCTTCCGATCT > demultiplexed/process_radtags.stdout 2>&1

