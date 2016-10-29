#!/bin/bash
#$ -cwd
#$ -V
#$ -N rad
#$ -pe shared 1
#$ -l highp,h_data=8G,time=25:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from the top level of the Plate_X directory
# using qsub.  i.e. qsub ~/genoscape-bioinformatics/plate-processing-scripts/02-rad-process.sh


source ~/genoscape-bioinformatics/program-defs.sh
source data-defs.sh # call this second so you can override program-defs.sh if you want to

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

