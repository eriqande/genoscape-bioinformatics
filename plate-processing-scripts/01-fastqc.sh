#!/bin/bash
#$ -cwd
#$ -V
#$ -N fastqc
#$ -o fastqc.log
#$ -e fastqc.error
#$ -pe shared 1
#$ -l highp,h_data=2G,time=1:20:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

# This should be run from the top level of the Plate_X directory
# using qsub.  i.e. qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh




source ~/genoscape-bioinformatics/program-defs.sh
source data-defs.sh # call this second so you can override program-defs.sh if you want to


# make necessary output directories
mkdir fastqc_reports
mkdir fastqc_reports/logs

# Do read 1
zcat $READ1 | $FASTQC stdin

# move output to output directory while renaming it
mv fastqc.error fastqc_reports/logs/$(basename $READ1)_fastqc.error
mv fastqc.log fastqc_reports/logs/$(basename $READ1)_fastqc.log
mv stdin_fastqc.html fastqc_reports/fastqc_report_$(basename $READ1).html
rm stdin_fastqc.zip

# Do read 2
zcat $READ2 | $FASTQC stdin

# move output to output directory while renaming it
mv fastqc.error fastqc_reports/logs/$(basename $READ2)_fastqc.error
mv fastqc.log fastqc_reports/logs/$(basename $READ2)_fastqc.log
mv stdin_fastqc.html fastqc_reports/fastqc_report_$(basename $READ2).html
rm stdin_fastqc.zip


