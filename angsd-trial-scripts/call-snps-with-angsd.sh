#!/bin/bash
#$ -cwd
#$ -V
#$ -N angsd-wifl
#$ -o angsd-wifl.log
#$ -e angsd-wifl.error
#$ -pe shared 8
#$ -l time=24:00:00
#$ -M eric.anderson@noaa.gov
#$ -m a

# run this in the directory you want the subsampled bams to go.

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load zlib/1.2.8



~/nobackup-klohmuel/bin/ngsTools/angsd/angsd  -P 8 -b bamlist_no_dupes.txt  -out Results/ALL  \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  \
        -minMapQ 10 -minQ 20 -minInd 108 -doCounts 1 \
        -GL 1 -doMajorMinor 1 -doMaf 2 -skipTriallelic 1  -minMaf 0.05 \
        -SNP_pval 1e-6\
        -doGeno 32 -doPost 2
