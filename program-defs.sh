
# this is the path to the file to source to get the module command
MODULE_SOURCE=/u/local/Modules/default/init/modules.sh

FASTQC=/u/nobackup/klohmuel/kruegg/bin/FastQC/fastqc
PROC_RADTAGS=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/process_radtags
CLONE_FILTER=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/clone_filter

# this can be obtained by, for example: 
#  wget https://github.com/broadinstitute/picard/releases/download/2.7.1/picard.jar
PICARD_JAR=/u/nobackup/klohmuel/kruegg/bin/picard.jar

# this can be obtained by downloading via a browser, then scp-ing it
# Old version 3.6
# GATK_JAR=/u/nobackup/klohmuel/kruegg/bin/GenomeAnalysisTK.jar

# This is the 3.8 version
GATK_JAR=/u/nobackup/klohmuel/kruegg/bin/GenomeAnalysisTK-3.8-0-ge9d806836/GenomeAnalysisTK.jar

NUCMER=/u/nobackup/klohmuel/kruegg/bin/MUMmer3.23/nucmer




# deal with MacOS's lameness re: zcat.  Users should not have to 
# alter anything here, but this must still be included!
ZCAT=zcat
if [ $(uname -a | awk '{print $1}') == "Darwin" ]; then
  ZCAT=gzcat
fi

