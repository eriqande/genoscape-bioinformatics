

FASTQC=/u/nobackup/klohmuel/kruegg/bin/FastQC/fastqc
PROC_RADTAGS=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/process_radtags
CLONE_FILTER=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/clone_filter

# this can be obtained by, for example: 
#  wget https://github.com/broadinstitute/picard/releases/download/2.7.1/picard.jar
PICARD_JAR=/u/home/k/kruegg/nobackup-klohmuel/bin/picard.jar


# deal with MacOS's lameness re: zcat.  Users should not have to 
# alter anything here, but this must still be included!
ZCAT=zcat
if [ $(uname -a | awk '{print $1}') == "Darwin" ]; then
  ZCAT=gzcat
fi

