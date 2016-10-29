

FASTQC=/u/nobackup/klohmuel/kruegg/bin/FastQC/fastqc
PROC_RADTAGS=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/process_radtags
CLONE_FILTER=/u/nobackup/klohmuel/kruegg/bin/stacks-1.32/clone_filter


# deal with MacOS's lameness re: zcat
ZCAT=zcat
if [ $(uname -a | awk '{print $1}') == "Darwin" ]; then
  ZCAT=gzcat
fi

