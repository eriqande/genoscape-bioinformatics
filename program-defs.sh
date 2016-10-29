

FASTQC=~/bin/FastQC/fastqc
PROC_RADTAGS=~/Documents/others_code/stacks-1.44/process_radtags
CLONE_FILTER=~/Documents/others_code/stacks-1.44/clone_filter


# deal with MacOS's lameness re: zcat
ZCAT=zcat
if [ $(uname -a | awk '{print $1}') == "Darwin" ]; then
  ZCAT=gzcat
fi

