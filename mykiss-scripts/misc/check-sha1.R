
# this should be run in the top dir of the repo.

# Kanaloa barfed out one day while I was transferring RAD seq files to 
# Hoffman2 and I am not sure what got corrupted.  So, I shasum-ed all the
# files and will check.

library(tidyverse)
library(readr)


hoff <- read_delim("./mykiss-scripts/other_inputs/shasum/hoffman_rad_seq_shasums.txt", delim = "  *")[,c(1,3)]
kana <- read_delim("./mykiss-scripts/other_inputs/shasum/kanaloa_rad_seq_shasums.txt", delim = "  *")[,c(1,3)]


full <- full_join(kana, hoff)

# note that there are no "extra" files on hoffman
full %>%
  filter(is.na(kana_hash))



# here is the one problem file
full %>%
  filter(kana_hash != hoff_hash)

# just copied that over on its own.

# and here are the ones that just didn't get there:
never_arrived <- full %>%
  filter(is.na(hoff_hash))

# i will put those file names into a file to put on a command line:
cat(never_arrived$file, sep = " ", file="/tmp/nev.txt")

# i put that on kanaloa and then did:
# scp $(cat /tmp/nev.txt)  kruegg@hoffman2.idre.ucla.edu:/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Prince_etal_raw/RAD_sequence/
  
  


