

# when you are getting a header from one of the bams and it has been
# turned into tab-delimited: @SQ\tsName\tLength
# this will assemble it into groups of about 1 million bp, and no bigger
# than 2 million bp, unless it is a single scaffold.


BEGIN {
  OFS = "\t";
  LIM=4e6   # this is the rough limit we want
}


{
  if($3 > LIM) {  # if the scaffold itself is bigger than the limit, just print it separately
    print $3, "-L " $2, $3
    next;
  }
  
  # otherwise we add stuff together.  
  sum +=$3
  comm = sprintf("%s -L %s", comm, $2);
  bps = sprintf("%s %s", bps, $3);
  if(sum > LIM) {
    print sum, comm, bps
    sum = 0
    comm = ""
    bps = ""
  }
}
