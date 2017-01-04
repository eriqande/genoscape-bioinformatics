Step-by-step Mykiss
================
04 January, 2017

-   [Introduction](#introduction)
-   [Build bowtie genome data base](#build-bowtie-genome-data-base)
-   [Mapping](#mapping)
-   [Dup-filtering](#dup-filtering)
-   [Merging the Bams](#merging-the-bams)

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

We have the data from the Prince et al. study. We have stuff that has already been demultiplexed, etc. Also, the names of the samples are not the same as the names of the files. So, we are going to use modified version of the genoscape-bioinformatics scripts. I am going to store those in the repo, but I will put them into `./mykiss_scripts/`

Build bowtie genome data base
-----------------------------

Simple, one script. Gonna call it omyV6. The first time I did it I think I might not have given it enough time, since bowtie2 never completed when aligning against it later. (And, now that I look over my emails from the cluster, that job was killed). So, this time I have given it 24 hours and we will give it a whirl.

``` sh
[kruegg@login1 Genome]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Genome
[kruegg@login1 Genome]$ qsub ~/genoscape-bioinformatics/mykiss-scripts/01-bowtie2-build-genome-database.sh  omyV6Chr.fasta  omyV6
JSV: PE=shared
Your job 1397081 ("bowtie2-build-db") has been submitted
[kruegg@login1 Genome]$ date
Tue Dec 27 10:02:16 PST 2016
```

Once that is done, we will try mapping some small bits first, then let it rip in a big job array.

Aha! Done now. And it only took a couple minutes over 1 hour.

Mapping
-------

This is a little different than before because we are going to do it as a job array and we want to only pick certain files out (because mykiss and chinook are all together in one directory), and we are going to want to name the sample with the DNA name, etc. I have all the IDs and the associated file names in the Excel file that Prince et al posted in the supplement on Biorxiv. I copied the text from the relevant sheet in there and then turned that into a text file with the stuff that I needed in it which we will use for the RG tags.

I did this:

``` sh
pbpaste | tr '\r' '\n' | awk 'BEGIN {print "0 ID PU SM PL LB"} NR > 1 {print ++n,$1, $3"."$4, $NF, "ILLUMINA", $4}' > mykiss_ids.txt 
```

and the first few lines of the file look like:

    0 ID PU SM PL LB
    1 SOMM024_NoIndex_AAGACGTGCAGG SOMM024.NoIndex DNAA004_A04 ILLUMINA NoIndex
    2 SOMM024_NoIndex_AAGCTATGCAGG SOMM024.NoIndex DNAA004_A05 ILLUMINA NoIndex
    3 SOMM024_NoIndex_AATATCTGCAGG SOMM024.NoIndex DNAA004_A06 ILLUMINA NoIndex
    4 SOMM024_NoIndex_AATGAGTGCAGG SOMM024.NoIndex DNAA004_A07 ILLUMINA NoIndex
    5 SOMM024_NoIndex_ACATACTGCAGG SOMM024.NoIndex DNAA001_D08 ILLUMINA NoIndex
    6 SOMM024_NoIndex_AGCGCATGCAGG SOMM024.NoIndex DNAA004_B04 ILLUMINA NoIndex
    7 SOMM024_NoIndex_AGGGTCTGCAGG SOMM024.NoIndex DNAA004_B05 ILLUMINA NoIndex

I am saving that in the repo at `./mykiss-scripts/mykiss_ids.txt`.

Now, I write a job-array script that picks out each line, then grabs the appropriate file and aligns it. That script is called `./mykiss-scripts/02-bowtie-map-job-array.sh`. That seems content to run on up to 20 nodes at 2 threads each at a time. So it should finish in just a few hours. The first time I ran it it turns out I hadn't transferred all the data yet, because Kanaloa had barfed while I was scp-ing. So, I finally got all the data, and then I ran the thing like this:

``` sh
[kruegg@login1 Mykiss_all_preps]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps
[kruegg@login1 Mykiss_all_preps]$ qsub ~/genoscape-bioinformatics/mykiss-scripts/02-bowtie-map-job-array.sh
JSV: PE=shared
Your job-array 1413857.1-242:1 ("radMap") has been submitted
[kruegg@login1 Mykiss_all_preps]$ date
Fri Dec 30 04:56:10 PST 2016
```

Holy Smokes! There must be very little load on hoffman at the moment---all 242 jobs started straight away! (And the whole thing took about half an hour, as far as I can tell.)

Dup-filtering
-------------

I am going to do this as a job array as well, using `samtools rmdup`. There is a little setup for that:

``` sh
[kruegg@n2194 bam]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/alignments/omyV6/bam
# break it into 25 chunks of 10 files each
[kruegg@n2194 bam]$ ls -l  *.bam | awk 'NR % 10 == 1 {printf("\n%d ",++n)} {printf(" %s", $NF);} END {printf("\n");}' > bamlist.txt 
[kruegg@n2194 bam]$ head bamlist.txt 

1  SOMM024_NoIndex_AAGACGTGCAGG.bam SOMM024_NoIndex_AAGCTATGCAGG.bam SOMM024_NoIndex_AATATCTGCAGG.bam SOMM024_NoIndex_AATGAGTGCAGG.bam SOMM024_NoIndex_ACATACTGCAGG.bam SOMM024_NoIndex_AGCGCATGCAGG.bam SOMM024_NoIndex_AGGGTCTGCAGG.bam SOMM024_NoIndex_AGGTGTTGCAGG.bam SOMM024_NoIndex_AGTAGGTGCAGG.bam SOMM024_NoIndex_ATCAAATGCAGG.bam
2  SOMM024_NoIndex_CAATCGTGCAGG.bam SOMM024_NoIndex_CACCTCTGCAGG.bam SOMM024_NoIndex_CAGGCATGCAGG.bam SOMM024_NoIndex_CATACTTGCAGG.bam SOMM024_NoIndex_CCATTTTGCAGG.bam SOMM024_NoIndex_CCGAGGTGCAGG.bam SOMM024_NoIndex_CGCGTGTGCAGG.bam SOMM024_NoIndex_CGGTCCTGCAGG.bam SOMM024_NoIndex_CGTCTATGCAGG.bam SOMM024_NoIndex_CGTGATTGCAGG.bam
3  SOMM024_NoIndex_CTACAGTGCAGG.bam SOMM024_NoIndex_CTGGTTTGCAGG.bam SOMM024_NoIndex_GACGACTGCAGG.bam SOMM024_NoIndex_GACTCTTGCAGG.bam SOMM024_NoIndex_GAGAGATGCAGG.bam SOMM024_NoIndex_GATCGTTGCAGG.bam SOMM024_NoIndex_GCAGATTGCAGG.bam SOMM024_NoIndex_GGGCGCTGCAGG.bam SOMM024_NoIndex_GGGGCGTGCAGG.bam SOMM024_NoIndex_GGTACATGCAGG.bam
4  SOMM024_NoIndex_GGTTTGTGCAGG.bam SOMM024_NoIndex_TACACATGCAGG.bam SOMM024_NoIndex_TACGGGTGCAGG.bam SOMM024_NoIndex_TAGTATTGCAGG.bam SOMM024_NoIndex_TATCACTGCAGG.bam SOMM024_NoIndex_TCGATTTGCAGG.bam SOMM024_NoIndex_TGACAATGCAGG.bam SOMM024_NoIndex_TGCCCGTGCAGG.bam SOMM024_NoIndex_TGCTTATGCAGG.bam SOMM024_NoIndex_TGGGGATGCAGG.bam
5  SOMM024_NoIndex_TTCTAGTGCAGG.bam SOMM025_NoIndex_ACCATGTGCAGG.bam SOMM025_NoIndex_ACCCCCTGCAGG.bam SOMM025_NoIndex_ATGCACTGCAGG.bam SOMM025_NoIndex_ATGTTGTGCAGG.bam SOMM025_NoIndex_CCGCATTGCAGG.bam SOMM025_NoIndex_CCTAACTGCAGG.bam SOMM025_NoIndex_CTTATGTGCAGG.bam SOMM025_NoIndex_CTTTGCTGCAGG.bam SOMM025_NoIndex_GCGCTGTGCAGG.bam
6  SOMM025_NoIndex_GCTCAATGCAGG.bam SOMM025_NoIndex_GTGTAATGCAGG.bam SOMM025_NoIndex_GTTGGATGCAGG.bam SOMM025_NoIndex_TCGGACTGCAGG.bam SOMM025_NoIndex_TCTCGGTGCAGG.bam SOMM025_NoIndex_TTTAATTGCAGG.bam SOMM025_NoIndex_TTTGTCTGCAGG.bam SOMM049_Index02_AATATCTGCAGG.bam SOMM049_Index02_AATGAGTGCAGG.bam SOMM049_Index02_AGGTGTTGCAGG.bam
7  SOMM049_Index02_AGTAGGTGCAGG.bam SOMM049_Index02_CATACTTGCAGG.bam SOMM049_Index02_CCATTTTGCAGG.bam SOMM049_Index02_CGTCTATGCAGG.bam SOMM049_Index02_CGTGATTGCAGG.bam SOMM049_Index02_CTACAGTGCAGG.bam SOMM049_Index02_GAGAGATGCAGG.bam SOMM049_Index02_GATCGTTGCAGG.bam SOMM049_Index02_GCAGATTGCAGG.bam SOMM049_Index02_GGTACATGCAGG.bam
8  SOMM049_Index02_GGTTTGTGCAGG.bam SOMM049_Index02_TAGTATTGCAGG.bam SOMM049_Index02_TATCACTGCAGG.bam SOMM049_Index02_TGCTTATGCAGG.bam SOMM049_Index02_TGGGGATGCAGG.bam SOMM050_Index03_AAACGGTGCAGG.bam SOMM050_Index03_AACGTTTGCAGG.bam SOMM050_Index03_AACTGATGCAGG.bam SOMM050_Index03_AAGACGTGCAGG.bam SOMM050_Index03_AAGCTATGCAGG.bam
9  SOMM050_Index03_AATATCTGCAGG.bam SOMM050_Index03_AATGAGTGCAGG.bam SOMM050_Index03_ACAAGATGCAGG.bam SOMM050_Index03_ACAGCGTGCAGG.bam SOMM050_Index03_ACATACTGCAGG.bam SOMM050_Index03_ACCATGTGCAGG.bam SOMM050_Index03_ACCCCCTGCAGG.bam SOMM050_Index03_ACTCTTTGCAGG.bam SOMM050_Index03_ACTGGCTGCAGG.bam SOMM050_Index03_AGCCATTGCAGG.bam
```

And with that done we can launch the job:

``` sh
[kruegg@n2194 bam]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/alignments/omyV6/bam
[kruegg@n2194 bam]$ qsub ~/genoscape-bioinformatics/mykiss-scripts/03-rmdup-alignments-job-array.sh 
JSV: PE=shared
Your job-array 1433016.1-25:1 ("rmdup") has been submitted
[kruegg@n2194 bam]$ myjobs
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
1432822 0.00208 QRLOGIN    kruegg       r     01/03/2017 13:28:38 inter_msa.q@n2194                  1        
1433016 0.00000 rmdup      kruegg       qw    01/03/2017 14:16:22                                    1 1-25:1
```

That takes about 10 minutes.

Merging the Bams
----------------

We need to merge the deduped bams to prepare for GATK SNP calling. I made a script to to this and I launched it like:

``` sh
[kruegg@n2194 omyV6]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/alignments/omyV6
[kruegg@n2194 omyV6]$ qsub ~/genoscape-bioinformatics/mykiss-scripts/04-merge-bams.sh 
JSV: PE=shared
Your job 1433096 ("merge-bams") has been submitted
```