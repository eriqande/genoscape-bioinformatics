Step-by-step Mykiss
================
05 January, 2017

-   [Introduction](#introduction)
-   [Build bowtie genome data base](#build-bowtie-genome-data-base)
-   [Mapping](#mapping)
-   [Dup-filtering](#dup-filtering)
-   [Merging the Bams](#merging-the-bams)
-   [Indexing the Genome Fasta and creating a dictionary](#indexing-the-genome-fasta-and-creating-a-dictionary)
-   [Calling SNPs with GATK](#calling-snps-with-gatk)
    -   [An experimental run](#an-experimental-run)
    -   [Running it as a job array](#running-it-as-a-job-array)
    -   [Merging SNPs](#merging-snps)
    -   [Filtering SNPs](#filtering-snps)
-   [Missing data visualize with genoscapeRtools](#missing-data-visualize-with-genoscapertools)

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

Indexing the Genome Fasta and creating a dictionary
---------------------------------------------------

This is just a small step (takes a couple minutes, tops) that needs to happen before using GATK to call SNPs.

``` sh
# login to a compute node
[kruegg@login3 ~]$ qrsh

# go to the genome
[kruegg@n2168 ~]$ cd nobackup-klohmuel/Mykiss/Genome/
[kruegg@n2168 Genome]$ ls
Annotation              bowtie2-omyV6                                 unmapped_scaffoldsV6.fasta
bowtie2-build-db.error  CIGENEomyV6-genes-longestTranscripts.gff3.gz
bowtie2-build-db.log    omyV6Chr.fasta
[kruegg@n2168 Genome]$ du -h omyV6Chr.fasta 
1.9G    omyV6Chr.fasta

# make a .fai with samtools faidx
[kruegg@n2168 Genome]$ module load samtools
[kruegg@n2168 Genome]$ samtools faidx omyV6Chr.fasta 

# make a dictionary with Picard tools
[kruegg@n2168 Genome]$ PICARD_JAR=/u/nobackup/klohmuel/kruegg/bin/picard.jar
[kruegg@n2168 Genome]$ module load java
[kruegg@n2168 Genome]$ java -jar $PICARD_JAR CreateSequenceDictionary R=omyV6Chr.fasta  O=omyV6Chr.dict
```

Calling SNPs with GATK
----------------------

Now we are ready to go. Because we have put the SM tags into each bam I am hoping that GATK will just call each individual from all of the reads it has across all libraries, seamlessly.

### An experimental run

At the moment, I am logged in to a compute node and I am going to try running this for one chromosome to see what it looks like...

``` sh
[kruegg@n2168 SNP-trial-one]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNP-trial-one

[kruegg@n2168 SNP-trial-one]$ module load java

[kruegg@n2168 SNP-trial-one]$ GATK_JAR=/u/nobackup/klohmuel/kruegg/bin/GATK-nightly-11-29-16/GenomeAnalysisTK.jar
[kruegg@n2168 SNP-trial-one]$ FASTA=../../Genome/omyV6Chr.fasta
[kruegg@n2168 SNP-trial-one]$ BAM=../alignments/omyV6/MergedBams/omyV6-merged.bam
[kruegg@n2168 SNP-trial-one]$ SEG=" -L omy01 "
[kruegg@n2168 SNP-trial-one]$ OUT=001

[kruegg@n2168 SNP-trial-one]$ java -Xmx4G -jar $GATK_JAR -T HaplotypeCaller \
 -R $FASTA \
 -I $BAM \
 -stand_call_conf 20.0  \
 -o $OUT.vcf --genotyping_mode DISCOVERY \
  $SEG  > $OUT.stdout 2>$OUT.stderr

# now have a look:
[kruegg@n2168 SNP-trial-one]$ cat 001.stderr
.
.
.
INFO  05:33:22,269 ProgressMeter -     omy01:47142              0.0    30.0 s           50.4 w        0.1%    15.0 h      15.0 h 
INFO  05:34:22,334 ProgressMeter -     omy01:47142              0.0    90.0 s          149.7 w        0.1%    45.0 h      45.0 h 
INFO  05:35:22,364 ProgressMeter -    omy01:149793              0.0     2.5 m          249.0 w        0.2%    23.6 h      23.6 h 
INFO  05:36:22,368 ProgressMeter -    omy01:149998              0.0     3.5 m          348.2 w        0.2%    33.0 h      33.0 h 
INFO  05:37:22,376 ProgressMeter -    omy01:195857              0.0     4.5 m          447.4 w        0.2%    32.5 h      32.4 h 
INFO  05:38:22,379 ProgressMeter -    omy01:195857              0.0     5.5 m          546.6 w        0.2%    39.7 h      39.6 h 
INFO  05:39:22,389 ProgressMeter -    omy01:286963              0.0     6.5 m          645.8 w        0.3%    32.0 h      31.9 h 
INFO  05:40:22,400 ProgressMeter -    omy01:340098              0.0     7.5 m          745.1 w        0.4%    31.2 h      31.1 h 
INFO  05:41:22,408 ProgressMeter -    omy01:407906              0.0     8.5 m          844.3 w        0.5%    29.5 h      29.3 h 
INFO  05:42:22,411 ProgressMeter -    omy01:465489              0.0     9.5 m          943.5 w        0.5%    28.9 h      28.7 h 
INFO  05:43:22,418 ProgressMeter -    omy01:520171              0.0    10.5 m         1042.7 w        0.6%    28.6 h      28.4 h 
```

So, that thing is chugging along. Barring any major troubles it looks like it would finish in between 1 and 2 days. So, I will consider running the each chromosome on 4 processors (-nct 4), and I will continue running it at -Xmx4G. I will run these in the default queues that have a 24 hour time limit and hope that it will work.

### Running it as a job array

**Note! I moved all the results of this section to `/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNPs/chromo_pieces` from `SNPs` **

First, I make a tab-delimited file that holds the chromosomes, so I can do it in a job array:

``` sh
[kruegg@n2168 SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNPs
[kruegg@n2168 SNPs]$ samtools view -H $BAM | awk 'BEGIN {OFS="\t"} /@SQ/ {print ++n, $2}' | sed 's/SN://g' > chromo_list.txt 
[kruegg@n2168 SNPs]$ cat chromo_list.txt 
1   omy01
2   omy02
3   omy03
4   omy04
5   omy05
6   omy06
7   omy07
8   omy08
9   omy09
10  omy10
11  omy11
12  omy12
13  omy13
14  omy14
15  omy15
16  omy16
17  omy17
18  omy18
19  omy19
20  omy20
21  omy21
22  omy22
23  omy23
24  omy24
25  omy25
26  omy26
27  omy27
28  omy28
29  omy29
```

Then we make a script called `mykiss-scripts/05-call-snps-array.sh`, and we call it like:

``` sh
[kruegg@n2168 SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNPs
[kruegg@n2168 SNPs]$ qsub ~/genoscape-bioinformatics/mykiss-scripts/05-call-snps-array.sh
JSV: PE=shared
Your job-array 1436581.1-29:1 ("snp-array") has been submitted
[kruegg@n2168 SNPs]$ date
Wed Jan  4 06:33:40 PST 2017
```

Holy Cow! The first 25 of those started immediately.

Some of these failed with the final statement being:

       WARN  06:34:12,610 PairHMMLikelihoodCalculationEngine$1 - Failed to load native library for VectorLoglessPairHMM - using Java implementation of LOGLESS_CACHING

Aha! Checking the abort messages that get emailed to me, it looks like these guys ended up using too much memory. They crumped out with just over 8Gb. So, I should up each thread to 3 Gb for 12 Gb total and try again... They were 1, 2, 7, 12, 14....so I will restart those.

``` sh
[kruegg@n2168 SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNPs
[kruegg@n2168 SNPs]$ rm 001.* 002.* 007.* 012.* 014.*

[kruegg@n2168 SNPs]$ awk 'BEGIN {OFS="\t"} $1==1 || $1==2 || $1==7 || $1==12 || $1==14 {print ++n, $1, $2}' chromo_list.txt > chromo_redos.txt 
[kruegg@n2168 SNPs]$ cat chromo_redos.txt 
1   1   omy01
2   2   omy02
3   7   omy07
4   12  omy12
5   14  omy14
```

Now, I make another script that picks out the line from the first column of chromo\_redos.txt and then creates the output number (like 001) from the second column, and does the -L statement from the 3rd. Call this \``mykiss-scripts/05-call-snps-array-redo.sh`. Try it:

``` sh
[kruegg@n2168 SNPs]$ qsub ~/genoscape-bioinformatics/mykiss-scripts/05-call-snps-array-redo.sh 
JSV: PE=shared
Your job-array 1436659.1-5:1 ("snp-array") has been submitted
[kruegg@n2168 SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNPs
[kruegg@n2168 SNPs]$ date
Wed Jan  4 07:26:26 PST 2017
```

That seems to be working now...

### Merging SNPs

That job array got done in less than a day. Now we put all the relevant files into a separate directory called `chromo_pieces` and then merge the files with Picard tools.

``` sh
[kruegg@n2008 chromo_pieces]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/Mykiss/Mykiss_all_preps/SNPs/chromo_pieces
[kruegg@n2008 chromo_pieces]$ PICARD_JAR=/u/nobackup/klohmuel/kruegg/bin/picard.jar
[kruegg@n2008 chromo_pieces]$ INPUTS=$(ls -l 0*.vcf | awk '{printf("I=%s ", $NF)}') 
[kruegg@n2008 chromo_pieces]$ module load java
[kruegg@n2008 chromo_pieces]$ java -jar $PICARD_JAR SortVcf $INPUTS O=full-omyV6.vcf
...
[kruegg@n2008 chromo_pieces]$ du -h full-omyV6.vcf
1.2G    full-omyV6.vcf
```

That takes less than two minutes. So, we only have 1.2G of VCF file.
That seems pretty small to me by comparison to the birds we have done. I moved that to: `/u/nobackup/klohmuel/kruegg/Mykiss/Mykiss_all_preps/SNPs/full-omyV6.vcf`.

### Filtering SNPs

Our filtering criteria:

1.  No indels
2.  Biallelic only
3.  Minor allele frequency &gt; 0.01
4.  minimum genotype quality = 30
5.  minimum depth = 8
6.  called in at least 10% of indivs

We only require it be called in 10% of individuals because when I required 50% we got &lt;1000 sites. So, let's do that:

``` sh
[kruegg@n2008 SNPs]$ module load vcftools
[kruegg@n2008 SNPs]$ vcftools --vcf full-omyV6.vcf --out full-omyV6-filtered  --remove-indels --min-alleles 2 --max-alleles 2 --maf 0.01  --minGQ 30 --minDP 8 --max-missing 0.1 --recode
```

This yields 83,314 sites.

No wonder they use ANGSD---their data quality if crapola!

So, let's make an 012 file

    [kruegg@n2008 SNPs]$ vcftools --vcf full-omyV6-filtered.recode.vcf --out full-omyV6-filtered --012

And I will bring that to my laptop:

    /Users/eriq/Documents/UnsyncedData/Mykiss/full-omyV6-filtered.012.gz
    /Users/eriq/Documents/UnsyncedData/Mykiss/full-omyV6-filtered.012.indv
    /Users/eriq/Documents/UnsyncedData/Mykiss/full-omyV6-filtered.012.pos

Missing data visualize with genoscapeRtools
-------------------------------------------

Let's have a look at how we are doing here:

``` r
library(genoscapeRtools)

mykiss <- read_012("~/Documents/UnsyncedData/Mykiss/full-omyV6-filtered", gz = TRUE)

indv <- miss_curves_indv(mykiss)
indv$plot
```

![](mykiss-step-by-step-figs/miss-viz-1.png) Wow! That is serious gargbage!

How about at the pos level?

``` r
loci <- miss_curves_locus(mykiss)
loci$plot
```

![](mykiss-step-by-step-figs/posey-1.png)

Gee, I have never seen data look so bad. If you want to see what things should look like with good data (Kristen's Zosterops data) check it out [here](https://github.com/eriqande/genoscapeRtools#doing-the-missing-data-calcs).

I think I shall have to run through everything using Stacks' `clone_filter` (what we used with Zosterops) instead of using `samtools rmdup` to see if that gives any different results.
