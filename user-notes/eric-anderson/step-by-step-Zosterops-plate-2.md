Step-by-step *Zosterops* Plate 2
================
18 November, 2016

-   [Introduction](#introduction)
-   [Download the data and move it where it needs to go (~ 1 hr)](#download-the-data-and-move-it-where-it-needs-to-go-1-hr)
    -   [Download with wget](#download-with-wget)
    -   [Move data to the proper directory](#move-data-to-the-proper-directory)
    -   [Remove the temporary download directory](#remove-the-temporary-download-directory)
    -   [A brief digression](#a-brief-digression)
-   [Set up the `data-defs.sh` and barcodes file](#set-up-the-data-defs.sh-and-barcodes-file)
-   [Run `fastqc` on it (~ 2 hr)](#run-fastqc-on-it-2-hr)
    -   [First run fastqc using the script provided](#first-run-fastqc-using-the-script-provided)
    -   [Copy the results back to your own machine to view](#copy-the-results-back-to-your-own-machine-to-view)
-   [Flip, Trim, and demultiplex, for the Best-RAD process (~ 1 day)](#flip-trim-and-demultiplex-for-the-best-rad-process-1-day)
    -   [Start the job](#start-the-job)
-   [Filter out PCR duplicates](#filter-out-pcr-duplicates)
    -   [Start the job](#start-the-job-1)
-   [Clean up stuff](#clean-up-stuff)
-   [Prepare and index the genome](#prepare-and-index-the-genome)
-   [Map the reads to the genome data base](#map-the-reads-to-the-genome-data-base)
-   [Merge the bams from Plate\_1 and Plate\_2](#merge-the-bams-from-plate_1-and-plate_2)
-   [Merge BAMs from a subset of individuals](#merge-bams-from-a-subset-of-individuals)
-   [Index the genome fasta and create a dictionary](#index-the-genome-fasta-and-create-a-dictionary)
-   [Small trial of SNPing](#small-trial-of-snping)
-   [Going for the full SNP calling at 20 individuals](#going-for-the-full-snp-calling-at-20-individuals)
-   [Trying to do SNP calling as a job array](#trying-to-do-snp-calling-as-a-job-array)
    -   [First, try running a scaffold that is about 100K bp](#first-try-running-a-scaffold-that-is-about-100k-bp)
    -   [Submitting as a job array](#submitting-as-a-job-array)
-   [Initial filtering of the VCF](#initial-filtering-of-the-vcf)
    -   [Getting a list of 1 Mb collections of scaffolds](#getting-a-list-of-1-mb-collections-of-scaffolds)

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

Here I will try to faithfully record the whole process I went through to analyse the *Zosterops* RAD data using the scripts in in [genoscape-bioinformatics](https://github.com/eriqande/genoscape-bioinformatics).

Let's see how it goes. For each step, I list the "prep" time, which is how long things take (if there are no problems) to actually interact with the terminal. The "cook" time is how long the computation takes once the job is submitted.

Download the data and move it where it needs to go (~ 1 hr)
-----------------------------------------------------------

*Prep Time: 10 minutes*
*Cook Time: 60 minutes*

### Download with wget

UC Davis sent an email with some links. I followed them and logged in an navigated my way to the directory that held the data I wanted in my web browser. The top of the page showed me that the path to these data is `/Data/873yxr7yvh/Unaligned/Project_TSKR_L3_Silvereye_Plate2`. So, that is what I am going to want to use with `wget` to download it to hoffman.

So, I log in to Hoffman and create a temporary directory on `nobackup-klohmuel` and then use wget with the domain and the path like this:

``` sh
# do this to show what directory I am in:
[kruegg@login4 tempDownload]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/tempDownload

# then start the download.
wget -r -nH -np -R index.html* "http://slims.bioinformatics.ucdavis.edu/Data/873yxr7yvh/Unaligned/Project_TSKR_L3_Silvereye_Plate2" &
```

Note that the `-R index.html*` part of the command just tells `wget` to not bother keeping all the little files that are there to allow the web-browing part of the filesystem. This took about 20 minutes for each of the two large files of data.

### Move data to the proper directory

I already have a directory `ZOLA` at `/u/home/k/kruegg/nobackup-klohmuel/ZOLA` for *Zosterops* work. This is the second plate, so I make the directory `Plate_2` inside that directory, and then I copy the relevant downloaded files to a `rawdata` directory inside that.

``` sh
# moving files.
[kruegg@login4 Project_TSKR_L3_Silvereye_Plate2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/tempDownload/Data/873yxr7yvh/Unaligned/Project_TSKR_L3_Silvereye_Plate2
[kruegg@login4 Project_TSKR_L3_Silvereye_Plate2]$ ls
laneBarcode.html  Silvereye-Plate2_S86_L003_R1_001.fastq.gz  Silvereye-Plate2_S86_L003_R2_001.fastq.gz
[kruegg@login4 Project_TSKR_L3_Silvereye_Plate2]$ mv * ~/nobackup-klohmuel/ZOLA/Plate_2/rawdata/

# once they are moved, here is what their destination looks like:
[kruegg@login1 rawdata]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2/rawdata
[kruegg@login1 rawdata]$ ls
laneBarcode.html  Silvereye-Plate2_S86_L003_R1_001.fastq.gz  Silvereye-Plate2_S86_L003_R2_001.fastq.gz
```

Note that I kept the `laneBarcode.html` file 'cuz it might have some information that we will want to look at some day. But it is not necessary for any analyses. And it is **not** related to the barcodes file we will use later in demultiplexing these data.

### Remove the temporary download directory

Once our data our safely in their final destination, we can delete the temporary download directory `/u/home/k/kruegg/nobackup-klohmuel/tempDownload`:

``` sh
[kruegg@login4 nobackup-klohmuel]$ pwd
/u/home/k/kruegg/nobackup-klohmuel
[kruegg@login4 nobackup-klohmuel]$ rm -rf tempDownload
```

### A brief digression

Here is a little digression (not something you have to do every time).
I ran this plate before and ran into some errors. I want to check if the files I just downloaded are exactly the same bit-for-bit as the ones I analyzed prevously. To do that, I compute the SHA-1 hash fingerprint for each of the original files. Here they are:

``` sh
[kruegg@login1 rawdata]$ sha1sum Silvereye-Plate2_S86_L003_R1_001.fastq.gz 
35696d76e18aa728320664e738cb5d6463e0021a  Silvereye-Plate2_S86_L003_R1_001.fastq.gz
[kruegg@login1 rawdata]$ sha1sum Silvereye-Plate2_S86_L003_R2_001.fastq.gz 
fb9c8cc7ff5d531f25e0469d8fb5c22d35471f54  Silvereye-Plate2_S86_L003_R2_001.fastq.gz
```

And, the new ones that I just downloaded look like this:

``` sh
[kruegg@login1 rawdata]$ sha1sum Silvereye-Plate2_S86_L003_R*
35696d76e18aa728320664e738cb5d6463e0021a  Silvereye-Plate2_S86_L003_R1_001.fastq.gz
fb9c8cc7ff5d531f25e0469d8fb5c22d35471f54  Silvereye-Plate2_S86_L003_R2_001.fastq.gz
```

So, it is exactly identical to what it was last time...

Set up the `data-defs.sh` and barcodes file
-------------------------------------------

*Prep Time: 25 minutes. Take time to double check everything and make sure that you have it all right. Check with the [directory structure specifications](https://github.com/eriqande/genoscape-bioinformatics#directory-structure) and the [progdefs specifications](https://github.com/eriqande/genoscape-bioinformatics#installation) to make sure that you are using the most current specs.*
*Cook Time: 0 minutes*

We are going to be working in `/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2` and that is where we need to make our file called `data-defs.sh`. This is what it looks like:

``` sh
PROGDEFS=/u/home/k/kruegg/genoscape-bioinformatics/program-defs.sh  # file that defines the paths
READ1=rawdata/Silvereye-Plate2_S86_L003_R1_001.fastq.gz  # fastq.gz file for Read one.
READ2=rawdata/Silvereye-Plate2_S86_L003_R2_001.fastq.gz # paths relative to the Plate_1 directory.
BC2COL=ZOLA-2-barcodes-2col.txt # The two column barcode file
SHORTNAME=zola-2   # species dash plate number is good.  This is used for output later  
GENOME_DB=/u/home/k/kruegg/nobackup-klohmuel/References/bowtie2-ZOLAv0/ZOLAv0   # must be the ABSOLUTE path AND the prefix
```

We see from the above that our barcodes file will be named `ZOLA-2-barcodes-2col.txt`. We make sure that is in the `Plate_2` directory. The first 10 of its 96 lines look like this:

    ACAAGCTA    SE_115
    AAACATCG    12-067
    ACATTGGC    D3
    ACCACTGT    SE_125
    AACGTGAT    SE_126
    CGCTGATC    CI13
    CAGATCTG    D20
    ATGCCTAA    SE_128
    AACGAACG    SE_118
    AGTACAAG    L536

Once I have everything in place, my files look like this:

``` sh
[kruegg@login4 Plate_2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2
[kruegg@login4 Plate_2]$ ls *
data-defs.sh  ZOLA-2-barcodes-2col.txt

rawdata:
laneBarcode.html  Silvereye-Plate2_S86_L003_R1_001.fastq.gz  Silvereye-Plate2_S86_L003_R2_001.fastq.gz
```

Run `fastqc` on it (~ 2 hr)
---------------------------

*Prep Time: 20 minutes (includes scanning the reports)* *Cook Time: 2 hours*

### First run fastqc using the script provided

This is pretty straightforward. Make sure you are in `Plate_2` and then run the script under `qsub`:

``` sh
[kruegg@login4 Plate_2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2
[kruegg@login4 Plate_2]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh 
JSV: PE=shared
Your job 1037389 ("fastqc") has been submitted
```

Note, when your job is queued, you can check on it with `myjobs`:

``` sh
[kruegg@login4 Plate_2]$ myjobs
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
1037389 0.00000 fastqc     kruegg       qw    11/06/2016 14:18:21                                    1        
```

The above shows that it is in stat "qw" which means it is queued and waiting. So we try it again later:

``` sh
[kruegg@login4 Plate_2]$ myjobs
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
1037389 0.00208 fastqc     kruegg       r     11/06/2016 14:20:26 eeskin_pod_16.q@n6130              1        
```

Now it is in state "r" for running.

### Copy the results back to your own machine to view

Once that has finished (in about 2 hours) the results are html files which are most easily viewed by copying them to your own computer and opening them in a browser.

``` sh
2016-11-07 02:13 /Desktop/--% scp -r  kruegg@hoffman2.idre.ucla.edu:/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2/fastqc_reports ./
kruegg@hoffman2.idre.ucla.edu's password: 
Silvereye-Plate2_S86_L003_R1_001.fastq.gz_fas 100%  166     0.2KB/s   00:00    
Silvereye-Plate2_S86_L003_R1_001.fastq.gz_fas 100%   56     0.1KB/s   00:01    
fastqc_report_Silvereye-Plate2_S86_L003_R1_00 100%  309KB 308.9KB/s   00:01    
fastqc_report_Silvereye-Plate2_S86_L003_R2_00 100%  305KB 305.4KB/s   00:01    
```

Then just double click them to open the reports in your browser.

Flip, Trim, and demultiplex, for the Best-RAD process (~ 1 day)
---------------------------------------------------------------

*Prep Time: 5 minutes*
*Cook Time: About 24 hours*

### Start the job

Just cd to the `Plate_2` directory and launch it.

``` sh
[kruegg@login4 Plate_2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2
[kruegg@login4 Plate_2]$ date
Mon Nov  7 02:34:24 PST 2016
[kruegg@login4 Plate_2]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/02-rad-flip-trim-and-demultiplex.sh 
JSV: PE=shared
Your job 1047166 ("demultiplex") has been submitted

# check to make sure it is in the queue
[kruegg@login4 Plate_2]$ myjobs
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
1047166 0.00000 demultiple kruegg       qw    11/07/2016 02:34:40                                    1        
```

And, then about a day later it is done...

Filter out PCR duplicates
-------------------------

*Prep Time: 2 minutes*
*Cook Time: 9 hours*

This is a quick script that relies on STACKS clone filter.

### Start the job

Just cd to the `Plate_2` directory and launch it

``` sh
[kruegg@login1 Plate_2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2
[kruegg@login1 Plate_2]$ date
Wed Nov  9 03:32:20 PST 2016
[kruegg@login1 Plate_2]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/03-rad-clone-filter.sh 
JSV: PE=shared
Your job 1067345 ("clone_filter") has been submitted
```

Clean up stuff
--------------

Just a quick script to move and rename a few things. Not compute intensive so it does not need to be run under `qusb`

``` sh
[kruegg@login1 Plate_2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2
[kruegg@login1 Plate_2]$ ~/genoscape-bioinformatics/plate-processing-scripts/04-clean-up-and-organize.sh 
```

Prepare and index the genome
----------------------------

I will discuss this later, as I already did it a while ago

Map the reads to the genome data base
-------------------------------------

*Prep Time: 5 minutes*
*Cook Time: 8 hours wallclock on Hoffman. 2 days of CPU time*

This just involves starting up the script:

``` sh
[kruegg@login1 Plate_2]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_2
[kruegg@login1 Plate_2]$ date
Wed Nov  9 14:30:17 PST 2016
[kruegg@login1 Plate_2]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/05-map-reads-to-genome-db.sh 
JSV: PE=shared
Your job 1075046 ("radMap") has been submitted
[kruegg@login1 Plate_2]$ 
[kruegg@login1 Plate_2]$ 
[kruegg@login1 Plate_2]$ myjobs
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
1075046 0.00000 radMap     kruegg       qw    11/09/2016 14:30:34                                    8        
```

Merge the bams from Plate\_1 and Plate\_2
-----------------------------------------

This is a "species-level" operation...

``` sh
# see that I am at the top level of the ZOLA directory
[kruegg@login2 ZOLA]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA

# here are the directories in here:
[kruegg@login2 ZOLA]$ ls
Genome  Plate_1  Plate_2

# note that each plate has an alignments directory with ZOLAv0 in it
[kruegg@login2 ZOLA]$ ls Plate_1/alignments Plate_2/alignments
Plate_1/alignments:
ZOLAv0

Plate_2/alignments:
ZOLAv0

# here is that help for 06:
[kruegg@login2 ZOLA]$~/Documents/git-repos/genoscape-bioinformatics/species-level-scripts.sh/06-merge-bams.sh -h 

Syntax:
  06-merge-bams.sh  AlignDir  PlateOne PlateTwo PlateThree ...

      Alignment:  The name of the alignment directories inside each of the Plate 
          directories, inside of which is are the bam directories that hold all 
          the bam files. For example, ZOLAv0.
      PlateOne ... The names of the directories holding all the aligned bams from 
          each of the plates.  For example, if you have done three plates it might
          be "Plate_1 Plate_2 PLate_3"


[kruegg@login2 ZOLA]$ date
Fri Nov 11 12:30:05 PST 2016

# let's give it a whirl
[kruegg@login4 ZOLA]$ qsub ~/genoscape-bioinformatics/species-level-scripts.sh/06-merge-bams.sh ZOLAv0 Plate_1 Plate_2
JSV: PE=shared
Your job 1090852 ("merge-bams") has been submitted
```

Merge BAMs from a subset of individuals
---------------------------------------

We are going to take 10 birds from the mainland and 10 from Heron Island and merge them to do some quick SNP calling as quickly as we can.

Here are the full paths of the bams we want to merge. The CDH's are mainland birds and the others are Heron Island. I took the 10 that had the most seqence (largest bams).

    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.CDH46.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH18.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH50.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH48.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH22.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH20.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH23.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH24.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH16.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.CDH49.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.12-052.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.12-064.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.12-068.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_1/alignments/ZOLAv0/bam/zola-1.12-061.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.12-066.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.12-053.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.12-051.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.12-049.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.12-060.bam
    /u/nobackup/klohmuel/kruegg/ZOLA/Plate_2/alignments/ZOLAv0/bam/zola-2.12-050.bam

I am going to put those in a file called `20_birds.txt` and then I will write script that will call it and merge the bam from it. That script is: `species-level-scripts.sh/06-merge-bams-subset.sh`

Let's fire it off:

``` sh
[kruegg@login1 ZOLA]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA
[kruegg@login1 ZOLA]$ date
Wed Nov 16 06:10:02 PST 2016
[kruegg@login1 ZOLA]$ qsub ~/genoscape-bioinformatics/species-level-scripts.sh/06-merge-bams-subset.sh twenty-birds-ZOLAv0 20-birds.txt 
JSV: PE=shared
Your job 1124228 ("merge-bam-subset") has been submitted
```

Index the genome fasta and create a dictionary
----------------------------------------------

Should be straightforward. Now that I can finally get an interactive shell, we will use that.

``` sh
# get a shell on a compute node (takes a few moments)
[kruegg@login3 ~]$ qrsh

# change directory to the one holding the genome
[kruegg@n9894 ~]$ cd nobackup-klohmuel/ZOLA/Genome/
[kruegg@n9894 Genome]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Genome
[kruegg@n9894 Genome]$ ls
GCA_001281735.1_ASM128173v1_genomic.fna

# make a .fai with samtools faidx
[kruegg@n9894 Genome]$ module load samtools
[kruegg@n9894 Genome]$ samtools faidx GCA_001281735.1_ASM128173v1_genomic.fna 

# see what that did:
[kruegg@n9894 Genome]$ ls
GCA_001281735.1_ASM128173v1_genomic.fna
GCA_001281735.1_ASM128173v1_genomic.fna.fai

# now make a dictionary with picard tools
# get the path to Picard tools:
[kruegg@n9894 Genome]$ source ~/genoscape-bioinformatics/program-defs.sh 
# see that we have it:
[kruegg@n9894 Genome]$ echo $PICARD_JAR 
/u/nobackup/klohmuel/kruegg/bin/picard.jar
# then run the command
[kruegg@n9894 Genome]$ module load java
[kruegg@n9894 Genome]$ java -jar $PICARD_JAR CreateSequenceDictionary R=GCA_001281735.1_ASM128173v1_genomic.fna O=GCA_001281735.1_ASM128173v1_genomic.dict

# now see what we have:
[kruegg@n9894 Genome]$ wc *
      2934      14668     453834 GCA_001281735.1_ASM128173v1_genomic.dict
  12954441   12977905 1049217754 GCA_001281735.1_ASM128173v1_genomic.fna
      2933      14665     107662 GCA_001281735.1_ASM128173v1_genomic.fna.fai
```

Small trial of SNPing
---------------------

Just running through a small piece:

``` sh
# get just a small piece
[kruegg@n9894 MergedBams]$ samtools view -bh -o small_seg.bam  twenty-birds-ZOLAv0-merged.bam LAII01000002.1

# index it
[kruegg@n9894 MergedBams]$ samtools index small_seg.bam 

# then try launching with no nct
[kruegg@n9894 MergedBams]$ module load java

[kruegg@n9894 SNPs]$ source ~/genoscape-bioinformatics/program-defs.sh 

# see if we have the address:
[kruegg@n9894 SNPs]$ echo $GATK_JAR 
/u/nobackup/klohmuel/kruegg/bin/GenomeAnalysisTK.jar

# then fire it off:
[kruegg@n9894 SNPs]$ java -Xmx8G -jar $GATK_JAR -T HaplotypeCaller -R ../Genome/GCA_001281735.1_ASM128173v1_genomic.fna -I ../MergedBams/small_seg.bam -stand_call_conf 20.0 -stand_emit_conf 20.0 -o small_seg.vcf --genotyping_mode DISCOVERY -nct 8 -L LAII01000002.1
Error: Unable to access jarfile HaplotypeCaller 
```

Going for the full SNP calling at 20 individuals
------------------------------------------------

Made a script for this and from my previous experiment, I think this should take about 100 hours. Gonna give it a whirl:

``` sh
# check out the syntax:
[kruegg@login1 SNPs]$ ~/genoscape-bioinformatics/species-level-scripts.sh/07-call-snps.sh 
Syntax:
  07-call-snps.sh  FASTA  BAM  OUTVCF

      FASTA:  The path to the fasta/fna file that holds the genome and the .fai and the .dict
          file.
      BAM  The path to the big merged bam file that is all indexed, etc.
      OUTVCF The path desired for the output file.  This should have a .vcf extension.
      
      This will likely take a long time....
      


# where are we:
[kruegg@login1 SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/SNPs
[kruegg@login1 SNPs]$ date
Wed Nov 16 12:35:17 PST 2016

# queue it:
[kruegg@login1 SNPs]$ qsub ~/genoscape-bioinformatics/species-level-scripts.sh/07-call-snps.sh ../Genome/GCA_001281735.1_ASM128173v1_genomic.fna ../MergedBams/twenty-birds-ZOLAv0-merged.bam twenty-birds-ZOLAv0.vcf
JSV: PE=shared
Your job 1125471 ("call-snps") has been submitted
```

Trying to do SNP calling as a job array
---------------------------------------

None of the scaffolds are very large in the ZOLA genome. The biggest is about 16 Mb, I think. That would take about 15 hours for 20 birds to do the SNP calling, it seems, from my trial runs. I am thinking that maybe I should just send this as a UGE job array for groups of scaffolds up to 1 Mb, or so, at a time. I have myself a shell, so I am going to just investigate a few things.

### First, try running a scaffold that is about 100K bp

``` sh
[kruegg@n2240 MergedBams]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/MergedBams
[kruegg@n2240 MergedBams]$ module load samtools

# here we sort the scaffolds, descending on length:
[kruegg@n2240 MergedBams]$ samtools view -H twenty-birds-ZOLAv0-merged.bam | awk '/SQ/' | sed 's/SN://g; s/LN://g;' | sort -n -b -r -k 3 | head 
@SQ LAII01000001.1  15146312
@SQ LAII01000002.1  11710427
@SQ LAII01000003.1  10855354
@SQ LAII01000004.1  10759109
@SQ LAII01000005.1  10442077
@SQ LAII01000006.1  10369376
@SQ LAII01000007.1  10101451
@SQ LAII01000008.1  9963541
@SQ LAII01000009.1  9896429
@SQ LAII01000010.1  9679279

# if we go down that a ways, we can pick out
# ont that is about 100K bp
@SQ LAII01000629.1  100961

# Do some standard steps to get some vars, etc:
[kruegg@n2240 MergedBams]$ source ~/genoscape-bioinformatics/program-defs.sh
[kruegg@n2240 MergedBams]$ source $MODULE_SOURCE
[kruegg@n2240 MergedBams]$ module load java

# and now change to a new directory 
[kruegg@n2240 try_array_SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/try_array_SNPs

# I want to launch the haplotyper and see if I can redirect stderr and stdout OK 
# on it.
[kruegg@n2240 try_array_SNPs]$ FASTA=../Genome/GCA_001281735.1_ASM128173v1_genomic.fna
[kruegg@n2240 try_array_SNPs]$ BAM=../MergedBams/twenty-birds-ZOLAv0-merged.bam
[kruegg@n2240 try_array_SNPs]$ SEG=LAII01000629.1
[kruegg@n2240 try_array_SNPs]$ OUT=trial1

# then this does it and captures output from stderr in a file:
java -Xmx2G -jar $GATK_JAR -T HaplotypeCaller \
-R $FASTA \
-I $BAM \
-stand_call_conf 20.0 -stand_emit_conf 20.0 \
-o $OUT.vcf --genotyping_mode DISCOVERY \
-L $SEG  > $OUT.stdout 2>$OUT.stderr
```

OK, that worked fine and took about 4 minutes for a 100 Kb segment.

Let's see if it would work for 1 million bp on lots of segments

``` sh
[kruegg@n2240 try_array_SNPs]$ OUT=trial2
[kruegg@n2240 try_array_SNPs]$ SEG="-L LAII01001698.1 -L LAII01001699.1 -L LAII01001700.1 -L LAII01001701.1 -L LAII01001702.1 -L LAII01001703.1 -L LAII01001704.1 -L LAII01001705.1 -L LAII01001706.1 -L LAII01001707.1 -L LAII01001708.1 -L LAII01001709.1 -L LAII01001710.1 -L LAII01001711.1 -L LAII01001712.1 -L LAII01001713.1 -L LAII01001714.1 -L LAII01001715.1 -L LAII01001716.1 -L LAII01001717.1 -L LAII01001718.1 -L LAII01001719.1 -L LAII01001720.1 -L LAII01001721.1 -L LAII01001722.1 -L LAII01001723.1 -L LAII01001724.1 -L LAII01001725.1 -L LAII01001726.1 -L LAII01001727.1 -L LAII01001728.1 -L LAII01001729.1 -L LAII01001730.1 -L LAII01001731.1 -L LAII01001732.1 -L LAII01001733.1 -L LAII01001734.1 -L LAII01001735.1 -L LAII01001736.1 -L LAII01001737.1 -L LAII01001738.1 -L LAII01001739.1 -L LAII01001740.1 -L LAII01001741.1 -L LAII01001742.1 -L LAII01001743.1 -L LAII01001744.1 -L LAII01001745.1 -L LAII01001746.1 -L LAII01001747.1 -L LAII01001748.1 -L LAII01001749.1 -L LAII01001750.1 -L LAII01001751.1 -L LAII01001752.1 -L LAII01001753.1 -L LAII01001754.1 -L LAII01001755.1 -L LAII01001756.1 -L LAII01001757.1 -L LAII01001758.1 -L LAII01001759.1 -L LAII01001760.1 -L LAII01001761.1 -L LAII01001762.1 -L LAII01001763.1 -L LAII01001764.1 -L LAII01001765.1 -L LAII01001766.1 -L LAII01001767.1 -L LAII01001768.1 -L LAII01001769.1 -L LAII01001770.1 -L LAII01001771.1 -L LAII01001772.1 -L LAII01001773.1 -L LAII01001774.1 -L LAII01001775.1 -L LAII01001776.1 -L LAII01001777.1 -L LAII01001778.1 -L LAII01001779.1 -L LAII01001780.1 -L LAII01001781.1 -L LAII01001782.1 -L LAII01001783.1 -L LAII01001784.1 -L LAII01001785.1 -L LAII01001786.1 -L LAII01001787.1 -L LAII01001788.1 -L LAII01001789.1 -L LAII01001790.1 -L LAII01001791.1 -L LAII01001792.1 -L LAII01001793.1 -L LAII01001794.1 -L LAII01001795.1 -L LAII01001796.1 -L LAII01001797.1 -L LAII01001798.1 -L LAII01001799.1 -L LAII01001800.1 -L LAII01001801.1 -L LAII01001802.1 -L LAII01001803.1 -L LAII01001804.1 -L LAII01001805.1 -L LAII01001806.1 -L LAII01001807.1 -L LAII01001808.1 -L LAII01001809.1 -L LAII01001810.1 -L LAII01001811.1 -L LAII01001812.1 -L LAII01001813.1 -L LAII01001814.1 -L LAII01001815.1 -L LAII01001816.1 -L LAII01001817.1 -L LAII01001818.1 -L LAII01001819.1 -L LAII01001820.1 -L LAII01001821.1 -L LAII01001822.1 -L LAII01001823.1 -L LAII01001824.1 -L LAII01001825.1 -L LAII01001826.1 -L LAII01001827.1 -L LAII01001828.1 -L LAII01001829.1 -L LAII01001830.1 -L LAII01001831.1 -L LAII01001832.1 -L LAII01001833.1 -L LAII01001834.1 -L LAII01001835.1 -L LAII01001836.1 -L LAII01001837.1 -L LAII01001838.1 -L LAII01001839.1 -L LAII01001840.1 -L LAII01001841.1 -L LAII01001842.1 -L LAII01001843.1 -L LAII01001844.1 -L LAII01001845.1 -L LAII01001846.1 -L LAII01001847.1 -L LAII01001848.1 -L LAII01001849.1 -L LAII01001850.1 -L LAII01001851.1 -L LAII01001852.1 -L LAII01001853.1 -L LAII01001854.1 -L LAII01001855.1 -L LAII01001856.1 -L LAII01001857.1 -L LAII01001858.1 -L LAII01001859.1 -L LAII01001860.1 -L LAII01001861.1 -L LAII01001862.1 -L LAII01001863.1 -L LAII01001864.1 -L LAII01001865.1 -L LAII01001866.1 -L LAII01001867.1 -L LAII01001868.1 -L LAII01001869.1 -L LAII01001870.1 -L LAII01001871.1 -L LAII01001872.1 -L LAII01001873.1 -L LAII01001874.1 -L LAII01001875.1 -L LAII01001876.1 -L LAII01001877.1 -L LAII01001878.1 -L LAII01001879.1 -L LAII01001880.1 -L LAII01001881.1 -L LAII01001882.1 -L LAII01001883.1 -L LAII01001884.1 -L LAII01001885.1 -L LAII01001886.1 -L LAII01001887.1 -L LAII01001888.1 -L LAII01001889.1 -L LAII01001890.1 -L LAII01001891.1 -L LAII01001892.1 -L LAII01001893.1 -L LAII01001894.1 -L LAII01001895.1 -L LAII01001896.1 -L LAII01001897.1 -L LAII01001898.1 -L LAII01001899.1 -L LAII01001900.1 -L LAII01001901.1 -L LAII01001902.1 -L LAII01001903.1 -L LAII01001904.1 -L LAII01001905.1 -L LAII01001906.1 -L LAII01001907.1 -L LAII01001908.1 -L LAII01001909.1 -L LAII01001910.1 -L LAII01001911.1 -L LAII01001912.1 -L LAII01001913.1 -L LAII01001914.1 -L LAII01001915.1 -L LAII01001916.1 -L LAII01001917.1 -L LAII01001918.1 -L LAII01001919.1 -L LAII01001920.1 -L LAII01001921.1 -L LAII01001922.1 -L LAII01001923.1 -L LAII01001924.1 -L LAII01001925.1 -L LAII01001926.1 -L LAII01001927.1 -L LAII01001928.1 -L LAII01001929.1 -L LAII01001930.1 -L LAII01001931.1 -L LAII01001932.1 -L LAII01001933.1 -L LAII01001934.1 -L LAII01001935.1 -L LAII01001936.1 -L LAII01001937.1 -L LAII01001938.1 -L LAII01001939.1 -L LAII01001940.1 -L LAII01001941.1 -L LAII01001942.1 -L LAII01001943.1 -L LAII01001944.1 -L LAII01001945.1 -L LAII01001946.1 -L LAII01001947.1 -L LAII01001948.1 -L LAII01001949.1 -L LAII01001950.1 -L LAII01001951.1 -L LAII01001952.1 -L LAII01001953.1 -L LAII01001954.1 -L LAII01001955.1 -L LAII01001956.1 -L LAII01001957.1 -L LAII01001958.1 -L LAII01001959.1 -L LAII01001960.1 -L LAII01001961.1 -L LAII01001962.1 -L LAII01001963.1 -L LAII01001964.1 -L LAII01001965.1 -L LAII01001966.1 -L LAII01001967.1 -L LAII01001968.1 -L LAII01001969.1 -L LAII01001970.1 -L LAII01001971.1 -L LAII01001972.1"

java -Xmx2G -jar $GATK_JAR -T HaplotypeCaller \
-R $FASTA \
-I $BAM \
-stand_call_conf 20.0 -stand_emit_conf 20.0 \
-o $OUT.vcf --genotyping_mode DISCOVERY \
 $SEG  > $OUT.stdout 2>$OUT.stderr
```

That actually works beautifully. So, now I just need to test the job array methodology. For that I will try 4 segments, each of about 200 Kbp.

So, let's say we have file called `comm_lines.txt` in the directory we are working in and it looks like this:

    1   207798   -L LAII01001011.1 -L LAII01001012.1 -L LAII01001013.1 -L LAII01001014.1 -L LAII01001015.1 -L LAII01001016.1 -L LAII01001017.1 -L LAII01001018.1 -L LAII01001019.1 -L LAII01001020.1 -L LAII01001021.1   19692 16201 19120 20605 16349 19576 19299 18801 19184 19986 18985
    2   204955   -L LAII01001022.1 -L LAII01001023.1 -L LAII01001024.1 -L LAII01001025.1 -L LAII01001026.1 -L LAII01001027.1 -L LAII01001028.1 -L LAII01001029.1 -L LAII01001030.1 -L LAII01001031.1 -L LAII01001032.1   19065 18681 18547 18817 19616 18835 18592 18325 18279 18171 18027
    3   202681   -L LAII01001033.1 -L LAII01001034.1 -L LAII01001035.1 -L LAII01001036.1 -L LAII01001037.1 -L LAII01001038.1 -L LAII01001039.1 -L LAII01001040.1 -L LAII01001041.1 -L LAII01001042.1 -L LAII01001043.1   18422 18486 18273 18313 24293 18123 18235 17937 18051 18187 14361
    4   209832   -L LAII01001044.1 -L LAII01001045.1 -L LAII01001046.1 -L LAII01001047.1 -L LAII01001048.1 -L LAII01001049.1 -L LAII01001050.1 -L LAII01001051.1 -L LAII01001052.1 -L LAII01001053.1 -L LAII01001054.1 -L LAII01001055.1     18029 15285 17730 17961 17841 21095 16944 18085 17573 17232 17225 14832
    5   206697   -L LAII01001095.1 -L LAII01001096.1 -L LAII01001097.1 -L LAII01001098.1 -L LAII01001099.1 -L LAII01001100.1 -L LAII01001101.1 -L LAII01001102.1 -L LAII01001103.1 -L LAII01001104.1 -L LAII01001105.1 -L LAII01001106.1 -L LAII01001107.1   15714 15783 16788 11951 15148 15774 14922 15693 15662 15367 23573 14772 15550
    6   209666   -L LAII01001108.1 -L LAII01001109.1 -L LAII01001110.1 -L LAII01001111.1 -L LAII01001112.1 -L LAII01001113.1 -L LAII01001114.1 -L LAII01001115.1 -L LAII01001116.1 -L LAII01001117.1 -L LAII01001118.1 -L LAII01001119.1 -L LAII01001120.1 -L LAII01001121.1     16653 15318 15285 13386 15237 14992 16778 15249 15074 15149 15090 14871 11478 15106
    7   200137   -L LAII01001122.1 -L LAII01001123.1 -L LAII01001124.1 -L LAII01001125.1 -L LAII01001126.1 -L LAII01001127.1 -L LAII01001128.1 -L LAII01001129.1 -L LAII01001130.1 -L LAII01001131.1 -L LAII01001132.1 -L LAII01001133.1 -L LAII01001134.1 -L LAII01001135.1     14983 14692 15045 11443 15078 15071 14862 15683 14712 14710 14708 13726 11197 14227
    8   207899   -L LAII01001136.1 -L LAII01001137.1 -L LAII01001138.1 -L LAII01001139.1 -L LAII01001140.1 -L LAII01001141.1 -L LAII01001142.1 -L LAII01001143.1 -L LAII01001144.1 -L LAII01001145.1 -L LAII01001146.1 -L LAII01001147.1 -L LAII01001148.1 -L LAII01001149.1 -L LAII01001150.1   15711 14820 14525 14229 14403 14284 13263 14511 14171 12651 14300 14092 14086 8881 13972
    9   202906   -L LAII01001151.1 -L LAII01001152.1 -L LAII01001153.1 -L LAII01001154.1 -L LAII01001155.1 -L LAII01001156.1 -L LAII01001157.1 -L LAII01001158.1 -L LAII01001159.1 -L LAII01001160.1 -L LAII01001161.1 -L LAII01001162.1 -L LAII01001163.1 -L LAII01001164.1 -L LAII01001165.1   8414 13869 17095 13709 9493 13679 12573 19121 13537 13564 13612 13525 13527 13402 13786
    10  204320   -L LAII01001166.1 -L LAII01001167.1 -L LAII01001168.1 -L LAII01001169.1 -L LAII01001170.1 -L LAII01001171.1 -L LAII01001172.1 -L LAII01001173.1 -L LAII01001174.1 -L LAII01001175.1 -L LAII01001176.1 -L LAII01001177.1 -L LAII01001178.1 -L LAII01001179.1 -L LAII01001180.1   13394 13312 13116 13087 13298 13346 13318 12917 17542 12729 12771 12973 16675 12890 12952
    11  205298   -L LAII01001181.1 -L LAII01001182.1 -L LAII01001183.1 -L LAII01001184.1 -L LAII01001185.1 -L LAII01001186.1 -L LAII01001187.1 -L LAII01001188.1 -L LAII01001189.1 -L LAII01001190.1 -L LAII01001191.1 -L LAII01001192.1 -L LAII01001193.1 -L LAII01001194.1 -L LAII01001195.1 -L LAII01001196.1     12847 13854 12733 12640 13726 12776 13419 12867 12679 12497 12640 12574 12481 12528 12579 12458

Each of those should take about five minutes to run, so we will give them each 15 minutes.

Now we need to make a UGE job array script to do it.

### Submitting as a job array

I made a test job array script and did this:

``` sh
[kruegg@login3 try_array_SNPs]$ qsub ~/genoscape-bioinformatics/species-level-scripts.sh/07-snps-via-jobarray.sh 
Your job-array 1128865.1-11:1 ("snp-array") has been submitted
```

It turns out that failed when it was run on any of the `pod` nodes. But did fine on the `c2` nodes. I tried forcing it to run only on the `c2` nodes:

``` sh
[kruegg@login3 try_array_SNPs]$ qsub  -q c2_smp.q ~/genoscape-bioinformatics/species-level-scripts.sh/07-snps-via-jobarray.sh 
You specified queue name = c2_smp.q : Usually this is not recommended; if the queue name conflicts with other job parameters, the job may not start. If you have questions, submit them to: support.idre.ucla.edu/helpdesk
Your job-array 1128941.1-11:1 ("snp-array") has been submitted
```

And it said that might not be recommended.

But it totally worked great. So I am going to change that to an array of 187 jobs that will finish off all the SNP calling, and see how it goes...

``` sh
[kruegg@login3 try_array_SNPs]$ qsub  -q c2_smp.q ~/genoscape-bioinformatics/species-level-scripts.sh/07-snps-via-jobarray.sh 
You specified queue name = c2_smp.q : Usually this is not recommended; if the queue name conflicts with other job parameters, the job may not start. If you have questions, submit them to: support.idre.ucla.edu/helpdesk
Your job-array 1128944.1-187:1 ("snp-array") has been submitted
[kruegg@login3 try_array_SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/try_array_SNPs
[kruegg@login3 try_array_SNPs]$ date
Thu Nov 17 04:59:33 PST 2016
```

Those are chunks of about 4 megabases each.

This is working. Woo-hoo! Once it is done I sort everything back together into a single VCF:

``` sh
[kruegg@n2239 try_array_SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/try_array_SNPs
[kruegg@n2239 try_array_SNPs]$ source ~/genoscape-bioinformatics/program-defs.sh 
[kruegg@n2239 try_array_SNPs]$ INPUTS=$(ls -l 0*.vcf | awk '{printf("I=%s ", $NF)}') 
[kruegg@n2239 try_array_SNPs]$ module load java
[kruegg@n2239 try_array_SNPs]$ java -jar $PICARD_JAR SortVcf $INPUTS O=twenty-birds-ZOLAv0.vcf

# that just took a couple of minutes.
```

Initial filtering of the VCF
----------------------------

It looks like there are 1,553,483 variants in that file. We are going to filter them like this:

1.  No indels
2.  Biallelic only
3.  Minor allele count &gt;= 2
4.  minimum genotype quality = 30
5.  minimum depth = 8
6.  called in at least 50% of indivs

We are dropping singletons amongst our 40 gene copies here.

``` sh
[kruegg@n2239 try_array_SNPs]$ module load vcftools
# move things to new directory:
[kruegg@n2239 try_array_SNPs]$ mv twenty-birds-ZOLAv0.vcf* ../SNPs/small_seg.vcf
small_seg.vcf      small_seg.vcf.idx  
[kruegg@n2239 try_array_SNPs]$ mv twenty-birds-ZOLAv0.vcf* ../SNPs/
[kruegg@n2239 try_array_SNPs]$ cd ../SNPs/
[kruegg@n2239 SNPs]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/SNPs
[kruegg@n2239 SNPs]$ vcftools --vcf twenty-birds-ZOLAv0.vcf --out zola-twenty-filtered --remove-indels --min-alleles 2 --max-alleles 2 --mac 2 --minGQ 30 --minDP 8 --max-missing 0.5 --recode 

VCFtools - 0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
    --vcf twenty-birds-ZOLAv0.vcf
    --mac 2
    --max-alleles 2
    --min-alleles 2
    --minDP 8
    --minGQ 30
    --max-missing 0.5
    --out zola-twenty-filtered
    --recode
    --remove-indels

After filtering, kept 20 out of 20 Individuals
Outputting VCF file...
After filtering, kept 364429 out of a possible 1553483 Sites
Run Time = 88.00 seconds

# now we also want to make the 012 file from that
[kruegg@n2239 SNPs]$ vcftools --vcf  zola-twenty-filtered.recode.vcf --out zola-twenty-filt --012

VCFtools - 0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
    --vcf zola-twenty-filtered.recode.vcf
    --012
    --out zola-twenty-filt

After filtering, kept 20 out of 20 Individuals
Writing 012 matrix files ... Done.
After filtering, kept 364429 out of a possible 364429 Sites
Run Time = 8.00 seconds
[kruegg@n2239 SNPs]$ gzip -9 zola-twenty-filt.012
```

### Getting a list of 1 Mb collections of scaffolds

Here is a little script that creates a tab-delimited file of `-L` commands to be used to pick out collections of scaffolds whose summed lengths try to be just a bit over 1 million base pairs:
