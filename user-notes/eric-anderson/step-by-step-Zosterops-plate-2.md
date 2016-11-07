Step-by-step *Zosterops* Plate 2
================
07 November, 2016

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
-   [Flip, Trim, and demultiplex, for the Best-RAD process](#flip-trim-and-demultiplex-for-the-best-rad-process)
    -   [Start the job](#start-the-job)

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

Flip, Trim, and demultiplex, for the Best-RAD process
-----------------------------------------------------

*Prep Time: 5 minutes*
*Cook Time: *

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
