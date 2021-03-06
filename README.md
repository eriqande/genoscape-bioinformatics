genoscape-bioinformatics
================
06 November, 2016

-   [TODO](#todo)
-   [Introduction](#introduction)
-   [Installation](#installation)
-   [Directory Structure](#directory-structure)
-   [Plate Processing Scripts](#plate-processing-scripts)
-   [Reference Genome DB Scripts](#reference-genome-db-scripts)
-   [Running the scripts](#running-the-scripts)
-   [Eric's notes to himself](#erics-notes-to-himself)

<!-- README.md is generated from README.Rmd. Please edit that file -->
TODO
----

Eric's list of things to do:

1.  Make it resilient to periods in the bird IDs.

Introduction
------------

This is Eric Anderson's attempt to pull together into a single, relatively-immutable set of scripts, the bioinformatic pipeline to be run on birds as part of the Genoscape Project. All the scripts were originally written and figured out by Rachael Bay, and Eric has just tried to organize things so that they are a little bit more reusable.

The basic idea is that you can clone this repository to your home directory in your computing environment (e.g. $HOME on the hoffman cluster at UCLA) and then have all the scripts that you will need. However, you do need to set some filepaths and options. Rather than changing those in the scripts themselves (which makes it really hard to version control them, etc.) we have set it up so that the scripts are designed to be launched with the current working directory being the directory that holds the `rawdata` directory for a particular Plate of data, and within that directory there must be a file called `data-defs.sh` that includes variable definitions (more on that later.)

Installation
------------

Simple. Clone the repository and then check the paths to the programs and modify them as necessary.

1.  Clone the repository into your home directory:

    ``` sh
    cd $HOME
    git clone https://github.com/eriqande/genoscape-bioinformatics
    ```

2.  Check out the file `program-defs.sh` in the top level of the directory. This is an example of how you would specify paths for all the different programs that are to be used (but which don't live in this repository, itself). (Of course, some, like `awk` and `perl` are just assumed to be on your $PATH already.). The file as it arrives with the repo looks like this:


        # this is the path to the file to source to get the module command
        MODULE_SOURCE=/u/local/Modules/default/init/modules.sh

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

    Inspect those paths (i.e. for `fastqc`, Stacks' `process_radtags`, and `clone_filter`, etc) and make sure that they are correct for your system. If you are working on Hoffman with Kristen, then these paths should be correct.

However, if you have working on a different system with different paths you will have to make a new file. You can put it anywhere you want to. In the repository is fine. Name it something else, and then change the paths as needed.

As you will see later in the `data-defs.sh` you have to specify which file should be read for the program definitions.

Directory Structure
-------------------

In order for this to work, we have to impose a specific file structure for each species. This seems like good practice anyway. Here is how we will do it:

1.  Each species gets its own directory on `nobackup-klohmuel`. For example: `ZOLA` at `/u/nobackup/klohmuel/kruegg/ZOLA` for *Zosterops lateralis*.
2.  Each **plate of data** from that species will get its own directory named `Plate_1` or `Plate_2`, etc., at the top level of the species directory.
3.  Inside the `Plate_X` directory must be:
    1.  a `rawdata` directory that holds the raw-data. These should be gzipped fastq files; one for read 1 and one for read 2. (Expect these to be about 30 to 40 Gb)
    2.  A two-column, tab-delimited, barcode file that can be named anything, but it ought to be something like `ZOLA-1-barcodes-2col.txt` which signifies it is for the species ZOLA, plate 1, and it has two columns. The first column is the barcode and the second is the ID of the sample. There are no headers on the columns. **Make sure that there are no spaces in any of the IDs** of the birds. This file should have Unix line endings, which means that you probably can't create it direcly from Excel (which gives you MacOS line endings.). As an example, here is what `ZOLA-1-barcodes-2col.txt` looks like. The first ten lines of the file are:

            ACAAGCTA    LH130
            AAACATCG    CDH19
            ACATTGGC    CDH49
            ACCACTGT    N110
            AACGTGAT    APB72
            CGCTGATC    CDH50
            CAGATCTG    LH104
            ATGCCTAA    CDH54
            AACGAACG    GT132
            AGTACAAG    L549

        It has 96 lines it. One line for each bird on the plate. \*\*Note: As currently written...the IDs can't have periods in them. I should fix that...
    3.  a file called `data-defs.sh` that holds all the variables that will be sent to the scripts. This file is merely a shell script that defines variables that will be used in the scripts. At the current time there are only 5 things that need to be set, and it is pretty self-explanatory. The definitions will grow, but for now they are as shown in the `data-defs.sh` file for Plate 1 of the Silvereye data. Note that any comments can be put in there following a `#`. Note that you **cannot have any spaces around the equals signs** when you make this file.

        ``` sh
        PROGDEFS=/u/home/k/kruegg/genoscape-bioinformatics/program-defs.sh  # file that defines the paths. 
        READ1=rawdata/Silvereye-Plate1_S36_L005_R1_001.fastq.gz # fastq.gz file for Read one.
        READ2=rawdata/Silvereye-Plate1_S36_L005_R2_001.fastq.gz # paths relative to the Plate_1 directory.
        BC2COL=ZOLA-1-barcodes-2col.txt # The two column barcode file
        SHORTNAME=zola-1   # species dash plate number is good.  This is used for output later  
        GENOME_DB=/u/home/k/kruegg/nobackup-klohmuel/References/bowtie2-ZOLAv0/ZOLAv0   # must be the ABSOLUTE path AND the prefix
        ```

        I want to say a little more about these:
        -   `PROGDEFS` should be given as an absolute path.
        -   `READ1` and `READ2` can be given as paths relative to the `Plate_X` directory.
        -   `BC2COL` should be the name of the file and it should be in the `Plate_X` directory.
        -   `SHORTNAME` must have no spaces in it, and no special characters that don't work well in file names.
        -   `GENOME_DB` gives the path and prefix for the indexed genome that you want to align to. This is something that you will create using, for example, `reference-genome-db-scripts/bowtie2-build-genome-database.sh`. I ran that in the directory `/u/home/k/kruegg/nobackup-klohmuel/References` and told it to name the data base `ZOLAv0`, so it put everything into `/u/home/k/kruegg/nobackup-klohmuel/References/bowtie2-ZOLAv0`. In particular, there are a lot of files like this:

                ZOLAv0.1.bt2
                ZOLAv0.4.bt2
                ZOLAv0.2.bt2
                ZOLAv0.rev.1.bt2
                ZOLAv0.3.bt2
                ZOLAv0.rev.2.bt2

            Notice that their prefix is `ZOLAv0`. That has to be specified as part of the `GENOME_DB` variable. That is what the `ZOLAv0` is doing at the end of the `/u/home/k/kruegg/nobackup-klohmuel/References/bowtie2-ZOLAv0/ZOLAv0`. Note also, `GENOME_DB` must be an absolute path!! Don't specify it like `../References/bowtie2-ZOLAv0/ZOLAv`.

Here is an example of what the directory structure looks like once we have downloaded two plates of ZOLA data:

    |-Plate_1:
    |----data-defs.sh
    |----ZOLA-1-barcodes-2col.txt
    |----rawdata:
    |-------Silvereye-Plate1_S36_L005_R1_001.fastq.gz
    |-------Silvereye-Plate1_S36_L005_R2_001.fastq.gz
    |-Plate_2:
    |----data-defs.sh
    |----ZOLA-2-barcodes-2col.txt
    |----rawdata:
    |-------Silvereye-Plate2_S86_L003_R1_001.fastq.gz
    |-------Silvereye-Plate2_S86_L003_R2_001.fastq.gz

Plate Processing Scripts
------------------------

The plate processing scripts in this repository (i.e. those in `plate-processing-scripts`) all operate on data at the Plate level which is why we have set the directory structure up as we have. Each script creates output in a new directory within the `Plate_X` directory in which it gets executed. *Now info about each script and its output here*

Reference Genome DB Scripts
---------------------------

If you have a genome. Then typically you are going to to need to index it or otherwise create a database from it in order to align sequence to it. This process is distinct from the process of actually assembling a genome for a species. Because you might want to align sequence from one species against the genome of a different species, it doesn't make a lot of sense to store the indexed genome data bases within the species directory itself. So, we are going to follow Rachael's lead here in recommending that the user maintain a directory called `References`, within which will be a series of directories that include the genome data bases.

We have a few different scripts that help the process of creating these genome data bases in a way that things are organized and it is clear what the original genome used was, etc.

Unlike the plate processing scripts, which use the `data-defs.sh` file to specify its input, these scripts take command line arguments.

These are all scripts that should be run inside the `References` directory that goes at the top level of the species directory.

Currently there is only one script: `reference-genome-db-scripts/bowtie2-build-genome-database.sh` It uses `bowtie2-build` which you can read about [here](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).

Here is an example of how it was used to create a DB for the *Zosterops* genome. Note that I had downloaded the *Zosterops* genome from NCBI, then I gunzipped it (that is important) and put it in the `Genome` directory inside the `ZOLA` directory. Note that you can use a relative path to the genome file, if desired. However, remember, always use **absolute paths** for things in the `data-defs.sh` file!

``` sh
[kruegg@login3 References]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/References
[kruegg@login3 References]$ qsub ~/genoscape-bioinformatics/reference-genome-db-scripts/bowtie2-build-genome-database.sh ../ZOLA/Genome/GCA_001281735.1_ASM128173v1_genomic.fna  ZOLAv0
```

This takes about half an hour and creates all the contents in the directory `/u/home/k/kruegg/nobackup-klohmuel/References/bowtie2-ZOLAv0`, including `README.txt` which tells us this:

    Created Thu Nov  3 03:56:37 PDT 2016,  by running  bowtie2-build-genome-database.sh script with arguments:
            FASTA: ../ZOLA/Genome/GCA_001281735.1_ASM128173v1_genomic.fna
            AbsolutePath: /u/nobackup/klohmuel/kruegg/ZOLA/Genome/GCA_001281735.1_ASM128173v1_genomic.fna
            OutputPrefix: ZOLAv0

Running the scripts
-------------------

Will put more in here later. But for now, it looks like:

``` sh
[kruegg@login2 Plate_1]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_1

[kruegg@login2 Plate_1]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh 
```

Eric's notes to himself
-----------------------

I can test this stuff locally on my laptop using data and files that I have put in: `/Users/eriq/Documents/UnsyncedData/ZOLA_100Kreads/Plate_1`. And making a symlink at `/u/nobackup/klohmuel/kruegg/bin` or at `/u/nobackup/klohmuel/kruegg/bin/stacks-1.32` that points where it needs to.

Here is a transcript of putting the fastQC and the rad\_processing into the job queue:

``` sh
[kruegg@login3 Plate_1]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh 
JSV: PE=shared
Your job 899364 ("fastqc") has been submitted
[kruegg@login3 Plate_1]$ cd ../Plate_2/
[kruegg@login3 Plate_2]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh 
JSV: PE=shared
Your job 899365 ("fastqc") has been submitted
[kruegg@login3 Plate_2]$ cd ../Plate_1/
[kruegg@login3 Plate_1]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/02-rad-process.sh 
JSV: PE=shared
Your job 899366 ("rad") has been submitted
[kruegg@login3 Plate_1]$ cd ../Plate_2/
[kruegg@login3 Plate_2]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/02-rad-process.sh 
JSV: PE=shared
Your job 899367 ("rad") has been submitted
[kruegg@login3 Plate_2]$ 
```
