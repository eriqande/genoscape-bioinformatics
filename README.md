genoscape-bioinformatics
================

-   [Introduction](#introduction)
-   [Example Directory Structure](#example-directory-structure)
-   [Examples of what the files look like](#examples-of-what-the-files-look-like)
    -   [ZOLA-1-barcodes-2col.txt](#zola-1-barcodes-2col.txt)
    -   [The `data-defs.sh` file for ZOLA Plate\_1](#the-data-defs.sh-file-for-zola-plate_1)
-   [Running the scripts](#running-the-scripts)

<!-- README.md is generated from README.Rmd. Please edit that file -->
Introduction
------------

This is Eric Anderson's attempt to pull together into a single, relatively-immutable set of scripts, the bioinformatic pipeline to be run on birds as part of the Genoscape Project. All the scripts were originally written and figured out by Rachael Bay, and Eric has just tried to organize things so that they are a little bit more reusable.

The basic idea is that you can clone this repository to your home directory in your computing environment (e.g. $HOME on the hoffman cluster at UCLA) and then have all the scripts that you will need. However, you do need to set some filepaths and options. Rather than changing those in the scripts themselves (which makes it really hard to version control them, etc.) we have set it up so that the scripts are designed to be launched with the current working directory being the directory that holds the `rawdata` directory for a particular Plate of data, and within that directory there must be a file called `data-defs.sh` that includes variable definitions (more on that later.)

So, in order for this to work, we have to impose a specific file structure for each species. This seems like good practice anyway. Here is how we will do it:

1.  Each species gets its own directory on `nobackup-klohmuel`. For example: `ZOLA` at `/u/nobackup/klohmuel/kruegg/ZOLA` for *Zosterops lateralis*.
2.  Each plate of data from that species will get its own directory named `Plate_1` or `Plate_2`, etc., at the top level of the species directory.
3.  Inside the `Plate_X` directory must be:
    1.  a `rawdata` directory that holds the raw-data. These should be gzipped fastq files; one for read 1 and one for read 2. (Expect these to be about 30 to 40 Gb)
    2.  A two-column, tab-delimited, barcode file that can be named anything, but it ought to be something like `ZOLA-1-barcodes-2col.txt` which signifies it is for the species ZOLA, plate 1, and it has two columns. The first column is the barcode and the second is the ID of the sample. There are no headers on the columns. **Make sure that there are no spaces in any of the IDs** of the birds. This file should have Unix line endings, which means that you probably can't create it direcly from Excel (which gives you MacOS line endings.)
    3.  a file called `data-defs.sh` that holds all the variables that will be sent to the scripts

The plate processing scripts in this repository (i.e. those in `plate-processing-scripts`) all operate on data at the Plate level which is why we have set the directory structure up as we have. Each script creates output in a new directory within the `Plate_X` directory in which it gets executed.

Example Directory Structure
---------------------------

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

Examples of what the files look like
------------------------------------

### ZOLA-1-barcodes-2col.txt

The first ten lines of the file look like this:

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

It has 96 lines it. One line for each bird on the plate.

### The `data-defs.sh` file for ZOLA Plate\_1

This is merely a shell script that defines variables that will be used in the scripts. The definitions will grow, but for now they are as follows. Note that any comments can be put in there following a `#`. Note that you **cannot have any spaces around the equals signs** when you make this file.

``` sh
READ1=rawdata/Silvereye-Plate1_S36_L005_R1_001.fastq.gz # fastq.gz file for Read one.
READ2=rawdata/Silvereye-Plate1_S36_L005_R2_001.fastq.gz # paths relative to the Plate_1 directory.
BC2COL=ZOLA-1-barcodes-2col.txt # The two column barcode file
```

Running the scripts
-------------------

Will put more in here later. But for now, it looks like:

``` sh
[kruegg@login2 Plate_1]$ pwd
/u/home/k/kruegg/nobackup-klohmuel/ZOLA/Plate_1

[kruegg@login2 Plate_1]$ qsub ~/genoscape-bioinformatics/plate-processing-scripts/01-fastqc.sh 
```
