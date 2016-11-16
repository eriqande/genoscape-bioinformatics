#!/bin/bash
#$ -cwd
#$ -V
#$ -N call-snps
#$ -o call-snps.log
#$ -e call-snps.error
#$ -pe shared 16
#$ -l highp,h_data=2G,time=200:00:00
#$ -M eric.anderson@noaa.gov
#$ -m bea

source ~/genoscape-bioinformatics/program-defs.sh
source $MODULE_SOURCE

module load java



function usage {
      echo Syntax:
      echo "  $(basename $0)  FASTA  BAM  OUTVCF"
      echo "
      FASTA:  The path to the fasta/fna file that holds the genome and the .fai and the .dict
          file.
      BAM  The path to the big merged bam file that is all indexed, etc.
      OUTVCF The path desired for the output file.  This should have a .vcf extension.
      
      This will likely take a long time....
      
"
}

if [ $# -eq 0 ]; then
    usage;
    exit 1;
fi;

while getopts ":h" opt; do
    case $opt in
	h    ) 
	    usage
	    exit 1
	    ;;
	#m    )  VAR=$OPTARG;
	#    ;;
	\?   )
	    usage
	    exit  1
	    ;;
    esac
done

shift $((OPTIND-1));


# uncomment to test for right number of required args
if [ $# -ne 3 ]; then
    usage;
    exit 1;
fi


java -Xmx8G -jar $GATK_JAR -T HaplotypeCaller \
  -R $FASTA \
  -I $BAM \
  -stand_call_conf 20.0 -stand_emit_conf 20.0 \
  -o $OUTVCF --genotyping_mode DISCOVERY \
  -nct 16
  
  
  