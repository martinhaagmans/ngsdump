#!/bin/sh
set -e

JAVA=/home/manager/programs/jre1.8.0_25/bin/java
TRIMMOMATIC=/home/manager/programs/Trimmomatic-0.36/trimmomatic-0.36.jar
HSQUTILS=/home/manager/programs/hsqutils_v1_1/hsqutils.jar
PROBES=/home/manager/Desktop/HEATSEQ/170928_HG19_CLL_AMCv2_HSQ_IRN200341732_probe_info.txt
REF=/home/manager/Documents/referentie/lifescope.hg19.fa

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r1|--read1)
    R1="$2"
    shift # past argument
    shift # past value
    ;;
    -r2|--read2)
    R2="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--sample)
    SAMPLE="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--probes)
    PROBES="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

$JAVA -jar $TRIMMOMATIC PE -phred33 $R1 $R2 $SAMPLE\_R1_quality_filtered.fastq $SAMPLE\_R1_unpaired.fastq $SAMPLE\_R2_quality_filtered.fastq $SAMPLE\_R2_unpaired.fastq TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50   

$JAVA -jar $HSQUTILS trim --r1 $SAMPLE\_R1_quality_filtered.fastq --r2 $SAMPLE\_R2_quality_filtered.fastq --probe $PROBES

bwa mem $REF -M trimmed_$SAMPLE\_R1_quality_filtered.fastq trimmed_$SAMPLE\_R2_quality_filtered.fastq | samtools view -Sb - > $SAMPLE\_initial_mapping.bam

$JAVA -jar $HSQUTILS dedup --r1 $SAMPLE\_R1_quality_filtered.fastq --r2 $SAMPLE\_R2_quality_filtered.fastq --probe $PROBES --inputBam $SAMPLE\_initial_mapping.bam --outputBamFileName $SAMPLE\_dedup.bam

mv HSQutils_dedup_summary.txt $SAMPLE\_dedup_summary.txt
