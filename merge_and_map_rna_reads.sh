#!/bin/sh

# $REF should be set
if [ -z ${REF+x} ]; then
  echo "REF is unset"
  REF=$HOME/Documents/referentie/hisat/grch37_tran/genome_tran
  echo "$REF used"  
fi

for i in $(ls *.gz | cut -d "_" -f 1 | sort | uniq) ; do
  R1=$(ls ${i}_*R1*.gz)
  R2=$(ls ${i}_*R2*.gz)
  echo $R1 $R2
  cope -a $R1 -b $R2 -o $i -m 0 -s 33 -L 0 -u 140
  
  cat $i.connect.fq $i.unConnect_1.fq > $i.connect_and_singleR1.fastq
  cat $i.connect.fq $i.unConnect_2.fq > $i.connect_and_singleR2.fastq  
  
  hisat2 -x $REF -U $i.connect_and_singleR1.fastq --threads 2 |\
  samtools view -Shu - |\
  samtools sort -T $i.tmp -O bam - > $i.connect_and_singleR1.bam
  samtools index $i.connect_and_singleR1.bam

  hisat2 -x $REF -U $i.connect_and_singleR2.fastq --threads 2 |\
  samtools view -Shu - |\
  samtools sort -T $i.tmp -O bam - > $i.connect_and_singleR2.bam
  samtools index $i.connect_and_singleR2.bam

done ;

