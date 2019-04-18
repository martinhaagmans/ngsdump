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

  echo $R1

  (hisat2 -x $REF -U $R1 --threads 2 |\
  samtools view -hu - |\
  samtools sort -o $i.R1only.sorted.bam) > $i.R1.hisat2.log 2>&1
  samtools index $i.R1only.sorted.bam
  
  echo $R2  

  (hisat2 -x $REF -U $R2 --threads 2 |\
  samtools view -hu - |\
  samtools sort -o $i.R2only.sorted.bam) > $i.R2.hisat2.log 2>&1
  samtools index $i.R2only.sorted.bam
  
done ;

