#!/bin/sh

# $REF should be set
if [ -z ${REF+x} ]; then
  echo "REF is unset"
  exit
fi

for i in $(ls *.gz | cut -d "_" -f 1 | sort | uniq) ; do
  R1=$(ls $i\_*R1*.gz)
  R2=$(ls $i\_*R2*.gz)
  echo $R1 $R2
  bwa mem  -R "@RG\\tID:$i\\tLB:$i\\tPL:ILLUMINA\\tPU:$i\\tSM:$i" -t 1 -M $REF $R1 $R2 |\
  samtools view -Shu - |\
  samtools sort -T $i.tmp -O bam - > $i.sorted.bam
  samtools index $i.sorted.bam
done ;
