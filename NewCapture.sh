#!/bin/sh

if [ $# -eq 0 ]; then
    echo "Usage: $0 Prefix NewCap.txt " ;
    exit 1;
fi

FILE=${BASH_ARGV[0]}
out=${BASH_ARGV[1]}


#remove carriage return
sed -i 's/\r//gi' $FILE

#create exonplus20 BED if input is txt
if [[ ${FILE: -4} != ".bed" ]] ; then
    awk -F '\t' '{print $1 "\t" $2-1 "\t" $3}' $FILE |  sort -k1,1V > $out\_exonplus20.bed
fi

#cat file if input is BED
if [ ${FILE: -4} == ".bed" ] ; then
    cat $FILE > $out\_exonplus20.bed
fi

echo "Created exonplus20.bed" 

#create Picard-target
cat PicardHeader.txt > $out\_target.picard.interval_list
awk -F '\t' '{print $0 "\t+\tExonplus20"}' $out\_exonplus20.bed  >> $out\_target.picard.interval_list

#create Picard-tiled-header (append probes later from CD)
cat PicardHeader.txt > $out\_tiled.picard.interval_list

echo "Created Picard target and tiled files" 
echo "Querying database for annotation" 

#annotate exonplus20 bed with transcription start to make generegions.bed
cat $out\_exonplus20.bed |\
awk -F '\t' '{printf("select \"%s\",\"%s\",\"%s\", G.name,G.name2,G.txStart,G.txEnd,G.strand from refGene as G  where chrom=\"%s\" and  not(%s>txEnd or %s<txStart);\n",$1,$2,$3,$1,$2,$3);}' |\
mysql --user=genome --host=genome-mysql.cse.ucsc.edu  -A  -N -D hg19  |\
awk '{print}' > UCSC.out

echo "Done." 

awk '{print $1"\t"$2"\t"$3"\t"$5"\t"$8}' UCSC.out | sort | uniq | sort -k 1.4,1n -k 2,2n -k 3,3n > $out\_exonplus20.annotated
awk '{print $1"\t"$6"\t"$7"\t"$5}' UCSC.out | sort | uniq | sort -k 1.4,1n -k 2,2n -k 3,3n > $out\_generegions.bed
 

